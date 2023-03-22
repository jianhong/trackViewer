/*
  The MIT License (MIT)

  Copyright (c) 2011-2016 Broad Institute, Aiden Lab

 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in
 all copies or substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 THE SOFTWARE.
*/
#include <cstring>
#include <fstream>
#include <sstream>
#include <map>
#include <cmath>
#include <set>
#include <vector>
#include <streambuf>
#include <curl/curl.h>
#include <Rcpp.h>
#include "zlib.h"
#include "straw.h"
using namespace std;

/*
  Straw: fast C++ implementation of dump. Not as fully featured as the
  Java version. Reads the .hic file, finds the appropriate matrix and slice
  of data, and outputs as text in sparse upper triangular format.

  Currently only supporting matrices.

  Usage: straw [observed/oe/expected] <NONE/VC/VC_SQRT/KR> <hicFile(s)> <chr1>[:x1:x2] <chr2>[:y1:y2] <BP/FRAG> <binsize>
 */
// this is for creating a stream from a byte array for ease of use
struct membuf : std::streambuf {
    membuf(char *begin, char *end) {
        this->setg(begin, begin, end);
    }
};

// for holding data from URL call
struct MemoryStruct {
    char *memory;
    size_t size;
};

// callback for libcurl. data written to this buffer
static size_t
WriteMemoryCallback(void *contents, size_t size, size_t nmemb, void *userp) {
    size_t realsize = size * nmemb;
    struct MemoryStruct *mem = (struct MemoryStruct *) userp;

    mem->memory = static_cast<char *>(realloc(mem->memory, mem->size + realsize + 1));
    if (mem->memory == nullptr) {
        /* out of memory! */
        Rcpp::stop("Not enough memory (realloc returned NULL).");
        return 0;
    }

    std::memcpy(&(mem->memory[mem->size]), contents, realsize);
    mem->size += realsize;
    mem->memory[mem->size] = 0;

    return realsize;
}

// get a buffer that can be used as an input stream from the URL
char *getData(CURL *curl, int64_t position, int64_t chunksize) {
    std::ostringstream oss;
    struct MemoryStruct chunk;

    chunk.memory = static_cast<char *>(malloc(1));
    chunk.size = 0;    /* no data at this point */
    oss << position << "-" << position + chunksize;
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, (void *) &chunk);
    curl_easy_setopt(curl, CURLOPT_RANGE, oss.str().c_str());
    CURLcode res = curl_easy_perform(curl);
    if (res != CURLE_OK) {
        Rcpp::stop("curl_easy_perform() failed: %s.",
                curl_easy_strerror(res));
    }
    //  printf("%lu bytes retrieved\n", (int64_t)chunk.size);

    return chunk.memory;
}

// returns whether or not this is valid HiC file
bool readMagicString(istream &fin) {
    string str;
    getline(fin, str, '\0');
    return str[0] == 'H' && str[1] == 'I' && str[2] == 'C';
}

char readCharFromFile(istream &fin) {
    char tempChar;
    fin.read(&tempChar, sizeof(char));
    return tempChar;
}

int16_t readInt16FromFile(istream &fin) {
    int16_t tempInt16;
    fin.read((char *) &tempInt16, sizeof(int16_t));
    return tempInt16;
}

int32_t readInt32FromFile(istream &fin) {
    int32_t tempInt32;
    fin.read((char *) &tempInt32, sizeof(int32_t));
    return tempInt32;
}

int64_t readInt64FromFile(istream &fin) {
    int64_t tempInt64;
    fin.read((char *) &tempInt64, sizeof(int64_t));
    return tempInt64;
}

float readFloatFromFile(istream &fin) {
    float tempFloat;
    fin.read((char *) &tempFloat, sizeof(float));
    return tempFloat;
}

double readDoubleFromFile(istream &fin) {
    double tempDouble;
    fin.read((char *) &tempDouble, sizeof(double));
    return tempDouble;
}

// reads the header, storing the positions of the normalization vectors and returning the masterIndexPosition pointer
map<string, chromosome> readHeader(istream &fin, int64_t &masterIndexPosition, string &genomeID, int32_t &numChromosomes,
                                   int32_t &version, int64_t &nviPosition, int64_t &nviLength, vector<int32_t> &bpResolutions) {
    map<string, chromosome> chromosomeMap;
    if (!readMagicString(fin)) {
        Rcpp::stop("Hi-C magic string is missing, does not appear to be a hic file.");
    }

    version = readInt32FromFile(fin);
    if (version < 6) {
        Rcpp::stop("Version %d no longer supported.", version);
    }
    masterIndexPosition = readInt64FromFile(fin);
    getline(fin, genomeID, '\0');

    if (version > 8) {
        nviPosition = readInt64FromFile(fin);
        nviLength = readInt64FromFile(fin);
    }

    int32_t nattributes = readInt32FromFile(fin);

    // reading and ignoring attribute-value dictionary
    for (int i = 0; i < nattributes; i++) {
        string key, value;
        getline(fin, key, '\0');
        getline(fin, value, '\0');
    }

    numChromosomes = readInt32FromFile(fin);
    // chromosome map for finding matrixType
    for (int i = 0; i < numChromosomes; i++) {
        string name;
        int64_t length;
        getline(fin, name, '\0');
        if (version > 8) {
            length = readInt64FromFile(fin);
        } else {
            length = (int64_t) readInt32FromFile(fin);
        }

        chromosome chr;
        chr.index = i;
        chr.name = name;
        chr.length = length;
        chromosomeMap[name] = chr;
    }

    int32_t nBpResolutions = readInt32FromFile(fin);
    for (int i=0; i<nBpResolutions; i++) {
      bpResolutions.push_back(readInt32FromFile(fin));
    }

    return chromosomeMap;
}

// reads the footer from the master pointer location. takes in the chromosomes,
// norm, unit (BP or FRAG) and resolution or binsize, and sets the file
// position of the matrix and the normalization vectors for those chromosomes
// at the given normalization and resolution
bool readFooter(istream &fin, int64_t master, int32_t version, int32_t c1, int32_t c2, const string &matrixType, const string &norm,
                const string &unit, int32_t resolution, int64_t &myFilePos, indexEntry &c1NormEntry, indexEntry &c2NormEntry,
                vector<double> &expectedValues) {
    if (version > 8) {
        /* int64_t nBytes = */ readInt64FromFile(fin);
    } else {
        /* int32_t nBytes = */ readInt32FromFile(fin);
    }

    stringstream ss;
    ss << c1 << "_" << c2;
    string key = ss.str();

    int32_t nEntries = readInt32FromFile(fin);
    bool found = false;
    for (int i = 0; i < nEntries; i++) {
        string str;
        getline(fin, str, '\0');
        int64_t fpos = readInt64FromFile(fin);
        /* int32_t sizeinbytes = */ readInt32FromFile(fin);
        if (str == key) {
            myFilePos = fpos;
            found = true;
        }
    }
    if (!found) {
        Rcpp::stop("File doesn't have the given chr_chr map %s.", key);
    }

    if ((matrixType == "observed" && norm == "NONE") || ((matrixType == "oe" || matrixType == "expected") && norm == "NONE" && c1 != c2))
        return true; // no need to read norm vector index

    // read in and ignore expected value maps; don't store; reading these to
    // get to norm vector index
    int32_t nExpectedValues = readInt32FromFile(fin);
    for (int i = 0; i < nExpectedValues; i++) {
        string unit0;
        getline(fin, unit0, '\0'); //unit
        int32_t binSize = readInt32FromFile(fin);

        int64_t nValues;
        if (version > 8) {
            nValues = readInt64FromFile(fin);
        } else {
            nValues = (int64_t) readInt32FromFile(fin);
        }

        bool store = c1 == c2 && (matrixType == "oe" || matrixType == "expected") && norm == "NONE" && unit0 == unit && binSize == resolution;

        if (version > 8) {
            for (int j = 0; j < nValues; j++) {
                double v = readFloatFromFile(fin);
                if (store) {
                    expectedValues.push_back(v);
                }
            }
        } else {
            for (int j = 0; j < nValues; j++) {
                double v = readDoubleFromFile(fin);
                if (store) {
                    expectedValues.push_back(v);
                }
            }
        }

        int32_t nNormalizationFactors = readInt32FromFile(fin);
        for (int j = 0; j < nNormalizationFactors; j++) {
            int32_t chrIdx = readInt32FromFile(fin);
            double v;
            if (version > 8) {
                v = readFloatFromFile(fin);
            } else {
                v = readDoubleFromFile(fin);
            }
            if (store && chrIdx == c1) {
                for (double &expectedValue : expectedValues) {
                    expectedValue = expectedValue / v;
                }
            }
        }
    }

    if (c1 == c2 && (matrixType == "oe" || matrixType == "expected") && norm == "NONE") {
        if (expectedValues.empty()) {
            Rcpp::stop("File did not contain expected values vectors at %d %s.", resolution, unit);
        }
        return true;
    }

    nExpectedValues = readInt32FromFile(fin);
    for (int i = 0; i < nExpectedValues; i++) {
        string type, unit0;
        getline(fin, type, '\0'); //typeString
        getline(fin, unit0, '\0'); //unit
        int32_t binSize = readInt32FromFile(fin);

        int64_t nValues;
        if (version > 8) {
            nValues = readInt64FromFile(fin);
        } else {
            nValues = (int64_t) readInt32FromFile(fin);
        }
        bool store = c1 == c2 && (matrixType == "oe" || matrixType == "expected") && type == norm && unit0 == unit && binSize == resolution;

        if (version > 8) {
            for (int j = 0; j < nValues; j++) {
                double v = readFloatFromFile(fin);
                if (store) {
                    expectedValues.push_back(v);
                }
            }
        } else {
            for (int j = 0; j < nValues; j++) {
                double v = readDoubleFromFile(fin);
                if (store) {
                    expectedValues.push_back(v);
                }
            }

        }

        int32_t nNormalizationFactors = readInt32FromFile(fin);
        for (int j = 0; j < nNormalizationFactors; j++) {
            int32_t chrIdx = readInt32FromFile(fin);
            double v;
            if (version > 8) {
                v = (double) readFloatFromFile(fin);
            } else {
                v = readDoubleFromFile(fin);
            }
            if (store && chrIdx == c1) {
                for (double &expectedValue : expectedValues) {
                    expectedValue = expectedValue / v;
                }
            }
        }
    }

    if (c1 == c2 && (matrixType == "oe" || matrixType == "expected") && norm != "NONE") {
        if (expectedValues.empty()) {
            Rcpp::stop("File did not contain normalized expected values vectors at %d %s.", resolution, unit);
        }
    }

    // Index of normalization vectors
    nEntries = readInt32FromFile(fin);
    bool found1 = false;
    bool found2 = false;
    for (int i = 0; i < nEntries; i++) {
        string normtype;
        getline(fin, normtype, '\0'); //normalization type
        int32_t chrIdx = readInt32FromFile(fin);
        string unit1;
        getline(fin, unit1, '\0'); //unit
        int32_t resolution1 = readInt32FromFile(fin);
        int64_t filePosition = readInt64FromFile(fin);
        int64_t sizeInBytes;
        if (version > 8) {
            sizeInBytes = readInt64FromFile(fin);
        } else {
            sizeInBytes = (int64_t) readInt32FromFile(fin);
        }

        if (chrIdx == c1 && normtype == norm && unit1 == unit && resolution1 == resolution) {
            c1NormEntry.position = filePosition;
            c1NormEntry.size = sizeInBytes;
            found1 = true;
        }
        if (chrIdx == c2 && normtype == norm && unit1 == unit && resolution1 == resolution) {
            c2NormEntry.position = filePosition;
            c2NormEntry.size = sizeInBytes;
            found2 = true;
        }
    }
    if (!found1 || !found2) {
        Rcpp::stop("File did not contain %s normalization vectors for one or both chromosomes at %d %s.", norm, resolution, unit);
    }
    return true;
}

// reads the raw binned contact matrix at specified resolution, setting the block bin count and block column count
map<int32_t, indexEntry> readMatrixZoomData(istream &fin, const string &myunit, int32_t mybinsize, float &mySumCounts,
                                        int32_t &myBlockBinCount, int32_t &myBlockColumnCount, bool &found) {

    map<int32_t, indexEntry> blockMap;
    string unit;
    getline(fin, unit, '\0'); // unit
    readInt32FromFile(fin); // Old "zoom" index -- not used
    float sumCounts = readFloatFromFile(fin); // sumCounts
    readFloatFromFile(fin); // occupiedCellCount
    readFloatFromFile(fin); // stdDev
    readFloatFromFile(fin); // percent95
    int32_t binSize = readInt32FromFile(fin);
    int32_t blockBinCount = readInt32FromFile(fin);
    int32_t blockColumnCount = readInt32FromFile(fin);

    found = false;
    if (myunit == unit && mybinsize == binSize) {
        mySumCounts = sumCounts;
        myBlockBinCount = blockBinCount;
        myBlockColumnCount = blockColumnCount;
        found = true;
    }

    int32_t nBlocks = readInt32FromFile(fin);

    for (int b = 0; b < nBlocks; b++) {
        int32_t blockNumber = readInt32FromFile(fin);
        int64_t filePosition = readInt64FromFile(fin);
        int32_t blockSizeInBytes = readInt32FromFile(fin);
        indexEntry entry = indexEntry();
        entry.size = (int64_t) blockSizeInBytes;
        entry.position = filePosition;
        if (found) blockMap[blockNumber] = entry;
    }
    return blockMap;
}

// reads the raw binned contact matrix at specified resolution, setting the block bin count and block column count
map<int32_t, indexEntry> readMatrixZoomDataHttp(CURL *curl, int64_t &myFilePosition, const string &myunit, int32_t mybinsize,
                                            float &mySumCounts, int32_t &myBlockBinCount, int32_t &myBlockColumnCount,
                                            bool &found) {

    map<int32_t, indexEntry> blockMap;
    char *buffer;
    int32_t header_size = 5 * sizeof(int32_t) + 4 * sizeof(float);
    char *first;
    first = getData(curl, myFilePosition, 1);
    if (first[0] == 'B') {
        header_size += 3;
    } else if (first[0] == 'F') {
        header_size += 5;
    } else {
        Rcpp::stop("Unit not understood.");
    }
    buffer = getData(curl, myFilePosition, header_size);
    membuf sbuf(buffer, buffer + header_size);
    istream fin(&sbuf);

    string unit;
    getline(fin, unit, '\0'); // unit
    readInt32FromFile(fin); // Old "zoom" index -- not used
    float sumCounts = readFloatFromFile(fin); // sumCounts
    readFloatFromFile(fin); // occupiedCellCount
    readFloatFromFile(fin); // stdDev
    readFloatFromFile(fin); // percent95
    int32_t binSize = readInt32FromFile(fin);
    int32_t blockBinCount = readInt32FromFile(fin);
    int32_t blockColumnCount = readInt32FromFile(fin);

    found = false;
    if (myunit == unit && mybinsize == binSize) {
        mySumCounts = sumCounts;
        myBlockBinCount = blockBinCount;
        myBlockColumnCount = blockColumnCount;
        found = true;
    }

    int32_t nBlocks = readInt32FromFile(fin);

    if (found) {
        int32_t chunkSize = nBlocks * (sizeof(int32_t) + sizeof(int64_t) + sizeof(int32_t));
        buffer = getData(curl, myFilePosition + header_size, chunkSize);
        membuf sbuf2(buffer, buffer + chunkSize);
        istream fin2(&sbuf2);
        for (int b = 0; b < nBlocks; b++) {
            int32_t blockNumber = readInt32FromFile(fin2);
            int64_t filePosition = readInt64FromFile(fin2);
            int32_t blockSizeInBytes = readInt32FromFile(fin2);
            indexEntry entry = indexEntry();
            entry.size = (int64_t) blockSizeInBytes;
            entry.position = filePosition;
            blockMap[blockNumber] = entry;
        }
    } else {
        myFilePosition = myFilePosition + header_size + (nBlocks * (sizeof(int32_t) + sizeof(int64_t) + sizeof(int32_t)));
    }
    delete buffer;
    return blockMap;
}

// goes to the specified file pointer in http and finds the raw contact matrixType at specified resolution, calling readMatrixZoomData.
// sets blockbincount and blockcolumncount
map<int32_t, indexEntry> readMatrixHttp(CURL *curl, int64_t myFilePosition, const string &unit, int32_t resolution,
                                    float &mySumCounts, int32_t &myBlockBinCount, int32_t &myBlockColumnCount) {
    char *buffer;
    int32_t size = sizeof(int32_t) * 3;
    buffer = getData(curl, myFilePosition, size);
    membuf sbuf(buffer, buffer + size);
    istream bufin(&sbuf);

    /* int32_t c1 = */ readInt32FromFile(bufin);
    /* int32_t c2 = */ readInt32FromFile(bufin);
    int32_t nRes = readInt32FromFile(bufin);
    int32_t i = 0;
    bool found = false;
    myFilePosition = myFilePosition + size;
    delete buffer;
    map<int32_t, indexEntry> blockMap;

    while (i < nRes && !found) {
        // myFilePosition gets updated within call
        blockMap = readMatrixZoomDataHttp(curl, myFilePosition, unit, resolution, mySumCounts, myBlockBinCount, myBlockColumnCount,
                                          found);
        i++;
    }
    if (!found) {
        Rcpp::stop("Error finding block data.");
    }
    return blockMap;
}

// goes to the specified file pointer and finds the raw contact matrixType at specified resolution, calling readMatrixZoomData.
// sets blockbincount and blockcolumncount
map<int32_t, indexEntry> readMatrix(istream &fin, int64_t myFilePosition, const string &unit, int32_t resolution,
                                float &mySumCounts, int32_t &myBlockBinCount, int32_t &myBlockColumnCount) {
    map<int32_t, indexEntry> blockMap;

    fin.seekg(myFilePosition, ios::beg);
    /* int32_t c1 = */ readInt32FromFile(fin);
    /* int32_t c2 = */ readInt32FromFile(fin);
    int32_t nRes = readInt32FromFile(fin);
    int32_t i = 0;
    bool found = false;
    while (i < nRes && !found) {
        blockMap = readMatrixZoomData(fin, unit, resolution, mySumCounts, myBlockBinCount, myBlockColumnCount, found);
        i++;
    }
    if (!found) {
        Rcpp::stop("Error finding block data.");
    }
    return blockMap;
}

// gets the blocks that need to be read for this slice of the data.  needs blockbincount, blockcolumncount, and whether
// or not this is intrachromosomal.
set<int32_t> getBlockNumbersForRegionFromBinPosition(const int64_t *regionIndices, int32_t blockBinCount, int32_t blockColumnCount,
                                                 bool intra) {
    int32_t col1 = static_cast<int32_t>(regionIndices[0] / blockBinCount);
    int32_t col2 = static_cast<int32_t>((regionIndices[1] + 1) / blockBinCount);
    int32_t row1 = static_cast<int32_t>(regionIndices[2] / blockBinCount);
    int32_t row2 = static_cast<int32_t>((regionIndices[3] + 1) / blockBinCount);

    set<int32_t> blocksSet;
    // first check the upper triangular matrixType
    for (int r = row1; r <= row2; r++) {
        for (int c = col1; c <= col2; c++) {
            int32_t blockNumber = r * blockColumnCount + c;
            blocksSet.insert(blockNumber);
        }
    }
    // check region part that overlaps with lower left triangle but only if intrachromosomal
    if (intra) {
        for (int r = col1; r <= col2; r++) {
            for (int c = row1; c <= row2; c++) {
                int32_t blockNumber = r * blockColumnCount + c;
                blocksSet.insert(blockNumber);
            }
        }
    }
    return blocksSet;
}

set<int32_t> getBlockNumbersForRegionFromBinPositionV9Intra(int64_t *regionIndices, int32_t blockBinCount, int32_t blockColumnCount) {
    // regionIndices is binX1 binX2 binY1 binY2
    set<int32_t> blocksSet;
    int32_t translatedLowerPAD = static_cast<int32_t>((regionIndices[0] + regionIndices[2]) / 2 / blockBinCount);
    int32_t translatedHigherPAD = static_cast<int32_t>((regionIndices[1] + regionIndices[3]) / 2 / blockBinCount + 1);
    int32_t translatedNearerDepth = static_cast<int32_t>(log2(
            1 + abs(regionIndices[0] - regionIndices[3]) / sqrt(2) / blockBinCount));
    int32_t translatedFurtherDepth = static_cast<int32_t>(log2(
            1 + abs(regionIndices[1] - regionIndices[2]) / sqrt(2) / blockBinCount));

    // because code above assume above diagonal; but we could be below diagonal
    int32_t nearerDepth = min(translatedNearerDepth, translatedFurtherDepth);
    if ((regionIndices[0] > regionIndices[3] && regionIndices[1] < regionIndices[2]) ||
        (regionIndices[1] > regionIndices[2] && regionIndices[0] < regionIndices[3])) {
        nearerDepth = 0;
    }
    int32_t furtherDepth = max(translatedNearerDepth, translatedFurtherDepth) + 1; // +1; integer divide rounds down

    for (int depth = nearerDepth; depth <= furtherDepth; depth++) {
        for (int pad = translatedLowerPAD; pad <= translatedHigherPAD; pad++) {
            int32_t blockNumber = depth * blockColumnCount + pad;
            blocksSet.insert(blockNumber);
        }
    }

    return blocksSet;
}

void appendRecord(vector<contactRecord> &vector, int32_t index, int32_t binX, int32_t binY, float counts) {
    contactRecord record = contactRecord();
    record.binX = binX;
    record.binY = binY;
    record.counts = counts;
    vector[index] = record;
}

// this is the meat of reading the data.  takes in the block number and returns the set of contact records corresponding to
// that block.  the block data is compressed and must be decompressed using the zlib library functions
vector<contactRecord> readBlock(istream &fin, CURL *curl, bool isHttp, indexEntry idx, int32_t version) {
    if (idx.size <= 0) {
        vector<contactRecord> v;
        return v;
    }
    char *compressedBytes = new char[idx.size];
    char *uncompressedBytes = new char[idx.size * 10]; //biggest seen so far is 3

    if (isHttp) {
        compressedBytes = getData(curl, idx.position, idx.size);
    } else {
        fin.seekg(idx.position, ios::beg);
        fin.read(compressedBytes, idx.size);
    }
    // Decompress the block
    // zlib struct
    z_stream infstream;
    infstream.zalloc = Z_NULL;
    infstream.zfree = Z_NULL;
    infstream.opaque = Z_NULL;
    infstream.avail_in = static_cast<uInt>(idx.size); // size of input
    infstream.next_in = (Bytef *) compressedBytes; // input char array
    infstream.avail_out = static_cast<uInt>(idx.size * 10); // size of output
    infstream.next_out = (Bytef *) uncompressedBytes; // output char array
    // the actual decompression work.
    inflateInit(&infstream);
    inflate(&infstream, Z_NO_FLUSH);
    inflateEnd(&infstream);
    int32_t uncompressedSize = static_cast<int32_t>(infstream.total_out);

    // create stream from buffer for ease of use
    membuf sbuf(uncompressedBytes, uncompressedBytes + uncompressedSize);
    istream bufferin(&sbuf);
    uint64_t nRecords = static_cast<uint64_t>(readInt32FromFile(bufferin));
    vector<contactRecord> v(nRecords);
    // different versions have different specific formats
    if (version < 7) {
        for (uInt i = 0; i < nRecords; i++) {
            int32_t binX = readInt32FromFile(bufferin);
            int32_t binY = readInt32FromFile(bufferin);
            float counts = readFloatFromFile(bufferin);
            appendRecord(v, i, binX, binY, counts);
        }
    } else {
        int32_t binXOffset = readInt32FromFile(bufferin);
        int32_t binYOffset = readInt32FromFile(bufferin);
        bool useShort = readCharFromFile(bufferin) == 0; // yes this is opposite of usual

        bool useShortBinX = true;
        bool useShortBinY = true;
        if (version > 8) {
            useShortBinX = readCharFromFile(bufferin) == 0;
            useShortBinY = readCharFromFile(bufferin) == 0;
        }

        char type = readCharFromFile(bufferin);
        int32_t index = 0;
        if (type == 1) {
            if (useShortBinX && useShortBinY) {
                int16_t rowCount = readInt16FromFile(bufferin);
                for (int i = 0; i < rowCount; i++) {
                    int32_t binY = binYOffset + readInt16FromFile(bufferin);
                    int16_t colCount = readInt16FromFile(bufferin);
                    for (int j = 0; j < colCount; j++) {
                        int32_t binX = binXOffset + readInt16FromFile(bufferin);
                        float counts;
                        if (useShort) {
                            counts = readInt16FromFile(bufferin);
                        } else {
                            counts = readFloatFromFile(bufferin);
                        }
                        appendRecord(v, index++, binX, binY, counts);
                    }
                }
            } else if (useShortBinX && !useShortBinY) {
                int32_t rowCount = readInt32FromFile(bufferin);
                for (int i = 0; i < rowCount; i++) {
                    int32_t binY = binYOffset + readInt32FromFile(bufferin);
                    int16_t colCount = readInt16FromFile(bufferin);
                    for (int j = 0; j < colCount; j++) {
                        int32_t binX = binXOffset + readInt16FromFile(bufferin);
                        float counts;
                        if (useShort) {
                            counts = readInt16FromFile(bufferin);
                        } else {
                            counts = readFloatFromFile(bufferin);
                        }
                        appendRecord(v, index++, binX, binY, counts);
                    }
                }
            } else if (!useShortBinX && useShortBinY) {
                int16_t rowCount = readInt16FromFile(bufferin);
                for (int i = 0; i < rowCount; i++) {
                    int32_t binY = binYOffset + readInt16FromFile(bufferin);
                    int32_t colCount = readInt32FromFile(bufferin);
                    for (int j = 0; j < colCount; j++) {
                        int32_t binX = binXOffset + readInt32FromFile(bufferin);
                        float counts;
                        if (useShort) {
                            counts = readInt16FromFile(bufferin);
                        } else {
                            counts = readFloatFromFile(bufferin);
                        }
                        appendRecord(v, index++, binX, binY, counts);
                    }
                }
            } else {
                int32_t rowCount = readInt32FromFile(bufferin);
                for (int i = 0; i < rowCount; i++) {
                    int32_t binY = binYOffset + readInt32FromFile(bufferin);
                    int32_t colCount = readInt32FromFile(bufferin);
                    for (int j = 0; j < colCount; j++) {
                        int32_t binX = binXOffset + readInt32FromFile(bufferin);
                        float counts;
                        if (useShort) {
                            counts = readInt16FromFile(bufferin);
                        } else {
                            counts = readFloatFromFile(bufferin);
                        }
                        appendRecord(v, index++, binX, binY, counts);
                    }
                }
            }
        } else if (type == 2) {
            int32_t nPts = readInt32FromFile(bufferin);
            int16_t w = readInt16FromFile(bufferin);

            for (int i = 0; i < nPts; i++) {
                //int32_t idx = (p.y - binOffset2) * w + (p.x - binOffset1);
                int32_t row = i / w;
                int32_t col = i - row * w;
                int32_t bin1 = binXOffset + col;
                int32_t bin2 = binYOffset + row;

                float counts;
                if (useShort) {
                    int16_t c = readInt16FromFile(bufferin);
                    if (c != -32768) {
                        appendRecord(v, index++, bin1, bin2, c);
                    }
                } else {
                    counts = readFloatFromFile(bufferin);
                    if (!isnan(counts)) {
                        appendRecord(v, index++, bin1, bin2, counts);
                    }
                }
            }
        }
    }
    delete[] compressedBytes;
    delete[] uncompressedBytes; // don't forget to delete your heap arrays in C++!
    return v;
}

// reads the normalization vector from the file at the specified location
vector<double> readNormalizationVector(istream &bufferin, int32_t version) {
    int64_t nValues;
    if (version > 8) {
        nValues = readInt64FromFile(bufferin);
    } else {
        nValues = (int64_t) readInt32FromFile(bufferin);
    }

    uint64_t numValues = static_cast<uint64_t>(nValues);
    vector<double> values(numValues);

    if (version > 8) {
        for (int i = 0; i < nValues; i++) {
            values[i] = (double) readFloatFromFile(bufferin);
        }
    } else {
        for (int i = 0; i < nValues; i++) {
            values[i] = readDoubleFromFile(bufferin);
        }
    }

    return values;
}

class FileReader {
public:
    string prefix = "http"; // HTTP code
    ifstream fin;
    CURL *curl;
    bool isHttp = false;

    static CURL *initCURL(const char *url) {
        CURL *curl = curl_easy_init();
        if (curl) {
            curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, WriteMemoryCallback);
            curl_easy_setopt(curl, CURLOPT_URL, url);
            curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1L);
            curl_easy_setopt(curl, CURLOPT_USERAGENT, "straw");
        }
        return curl;
    }

    explicit FileReader(const string &fname) {

        // read header into buffer; 100K should be sufficient
        if (std::strncmp(fname.c_str(), prefix.c_str(), prefix.size()) == 0) {
            isHttp = true;
            curl = initCURL(fname.c_str());
            if (!curl) {
                Rcpp::stop("URL %s cannot be opened for reading.", fname);
            }
        } else {
            fin.open(fname, fstream::in | fstream::binary);
            if (!fin) {
                Rcpp::stop("File %s cannot be opened for reading.", fname);
            }
        }
    }

    void close(){
        if(isHttp){
            curl_easy_cleanup(curl);
        } else {
            fin.close();
        }
    }
};

class HiCFile {
public:
    string prefix = "http"; // HTTP code
    bool isHttp = false;
    ifstream fin;
    CURL *curl;
    int64_t master = 0LL;
    map<string, chromosome> chromosomeMap;
    vector<int32_t> bpResolutions;
    string genomeID;
    int32_t numChromosomes = 0;
    int32_t version = 0;
    int64_t nviPosition = 0LL;
    int64_t nviLength = 0LL;
    static int64_t totalFileSize;

    static size_t hdf(char *b, size_t size, size_t nitems, void *userdata) {
        size_t numbytes = size * nitems;
        b[numbytes + 1] = '\0';
        string s(b);
        int32_t found = static_cast<int32_t>(s.find("content-range"));
        if ((size_t)found == string::npos) {
          found = static_cast<int32_t>(s.find("Content-Range"));
        }
        if ((size_t)found != string::npos) {
            int32_t found2 = static_cast<int32_t>(s.find("/"));
            //content-range: bytes 0-100000/891471462
            if ((size_t)found2 != string::npos) {
                string total = s.substr(found2 + 1);
                totalFileSize = stol(total);
            }
        }

        return numbytes;
    }

    static CURL *initCURL(const char *url) {
        CURL *curl = curl_easy_init();
        if (curl) {
            curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, WriteMemoryCallback);
            curl_easy_setopt(curl, CURLOPT_URL, url);
            curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1L);
            curl_easy_setopt(curl, CURLOPT_HEADERFUNCTION, hdf);
            curl_easy_setopt(curl, CURLOPT_USERAGENT, "straw");
        }
        return curl;
    }

    explicit HiCFile(const string &fname) {

        // read header into buffer; 100K should be sufficient
        if (std::strncmp(fname.c_str(), prefix.c_str(), prefix.size()) == 0) {
            isHttp = true;
            char *buffer;
            curl = initCURL(fname.c_str());
            if (curl) {
                buffer = getData(curl, 0, 100000);
            } else {
                Rcpp::stop("URL %s cannot be opened for reading.", fname);
            }
            membuf sbuf(buffer, buffer + 100000);
            istream bufin(&sbuf);
            chromosomeMap = readHeader(bufin, master, genomeID, numChromosomes,
                                       version, nviPosition, nviLength, bpResolutions);
            delete buffer;
        } else {
            fin.open(fname, fstream::in | fstream::binary);
            if (!fin) {
                Rcpp::stop("File %s cannot be opened for reading.", fname);
            }
            chromosomeMap = readHeader(fin, master, genomeID, numChromosomes,
                                       version, nviPosition, nviLength, bpResolutions);
        }
    }

    void close(){
        if(isHttp){
            curl_easy_cleanup(curl);
        } else {
            fin.close();
        }
    }

    vector<double> readNormalizationVectorFromFooter(indexEntry cNormEntry) {
        char *buffer;
        if (isHttp) {
            buffer = getData(curl, cNormEntry.position, cNormEntry.size);
        } else {
            buffer = new char[cNormEntry.size];
            fin.seekg(cNormEntry.position, ios::beg);
            fin.read(buffer, cNormEntry.size);
        }
        membuf sbuf3(buffer, buffer + cNormEntry.size);
        istream bufferin(&sbuf3);
        vector<double> cNorm = readNormalizationVector(bufferin, version);
        delete buffer;
        return cNorm;
    }
};

int64_t HiCFile::totalFileSize = 0LL;

class MatrixZoomData {
public:
    indexEntry c1NormEntry, c2NormEntry;
    int64_t myFilePos = 0LL;
    vector<double> expectedValues;
    bool foundFooter = false;
    vector<double> c1Norm;
    vector<double> c2Norm;
    int32_t c1 = 0;
    int32_t c2 = 0;
    string matrixType;
    string norm;
    string unit;
    int32_t resolution = 0;
    int32_t numBins1 = 0;
    int32_t numBins2 = 0;

    MatrixZoomData(HiCFile *hiCFile, const chromosome &chrom1, const chromosome &chrom2, const string &matrixType,
                   const string &norm, const string &unit, int32_t resolution) {

        int32_t c01 = chrom1.index;
        int32_t c02 = chrom2.index;
        if (c01 <= c02) { // default is ok
            this->c1 = c01;
            this->c2 = c02;
            this->numBins1 = static_cast<int32_t>(chrom1.length / resolution);
            this->numBins2 = static_cast<int32_t>(chrom2.length / resolution);
        } else { // flip
            this->c1 = c02;
            this->c2 = c01;
            this->numBins1 = static_cast<int32_t>(chrom2.length / resolution);
            this->numBins2 = static_cast<int32_t>(chrom1.length / resolution);
        }

        this->matrixType = matrixType;
        this->norm = norm;
        this->unit = unit;
        this->resolution = resolution;

        if (hiCFile->isHttp) {
            char *buffer2;
            int64_t bytes_to_read = hiCFile->totalFileSize - hiCFile->master;
            buffer2 = getData(hiCFile->curl, hiCFile->master, bytes_to_read);
            membuf sbuf2(buffer2, buffer2 + bytes_to_read);
            istream bufin2(&sbuf2);
            foundFooter = readFooter(bufin2, hiCFile->master, hiCFile->version, c1, c2, matrixType, norm, unit,
                                     resolution,
                                     myFilePos,
                                     c1NormEntry, c2NormEntry, expectedValues);
            delete buffer2;
        } else {
            hiCFile->fin.seekg(hiCFile->master, ios::beg);
            foundFooter = readFooter(hiCFile->fin, hiCFile->master, hiCFile->version, c1, c2, matrixType, norm,
                                     unit,
                                     resolution, myFilePos,
                                     c1NormEntry, c2NormEntry, expectedValues);
        }

        if (!foundFooter) {
            return;
        }

        if (norm != "NONE") {
            c1Norm = hiCFile->readNormalizationVectorFromFooter(c1NormEntry);
            if (c1 == c2) {
                c2Norm = c1Norm;
            } else {
                c2Norm = hiCFile->readNormalizationVectorFromFooter(c2NormEntry);
            }
        }
    }
};

MatrixZoomData *
getMatrixZoomData(HiCFile *hiCFile, const string &chr1, const string &chr2, string matrixType, string norm,
                  string unit, int32_t resolution) {

    chromosome chrom1 = hiCFile->chromosomeMap[chr1];
    chromosome chrom2 = hiCFile->chromosomeMap[chr2];
    return new MatrixZoomData(hiCFile, chrom1, chrom2, std::move(matrixType), std::move(norm), std::move(unit),
                              resolution);
}

void parsePositions(const string &chrLoc, string &chrom, int64_t &pos1, int64_t &pos2, map<string, chromosome> map) {
    string x, y;
    stringstream ss(chrLoc);
    getline(ss, chrom, ':');
    if (map.count(chrom) == 0) {
        Rcpp::stop("%s not found in the file.", chrom);
    }

    if (getline(ss, x, ':') && getline(ss, y, ':')) {
        pos1 = stol(x);
        pos2 = stol(y);
    } else {
        pos1 = 0LL;
        pos2 = map[chrom].length;
    }
}

class BlocksRecords {
public:
    float sumCounts;
    int32_t blockBinCount, blockColumnCount;
    map<int32_t, indexEntry> blockMap;
    double avgCount;
    bool isIntra;

    BlocksRecords(FileReader *fileReader, const footerInfo &footer) {

        isIntra = footer.c1 == footer.c2;

        if (fileReader->isHttp) {
            // readMatrix will assign blockBinCount and blockColumnCount
            blockMap = readMatrixHttp(fileReader->curl, footer.myFilePos, footer.unit, footer.resolution, sumCounts,
                                      blockBinCount,
                                      blockColumnCount);
        } else {
            // readMatrix will assign blockBinCount and blockColumnCount
            blockMap = readMatrix(fileReader->fin, footer.myFilePos, footer.unit, footer.resolution, sumCounts,
                                  blockBinCount,
                                  blockColumnCount);
        }

        if (!isIntra) {
            avgCount = (sumCounts / footer.numBins1) / footer.numBins2;   // <= trying to avoid overflows
        }
    }

    set<int32_t> getBlockNumbers(int32_t version, bool isIntra, int64_t *regionIndices, int32_t blockBinCount, int32_t blockColumnCount) {
        set<int32_t> blockNumbers;
        if (version > 8 && isIntra) {
            blockNumbers = getBlockNumbersForRegionFromBinPositionV9Intra(regionIndices, blockBinCount,
                                                                          blockColumnCount);
        } else {
            blockNumbers = getBlockNumbersForRegionFromBinPosition(regionIndices, blockBinCount, blockColumnCount,
                                                                   isIntra);
        }
        return blockNumbers;
    }

    vector<contactRecord>
    getRecords(FileReader *fileReader, int64_t regionIndices[4],
               const int64_t origRegionIndices[4], const footerInfo &footer) {

        set<int32_t> blockNumbers = getBlockNumbers(footer.version, isIntra, regionIndices, blockBinCount,
                                                blockColumnCount);

        vector<contactRecord> records;
        for (int32_t blockNumber : blockNumbers) {
            // get contacts in this block
            //cout << *it << " -- " << blockMap.size() << endl;
            //cout << blockMap[*it].size << " " <<  blockMap[*it].position << endl;
            vector<contactRecord> tmp_records = readBlock(fileReader->fin, fileReader->curl, fileReader->isHttp,
                                                          blockMap[blockNumber], footer.version);
            for (contactRecord rec : tmp_records) {
                int64_t x = rec.binX * footer.resolution;
                int64_t y = rec.binY * footer.resolution;

                if ((x >= origRegionIndices[0] && x <= origRegionIndices[1] &&
                     y >= origRegionIndices[2] && y <= origRegionIndices[3]) ||
                    // or check regions that overlap with lower left
                    (isIntra && y >= origRegionIndices[0] && y <= origRegionIndices[1] && x >= origRegionIndices[2] &&
                     x <= origRegionIndices[3])) {

                    float c = rec.counts;
                    if (footer.norm != "NONE") {
                        c = static_cast<float>(c / (footer.c1Norm[rec.binX] * footer.c2Norm[rec.binY]));
                    }
                    if (footer.matrixType == "oe") {
                        if (isIntra) {
                            c = static_cast<float>(c / footer.expectedValues[min(footer.expectedValues.size() - 1,
                                                                                 (size_t) floor(abs(y - x) /
                                                                                                footer.resolution))]);
                        } else {
                            c = static_cast<float>(c / avgCount);
                        }
                    } else if (footer.matrixType == "expected") {
                        if (isIntra) {
                            c = static_cast<float>(footer.expectedValues[min(footer.expectedValues.size() - 1,
                                                                                 (size_t) floor(abs(y - x) /
                                                                                                footer.resolution))]);
                        } else {
                            c = static_cast<float>(avgCount);
                        }
                    }

                    contactRecord record = contactRecord();
                    record.binX = static_cast<int32_t>(x);
                    record.binY = static_cast<int32_t>(y);
                    record.counts = c;
                    records.push_back(record);
                }
            }
        }
        return records;
    }
};

vector<contactRecord> getBlockRecords(FileReader *fileReader, int64_t origRegionIndices[4], const footerInfo &footer) {
    if (!footer.foundFooter) {
        vector<contactRecord> v;
        return v;
    }

    int64_t regionIndices[4]; // used to find the blocks we need to access
    regionIndices[0] = origRegionIndices[0] / footer.resolution;
    regionIndices[1] = origRegionIndices[1] / footer.resolution;
    regionIndices[2] = origRegionIndices[2] / footer.resolution;
    regionIndices[3] = origRegionIndices[3] / footer.resolution;

    BlocksRecords *blocksRecords = new BlocksRecords(fileReader, footer);
    return blocksRecords->getRecords(fileReader, regionIndices, origRegionIndices, footer);
}

footerInfo getNormalizationInfoForRegion(string fname, string chr1, string chr2,
                                         const string &matrixType, const string &norm,
                                         const string &unit, int32_t binsize) {

    HiCFile *hiCFile = new HiCFile(std::move(fname));
    MatrixZoomData *mzd = getMatrixZoomData(hiCFile, chr1, chr2, std::move(matrixType), std::move(norm), unit,
                                            binsize);
    footerInfo footer = footerInfo();
    footer.resolution = mzd->resolution;
    footer.foundFooter = mzd->foundFooter;
    footer.version = hiCFile->version;
    footer.c1 = mzd->c1;
    footer.c2 = mzd->c2;
    footer.numBins1 = mzd->numBins1;
    footer.numBins2 = mzd->numBins2;
    footer.myFilePos = mzd->myFilePos;
    footer.unit = mzd->unit;
    footer.norm = mzd->norm;
    footer.matrixType = mzd->matrixType;
    footer.c1Norm = mzd->c1Norm;
    footer.c2Norm = mzd->c2Norm;
    footer.expectedValues = mzd->expectedValues;
    hiCFile->close();
    return footer;
}

vector<contactRecord>
getBlockRecordsWithNormalization(string fname,
                                 int64_t c1pos1, int64_t c1pos2, int64_t c2pos1, int64_t c2pos2,
                                 int32_t resolution, bool foundFooter, int32_t version, int32_t c1, int32_t c2,
                                 int32_t numBins1, int32_t numBins2, int64_t myFilePos, string unit, string norm,
                                 string matrixType,
                                 vector<double> c1Norm, vector<double> c2Norm, vector<double> expectedValues) {
    int64_t origRegionIndices[4]; // as given by user
    origRegionIndices[0] = c1pos1;
    origRegionIndices[1] = c1pos2;
    origRegionIndices[2] = c2pos1;
    origRegionIndices[3] = c2pos2;

    FileReader *fileReader = new FileReader(std::move(fname));
    footerInfo footer = footerInfo();
    footer.resolution = resolution;
    footer.foundFooter = foundFooter;
    footer.version = version;
    footer.c1 = c1;
    footer.c2 = c2;
    footer.numBins1 = numBins1;
    footer.numBins2 = numBins2;
    footer.myFilePos = myFilePos;
    footer.unit = unit;
    footer.norm = norm;
    footer.matrixType = matrixType;
    footer.c1Norm = c1Norm;
    footer.c2Norm = c2Norm;
    footer.expectedValues = expectedValues;
    vector<contactRecord> v = getBlockRecords(fileReader, origRegionIndices, footer);
    fileReader->close();
    return v;
}

//' Straw Quick Dump
//'
//' fast C++ implementation of dump. Not as fully featured as the
//' Java version. Reads the .hic file, finds the appropriate matrix and slice
//' of data, and outputs as data.frame in sparse upper triangular format.
//' Currently only supporting matrices.
//'
//' Usage: straw <NONE/VC/VC_SQRT/KR> <hicFile(s)> <chr1>[:x1:x2] <chr2>[:y1:y2] <BP/FRAG> <binsize> [observed/oe/expected]
//'
//' @param norm Normalization to apply. Must be one of NONE/VC/VC_SQRT/KR.
//'     VC is vanilla coverage, VC_SQRT is square root of vanilla coverage, and KR is Knight-Ruiz or
//'     Balanced normalization.
//' @param fname path to .hic file
//' @param chr1loc first chromosome location
//' @param chr2loc second chromosome location
//' @param unit BP (BasePair) or FRAG (FRAGment)
//' @param binsize The bin size. By default, for BP, this is one of <2500000, 1000000, 500000,
//'     250000, 100000, 50000, 25000, 10000, 5000> and for FRAG this is one of <500, 200,
//'     100, 50, 20, 5, 2, 1>.
//' @param matrix Type of matrix to output. Must be one of observed/oe/expected.
//'     observed is observed counts, oe is observed/expected counts, expected is expected counts.
//' @return Data.frame of a sparse matrix of data from hic file. x,y,counts
//' @examples
//' straw("NONE",
//'   system.file("extdata", "test_chr22.hic", package = "trackViewer"),
//'  "22", "22", "BP", 2500000)
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame
straw(std::string norm, std::string fname, std::string chr1loc, std::string chr2loc, const std::string &unit, int32_t binsize, std::string matrix = "observed") {
    if (!(unit == "BP" || unit == "FRAG")) {
        Rcpp::stop("Norm specified incorrectly, must be one of <BP/FRAG>.\nUsage: straw <NONE/VC/VC_SQRT/KR> <hicFile(s)> <chr1>[:x1:x2] <chr2>[:y1:y2] <BP/FRAG> <binsize> [observed/oe/expected].");
    }

    HiCFile *hiCFile = new HiCFile(std::move(fname));

    string chr1, chr2;
    int64_t c1pos1 = -100LL, c1pos2 = -100LL, c2pos1 = -100LL, c2pos2 = -100LL;
    parsePositions(std::move(chr1loc), chr1, c1pos1, c1pos2, hiCFile->chromosomeMap);
    parsePositions(std::move(chr2loc), chr2, c2pos1, c2pos2, hiCFile->chromosomeMap);

    // from header have size of chromosomes, set region to read

    int64_t origRegionIndices[4]; // as given by user
    // reverse order if necessary
    if (hiCFile->chromosomeMap[chr1].index > hiCFile->chromosomeMap[chr2].index) {
        origRegionIndices[0] = c2pos1;
        origRegionIndices[1] = c2pos2;
        origRegionIndices[2] = c1pos1;
        origRegionIndices[3] = c1pos2;
    } else {
        origRegionIndices[0] = c1pos1;
        origRegionIndices[1] = c1pos2;
        origRegionIndices[2] = c2pos1;
        origRegionIndices[3] = c2pos2;
    }
    hiCFile->close();

    footerInfo footer = getNormalizationInfoForRegion(fname, chr1, chr2, matrix, norm, unit, binsize);

    vector<contactRecord> records = getBlockRecordsWithNormalization(fname,
                                            origRegionIndices[0], origRegionIndices[1],
                                            origRegionIndices[2], origRegionIndices[3],
                                            footer.resolution, footer.foundFooter, footer.version,
                                            footer.c1, footer.c2, footer.numBins1, footer.numBins2,
                                            footer.myFilePos, footer.unit, footer.norm, footer.matrixType,
                                            footer.c1Norm, footer.c2Norm, footer.expectedValues);
    vector<int32_t> xActual_vec, yActual_vec;
    vector<float> counts_vec;
    for (vector<contactRecord>::iterator it=records.begin(); it!=records.end(); ++it) {
      xActual_vec.push_back(it->binX);
      yActual_vec.push_back(it->binY);
      counts_vec.push_back(it->counts);
    }
    return Rcpp::DataFrame::create(Rcpp::Named("x") = xActual_vec, Rcpp::Named("y") = yActual_vec, Rcpp::Named("counts") = counts_vec);
}

vector<chromosome> getChromosomes(string fname){
    HiCFile *hiCFile = new HiCFile(std::move(fname));
    vector<chromosome> chromosomes;
    std::map<std::string, chromosome>::iterator iter = hiCFile->chromosomeMap.begin();
    while (iter != hiCFile->chromosomeMap.end()) {
        chromosomes.push_back(static_cast<chromosome>(iter->second));
        iter++;
    }
    hiCFile->close();
    return chromosomes;
}

//' Function for reading chromosomes from .hic file
//'
//' @param fname path to .hic file
//' @return Data frame of chromosome indices, names and lengths
//' @examples
//' readHicChroms(system.file("extdata", "test_chr22.hic", package = "trackViewer"))
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame readHicChroms(std::string fname)
{
  vector<chromosome> chroms = getChromosomes(std::move(fname));
  Rcpp::NumericVector indices;
  Rcpp::StringVector names;
  Rcpp::NumericVector lengths;
  for (std::vector<chromosome>::iterator it = chroms.begin(); it != chroms.end(); ++it) {
    indices.push_back(it->index);
    names.push_back(it->name);
    lengths.push_back(it->length);
  }
  return Rcpp::DataFrame::create(Rcpp::Named("index") = indices, Rcpp::Named("name") = names, Rcpp::Named("length") = lengths);
}

//' Function for reading basepair resolutions from .hic file
//'
//' @param fname path to .hic file
//' @return Vector of basepair resolutions
//' @examples
//' readHicBpResolutions(system.file("extdata", "test_chr22.hic", package = "trackViewer"))
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector readHicBpResolutions(std::string fname)
{
  HiCFile *hiCFile = new HiCFile(std::move(fname));
  Rcpp::NumericVector bpResolutions;
  for (std::vector<int32_t>::iterator it = hiCFile->bpResolutions.begin(); it != hiCFile->bpResolutions.end(); ++it) {
    bpResolutions.push_back(*it);
  }
  hiCFile->close();
  return bpResolutions;
}

// Reads all normalizations from the footer
Rcpp::CharacterVector readNormsFromFooter(istream &fin, int64_t master, int32_t version) {
    
    // Initialize variable to store norm types
    Rcpp::CharacterVector normTypes;
    
    // Read through the footer section
    //nBytes
    if (version > 8) {
        readInt64FromFile(fin);
    } else {
        readInt32FromFile(fin);
    }
    
    // nEntries
    int32_t nEntries = readInt32FromFile(fin);
    for (int i = 0; i < nEntries; i++) {
        string str;
        getline(fin, str, '\0');
        readInt64FromFile(fin); //fpos
        readInt32FromFile(fin); //sizeInBytes
    }

    // nExpectedValues
    int32_t nExpectedValues = readInt32FromFile(fin);
    for (int i = 0; i < nExpectedValues; i++) {
        string unit0;
        getline(fin, unit0, '\0'); //unit
        readInt32FromFile(fin);

        int64_t nValues;
        if (version > 8) {
            nValues = readInt64FromFile(fin);
            for (int j = 0; j < nValues; j++) {
                readFloatFromFile(fin);
            }
        } else {
            nValues = (int64_t) readInt32FromFile(fin);
            for (int j = 0; j < nValues; j++) {
                readDoubleFromFile(fin);
            }
        }

        int32_t nNormalizationFactors = readInt32FromFile(fin);
        for (int j = 0; j < nNormalizationFactors; j++) {
            readInt32FromFile(fin); //chrIdx
            if (version > 8) {
                readFloatFromFile(fin); //v
            } else {
                readDoubleFromFile(fin);//v
            }
        }
    }

    // Needs to be read like this (readInt32FromFile doesn't work)
    fin.read((char*)&nExpectedValues, sizeof(int32_t));
    for (int i = 0; i < nExpectedValues; i++) {
        //Record available norm types (handling empty strings as NONE)
        string type;
        getline(fin, type, '\0'); //typeString
        if (type == "") {
            type = "NONE";
        }
        normTypes.push_back(type);

        string unit0;
        getline(fin, unit0, '\0'); //unit
        readInt32FromFile(fin);

        int64_t nValues;
        if (version > 8) {
            nValues = readInt64FromFile(fin);
            for (int j = 0; j < nValues; j++) {
                readFloatFromFile(fin); //v
            }
        } else {
            nValues = (int64_t) readInt32FromFile(fin);
            for (int j = 0; j < nValues; j++) {
                readDoubleFromFile(fin); //v
            }
        }

        int32_t nNormalizationFactors = readInt32FromFile(fin);
        for (int j = 0; j < nNormalizationFactors; j++) {
            readInt32FromFile(fin); //chrIdx
            if (version > 8) {
                readFloatFromFile(fin); //v
            } else {
                readDoubleFromFile(fin); //v
            }
        }
    }
    
    // Include "NONE"
    normTypes.push_back("NONE");
    
    // Return unique norms
    return unique(normTypes);
}

//' Function for reading available normalizations from .hic file
//' 
//' @param fname path to .hic file
//' @return Vector of available normalizations
//' @examples
//' readHicNormTypes(system.file("extdata", "test_chr22.hic", package = "trackViewer"))
//' @export
// [[Rcpp::export]]
Rcpp::CharacterVector readHicNormTypes(std::string fname)
{
    HiCFile *hiCFile = new HiCFile(std::move(fname));
    Rcpp::CharacterVector normTypes;
    hiCFile->fin.seekg(hiCFile->master, ios::beg);
    normTypes = readNormsFromFooter(hiCFile->fin,
                                    hiCFile->master,
                                    hiCFile->version);
    hiCFile->close();
    return normTypes;
}
