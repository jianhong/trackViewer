#include <Rcpp.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <map>
#include <zlib.h>
using namespace Rcpp;

#define CHUNK 16384

bool debug = false;
long footerPosition;
long bodyPosition;
int version;

std::string readString(std::ifstream& ifile){
  char s;
  char* out = new char[1024];
  ifile.read(&s, sizeof(s));
  int count=0;
  while(s != '\0'){
    out[count++]= s;
    ifile.read(&s, sizeof(s));
  }
  std::string str(out, count);
  return(str);
}
int readInt(std::ifstream& ifile){
  int i;
  ifile.read(reinterpret_cast<char *>(&i), sizeof(int));
  return(i);
}
long readLong(std::ifstream& ifile){
  long i;
  ifile.read(reinterpret_cast<char *>(&i), sizeof(long));
  return(i);
}
float readFloat(std::ifstream& ifile){
  float i;
  ifile.read(reinterpret_cast<char *>(&i), sizeof(float));
  return(i);
}
double readDouble(std::ifstream& ifile){
  double i;
  ifile.read(reinterpret_cast<char *>(&i), sizeof(double));
  return(i);
}
int getInt(std::istringstream& ifile){
  int i;
  ifile.read(reinterpret_cast<char*>(&i), sizeof(int));
  return(i);
}
short getShort(std::istringstream& ifile){
  short i;
  ifile.read(reinterpret_cast<char*>(&i), sizeof(short));
  return(i);
}
long getLong(std::istringstream& ifile){
  long i;
  ifile.read(reinterpret_cast<char*>(&i), sizeof(long));
  return(i);
}
float getFloat(std::istringstream& ifile){
  float i;
  ifile.read(reinterpret_cast<char*>(&i), sizeof(float));
  return(i);
}
Byte getByte(std::istringstream& ifile){
  Byte i;
  ifile.read(reinterpret_cast<char*>(&i), sizeof(Byte));
  return(i);
}

class intValues {
public: 
  int nValues;
  int* values;
  
  void read(std::ifstream& ifile){
    nValues = readInt(ifile);
    values = new int[nValues];
    for(int i=0; i<nValues; i++){
      values[i] = readInt(ifile);
    }
  }
};

class doubleValues {
public: 
  int nValues;
  double* values;
  
  void read(std::ifstream& ifile){
    nValues = readInt(ifile);
    values = new double[nValues];
    for(int i=0; i<nValues; i++){
      values[i] = readDouble(ifile);
    }
  }
};

class intDoubleValues {
public:
  int nValues;
  std::map<int, double> values;
  
  void read(std::ifstream& ifile){
    nValues = readInt(ifile);
    for(int i=0; i<nValues; i++){
      int index = readInt(ifile);
      double val = readDouble(ifile);
      values.insert( std::pair<int, double>(index, val));
    }
  }
  
  void printValues(){
    Rcout << "nValues " << nValues << std::endl;
    
    for(std::map<int, double>::iterator it=values.begin();
        it!=values.end(); ++it){
      Rcout << it->first << "\t" << it->second << std::endl;
    }
  }
};

class chromosome{
public:
  std::string chrName;
  long chrLength;
  
  void read(std::ifstream& ifile){
    chrName = readString(ifile);
    if(version < 9){
      chrLength = (long)readInt(ifile);
    }else{
      chrLength = readLong(ifile);
    }
  }
};
class chromosomes{
public:
  int nChrs;
  chromosome* chromosomeIndexMap;
  
  void printChrIndexMap(){
    for(int i=0; i<nChrs; i++){
      Rcout << chromosomeIndexMap[i].chrName << "\t" << chromosomeIndexMap[i].chrLength << std::endl;
    }
  }
  void read(std::ifstream& ifile){
    nChrs = readInt(ifile);
    chromosomeIndexMap = new chromosome[nChrs];
    for(int i=0; i<nChrs; i++){
      chromosomeIndexMap[i].read(ifile);
    }
  }
  int find(std::string chr){
    for(int i=0; i<nChrs; i++){
      if(chromosomeIndexMap[i].chrName==chr){
        return(i);
      }
    }
    for(int i=0; i<nChrs; i++){
      if(chromosomeIndexMap[i].chrName=="chr"+chr){
        return(i);
      }
    }
    for(int i=0; i<nChrs; i++){
      if("chr"+chromosomeIndexMap[i].chrName==chr){
        return(i);
      }
    }
    stop("seqname is not available in the file.");
  }
};

class HicFileHead {
public:
  std::string magic;
  std::string genomeId;
  int nAttributes;
  std::map<std::string, std::string> attributes;
  chromosomes chromosomeIndex;
  intValues bpResolutions;
  intValues fragResolutions;
  intValues sites;
  
  void printAttributes(){
    for(std::map<std::string, std::string>::iterator it=attributes.begin();
        it!=attributes.end(); ++it){
      Rcout << it->first << "\t" << it->second << std::endl;
    }
  }
  void readHeader(std::ifstream& ifile){
    if(ifile.is_open()){
      ifile.seekg(0, std::ios::end);
      long length = ifile.tellg();
      ifile.seekg(0, std::ios::beg);
      if(length>0){
        magic = readString(ifile);
        if(debug) {
          Rcout << "magic '" << magic << "'" << ifile.tellg() << std::endl;
        }
        version = readInt(ifile);
        if(debug) {
          Rcout << "version " << version << "\t" << ifile.tellg() << std::endl;
        }
        footerPosition = readLong(ifile);
        if(debug) {
          Rcout << "footerPosition " << footerPosition << "\t" << ifile.tellg() << std::endl;
        }
        genomeId = readString(ifile);
        if(debug) {
          Rcout << "genomeId " << genomeId << "\t" << ifile.tellg() << std::endl;
        }
        nAttributes = readInt(ifile);
        if(debug) {
          Rcout << "nAttributes " << nAttributes << "\t" << ifile.tellg() << std::endl;
        }
        for(int i=0; i<nAttributes; i++){
          std::string attr_key = readString(ifile);
          std::string attr_value = readString(ifile);
          attributes[attr_key] = attr_value;
        }
        if(debug) {
          printAttributes();
        }
        chromosomeIndex.read(ifile);
        if(debug) {
          Rcout << "nChrs " << chromosomeIndex.nChrs << std::endl;
          chromosomeIndex.printChrIndexMap();
        }
        bpResolutions.read(ifile);
        if(debug){
          Rcout << "nBPResolutions " << bpResolutions.nValues << "\t" << ifile.tellg() << std::endl;
        }
        fragResolutions.read(ifile);
        if(debug){
          Rcout << "nFragResolutions " << fragResolutions.nValues << "\t" << ifile.tellg() << std::endl;
        }
        sites.read(ifile);
        if(debug){
          Rcout << "nSites " << sites.nValues << "\t" << ifile.tellg() << std::endl;
        }
        bodyPosition = ifile.tellg();
      }
    }
  }
};

struct mIdx{
  long start;
  int size;
};
class entry{
public:
  int nEntries;
  std::map<std::string, mIdx> masterIndex;
  
  void read(std::ifstream& ifile){
    nEntries = readInt(ifile);
    for(int i=0; i<nEntries; i++){
      std::string key = readString(ifile);
      mIdx idx;
      idx.start = readLong(ifile);
      idx.size = readInt(ifile);
      masterIndex[key] = idx; 
    }
  }
  
  void printEntries(){
    Rcout << "nEntries " << nEntries << std::endl;
    
    for(std::map<std::string, mIdx>::iterator it=masterIndex.begin();
        it!=masterIndex.end(); ++it){
      Rcout << it->first << "\t" << it->second.start << "\t" << it->second.size << std::endl;
    }
  }
};

class expectedValues{
public:
  int binSize;
  std::string unit;
  doubleValues expectedValues;
  intDoubleValues chrScaleFactors;
  
  void read(std::ifstream& ifile){
    binSize = readInt(ifile);
    unit = readString(ifile);
    expectedValues.read(ifile);
    chrScaleFactors.read(ifile);
  }
};

class expectedValueVectors{
public:
  int nExpectedValueVectors;
  expectedValues* expValues;
  
  void read(std::ifstream& ifile){
    nExpectedValueVectors = readInt(ifile);
    expValues = new expectedValues[nExpectedValueVectors];
    for(int i=0; i<nExpectedValueVectors; i++){
      expValues[i].read(ifile);
    }
  }
};

class normExpectedValueVector{
public:
  std::string type;
  int binSize;
  std::string unit;
  doubleValues expectedValues;
  intDoubleValues chrNormalizationFactors;
  
  void read(std::ifstream& ifile){
    type= readString(ifile);
    binSize = readInt(ifile);
    unit= readString(ifile);
    expectedValues.read(ifile);
    chrNormalizationFactors.read(ifile);
  }
};

class normExpectedValueVectors{
public:
  int nNormExpectedValueVectors;
  normExpectedValueVector* normExpValVecs;
  
  void read(std::ifstream& ifile){
    nNormExpectedValueVectors = readInt(ifile);
    normExpValVecs = new normExpectedValueVector[nNormExpectedValueVectors];
    for(int i=0; i<nNormExpectedValueVectors; i++){
      normExpValVecs[i].read(ifile);
    }
  }
};

class normVectorEntry{
public:
  std::string type;
  int chrIdx;
  std::string unit;
  int binSize;
  long position;
  long nBytes;
  std::string key;
  
  void read(std::ifstream& ifile){
    type= readString(ifile);
    chrIdx = readInt(ifile);
    unit= readString(ifile);
    binSize = readInt(ifile);
    position = readLong(ifile);
    nBytes = version<9?(long)readInt(ifile):readLong(ifile);
    key= type+"_"+std::to_string(chrIdx)+"_"+unit+"_"+std::to_string(binSize);
  }
  
  doubleValues getValues(std::ifstream& ifile){
    doubleValues values;
    if(ifile.is_open()){
      ifile.clear();
      ifile.seekg(position, std::ios::beg);
      values.read(ifile);
    }
    return(values);
  }
};

class normVectorEntries{//normalization vector arrays
public:
  long nNormVectors;
  normVectorEntry* normVecs;
  
  void read(std::ifstream& ifile){
    if(ifile.is_open()){
      nNormVectors = version<9?(long)readInt(ifile):readLong(ifile);
      normVecs = new normVectorEntry[nNormVectors];
      for(int i=0; i<nNormVectors; i++){
        normVecs[i].read(ifile);
      }
    }
  }
  
  int find(std::string type, int chrIdx, std::string unit, int binSize){
    std::string key = type+"_"+std::to_string(chrIdx)+"_"+unit+"_"+std::to_string(binSize);
    for(int i=0; i<nNormVectors; i++){
      if(key==normVecs[i].key){
        return(i);
      }
    }
    return(-1);
  }
  
  doubleValues getValues(std::ifstream& ifile, int id){
    normVectorEntry nVE = normVecs[id];
    doubleValues values = nVE.getValues(ifile);
    //int adjustedStart = max(0, start - 1000);
    //int adjustedEnd = min(values.nValues, end+1000);
    //int startPosition = nVE.position + adjustedStart*sizeof(double);
    return(values);
  }
};

class HicFileFoot {
public:
  long nBytes;
  entry entries;
  expectedValueVectors expValVecs;
  normExpectedValueVectors normExpValVecs;
  normVectorEntries normVecs;
  
  void read(std::ifstream& ifile){
    ifile.clear() ;
    ifile.seekg(footerPosition, std::ios::beg);
    if(version<9){
      nBytes = (long) readInt(ifile);
    }else{
      nBytes = readLong(ifile);
    }
    if(debug){
      Rcout << "nBytes " << nBytes << std::endl;
    }
    entries.read(ifile);
    if(debug){
      entries.printEntries();
    }
    expValVecs.read(ifile);
    if(debug){
      Rcout << "nExpectedValueVectors = " << expValVecs.nExpectedValueVectors << std::endl;
    }
    normExpValVecs.read(ifile);
    if(debug){
      Rcout << "nNormExpectedValueVectors = " << normExpValVecs.nNormExpectedValueVectors << std::endl;
    }
    normVecs.read(ifile);
    if(debug){
      Rcout << "nNormVectors = " << normVecs.nNormVectors << std::endl;
    }
  }
};

class cellData{
public:
  int binX;
  long binY;
  float value;
  
  void read(std::istringstream& iss){
    binX = getInt(iss);
    binY = getLong(iss);
    value = getFloat(iss);
  }
  void print(){
    Rcout << "x:" << binX << "\ty:" << binY << "\t=" << value << std::endl;
  }
};

// decompress by zlib
int inflate(char *compr, int comprLen, std::string& uncompr){
  int ret;
  unsigned have;
  z_stream strm; /* decompression stream */
  unsigned char out[CHUNK];
  
  strm.zalloc = Z_NULL;
  strm.zfree = Z_NULL;
  strm.opaque = Z_NULL;
  strm.avail_in = 0;
  strm.next_in = Z_NULL;
  ret = inflateInit2(&strm, MAX_WBITS|32);
  if (ret != Z_OK)
    return ret;
  
  strm.next_in  = (unsigned char *)compr;
  strm.avail_in = comprLen;
  do{
    strm.next_out = out;
    strm.avail_out = CHUNK;
    
    ret = inflate(&strm, Z_NO_FLUSH);
    assert(ret != Z_STREAM_ERROR);  /* state not clobbered */
    switch (ret) {
    case Z_NEED_DICT:
      ret = Z_DATA_ERROR;     /* and fall through */
    case Z_DATA_ERROR:
    case Z_MEM_ERROR:
      (void)inflateEnd(&strm);
      return ret;
    }
    
    // add data to uncompr
    have = CHUNK - strm.avail_out;
    for(int i=0; i<have; i++){
      uncompr += out[i];
    }
  }while(strm.avail_out == 0);
  
  (void)inflateEnd(&strm);
  return ret == Z_STREAM_END ? Z_OK : Z_DATA_ERROR;
}

class cellDatas{
public:
  int cellCount;
  cellData* cells;
  
  void read(std::ifstream& ifile, long blockPosition, int blockSize){
    if(ifile.is_open()){
      ifile.clear();
      ifile.seekg(blockPosition, std::ios::beg);
      char* compr = new char[blockSize];
      std::string uncompr;
      ifile.read(compr, blockSize);
      int decompression_result = inflate(compr, blockSize, uncompr);
      if(decompression_result!=F_OK){
        stop("Something wrong when inflate the .hic file.");
      }
      std::istringstream ifs(uncompr);
      cellCount = getInt(ifs);
      cells = new cellData[cellCount];
      if(version<7){
        for(int i=0; i<cellCount; i++){
          cells[i].read(ifs);
        }
      }else{
        int binXoffset = getInt(ifs);
        int binYoffset = getInt(ifs);
        
        bool useFloatContact = getByte(ifs) == 1;
        bool useIntXPos = version<9 ? false: getByte(ifs) == 1;
        bool useIntYPos = version<9 ? false: getByte(ifs) == 1;
        Byte type = getByte(ifs);
        
        if(type==1){
          int rowCount = useIntYPos ? getInt(ifs) : (int)getShort(ifs);
          int k=0;
          for(int i=0; i<rowCount; i++){
            int dy = useIntYPos ? getInt(ifs) : (int)getShort(ifs);
            long binY = (long)(binYoffset + dy);
            int colCount = useIntXPos ? getInt(ifs) : (int)getShort(ifs);
            for(int j=0; j<colCount; j++){
              int dx = useIntXPos ? getInt(ifs) : (int)getShort(ifs);
              int binX = (int)(binXoffset + dx);
              float counts = useFloatContact ? getFloat(ifs) : (float)getShort(ifs);
              cells[k].binX = binX;
              cells[k].binY = binY;
              cells[k].value = counts;
              k++;
            }
          }
        }else{
          if(type==2){
            int nPts = getInt(ifs);
            short w = getShort(ifs);
            int k = 0;
            for(int i=0; i<nPts; i++){
              //int idx = (p.y - binOffset2) * w + (p.x - binOffset1);
              int row = floor(i/w);
              int col = i - row*w;
              int bin1 = binXoffset + col;
              int bin2 = binYoffset + row;
              if(useFloatContact){
                float counts = getFloat(ifs);
                if(counts==counts){//not NAN
                  cells[k].binX = bin1;
                  cells[k].binY = (long) bin2;
                  cells[k].value = counts;
                  k++;
                }
              }else{
                short counts = getShort(ifs);
                if(counts!=-32768){
                  cells[k].binX = bin1;
                  cells[k].binY = (long) bin2;
                  cells[k].value = (float)counts;
                  k++;
                }
              }
            }
          }else{
            stop("invalid block type.");
          }
        }
      }
      if(debug){
        this->print();
      }
    }
  }
  void print(){
    Rcout << "cellCount= " << cellCount << std::endl;
    for(int i=0; i<std::min(cellCount, 10); i++){
      cells[i].print();
    }
  }
};

class blockIndex{
public:
  int blockID;
  long blockPosition;
  int blockSize;
  
  void read(std::ifstream& ifile){
    blockID = readInt(ifile);
    blockPosition = readLong(ifile);
    blockSize = readInt(ifile);
  }
  void print(){
    Rcout << "blockID=" << blockID << "\tblockPosition=" << blockPosition;
    Rcout << "\tblockSize=" << blockSize << std::endl;
  }
};

class blockIndexs{
public:
  int blockCount;
  blockIndex* blockIdx;
  
  void read(std::ifstream& ifile){
    blockCount = readInt(ifile);
    blockIdx = new blockIndex[blockCount];
    for(int i=0; i<blockCount; i++){
      blockIdx[i].read(ifile);
    }
  }
  
  void print(){
    Rcout << "blockCount = " << blockCount << std::endl;
    for(int i=0; i<blockCount; i++){
      blockIdx[i].print();
    }
  }
};

class resolution{
public:
  std::string unit;
  int resIdx;
  float sumCounts;
  float occupiedCellCount;
  float stdDev;
  float percent95;
  int binSize;
  int blockBinCount;
  int blockColumnCount;
  blockIndexs blocks;
  
  void read(std::ifstream& ifile){
    unit= readString(ifile);
    resIdx = readInt(ifile);
    sumCounts = readFloat(ifile);
    occupiedCellCount = readFloat(ifile);
    stdDev = readFloat(ifile);
    percent95 = readFloat(ifile);
    binSize = readInt(ifile);
    blockBinCount = readInt(ifile);
    blockColumnCount = readInt(ifile);
    blocks.read(ifile);
  }
  
  void print(){
    Rcout << "unit = " << unit << "\tresIdx = " << resIdx;
    Rcout << "\tbinSize = " << binSize << "\tblockBinCount = " << blockBinCount;
    Rcout << "\tblockColumnCount = " << blockColumnCount << std::endl;
    blocks.print();
  }
};

class resolutions{
public:
  int nResolutions;
  resolution* res;
  
  void read(std::ifstream& ifile){
    nResolutions = readInt(ifile);
    res = new resolution[nResolutions];
    for(int i=0; i<nResolutions; i++){
      res[i].read(ifile);
    }
  }
  
  void print(){
    Rcout << "nResolutions = " << nResolutions << std::endl;
    for(int i=0; i<nResolutions; i++){
      res[i].print();
    }
  }
};

class body{
public:
  int chr1Idx;
  int chr2Idx;
  resolutions data;
  
  void read(std::ifstream& ifile, long position, int size){
    if(ifile.is_open()){
      ifile.clear();
      ifile.seekg(position, std::ios::beg);
      chr1Idx = readInt(ifile);
      chr2Idx = readInt(ifile);
      data.read(ifile);
    }
  }
  
  void print(){
    Rcout << "id1=" << chr1Idx << "\tid2="<< chr2Idx << std::endl;
    data.print();
  }
  
  int getResolution(int binSize){
    for(int i=0; i<data.nResolutions; i++){
      if(data.res[i].binSize==binSize){
        return(i);
      }
    }
    return(-1);
  }
};

std::vector<cellData> getByCoor(std::string seq1, int start1, int end1, 
                           std::string seq2, int start2, int end2,
                           int binSize, std::string normalization, std::string unit,
                           HicFileHead& hicHead,
                           HicFileFoot& hicFoot, 
                           std::ifstream& hicfile){
  int chr1 = hicHead.chromosomeIndex.find(seq1);
  int chr2 = hicHead.chromosomeIndex.find(seq2);
  //chr1 must no greater than chr2
  //if(chr1==chr2) region1 must be upstream of region2
  bool transpose = (chr1 > chr2) || (chr1 == chr2 && start1 >= end2);
  if(transpose){
    int tmp = start1;
    start1 = start2;
    start2 = tmp;
    tmp = end1;
    end1 = end2;
    end2 = tmp;
  }
  
  int x1 = start1 / binSize;
  int x2 = end1 / binSize;
  int y1 = start2 / binSize;
  int y2 = end2 / binSize;
  
  if(debug){
    Rcout << x1 << "\t" << x2 << "\t" << y1 << "\t" << y2 << std::endl;
  }
  
  std::string id = std::to_string(chr1)+"_"+std::to_string(chr2);
  std::map<std::string, mIdx>::iterator it = hicFoot.entries.masterIndex.find(id);
  std::vector<cellData> contactRecords;
  if(it != hicFoot.entries.masterIndex.end()){
    if(debug){
      Rcout << it->first << "\tstart=" << it->second.start;
      Rcout << "\tsize=" << it->second.size << std::endl;
    }
    body bd;
    bd.read(hicfile, it->second.start, it->second.size);
    int resId = bd.getResolution(binSize);
    if(resId==-1){
      stop("The resolution is not available");
    }
    resolution res = bd.data.res[resId];
    // determine blockIDs to be read;
    for(int i=0; i<res.blocks.blockCount ; i++){
      //res.blocks.blockIdx[i].print();
      cellDatas cd;
      cd.read(hicfile, 
              res.blocks.blockIdx[i].blockPosition,
              res.blocks.blockIdx[i].blockSize);
      bool isNorm = normalization!="" && normalization!="NONE";
      doubleValues normVector1;
      doubleValues normVector2;
      if(isNorm){
        //getNormalizationVector(type, chr, unit, binSize)
        int nv1 = hicFoot.normVecs.find(normalization, chr1, unit, binSize);
        int nv2 = hicFoot.normVecs.find(normalization, chr2, unit, binSize);
        
        if(nv1!=-1 && nv2!=-1){
          normVector1 = hicFoot.normVecs.getValues(hicfile, nv1);
          normVector2 = hicFoot.normVecs.getValues(hicfile, nv2);
        }else{
          isNorm = false;
        }
      }
      for(int i=0; i<cd.cellCount; i++){
        cellData rec = cd.cells[i];
        if(rec.binX>=x1 && rec.binX<x2 && rec.binY>=y1 && rec.binY<y2){
          if(isNorm){
            cellData rec_norm;
            rec_norm.binX = rec.binX;
            rec_norm.binY = rec.binY;
            double nvnv = normVector1.values[rec_norm.binX] *
              normVector2.values[rec_norm.binY];
            if(nvnv!=0 && nvnv==nvnv){
              rec_norm.value = rec.value/nvnv;
              contactRecords.push_back(rec_norm);
            }
          }else{
            contactRecords.push_back(rec);
          }
        }
      }
    }
  }
  return(contactRecords);
}

// [[Rcpp::export]]
IntegerVector listResolutions(CharacterVector hicfilename,
                              CharacterVector unit){
  IntegerVector resolutions;
  std::string filename = as<std::string>(hicfilename);
  std::ifstream hicfile(filename, std::ios::binary);
  std::string type = as<std::string>(unit);
  HicFileHead hicHead;
  hicHead.readHeader(hicfile);
  if(type=="BP"){
    for(int i=0; i<hicHead.bpResolutions.nValues; i++){
      resolutions.push_back(hicHead.bpResolutions.values[i]);
    }
  }else{
    if(type=="FRAG"){
      for(int i=0; i<hicHead.fragResolutions.nValues; i++){
        resolutions.push_back(hicHead.fragResolutions.values[i]);
      }
    }
  }
  return(resolutions);
}
// [[Rcpp::export]]
CharacterVector listUnits(CharacterVector hicfilename){
  CharacterVector units;
  std::string filename = as<std::string>(hicfilename);
  std::ifstream hicfile(filename, std::ios::binary);
  HicFileHead hicHead;
  hicHead.readHeader(hicfile);
  if(hicHead.bpResolutions.nValues>0){
    units.push_back("BP");
  }
  if(hicHead.fragResolutions.nValues>0){
    units.push_back("FRAG");
  }
  return(units);
}
// [[Rcpp::export]]
CharacterVector listNormalizations(CharacterVector hicfilename){
  CharacterVector norm;
  std::string filename = as<std::string>(hicfilename);
  std::ifstream hicfile(filename, std::ios::binary);
  HicFileHead hicHead;
  hicHead.readHeader(hicfile);
  HicFileFoot hicFoot;
  hicFoot.read(hicfile);
  for(int i=0; i<hicFoot.normVecs.nNormVectors; i++){
    norm.push_back(hicFoot.normVecs.normVecs[i].type);
  }
  norm=unique(norm);
  return(norm);
}
// [[Rcpp::export]]
CharacterVector listChroms(CharacterVector hicfilename){
  CharacterVector chrs;
  std::string filename = as<std::string>(hicfilename);
  std::ifstream hicfile(filename, std::ios::binary);
  HicFileHead hicHead;
  hicHead.readHeader(hicfile);
  for(int i=0; i<hicHead.chromosomeIndex.nChrs; i++){
    chrs.push_back(hicHead.chromosomeIndex.chromosomeIndexMap[i].chrName);
  }
  chrs=unique(chrs);
  return(chrs);
}
// [[Rcpp::export]]
DataFrame getContactRecords(CharacterVector hicfilename,
                            CharacterVector qname1,
                            IntegerVector start1,
                            IntegerVector end1,
                            CharacterVector qname2,
                            IntegerVector start2, 
                            IntegerVector end2,
                            IntegerVector binSize,
                            CharacterVector normalization,
                            CharacterVector unit) {
  std::string filename = as<std::string>(hicfilename);
  std::string seq1 = as<std::string>(qname1);
  int start_pos1 = as<int>(start1);
  int end_pos1 = as<int>(end1);
  std::string seq2 = as<std::string>(qname2);
  int start_pos2 = as<int>(start2);
  int end_pos2 = as<int>(end2);
  int bin_size = as<int>(binSize);
  std::string norm = as<std::string>(normalization);
  std::string type = as<std::string>(unit);
  std::ifstream hicfile(filename, std::ios::binary);
  HicFileHead hicHead;
  hicHead.readHeader(hicfile);
  HicFileFoot hicFoot;
  hicFoot.read(hicfile);
  std::vector<cellData> rec = getByCoor(seq1, start_pos1, end_pos1, 
                                   seq2, start_pos2, end_pos2,
                                   bin_size, norm, type,
                                   hicHead,
                                   hicFoot, 
                                   hicfile);
  if(debug){
    for(std::vector<cellData>::iterator it=rec.begin(); it!=rec.end(); ++it){
      it->print();
    }
  }
  
  hicfile.close();
  if(rec.size()<=1){
    return(DataFrame::create());
  }
  CharacterVector seqn1;
  IntegerVector a1;
  IntegerVector b1;
  CharacterVector seqn2;
  IntegerVector a2;
  IntegerVector b2;
  NumericVector score;
  for(std::vector<cellData>::iterator it=rec.begin(); it!=rec.end(); ++it){
    seqn1.push_back(seq1);
    a1.push_back(it->binX*bin_size);
    b1.push_back((it->binX+1)*bin_size);
    seqn2.push_back(seq2);
    a2.push_back(it->binY*bin_size);
    b2.push_back((it->binY+1)*bin_size);
    score.push_back(it->value);
  }
  return(DataFrame::create(Named("seq1")=seqn1, Named("start1")=a1, Named("end1")=b1,
                           Named("seq2")=seqn2, Named("start2")=a2, Named("end1")=b2,
                           Named("score")=score));
}

