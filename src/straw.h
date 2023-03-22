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
#ifndef STRAW_H
#define STRAW_H

#include <fstream>
#include <set>
#include <vector>
#include <map>

// pointer structure for reading blocks or matrices, holds the size and position
struct indexEntry {
    int64_t size;
    int64_t position;
};

// sparse matrixType entry
struct contactRecord {
  int32_t binX;
  int32_t binY;
  float counts;
};

struct footerInfo {
    int32_t resolution;
    bool foundFooter;
    int32_t version;
    int32_t c1;
    int32_t c2;
    int32_t numBins1;
    int32_t numBins2;
    int64_t myFilePos;
    std::string unit;
    std::string norm;
    std::string matrixType;
    std::vector<double> c1Norm;
    std::vector<double> c2Norm;
    std::vector<double> expectedValues;
};

// chromosome
struct chromosome {
    std::string name;
    int32_t index;
    int64_t length;
};

#endif
