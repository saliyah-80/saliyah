//
//  mvm.hpp
//  
//
//  Created by Saliyah on 25/03/1444 AH.
//

#ifndef mvm_hpp
#define mvm_hpp

//#include "SR.hpp"
#include <stdio.h>
#include <omp.h>
#include <cmath>
#include <iostream>
#include "mtxReader.hpp"
#include "pcl_stack.h"
#include <immintrin.h>
#include <getopt.h>
#include <float.h>
#include <stdexcept>
#include <limits>
#include "mic_utils.h"

using namespace std;

#define CHUNK 256
#define DYN
#define BSIZE 256
#define LBUF
#define MAX_VAL FLT_MAX/2.0





struct Param
{
   int cstart;
   int cend;
   int ostart;
   int oend;
};

struct Info
{
   int id;
   float weight;
};

class Node
{
   public:
   int maxSize;
   int curSize;
   //Info* heap;
   Info minEntry;


   void print();


   void Add(float wt, int id);

   
};

struct Walker
{
    int midx1;
    int midx2;
    float W1;
    float W2;
    float* M1;
    float* M2;
};





//void mvmHalfPar(DCSR* G, long long* mate, double* vtxWght, bool verbose);
void mvmHalfPar(DCSR* G, long long* mate, double* vtxWght, int *nlocks, Node* M, long long* start, long long* end, bool verbose,char* mark, int type);
//void mvmHalfPar(DCSR* G, long long* mate,long long *b, double* vtxWght, int *nlocks, Node* S, long long* start, long long* end, bool verbose);
//void bSuitorBPQD(DCSR* G, long long *b, int *nlocks, Node* S, long long* start,
        //long long* end, long long stepM, long long* sVer, long long npart, bool nsort,  bool //verbose);


//void suitorPar(DCSR* G, long long* mate, bool verbose);
void suitorPar(DCSR* G,long long *b, long long* mate,double* vtxWght, int *nlocks, Node* S, long long* start, long long* end, long long stepM, bool verbose,char* mark);

bool comparator(Edge left, Edge right);
bool predicate(Edge left);

int custom(Edge* verInd,  long long start, long long end, long long step);


bool verifyMatching(DCSR* G, long long* mate, double* vtxWght);
bool verifyMatchingEdgeWght(DCSR* G, long long* mate);
bool vtxWghtToEdgeWght(DCSR* G, double* vtxWght);
void printMatching( long long n);

bool comparator(Edge left, Edge right);

const long long cNullItm = std::numeric_limits<long long>::max() - 1;

struct matching_parameters
{
   char* problemname;
   char* wghtfilename;
   //char* bFileName;
   int algorithm;
  // int b;
   int step;
   bool verbose;
   bool bipartite;

   matching_parameters();
   void usage();
   bool parse(int argc, char** argv);
};
#endif /* mvm_hpp */
