 
#ifndef MTXREADER_H
#define MTXREADER_H

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>
#include <cassert>
#include <cstring>

//#include "mvm.hpp"
//#include"mtxReader.hpp"
//#include "MPI.hpp"
#include <stdlib.h>
#include <time.h>
#include <getopt.h>
#include <float.h>
#include <math.h>
#include <mm_malloc.h>

using namespace std;

struct EdgeE
{
    int head;
    int id;         // Edge tail
    float weight;  // Edge weight
};
struct Edge
{
    long long id;         // Edge tail
    float weight;  // Edge weight
};

class CSR
{
    public:
    int nVer;       // number of vertices 
    int nEdge;      // number of edges
    int maxDeg;
    long long* verPtr;    // vertex pointer array of size nVer+1
    Edge* verInd;   // Edge array
    bool bipartite; // general graph or bipartite graph
    int rVer;       // The number of vertices on Right for bipartite graph;
    int lVer;       // The number of vertices on left for bipartite graph;
    float* verWt;
    
    bool readMtxG(char * filename); // reading as a general graph
    bool readMtxB(char * filename); // reading as a bipartite graph
    bool mtxG2csrbin(char* filename, char* outfile); //converting mtx to binary file
    bool mtxB2csrbin(char* filename, char* outfile); //converting mtx to binary file
    bool readCSRbin(char* filename, int opt); // reading binary file of csr format
    bool readCSRbinBipartite(char* filename, int opt); // reading binary file of csr format
    
    CSR():nVer(0),nEdge(0),verPtr(NULL),verInd(NULL),bipartite(false){}
    CSR(long long n, long long m)
    {
        nVer=n;
        nEdge=m;
        bipartite=false;
        verPtr=(long long*)_mm_malloc((n+1)*sizeof(long long),64);
        verInd=(Edge*)_mm_malloc(m*sizeof(Edge),64);
    }
    ~CSR()
    {
        if(verPtr!=NULL)
            delete verPtr;

        if(verInd!=NULL)
            delete verInd;
    }

};

class TCSR
{
    public:
    long long nVer;       // number of vertices
    long long nEdge;      // number of edges
    //long long nEdge;
    long long maxDeg;
    long long* verPtr;    // vertex pointer array of size nVer+1
    long long* verInd;    // Edge Indices
    float* verWt;         // Edge weights
    bool bipartite;       // general graph or bipartite graph
    //long long sVer;       // The number of vertices on left for bipartite graph;
    
    bool readMtxG(char * filename); // reading as a general graph
    bool readMtxB(char * filename); // reading as a bipartite graph
    bool mtxG2csrbin(char* filename, char* outfile); //converting mtx to binary file
    bool mtxB2csrbin(char* filename, char* outfile); //converting mtx to binary file
    bool readCSRbin(char* filename); // reading binary file of csr format

    TCSR():nVer(0),nEdge(0),verPtr(NULL),verInd(NULL),verWt(NULL),bipartite(false){}

};


class DCSR
{
public:
long long gnVer;
long long gnEdge;
long long nVer;       // number of vertices
long long nEdge;      // number of edges
long long maxDeg;                          //max degree
long long nProc;      // number of process
//long long rank;       //rank
long long startVer;  // start Node
long long endVer;    // end Node
long long VerPerProc;   // number per process
long long* verPtr;    // vertex pointer array of size nVer+1
long long* verWt; 
Edge* verInd;         // Edge array

// long long* bVer;      // b Values                
// long long* sVer;      // s values
};

class MCSR
{
    public:
    long long nVer;       // number of vertices
    long long nEdge;      // number of edges
    //long long nEdge;
    long long* verPtr;    // vertex pointer array of size nVer+1
    long long* verInd;    // Edge indices
    float* verWt;   // Edge Weight
};

#endif //MTXREADER_H
