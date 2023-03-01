//
//  mvm.cpp
//  
//
//  Created by Saliyah on 25/03/1444 AH.
//
/********************************/ /* Include Statements */ /********************************/

#include "mvm.hpp"
#include"mtxReader.hpp"
//#include "MPI.hpp"
#include <stdlib.h>
#include <time.h>
#include <getopt.h>
#include <float.h>
#include <math.h>
#include <mm_malloc.h>
#include <cstring>
#include <xmmintrin.h> 
#include <mm_malloc.h>



//#include <parallel/algorithm>

using namespace std;
//Reading graph

//bool GraphReading(DCSR* G, char* pf)
bool distGraphReading(DCSR* DG, char* pf, char* bf,unsigned long np, MPI_Comm comm)

{
    int rank;
    MPI_Status stat;
    MPI_Comm_rank(comm, &rank);
    //MPI_Request request;
    long long start,offset,count;
    long long nodePerProcess;
    long long gnver,gnedge,maxdeg;
    ifstream inf;
    inf.open(pf,ios::in|ios::binary);
    if(inf.is_open())
    {
        inf.read((char*)&gnver,sizeof(long long));
        inf.read((char*)&gnedge,sizeof(long long));
        inf.read((char*)&maxdeg,sizeof(long long));
        nodePerProcess=ceil(gnver/(np*1.0));
        
        if(rank==0) // Extra adjustment for the master node because it will have the excess nodes
        {
            DG->gnVer=gnver;
            DG->gnEdge=gnedge;
            DG->nVer=gnver-(nodePerProcess*(np-1));

            //Copy vertex Pointers
            DG->verPtr=(long long*)_mm_malloc((DG->nVer+1)*sizeof(long long),64);
            inf.read((char*)&DG->verPtr[0],sizeof(long long)*(DG->nVer+1));
            inf.seekg((3+gnver+1)*sizeof(long long),inf.beg);
           

            DG->nEdge=DG->verPtr[DG->nVer];

            long long* verInd=(long long*)_mm_malloc(DG->nEdge*sizeof(long long),64);
            float* verWt=(float*)_mm_malloc(DG->nEdge*sizeof(float),64);

            inf.read((char*)&verInd[0],sizeof(long long)*DG->nEdge);
            inf.seekg((3+gnver+1+gnedge)*sizeof(long long),inf.beg);
            //cout<<"Rank 0, after verInd:"<<DG->nEdge<<" "<<inf.tellg()<<endl;
            inf.read((char*)&verWt[0],sizeof(float)*DG->nEdge);

            DG->nProc=np;
            DG->rank=rank;
            DG->startNode=0;
            DG->endNode=gnver-(nodePerProcess*(np-1))-1;
            DG->nPerProc=nodePerProcess;

            //Copy Edge Indices and Edge Weights
            DG->verInd=(Edge*)_mm_malloc(DG->nEdge*sizeof(Edge),64);
            
            #pragma omp parallel for
            for(long long i=0;i<DG->nEdge;i++)
            {
                DG->verInd[i].id=verInd[i];
                DG->verInd[i].weight=(float)verWt[i];
            }
            
            delete verInd;
            delete verWt;
           // _mm_free(verInd);
           // _mm_free(verWt);

/*
            //Copy lVer value
            DG->lVer=(long long*)_mm_malloc(DG->nVer*sizeof(long long),64);
            if(b>0)
            {
                #pragma omp parallel for
                for(long long i=0;i<DG->nVer;i++)
                {
                    long long deg=DG->verPtr[i+1]-DG->verPtr[i];
                    if(deg>b)
                        //DG->bVer[i]=b; //Matching
                        DG->lVer[i]=deg-b; // Edge Cover
                    else
                    {
                        /*if(deg>0)
                            DG->bVer[i]=deg;
                        else
                            DG->bVer[i]=1; */ //Matching
/*
                        DG->lVer[i]=1; // Edge Cover
                    }
                }
            }
            else
            {
                #pragma omp parallel for
                for(long long i=0;i<DG->nVer;i++)
                {
                    long long deg=DG->verPtr[i+1]-DG->verPtr[i];
                    if(deg<=2)
                        DG->lVer[i]=1;
                    else
                        DG->lVer[i]=(long long)(deg/(log2(deg)*1.0));
                }
            }
            //Copy sV value
            DG->rVer=(long long*)_mm_malloc(DG->nVer*sizeof(long long),64);
            
            if(s>0)
            {
                #pragma omp parallel for
                for(long long i=0;i<DG->nVer;i++)
                    DG->rVer[i]=s;
            }
            else
            {
                #pragma omp parallel for
                for(long long i=0;i<DG->nVer;i++)
                {
                    if(DG->lVer[i]>=1)
                        DG->rVer[i]=DG->lVer[i]/2;
                    else
                        DG->rVer[i]=1;
                }
            }
        }
        else   // for other nodes
        {
            start=gnver-(nodePerProcess*(np-1));
            start=start+(rank-1)*nodePerProcess;
                
            DG->gnVer=gnver;
            DG->gnEdge=gnedge;
            DG->nVer=nodePerProcess;
            DG->nProc=np;
            DG->rank=rank;
            DG->nPerProc=nodePerProcess;
            
*/
            //Copy vertex Pointers
            DG->verPtr=(long long*)_mm_malloc((DG->nVer+1)*sizeof(long long),64);
            inf.seekg((3+start)*sizeof(long long),inf.beg);
            //if(rank==5)
                //cout<<"Rank 5 initially :"<<start<<" "<<inf.tellg()<<endl;
            inf.read((char*)&DG->verPtr[0],sizeof(long long)*(DG->nVer+1));

            offset=DG->verPtr[0];
            
            #pragma omp parallel for
            for(long long i=0;i<DG->nVer+1;i++)
                DG->verPtr[i]=DG->verPtr[i]-offset;
            
            DG->startNode=start;
            DG->endNode=start+nodePerProcess-1;
            cout<<"R 10 "<<endl;
            //Copy Edge Indices and Edge Weights
            DG->nEdge=DG->verPtr[DG->nVer];
            DG->verInd=(Edge*)_mm_malloc(DG->nEdge*sizeof(Edge),64);
            
            long long* verInd=(long long*)_mm_malloc(DG->nEdge*sizeof(long long),64);
            float* verWt=(float*)_mm_malloc(DG->nEdge*sizeof(float),64);
            
            inf.seekg((3+gnver+1+offset)*sizeof(long long),inf.beg);

            //if(rank==5)
                //cout<<"Rank 5 after verptr :"<<offset<<" "<<inf.tellg()<<endl;

            inf.read((char*)&verInd[0],sizeof(long long)*DG->nEdge);
            inf.seekg((3+gnver+1+gnedge)*sizeof(long long),inf.beg);
            inf.seekg(offset*sizeof(float),inf.cur);

            //if(rank==5)
                //cout<<"Rank 5 after verind :"<<DG->nEdge<<" "<<inf.tellg()<<endl;
            inf.read((char*)&verWt[0],sizeof(float)*DG->nEdge);

            #pragma omp parallel for
            for(long long i=0;i<DG->nEdge;i++)
            {
                DG->verInd[i].id=verInd[i];
                DG->verInd[i].weight=verWt[i];
            }
            delete verInd;
            delete verWt;
           // _mm_free(verInd);
           // _mm_free(verWt);
            
/*                /// recieve the b values from the master..
            DG->lVer=(long long*)_mm_malloc(DG->nVer*sizeof(long long),64);
            if(b>0)
            {
                #pragma omp parallel for
                for(long long i=0;i<DG->nVer;i++)
                {
                    long long deg=DG->verPtr[i+1]-DG->verPtr[i];
                    if(deg>b)
                        //DG->bVer[i]=b; //Matching
                        DG->lVer[i]=deg-b; // Edge Cover
                    else
                    {
                        if(deg>0)
                            DG->lVer[i]=deg;
                        else
                            DG->lVer[i]=1; // Matching
                        
                        //DG->bVer[i]=1; // Edge Cover
                    }
                }
            }
            else
            {
                cout<<"R 12 "<<endl;
                #pragma omp parallel for
                for(long long i=0;i<DG->nVer;i++)
                {
                    long long deg=DG->verPtr[i+1]-DG->verPtr[i];
                    if(deg<=2)
                        DG->lVer[i]=1;
                    else
                        DG->lVer[i]=(long long)(deg/(log2(deg)*1.0));
                }
            }
            
            //Copy sV value
            DG->rVer=(long long*)_mm_malloc(DG->nVer*sizeof(long long),64);
            
            if(s>0)
            {
                #pragma omp parallel for
                for(long long i=0;i<DG->nVer;i++)
                    DG->rVer[i]=s;
            }
            else
            {
                #pragma omp parallel for
                for(long long i=0;i<DG->nVer;i++)
                {
                    if(DG->lVer[i]>=1)
                        DG->rVer[i]=DG->lVer[i]/2;
                    else
                        DG->rVer[i]=1;
                }
            }
         }
    
    
    }
    inf.close();
    return 1;

}


*/


{
 
    return G->readMtxG(pf);    // return all data of graph # ver , edges , verind , verptr
}
  
bool VtxWghtReading(double* vtxWght, char* filename)
{
    
    int count,i,n;
    ifstream inf;

    inf.open(filename, ios::in);

    if(inf.is_open())
    {
        
        inf>>inp;    // # of Ver
        count=inp;
		i=0;
        while(count>0) // count # of Ver
        {     
            inf>>f;        // red value then store it in vtxWght array
			vtxWght[i]=f;  // store
			i++;
            count--; 
        }     
        inf.close(); 
    }
    else return false;
   
   if(i == inp)    // read all Ver
   return true;
   else 
	return false;
    
  }


/*////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/

bool verifyMatching(DCSR* G, long long* mate, double* vtxWght)
{
  
  bool flag =false;
  long long n=G->nVer;
  double weight=0.0;     //wight
  long long card =0;     // cardinalty
  for (long long s = 0; s < n; ++s) {
    long long t = mate[s];
    if (t != cNullItm) {
      assert(t < n);
      assert(s != t);
      /*if(mate[t] != s)
      {
      	//cout << "the problem is the mate of " <<s<< " is  "<<mateArr[s]<< " but the mate of "<<t<< " is "<< mateArr[t]<<endl;
		
	  }*/
      assert(mate[t] == s);
      // check if there is an edge between s and t.
      flag=false;
	  for(long long j=G->verPtr[s];j<G->verPtr[s+1];j++)
        {
           long long a=G->verInd[j].id;
		   if(a == t)
		   {
			   flag =true;
			   weight+=vtxWght[s]; 
			   card++;
			   break;
		   }
		}
		if(!flag)
			return false;
	  
    }
  }
   cout<<"Matching weight: "<<weight<<endl;
   cout<<"Cardinality: "<<card/2<<endl;
    return true;
    

}

bool verifyMatchingEdgeWght(DCSR* G, long long* mate)
{
  
  bool flag =false;
  long long n=G->nVer;
  //Edge* verInd=G->verInd;
  double weight=0.0;
  long long card=0;
  for (long long s = 0; s < n; ++s) {
    long long t = mate[s];
    if (t != cNullItm) {
      assert(t < n);
      assert(s != t);
	  
      if(mate[t] != s)
      {
      	cout << "the problem is the mate of " <<s<< " is  "<<mate[s]<< " but the mate of "<<t<< " is "<< mate[t]<<endl;
		for(long long j=G->verPtr[s];j<G->verPtr[s+1];j++)
        {
           long long a=G->verInd[j].id;
		   if(a == t)
		   {
			   
			   cout << "weight of s, t " <<G->verInd[j].weight<<endl;
			   break;
		   }
		}
		long long mt =mate[t];
		for(long long j=G->verPtr[t];j<G->verPtr[t+1];j++)
        {
           long long a=G->verInd[j].id;
		   if(a == mt)
		   {
			   
			   cout << "weight of mt, t " <<G->verInd[j].weight<<endl;
			   break;
		   }
		}
	  }
	  
      assert(mate[t] == s);
      // check if there is an edge between s and t.
	  flag=false;
	  for(long long j=G->verPtr[s];j<G->verPtr[s+1];j++)
        {
           long long a=G->verInd[j].id;
		   if(a == t)
		   {
			   flag =true;
			   weight+=G->verInd[j].weight;
			   card++;
			   break;
		   }
		}
		
		if(!flag)
			return false;
        
    }
  }
   cout<<"Matching Weight: "<<weight/2.0<<endl;
    cout<<"Cardinality: "<<card/2<<endl;
    return true;
    

}

bool vtxWghtToEdgeWght(DCSR* G, double* vtxWght)
{
  
  long long n=G->nVer;
  //Edge* verInd=G->verInd;
  long long count =0;
  for (long long s = 0; s < n; ++s) {
    for(long long j=G->verPtr[s];j<G->verPtr[s+1];j++)
    {
		long long t=G->verInd[j].id;
		//if (s<t) 
		//{
			G->verInd[j].weight=vtxWght[s]+vtxWght[t];
			count++;
		//}
	}
  }
  if(count/2==G->nEdge)
    return true;
  else
	  return false;
}

//bmatching_parameters::bmatching_parameters()
  //  :problemname(NULL),bFileName(NULL),b(10),algorithm(8),step(3),sstep(1),nsort(false),verbose(false),hybrid(true),nstep(1)
matching_parameters::matching_parameters()
    :problemname(NULL),algorithm(2),verbose(false),wghtfilename(NULL)
{}

void matching_parameters::usage()
{
    const char *params =
    "\n\n"
    "Usage: %s -f <problemname>  -w <wghtfilename> -t <type>  -v \n\n"
    "    -f problemname  : file containing graph. Currently inputs .mtx files\n"
	"	 -w wghtfilename  : file containing weights of vertices\n"
    "    -t type         : Which algorithm to use\n"
    "                           1= 1/2-Approx MVM Algorithm  \n"
    "                           2= Suitor Algorithm  \n";
    fprintf(stderr, params);
}

bool matching_parameters::parse(int argc, char** argv)
{
    static struct option long_options[]=
    {
        // These options don't take extra arguments
        {"verbose", no_argument, NULL, 'v'},
        {"help", no_argument, NULL, 'h'},
        
        // These do
        {"problem", required_argument, NULL, 'f'},
		{"weight", required_argument, NULL, 'w'},
        {"type", required_argument, NULL, 't'},
        {NULL, no_argument, NULL, 0}
    };

    static const char *opt_string="vhmf:w:t:";
    int opt, longindex;
    opt=getopt_long(argc,argv,opt_string,long_options,&longindex);
    while(opt != -1)
    {
        switch(opt)
        {
            case 'v': verbose=true;
                      break;

            case 'h': usage();
                      return false;
                      break;
                
            case 'm': hybrid=false;
                      break;
                
            case 'f': problemname=optarg;
                      if(problemname==NULL)
                      {
                        cerr<<"Problem file is not speficied"<<endl;
                        return false;
                      }
                      break;
					  
		    case 'w': wghtfilename=optarg;
                      if(wghtfilename==NULL)
                      {
                        cerr<<"weight file is not speficied"<<endl;
                        return false;
                      }
                      break;

            case 't': algorithm=atoi(optarg);
                      if(algorithm < 1 || algorithm > 2)
                      {
                        cerr<<"The algorithm type should be between 1 and 2, please type -h to see the details"<<endl;
                        return false;
                      }
                      break;

        
        }
        opt=getopt_long(argc,argv,opt_string,long_options,&longindex);
    }

    return true;
}
int main(int argc, char** argv)
{
    MPI_Request* ssreq;
    MPI_Request* rsreq;
    MPI_Request* sdreq;
    MPI_Request* rdreq;
    MPI_Status* status;
    
   /* MPI_Request* ssreq=(MPI_Request*)_mm_malloc(G->nProc*sizeof(MPI_Request),64);
    MPI_Request* rsreq=(MPI_Request*)_mm_malloc(G->nProc*sizeof(MPI_Request),64);
    MPI_Request* sdreq=(MPI_Request*)_mm_malloc(G->nProc*sizeof(MPI_Request),64);
    MPI_Request* rdreq=(MPI_Request*)_mm_malloc(G->nProc*sizeof(MPI_Request),64);
    MPI_Status* status=(MPI_Status*)_mm_malloc(G->nProc*sizeof(MPI_Status),64);
    */
    
    cout<<"Done 1..!!"<<endl;
    matching_parameters opts;
    if(!opts.parse(argc,argv))
    {
        cout<<"Argument Parsing Is Not Correctly Done..!!"<<endl;
        return -1;
    }
	cout<<"Argument Parsing Is Done..!!"<<endl;

    /*********MPI INIT ****/
    int rank, len, rc, numProc;
    char hostname[MPI_MAX_PROCESSOR_NAME];
  
     rc = MPI_Init(&argc,&argv);
     if (rc != MPI_SUCCESS)
     {
        printf ("Error starting MPI program. Terminating.\n");
        MPI_Abort(MPI_COMM_WORLD, rc);
     }
     cout<<"Done 2..!!"<<endl;
     MPI_Comm_size(MPI_COMM_WORLD,&numProc);
     MPI_Comm_rank(MPI_COMM_WORLD,&rank);
     MPI_Get_processor_name(hostname, &len);
    
    /********Data distribute ********/

    double t1;
    int numThreads,best;
    double rt_start = omp_get_wtime();

    /********* Reading the INPUT ******/
  DCSR G;                //return reading graph
    //distNodePartition(&G, opts.problemname, opts.bFileName,numProc, MPI_COMM_WORLD);
    distGraphReading(&G, opts.problemname, opts.bFileName,opts.b, opts.sstep,numProc, MPI_COMM_WORLD);
  if(!GraphReading(&G, opts.problemname))
  {
	  cout<<"Problem in reading graph "<<endl;
	  return -1;
  }
  
  
    //MPI_Barrier(MPI_COMM_WORLD);
    cout<<"Done 4..!!"<<endl;
   // if(rank==0)
   // {

    double rt_end = omp_get_wtime();
    cout<<"Graph (" << G.nVer << ", " << G.nEdge/2 <<
            ") Reading Done....!! took " << rt_end - rt_start <<endl;
    
    /*********** Memory Allocation *************/
    
  
    #pragma omp parallel
    numThreads=omp_get_num_threads();
    
  double* vtxWght=(double*)_mm_malloc(G.nVer*sizeof(double),64);   //do array to storing the wight of Ver with size nVer
  
  if(!VtxWghtReading(vtxWght,opts.wghtfilename))
  {
	   cout<<"Problem in reading vertex weights "<<endl;
	   return -1;
  }
   
  long long* mate=(long long*)_mm_malloc(G.nVer*sizeof(long long),64);  // mate store each Ver connect with which Ver
  for(long long i=0;i<G.nVer;i++)       
    {
        // في البداية كل الفيرتكس راح نبداها بحيث تكون مومتصله بعدها بيخزن باقي القيم
       mate[i] = cNullItm;
    }   

    
    int type=opts.algorithm;

  double t2,t3,t4,t5;
 switch(type)
    {
        case 1: {
				t2=omp_get_wtime();
                mvmHalfPar(&G,mate,vtxWght,opts.verbose);
				t3=omp_get_wtime()-t2;                     // record time
				cout<<"Matching done in :"<<t3<<" sec"<<endl;
				if(!verifyMatching(&G, mate, vtxWght))     //
					cout<<"The vertex weighted matching is not correct"<<endl;
				break;
				}

        case 2: {
				vtxWghtToEdgeWght(&G, vtxWght);
				t4=omp_get_wtime();
				suitorPar(&G,mate,opts.verbose);
				t5=omp_get_wtime()-t4;
				cout<<"Matching done in :"<<t5<<" sec"<<endl;
				if(!verifyMatchingEdgeWght(&G, mate))
					cout<<"The edge weighted matching is not correct"<<endl;
                break;
				}

    }
	 
    cout<<"Done 8..!!"<<endl;


    cout<<"Done 10..!!"<<endl;
    return 0;

}
