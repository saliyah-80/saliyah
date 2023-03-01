//
//  SRD.cpp
//
//
//  Created by Saliyah on 28/03/1444 AH.
//


#include "mvm.hpp"
#include <cstring>
#include <parallel/algorithm>
using namespace std;

#define MBUFLENGTH 1024*4       //size memoy chip
#define TBUFLENGTH 256*4    

long long plock;                // process lock
long long* pcount;
double** ptemp;                 
long long* endT;                //end thread
long long* startT;              //start thread
//long long* bT;                // buffer thread
//Edge* eT;                     //edge thread


vector<vector<double> >psend;    //preocess send
vector<vector<double> >precv;    //process recive


 /////////////////////////////////////////////////////////////////////////////////////
bool insertMessage(DCSR* G, long long sid, long long vid, float weight, long long type, long long tid)


// DCSR* G,graph
// sid,    sender id
// vid,    vertices id
// weight, wight
// type
// tid ,   thread id

{
    //assert(tid==0);
    assert(vid>=0 && vid<G->gnVer);
    long long count=pcount[tid];
    long long pid,id;
    if(count==TBUFLENGTH)
    {//type __sync_lock_test_and_set (type *ptr, type value, ...)
        while(__sync_lock_test_and_set(&plock,1)==1);

        for(long long i=0;i<TBUFLENGTH;i=i+4)
        {
            pid=G->gnVer-G->VerPerProc*(G->nProc-1);
            id=(long long)ptemp[tid][i+1];
            assert(pid>=0);

            if(id<pid)
                pid=0;
            else
                pid=((id-pid)/G->VerPerProc)+1;
            
            assert(pid>=0);
            if(pid>=G->nProc)
            {
                cout<<pid<<" "<<id<<" "<<G->VerPerProc<<" "<<G->gnVer<<" "
                    <<(long long)ptemp[tid][i+1]<<" "<<(i+1)%4<<endl;
            }
            assert(pid<G->nProc);
            psend[pid].push_back(ptemp[tid][i]);
            psend[pid].push_back(ptemp[tid][i+1]);
            psend[pid].push_back(ptemp[tid][i+2]);
            psend[pid].push_back(ptemp[tid][i+3]);
        }
        //  __sync_lock_release (type *ptr, ...)
        //  This builtin releases the lock acquired by __sync_lock_test_and_set. Normally this means writing the constant 0 to *ptr.
        __sync_lock_release(&plock);
        count=0;
    }
    assert(count>=0 && count+3<TBUFLENGTH);
    ptemp[tid][count]=sid;
    ptemp[tid][count+1]=vid;
    ptemp[tid][count+2]=weight;
    ptemp[tid][count+3]=type;
    count=count+4;
    pcount[tid]=count;

    return 1;
}
 /////////////////////////////////////////////////////////////////////////////////////
 
void mvmHalfPar(DCSR* G,long long *b, long long* mate,double* vtxWght, int *nlocks, Node* M, long long* start, long long* end, long long stepM, bool verbose,char* mark)
{

    long long n=G->nVer;
    long long m=G->nEdge;
    long long* ver=G->verPtr;
    Edge* verInd=G->verInd;
    long long* Q1=(long long*)_mm_malloc(n*sizeof(long long),64);
    long long* Q2=(long long*)_mm_malloc(n*sizeof(long long),64);
    long long* Qb1=(long long*)_mm_malloc(n*sizeof(long long),64);
    long long* Qb2=(long long*)_mm_malloc(n*sizeof(long long),64);
    
    
    long long tQsize;
    long long npart;
    long long pQsize=n/npart;   //number of the Ver that we split over the threads
    long long startQ,endQ;
    long long addToGlobal=G->startVer;
    long long state=0,donecount=0;
    long long terminate=0,msize,count,vid;
    int pid;
    
  ///////all about thread///////////
    long long numThreads;
    #pragma omp parallel
    numThreads=omp_get_num_threads();
    //cout<<"Processes: "<<G->nProc<<endl;
    //cout<<"Threads: "<<numThreads<<endl;
    
    long long* Qtemp;
    long long Qsize=n, Qindx=0, iter=1;
    double t1,t2,t3,t4;
    long long candC=0,rejC=0,kickC=0,msgC=0,conC=0;
    
    
    #ifdef LBUF
    long long** TQ=(long long**)_mm_malloc(numThreads*sizeof(long long*),64);
    long long* tindx=(long long*)_mm_malloc(numThreads*sizeof(long long),64);

    for(long long i=0;i<numThreads;i++)
    {
        tindx[i]=0;
        TQ[i]=(long long*)_mm_malloc(BSIZE*sizeof(long long),64);
    }
    #endif
 /////////////////////////////////////////////////////////////////////////////////////
    
    
    ///////////// memory allocation for threads////////////
    ///////////// Data structures for messages/////////////
    plock=0;                                                        // initialize process lock
    pcount=(long long*)_mm_malloc(numThreads*sizeof(long long),64); //pcount
    
    ptemp=(double**)_mm_malloc(numThreads*sizeof(double*),64);
    for(long long i=0;i<numThreads;i++)
    {
        ptemp[i]=(double*)_mm_malloc(TBUFLENGTH*sizeof(double),64);  //ptemp[i]thread buffer length
        pcount[i]=0;
    }
    
   /////////////////////////////////////////////////////////////////////////////////////
    
    //////////////send and receive//////////////
    vector<double> temp;         //assign data and get data from temp
    for(long long i=0;i<G->nProc;i++)
    {
        psend.push_back(temp);   //used to push elements into a vector from the back.
        precv.push_back(temp);
        
        
    }
    for(long long i=0;i<G->nProc;i++)
    {
        psend[i].clear();      //clear the queue
        precv[i].clear();
    }

    double** ssize=(double**)_mm_malloc(G->nProc*sizeof(double*),64);  //send size
    double** rsize=(double**)_mm_malloc(G->nProc*sizeof(double*),64);  //recive size
    for(long long i=0;i<G->nProc;i++)
    {
        ssize[i]=(double*)_mm_malloc(4*sizeof(double),64);
        rsize[i]=(double*)_mm_malloc(4*sizeof(double),64);
    }
  
    cout<<"Memory Allocation done..!!"<<endl;
 /////////////////////////////////////////////////////////////////////////////////////
 
    
    if(verbose)
        t1=omp_get_wtime();

    #ifdef DYN
        #pragma omp parallel for schedule(dynamic,CHUNK)
    #else
        #pragma omp parallel for schedule(static, CHUNK)
    #endif
    
    for(long long i=0;i<n;i++)
    {
        nlocks[i]=0;            //Ver locks
        start[i]=ver[i];        // Adj List start pointer
        //int custom(Edge* verInd,  long long start, long long end, long long step);
        end[i]=custom(verInd,ver[i],ver[i+1],stepM*b[i]);
        Q1[i]=i;
        Qb1[i]=mate[i];   //qbuffer=all vertices
        Qb2[i]=0;
    }
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      ///////////// check matching over all threads////////////
    #ifdef DYN
        #pragma omp parallel for schedule(dynamic,CHUNK)
    #else
        #pragma omp parallel for schedule(static, CHUNK)
    #endif
    for(long long i=0;i<G->gnVer;i++)
    {
        if(i>=G->startVer && i<=G->endVer)
        {
            M[i].curSize=0;                    // current size for node M=0
            M[i].maxSize=mate[i-G->startVer]; // maximum size for node M=mate
            M[i].minEntry.id=-1;
            M[i].minEntry.weight=0.0;
        }
        else
        {
            M[i].curSize=0;         // current size for node M=0
            M[i].maxSize=1;         // current size for node M=1
            M[i].minEntry.id=-1;
            M[i].minEntry.weight=0.0;
        }
      
    }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    if(verbose)
        t2=omp_get_wtime();

    
    endT=end;      //end thread
    startT=start;  //start thread
   // bT=mate;       //buffer thread
   // eT=verInd;     //edge thread
    
 while(Qsize>0 || donecount !=(G->nProc-1))
        {
            cout<<"Qsize: "<<Qsize<<endl;
           // if(nsort)
             // Explicitly force a call to parallel sort.
             // __gnu_parallel::sort(v.begin(), v.end());
                __gnu_parallel::sort(&Q1[0],&Q1[Qsize]);
                
            
            startQ=0;
            endQ=0;
            while(true)
            {
                startQ=endQ;
                endQ=startQ+pQsize;
                if(endQ>=Qsize)
                    endQ=Qsize;
                
                    #ifdef DYN
                        #pragma omp parallel for schedule(guided,CHUNK)
                    #else
                        #pragma omp parallel for schedule(static, CHUNK)
                    #endif
                    for(int i=0;i<m;i++)
                        mark[i]=0;
                
               
                cout<<"Initialization Done...!!"<<endl;
             
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                if(verbose)
                    t2=omp_get_wtime();
                
                cout << "Start Matching" << Qsize << endl;
                

     #pragma omp parallel for schedule(guided,CHUNK)
     for(long long i=startQ;i<endQ;i++)
                {
                    float min_w,heaviest,vtxWght;   //minmum wight, heaviest , vtxWght
                    long long  min_idx,Ver,current,gcurrent,next_vertex,sold,eold,skip;
                    long long  Ver_XIdx,Ver_X,lid,y,j,done,tid=0;
                    long long excess;
                    
                    tid=omp_get_thread_num();     // returns the thread number of the calling thread.0 to
                    //assert(tid==0);
                    
                    #ifdef LBUF
                    long long* LQ=TQ[tid];
                    long long* lindx=&tindx[tid];
                    #endif
                    
                    current=Q1[i];
                    gcurrent=current+addToGlobal;
                    assert(current>=0 && current<n);
                    Ver=Qb1[current];
                    Qb1[current]=0;
                    
                    
                 /*   /// Super stepping logic
                    if(Ver>stepM)
                    {
                        excess=Ver-stepM;
                         Ver=stepM;
                 */
                    if(end[current]!=-1)
                        {
                            #ifdef LBUF
                            //T __sync_fetch_and_add (T* __p, U __v, ...);
                            //__sync_fetch_and_add (type *ptr, type value, ...)
                            //This function atomically adds the value of __v to the variable that __p points to. The result is stored in the address that is specified by __p. A full memory barrier is created when this function is invoked.
                            
                            if(__sync_fetch_and_add(&Qb2[current],excess)==0) //if excess val that added to Qeue buffer 2 is equal 0
                            {
                                LQ[(*lindx)]=current;
                                (*lindx)++;
                                if((*lindx)==BSIZE)
                                {
                                    long long start=__sync_fetch_and_add(&Qindx,BSIZE);//start is add buffer size value to queue indx
                                    long long count=0;
                                    for(long long k=start;k<start+BSIZE;k++)
                                    {
                                        Q2[k]=LQ[count];
                                        count++;
                                    }
                                    
                                    (*lindx)=0;
                                }
                            }
                            #else
                            if(__sync_fetch_and_add(&Qb2[current],excess)==0)
                               Q2[__sync_fetch_and_add(&Qindx,1)]=current;
                            #endif
                        }

                 

                    while(Ver>0)
                    { // For Each vertex we want to find match
                        cout<<i<<" "<<Ver<<endl;
                        Ver=-1;
                        heaviest=0.0;
                        next_vertex=-1;
                        done=1;

                        sold=start[current];
                        eold=end[current];
                        j=sold;
                        assert(j>=0);
                        //assert(j<m);
                        while(j<end[current]) // Loop over neighbors of the current vertex
                        {
                            //cout<<"neighbor: "<<Ver<<" "<<j<<" "<<end[current]<<endl;
                            y = verInd[j].id;               // y is the neighbor of the current vertex
                            assert(y>=0 && y<G->gnVer);
                           vtxWght=verInd[j].weight;        // weight is w(current,y)
                            if(vtxWght<=0)
                            {
                                heaviest=-1.0;
                                break;
                            }
                            min_w=M[y].minEntry.weight;

                            if(vtxWght <=0 || min_w >vtxWght)
                            {
                                j++;
                                continue;
                            }
                         else
                         {
              // KICK Message
                                    lid=y-addToGlobal;
                                       while(__sync_lock_test_and_set(&nlocks[y],1)==1);
                                       insertMessage(G,y,min_idx,
                                                M[y].minEntry.weight,M[y].minEntry.id,tid);
                                       __sync_lock_release(&nlocks[lid]);
                                        
                                       if(verbose)
                                        __sync_fetch_and_add(&kickC,1);
                                  //  }
                                    
                                    break; //breaking from neighbor search
                                }
                          }
                         //   else
                            {
                                if(min_w > vtxWght)
                                    j++;
                                else
                                 
                                {
                                
                                    while(__sync_lock_test_and_set(&M[y].curSize,1)==1);
                                    M[y].minEntry.weight=vtxWght;
                                    M[y].minEntry.id=gcurrent;
                                    __sync_lock_release(&M[y].curSize);
                                    
                                    assert(y<G->gnVer);
                                    insertMessage(G,gcurrent,y,vtxWght,-2,tid);
                                    // ADD To The proposal queue
                                    heaviest=vtxWght;
                                    start[current]=j+1;
                                    
                                    if(verbose)
                                        __sync_fetch_and_add(&candC,1);
                                    
                                    break;
                                }
                          //  }

                        }   // while(j<end[current])

                        //cout<<"More neighbor"<<endl;
                        if(heaviest<=0)
                        {
                            start[current]=j;
                            if(end[current]<ver[current+1] && heaviest==0)
                            {
                            
                            end[current]=custom(verInd,start[current],ver[current+1],stepM*b[current]);
                                done=0;
                            }
                            else
                                end[current]=-1;
                        }
                      
                        if(end[current]==-1)  // if the vertex is alive
                            break;            // Do not forget to decrease Ver ...!!
                        else
                            if(next_vertex!=current && done==1)
                                Ver--;
                    } 
              //  }
/////////////////////////////////////////////////////////////////////////////////////////
                 j=ver[current];
                
                    for(;j<ver[current+1] && j!=-1;j++)
                    { // Loop over neighbors of the current vertex
                
                      if(mark[j]==1)
                            continue;
                
                        y = verInd[j].id;       // y is the neighbor of the current vertex
                        vtxWght=verInd[j].weight;        // weight is w(current,y)
                        min_w=M[y].minEntry.weight;

                        if((vtxWght<heaviest)|| (vtxWght == heaviest && y < Ver_X)||min_w>vtxWght)
                            continue;

                        if(min_w==vtxWght)
                        {
                            skip=0;
                    
                            while(__sync_lock_test_and_set(&nlocks[y],1)==1);
                    
                            min_w=M[y].minEntry.weight;
                            min_idx=M[y].minEntry.id;
                            if((min_w >vtxWght)||(vtxWght == min_w && current < min_idx))
                            skip=1;
                            __sync_lock_release(&nlocks[y]);
                       
                            if(skip==1)
                                continue;
                        }
                       
                        heaviest = vtxWght;  // Store the weight of the heaviest edge found so far
                        Ver_X = y;
                        Ver_XIdx=j;
                    
                    
                    } // loop over neighbors
                    
                    if (heaviest > 0) // True if there is a new Ver
                    {
                        if(__sync_lock_test_and_set(&nlocks[Ver_X],1)==0) //Locking Ver
                        {
                            min_w=M[Ver_X].minEntry.weight;
                            min_idx=M[Ver_X].minEntry.id;

                          
                               // the state so search Ver all over again
                             next_vertex=current;
                         
                              mark[Ver_XIdx]=1;
                            __sync_lock_release(&nlocks[Ver_X]); // Unlocking Ver
           
                        }
                        else
                            // Missed the lock, so someone else may or
                            // may not changed the state. Hence, RE DO the next_vertex search
                            next_vertex=current;
                  
                        if (next_vertex != -1 && next_vertex!=current)
                            // True if current vertex just kicked another vertex which is alive
                        {
                            #ifdef LBUF

                            if(__sync_fetch_and_add(&Qb2[next_vertex],1)==0)
                            {
                                LQ[(*lindx)]=next_vertex;
                                (*lindx)++;
                                if((*lindx)==BSIZE)
                                {
                                    int start=__sync_fetch_and_add(&Qindx,BSIZE);
                                    int count=0;
                                    for(int k=start;k<start+BSIZE;k++)
                                    {
                                        Q2[k]=LQ[count];
                                        count++;
                                    }
                                    (*lindx)=0;
                                }
                            }
                            #else
                            if(__sync_fetch_and_add(&Qb2[next_vertex],1)==0)
                                Q2[__sync_fetch_and_add(&Qindx,1)]=next_vertex;
                            #endif
                        }
                    
                    else
                    end[current]=-1;
                    //This means the neighbors list is exhausted,
                    //this vertex will never be considered again .!!
                }
            
          
            } // while
        } // loop over vertices
    }
//}

//////////////////////////////////////////////////
     //  int countm=0;
       for(long long i =0; i< n;i++)
        {
        
            int tid=omp_get_thread_num();
            nlocks[i]=0;            // Initialize locks
            M[i].curSize=0;         // current matching=0
            M[i].maxSize=mate[i];         // maximum matching=mate
            M[i].minEntry.id=-1;
            M[i].minEntry.weight=0.0;
           
           
        if(mate[i]!=cNullItm)
            continue;
            start[i]=ver[i];    // Adj List start pointer
            if(mate[i]>0){
               // end[i]=custom_sort(verInd,ver[i],ver[i+1],stepM*b[i],NULL,type);
                end[i]=custom(verInd,ver[i],ver[i+1],stepM*b[i]);}
            else
            {
                end[i]=-1;
                M[i].minEntry.id=G->nEdge+1;
                M[i].minEntry.weight=MAX_VAL;
            }
               
    //    }
        
        Q1[i]=i;
        Qb1[i]=mate[i];
        Qb2[i]=0;
      }
    }
   
//match all ver without condition
         
            while(j < ver[i+1])
            {
                long long j = ver[i];
                long long y = verInd[j].id;
                if(mate[y]==cNullItm)
                do {
#ifdef DYN
    #pragma omp parallel for schedule(guided,CHUNK)
#else
    #pragma omp parallel for schedule(static, CHUNK)
#endif
for(int i=0;i<m;i++)
    mark[i]=0;
                    mate[i]=y;
                    mate[y]=i;  // match all ver
                    break;
                }
              j++;
            
        }
       bool change = false;
        do{
            
        
                if(mate[y]==cNullItm)
                {
#ifdef DYN
    #pragma omp parallel for schedule(guided,CHUNK)
#else
    #pragma omp parallel for schedule(static, CHUNK)
#endif
for(int i=0;i<m;i++)
    mark[i]=0;
                    mate[i]=y;
                    mate[y]=i;  // match all ver
                    
                    
                     bool change = false;
            for(long long i =0; i< n;i++)
            {
                if(mate[i]!=cNullItm)
                continue;
                
                long long cand1=cNullItm;
                long long cand2=cNullItm;
                long long j = ver[i];
                    double lightest = vtxWght[i];
                
 ///////////////////////////////////////////////////////////////////////////////////////

// to find unmatchd Ver and matched edge with lightest wight z
#pragma omp parallel num_threads(numThreads)

                while(j < ver[i+1])
                {
                    long long y = verInd[j].id;
                    if(mate[y]==cNullItm )   //if y unmatched we match it
                    {
                        cand1=y;
                        cand2=cNullItm;
                        break;
                        
                    }
                    else
                    {
                        long long z = mate[y];   //if y matched we take wight of z in cand2
                        if(vtxWght[z] < lightest)
                        {
                            lightest = vtxWght[z];
                            cand1 = y;
                            cand2 = z;
                        }
                    }
                  j++;
                }
                if(cand1!=cNullItm && cand2==cNullItm)
                {
                    mate[i]=cand1;
                    mate[cand1]=i;
                    countm++;
                    change = true;
                }

                else if(cand1!=cNullItm && cand2 !=cNullItm)
                {
                    mate[i]=cand1;
                    mate[cand1]=i;
                    mate[cand2]=cNullItm;
                    change =true;
                }
                
            }
           }
        while(change){

    
cout<<"Matching Done....number of mathing!! " <<countm<<endl;}

////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////// from here
                /// Super stepping logic : is show to load-balance the communications among the compute nodes by making a trade-off between the freshness of the lastinformation and the volume of communication.

                if(end[current]!=-1)
                {
                    #ifdef LBUF
                    if(__sync_fetch_and_add(&Qb2[current],Ver)==0)
                    {
                        LQ[(*lindx)]=current;
                        (*lindx)++;
                        if((*lindx)==BSIZE)
                        {
                            long long start=__sync_fetch_and_add(&Qindx,BSIZE);
                            long long count=0;
                            for(long long k=start;k<start+BSIZE;k++)
                            {
                                Q2[k]=LQ[count];
                                count++;
                            }
                            
                            (*lindx)=0;
                        }
                    }
                    #else
                    if(__sync_fetch_and_add(&Qb2[current],excess)==0)
                       Q2[__sync_fetch_and_add(&Qindx,1)]=current;
                    #endif
                }
         //   }
            
#pragma omp parallel num_threads(numThreads)
//Add the rest messages to process message queues
        for(long long i=0;i<numThreads;i++)
        {
            count=pcount[i];
            long long tpid;
            for(long long j=0;j<count;j=j+4)
            {
                tpid=G->gnVer-G->VerPerProc*(G->nProc-1);
                vid=(long long)ptemp[i][j+1];

                if(vid<tpid)
                    pid=0;
                else
                    pid=((vid-tpid)/G->VerPerProc)+1;
                
                psend[pid].push_back(ptemp[i][j]);
                psend[pid].push_back(ptemp[i][j+1]);
                psend[pid].push_back(ptemp[i][j+2]);
                psend[pid].push_back(ptemp[i][j+3]);
            }
            pcount[i]=0;
        }
        
        //cout<<"(Process "<<G->rank<<"): "<<"Sending size messages"<<endl;
        for(int i=0;i<G->nProc;i++)
        {
            if(i==G->rank)
            {
                // Let's wait for the MPI_Isend to complete before progressing further and check the request has been set to MPI_REQUEST_NULL.
                ssreq[i]=MPI_REQUEST_NULL;
                continue;
            }
// send size messages

            if(psend[i].size()>0)
                msize=psend[i].size();
            else
                msize=4;
            //cout<<"Sending size: "<<msize<<endl;
            
            
            ssize[i][2]=msize;
            ssize[i][3]=-1;
           // int MPI_Isend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request)
            MPI_Isend(&ssize[i][0],4,MPI_DOUBLE,i,10,MPI_COMM_WORLD,&ssreq[i]);
            
            if(verbose)
            {
                msize=msize/4;
                __sync_fetch_and_add(&msgC,msize);
              //  __sync_fetch_and_add(&conC,1);
            }
        }

        //cout<<"(Process "<<G->rank<<"): "<<"Sending Data messages"<<endl;
        
        
// send data message
        for(int i=0;i<G->nProc;i++)
        {
            if(i==G->rank)
            {

                sdreq[i]=MPI_REQUEST_NULL;
                continue;
            }

//int MPI_Waitany(int count, MPI_Request array_of_requests[], int *indx, MPI_Status * status)
            MPI_Waitany(G->nProc,ssreq,&pid,status);
            
            

            if(psend[pid].size()>0)
            {
                msize=psend[pid].size();
                //cout<<"Sending DATA to: "<<pid<<endl;
                MPI_Isend(&psend[pid][0],msize,MPI_DOUBLE,pid,11,MPI_COMM_WORLD,&sdreq[pid]);
                
              //  if(verbose)
              //      __sync_fetch_and_add(&conC,1);

           
            else // Termination message
           
                //cout<<"Sending TERM to:"<<pid<<endl;
                terminate++;
                if(state==0)
                    ssize[pid][3]=-4;
                else
                    ssize[pid][3]=-5;
                MPI_Isend(&ssize[pid][0],4,MPI_DOUBLE,pid,11,MPI_COMM_WORLD,&sdreq[pid]);
            }
        }
 //Recieving size messages       
        //cout<<"(Process "<<G->rank<<"): "<<"Recieving size messages"<<endl;
        for(int i=0;i<G->nProc;i++)
        {
            if(i==G->rank)
            {
                rsreq[i]=MPI_REQUEST_NULL;
                continue;
            }
            MPI_Irecv(&rsize[i][0],4,MPI_DOUBLE,i,10,MPI_COMM_WORLD,&rsreq[i]);
        }
//Recieving data messages
        //cout<<"(Process "<<G->rank<<"): "<<"Recieving data messages"<<endl;
        for(int i=0;i<G->nProc;i++)
        {
            // Do waitany on size request and call data message for that process
            if(i==G->rank)
            {
                rdreq[i]=MPI_REQUEST_NULL;
                continue;
            }
            
            MPI_Waitany(G->nProc,rsreq,&pid,status);
            msize=rsize[pid][2];
            precv[pid].resize(msize);
            
            MPI_Irecv(&precv[pid][0],msize,MPI_DOUBLE,pid,11,MPI_COMM_WORLD,&rdreq[pid]);
        }
        
        /// We need to clear all the send buffer
        /// We can avoid this Waitall() by using another set of send buffers
        
        MPI_Waitall(G->nProc,sdreq,status);
        
        
        for(int i=0;i<G->nProc;i++)
            psend[i].clear();

        //cout<<"(Process "<<G->rank<<"): "<<"Processing messages.."<<endl;
        for(int i=0;i<G->nProc;i++)
        {
            // Do waitany on data request and process parallely
            if(i==G->rank)
                continue;
            MPI_Waitany(G->nProc,rdreq,&pid,status);
            //pid=i;
            
            #pragma omp parallel for
            for(long long k=0;k<(long long)rsize[pid][2];k=k+4)
            {
                
                long long gcurrent=precv[pid][k];
                long long current=gcurrent-addToGlobal;
                long long y=precv[pid][k+1];
                float weight=precv[pid][k+2];
                long long mtype=precv[pid][k+3];
                long long min_idx=-1;
                float heaviest=0.0,min_w;
                long long tid=omp_get_thread_num();
                #ifdef LBUF
                long long* LQ=TQ[tid];
                long long* lindx=&tindx[tid];
                #endif
                if(mtype<0)
                {
                    switch(mtype)
                    {
                        case -1: // SIZE
                          //  cout<<"Size message should not show up here..!!"<<endl;
                            break;
                        case -2: // PROPOSAL
                             
                            while(__sync_lock_test_and_set(&nlocks[lid],1)==1);
                                           
                            min_w=M[y].minEntry.weight;
                            min_idx=M[y].minEntry.id;

                            if((min_w > vtxWght) || (vtxWght == min_w && gcurrent < min_idx))
                            {
                               
                                insertMessage(G,y,gcurrent,min_w,min_idx,tid);
                                // ADD REJECT MESSAGE
                               
                                __sync_lock_release(&nlocks[lid]);
                                
                                if(verbose)
                                __sync_fetch_and_add(&rejC,1);
                            }
                            else
                            {
                                heaviest=weight;
                                            
                                __sync_lock_release(&nlocks[lid]);

                                if(min_idx==-1)
                                    break;
                                
                                /// Check if dislodged is local
                                if(min_idx>=G->startVer && min_idx<=G->endVer)
                                {
                                    //gid=min_idx;
                                    lid=min_idx-addToGlobal;
                                    
                                    if(end[lid]!=-1)
                                    {
                                        #ifdef LBUF
                                        if(__sync_fetch_and_add(&Qb2[lid],1)==0)
                                        {
                                            LQ[(*lindx)]=lid;
                                            (*lindx)++;
                                            if((*lindx)==BSIZE)
                                            {
                                                long long start=__sync_fetch_and_add(&Qindx,BSIZE);
                                                long long count=0;
                                                for(long long k=start;k<start+BSIZE;k++)
                                                {
                                                    Q2[k]=LQ[count];
                                                    count++;
                                                }
                                                
                                                (*lindx)=0;
                                            }
                                        }
                                        #else
                                        if(__sync_fetch_and_add(&Qb2[lid],1)==0)
                                           Q2[__sync_fetch_and_add(&Qindx,1)]=lid;
                                        #endif
                                    }
                                }
                                else
                                {
                                   
                                   lid=y-addToGlobal;
                                   while(__sync_lock_test_and_set(&nlocks[lid],1)==1);
                                   insertMessage(G,y,min_idx,
                                            M[y].minEntry.weight,M[y].minEntry.id,tid);
                                   __sync_lock_release(&nlocks[lid]);

                                    if(verbose)
                                        __sync_fetch_and_add(&kickC,1);
                                }
                            }
                            break;

                        case -3: // UPDATE
                            break;
                        case -4: // TERMINATE
                            __sync_fetch_and_add(&terminate,1);
                            break;
                        case -5: // DONE
                            __sync_fetch_and_add(&terminate,1);
                            __sync_fetch_and_add(&donecount,1);
                            break;

                    }
                }
                else // REJECT
                {
                    
                    while(__sync_lock_test_and_set(&M[gcurrent].curSize,1)==1);
                    M[gcurrent].minEntry.weight=weight;
                    M[gcurrent].minEntry.id=mtype;
                    __sync_lock_release(&M[gcurrent].curSize);

             
                    if(end[lid]!=-1)
                    {
                        #ifdef LBUF
                        if(__sync_fetch_and_add(&Qb2[lid],1)==0)
                        {
                            LQ[(*lindx)]=lid;
                            (*lindx)++;
                            if((*lindx)==BSIZE)
                            {
                                long long start=__sync_fetch_and_add(&Qindx,BSIZE);
                                long long count=0;
                                for(long long k=start;k<start+BSIZE;k++)
                                {
                                    Q2[k]=LQ[count];
                                    count++;
                                }
                                
                                (*lindx)=0;
                            }
                        }
                        #else
                        if(__sync_fetch_and_add(&Qb2[lid],1)==0)
                           Q2[__sync_fetch_and_add(&Qindx,1)]=lid;
                        #endif
                    }
                }
            } /// END FOR current process
        }/// END FOR all process
        #ifdef LBUF
        #pragma omp parallel for
        for(long long i=0;i<numThreads;i++)
        {
            long long* LQ=TQ[i];
            long long* lindx=&tindx[i];
            if((*lindx)>0)
            {
                long long start=__sync_fetch_and_add(&Qindx,(*lindx));
                long long count=0;
                for(long long k=start;k<start+(*lindx);k++)
                {
                    Q2[k]=TQ[i][count];
                    count++;
                }
                (*lindx)=0;
            }
       // }
    
    
        if(endQ==Qsize)
            break;
    }//// END WHILE NODE STEPPING
    
    #endif
    Qtemp=Q1;
    Q1=Q2;
    Q2=Qtemp;
    Qtemp=Qb1;
    Qb1=Qb2;
    Qb2=Qtemp;
    Qsize=Qindx;
    Qindx=0;


endT=end;
startT=start;
//bT=mate;
//eT=verInd;

    
    /////////// Termination criteria... a bit critical
    if(terminate==(2*G->nProc-2))
    {
        if(state==0)
        {
            //cout<<"Terminate: "<<G->rank<<" "<<terminate<<endl;
            state=1;
            terminate=0;
            donecount=0;
        
        elseif
        
            if(donecount==(G->nProc-1))
                break;
                
            else
            
                terminate=0;
                donecount=0;
            }
       

    else
   
        state=0;
        donecount=0;
        terminate=0;
    }


MPI_Barrier(MPI_COMM_WORLD);
MPI_Waitall(G->nProc,sdreq,status);
//}// END WHILE all done

if(G->rank==0)
{
cout<<"Matching Done....!!"<<endl;
}
if(verbose)
{
double* lcounts=(double*)_mm_malloc(5*sizeof(double),64);
    
lcounts[0]=msgC;
lcounts[1]=candC;
lcounts[2]=rejC;
lcounts[3]=kickC;
//lcounts[4]=conC;

MPI_Reduce(&lcounts[0],&gcounts[0],5,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

long long giga=1024*1024*1024;
if(G->rank==0)
{
    long long data=(long long)((gcounts[0]*32)/giga);
    double bandwidth=data/(23.7*1024);

    cout<<"# of communications: "<<(long long)gcounts[4]<<endl;
    cout<<"# of data xfer: "<<data<<" GB in "<<bandwidth<<" sec"<<endl;
    cout<<"# of messages: "<<(long long)gcounts[0]<<endl;
    cout<<"# of candidate: "<<(long long)gcounts[1]<<endl;
    cout<<"# of rejections: "<<(long long)gcounts[2]<<endl;
    cout<<"# of annullements: "<<(long long)gcounts[3]<<endl;
 }
        _mm_free(lcounts);
   // }
    
    _mm_free(Q1);
    _mm_free(Q2);
    _mm_free(Qb1);
    _mm_free(Qb2);
    _mm_free(tindx);
    _mm_free(pcount);
    _mm_free(ssreq);
    _mm_free(rsreq);
    _mm_free(sdreq);
    _mm_free(rdreq);
    _mm_free(status);

    for(long long i=0;i<G->nProc;i++)
    {
        _mm_free(ssize[i]);
        _mm_free(rsize[i]);
    }
    for(long long i=0;i<numThreads;i++)
    {
        _mm_free(TQ[i]);
        _mm_free(ptemp[i]);
    }
    
    _mm_free(ptemp);
    _mm_free(TQ);

}


