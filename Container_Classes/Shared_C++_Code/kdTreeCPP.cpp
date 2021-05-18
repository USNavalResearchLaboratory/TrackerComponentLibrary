/*
 * December 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
 *(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#include <limits>
#include <algorithm>
#include <cstring>
#include "mathFuncs.hpp"
#include "kdTreeCPP.hpp"

#include "mex.h"


using namespace std;

//Prototypes for functions not declared in external headers.
double dist(const double *a, const double *b,const size_t numEl);
bool rectsIntersect(const double *rectMin1,const double *rectMax1,const double *rectMin2,const double *rectMax2,const size_t numEl);
bool rectContained(const double *rectMin1,const double *rectMax1,const double *rectMin2,const double *rectMax2,const size_t numEl);
bool inHyperrect(const double *P,const double *rectMin,const double *rectMax,const size_t numEl);
bool boundsIntersectBall(const double *point, const double rSquared,const double *rectMin,const double *rectMax, const size_t numEl);

/*This structure is used with the sort function to sort one array according
 *to the values in another.*/
struct CompVal {
    const double *compList;
    
    CompVal(const double *compList1): compList(compList1)
    {}

    bool operator()(const size_t Lhs, const size_t Rhs)const
    {
        return compList[Lhs] < compList[Rhs];
    }
};


double dist(const double *a, const double *b,const size_t numEl) {
//DIST Compute the squared Euclidean distance between two vectors.
    double temp,distVal=0;
    size_t i;
    
    for(i=0;i<numEl;i++) {
        temp=a[i]-b[i];
        distVal+=temp*temp;
    }
    
    return distVal;
}


bool rectsIntersect(const double *rectMin1,const double *rectMax1,const double *rectMin2,const double *rectMax2,const size_t numEl) {
//RECTSINTERSECT  Determine whether two hyperrectangles intersect.
    size_t i;
    
    for(i=0;i<numEl;i++) {
        if(rectMax1[i]<rectMin2[i]||rectMax2[i]<rectMin1[i])
            return false;
    }
    
    return true;
}

bool rectContained(const double *rectMin1,const double *rectMax1,const double *rectMin2,const double *rectMax2,const size_t numEl) {
//RECTCONTAINED Determine whether a hyperrectangle is completely contained
//              in another hyperrectangle.
    size_t i;

    for (i=0;i<numEl;i++) {
        if(rectMin1[i]<rectMin2[i]|| rectMax1[i]>rectMax2[i])
            return false;
    }
    
    return true;
}

bool inHyperrect(const double *P,const double *rectMin,const double *rectMax,const size_t numEl) {
//INHYPERRECT Determine whether a point is within a given hyperectangular
//            region.
    size_t i;
    
    for(i=0;i<numEl;i++) {
        if(P[i]<rectMin[i]||P[i]>rectMax[i])
            return false;
    }
            
    return true;
}

bool boundsIntersectBall(const double *point, const double rSquared,const double *rectMin,const double *rectMax, const size_t numEl) {
//BOUNDSINTERSECTBALL Determines whether a sphere of a given squared radius
//                    centered at the given point intersects a
//                    hyperrectangular region.
    size_t i;
    double cumDist=0;

    for(i=0;i<numEl;i++) {
        double dist1, dist2;

        dist1=point[i]-rectMin[i];
        dist1=dist1*dist1;

        dist2=point[i]-rectMax[i];
        dist2=dist2*dist2;

        if(dist1<rSquared && dist2 < rSquared)
            continue;

        cumDist+=min(dist1,dist2);

        if(cumDist>rSquared)
            return false;
    }

    return true;
}

kdTreeCPP::kdTreeCPP() {
    buffer=NULL;
    N=0;
    k=0;
}

kdTreeCPP::kdTreeCPP(const size_t kDes, const size_t NDes) {
    char *basePtr;
    N=NDes;
    k=kDes;
    
/*To minimize the number of calls to memory allocation and deallocation
 * routines, a big chunk of memory is allocated at once and pointers
 * to parts of it for the different variables are saved.*/
    buffer=new char[sizeof(ptrdiff_t)*2*N+sizeof(size_t)*3*N+sizeof(double)*k*N*3];
    basePtr=buffer;

    LOSON=reinterpret_cast<ptrdiff_t*>(basePtr);
    basePtr+=sizeof(ptrdiff_t)*N;
    HISON=reinterpret_cast<ptrdiff_t*>(basePtr);
    basePtr+=sizeof(ptrdiff_t)*N;
    DATAIDX=reinterpret_cast<size_t*>(basePtr);
    basePtr+=sizeof(size_t)*N;
    DISC=reinterpret_cast<size_t*>(basePtr);
    basePtr+=sizeof(size_t)*N;
    subtreeSizes=reinterpret_cast<size_t*>(basePtr);
    basePtr+=sizeof(size_t)*N;
    BMin=reinterpret_cast<double*>(basePtr);
    basePtr+=sizeof(double)*k*N;
    BMax=reinterpret_cast<double*>(basePtr);
    basePtr+=sizeof(double)*k*N;
    data=reinterpret_cast<double*>(basePtr);
}

void kdTreeCPP::buildTreeFromBatch(const double *dataBatch){
/* This builds a balanced kd tree from a batch of data. Much of the
 * complexity comes from having to repeatedly sort the data to
 * determine the median element to use.
 * 
 * This sets up a variety of data structures to avoid having to
 * call any allocation routines during the recursion to add all of the
 * data.
 * 
 * Unlike the Matlab implementation, a subarray is not explicitly sorted
 * during each recursion. Rather, an array of indices, idx, is sorted
 * according to the values in the current dimension  of the data. To be
 * able to scan the data across all points for a particular dimension of
 * the data, the temporary transpose matrix dataTrans is used.
 *
 */

    /*To hold a transposed version of the data batch as well as the other
     *buffers needed while creating the tree.*/
    char *buffLoc;
    char *tempPtr;
    double *dataTrans;
    //This will also be recast and used to temporarily store indices.
    double *tempSortRow;
    size_t *idx;
    size_t i, j;

    //Allocate space for the array of indices and the transposed data
    //batch. The k*N is for the transposed data matrix. The N is for one
    //rows for sorting data and the N size_t is for an array of indices.
    buffLoc = new char[sizeof(double)*k*N+sizeof(double)*N+sizeof(size_t)*N];

    tempPtr=buffLoc;
    //Store a transposed temporary copy of the dataBatch.
    dataTrans= reinterpret_cast<double*>(tempPtr);
    for(i=0;i<k;i++){
        for(j=0;j<N;j++){
            dataTrans[j+i*N]=dataBatch[i+j*k];
        }
    }

    //Allocate space for a temporary row for sorting.
    tempPtr+=sizeof(double)*k*N;
    tempSortRow=reinterpret_cast<double*>(tempPtr);

    //Allocate space for the indices.
    tempPtr+=sizeof(double)*N;
    idx=reinterpret_cast<size_t*>(tempPtr);
    //Fill in the indices from 0 to N.
    for(i=0;i<N;i++){
        idx[i]=i;
    }
     
    //Copy the batch of data into the  memory for the class.
    memcpy(data,dataBatch,k*N*sizeof(double));

    //Initialize the BMin and BMax arrays.
    for(i=0;i<k*N;i++) {
        BMin[i]=numeric_limits<double>::infinity();
        BMax[i]=-numeric_limits<double>::infinity();
    }
    
    treeGrow(dataTrans,idx,tempSortRow,0,0,N);

    delete[] buffLoc;
}

size_t kdTreeCPP::treeGrow(double *dataTrans, size_t *idx, double *tempSortRow,const size_t level,const size_t curNode,const size_t NSubTree) {
    size_t nextFreeNode, nextNode;
    double *BMinBase=BMin+curNode*k;
    double *BMaxBase=BMax+curNode*k; 
    size_t curRow,i, midIdx, nextLevel;

    //Record the number of children of this node, this code included.
    subtreeSizes[curNode]=NSubTree;

    //Record the discriminating index at this level.
    DISC[curNode]=level;
    nextFreeNode=curNode+1;
    
    //First, check whether this is a leaf node. If so, then it has no
    //children and the index of the discriminator is just idx.
    if(NSubTree==1) {
        LOSON[curNode]=-1;
        HISON[curNode]=-1;
        DATAIDX[curNode]=*idx;
        
        //Copy the data vector of the current point into the min and
        //max arrays
        memcpy(BMinBase,data+k*(*idx),k*sizeof(double));
        memcpy(BMaxBase,data+k*(*idx),k*sizeof(double));
        return nextFreeNode;
    }

    //Sort the indicies according to the entries in the current
    //level of sortBatch. That is, according to the level-th dimension
    //of the data.    
    sort(idx,idx+NSubTree,CompVal(dataTrans+level*N));

    //Now, use the indicies to make a copy of the sorted data elements that
    //can be searched.
    for(i=0;i<NSubTree;i++) {
        tempSortRow[i]=dataTrans[idx[i]+level*N];
    }
    
    //Next, find the first occurence of the median element.
    midIdx=findFirstMaxCPP(tempSortRow,(NSubTree+1)/2);
    
    //Save the median (split) element.
    DATAIDX[curNode]=idx[midIdx];

    //The next level's index.
    nextLevel=(level+1)%k;        

    //If duplicate values are present, there might be nothing before the
    //current node and LOSON will be empty.
    if(midIdx==0){
        LOSON[curNode]=-1;//There is no splitting.
    } else {
        LOSON[curNode]=static_cast<ptrdiff_t>(nextFreeNode);
        nextNode=this->treeGrow(dataTrans,idx,tempSortRow,nextLevel,nextFreeNode,midIdx);

        //Record the minimum and maximum values.
        memcpy(BMinBase,BMin+nextFreeNode*k,k*sizeof(double));
        memcpy(BMaxBase,BMax+nextFreeNode*k,k*sizeof(double));
        nextFreeNode=nextNode;
    }

    //The HISON node will never be empty, except at a leaf node, since
    //the median is always taken using the floor function.
    HISON[curNode]=static_cast<ptrdiff_t>(nextFreeNode);
    nextNode=this->treeGrow(dataTrans,idx+midIdx+1,tempSortRow,nextLevel,nextFreeNode,NSubTree-midIdx-1);
    
    //Record the minimum and maximum values due to the child node and
    //due to the contribution of the current node.
    {
        const double *BMinNextBase=BMin+nextFreeNode*k;
        const double *BMaxNextBase=BMax+nextFreeNode*k;
        const double *dataBase=data+k*idx[midIdx];

        for(curRow=0;curRow<k;curRow++) {
            BMinBase[curRow]=min(BMinBase[curRow],BMinNextBase[curRow]);
            BMinBase[curRow]=min(BMinBase[curRow],dataBase[curRow]);
            BMaxBase[curRow]=max(BMaxBase[curRow],BMaxNextBase[curRow]);
            BMaxBase[curRow]=max(BMaxBase[curRow],dataBase[curRow]);
        }
    }

    nextFreeNode=nextNode;
    return nextFreeNode;
}

size_t *kdTreeCPP::rangeCount(const double *rectMin,const double *rectMax,const size_t numRanges) const {
    size_t *rangeCounts=new size_t[numRanges];
    size_t i;

    for(i=0;i<numRanges;i++) {
        rangeCounts[i]=this->rangeCountRecur(0,rectMin+i*k,rectMax+i*k);
    }
    
    return rangeCounts;
}

void kdTreeCPP::rangeQuery(ClusterSetCPP<size_t> &rangeClust,const double *rectMin,const double *rectMax,const  size_t numRanges) const {
    size_t *clusterSizes;
    size_t i;
    
    //Calculate the precise number of values to return by traversing the
    //tree once for each rectangle
    clusterSizes=this->rangeCount(rectMin,rectMax,numRanges);
    rangeClust.initWithClusterSizes(clusterSizes,numRanges);
    delete[] clusterSizes;
    
    for(i=0;i<numRanges;i++) {
        size_t numFound=0;
        
        this->rangeQueryRecur(0,rectMin+i*k,rectMax+i*k,rangeClust[i],numFound,clusterSizes[i]);
    }
}

void kdTreeCPP::rangeQueryRecur(const size_t curNode, const double *rectMin, const double *rectMax, size_t *idxRange, size_t &numFound, const size_t numInRange) const {
//The numFound parameter keeps track of how many solutions have been found for a particular range query and the search is broken off early if numFound=numInRange.
    
    double *P;
    ptrdiff_t childNode;

    //If curNode and all of its children are in the box, then return the
    //tree and all of its children.
    if(rectContained(BMin+curNode*k,BMax+curNode*k,rectMin,rectMax,k)) {
        getSubtreeIdx(curNode, idxRange,numFound);
        return;
    }

    //Otherwise, return the point, if it is in the range, and all of the
    //points from subtrees that overlap.
    P=data+k*DATAIDX[curNode];
    if(inHyperrect(P,rectMin,rectMax,k)) {
        idxRange[numFound]=DATAIDX[curNode];
        numFound++;
        if(numFound==numInRange) {
            return;
        }
    }
    
    childNode=HISON[curNode];
    if(childNode!=-1) {
        //If points might be in the child nodes.
        if(rectsIntersect(BMin+k*static_cast<size_t>(childNode),BMax+k*static_cast<size_t>(childNode),rectMin,rectMax,k)) {
            this->rangeQueryRecur(static_cast<size_t>(childNode), rectMin, rectMax,idxRange, numFound,numInRange);
            if(numFound==numInRange) {
                return;
            }
        }
    }
    
    childNode=LOSON[curNode];
    if(childNode!=-1) {
        //If points might be in the child nodes.
        if(rectsIntersect(BMin+k*static_cast<size_t>(childNode),BMax+k*static_cast<size_t>(childNode),rectMin,rectMax,k)) {
            this->rangeQueryRecur(static_cast<size_t>(childNode), rectMin, rectMax,idxRange, numFound,numInRange);
            if(numFound==numInRange) {
                return;
            }
        }
    }
}

size_t kdTreeCPP::rangeCountRecur(const size_t curNode,const double *rectMin,const double *rectMax) const {
    size_t numInRange=0;
    double *P;
    ptrdiff_t childNode;
    
    //If curNode and all of its children are in the box, then return the
    //number of items in the whole subtree.    
    if(rectContained(BMin+curNode*k,BMax+curNode*k,rectMin,rectMax,k)) {
        return subtreeSizes[curNode];
    }
    
    //Otherwise, increment the solution by one, if the point is in the
    //range, and all of the points from subtrees that overlap.
    P=data+k*DATAIDX[curNode];
    if(inHyperrect(P,rectMin,rectMax,k)) {
        numInRange++;
    }
    
    childNode=HISON[curNode];
    if(childNode!=-1) {
        //If points might be in the child nodes.
        if(rectsIntersect(BMin+k*static_cast<size_t>(childNode),BMax+k*static_cast<size_t>(childNode),rectMin,rectMax,k)) {
            numInRange+=this->rangeCountRecur(static_cast<size_t>(childNode),rectMin,rectMax);
        }
    }

    childNode=LOSON[curNode];
    if(childNode!=-1) {
        //If points might be in the child nodes.
        if(rectsIntersect(BMin+k*static_cast<size_t>(childNode),BMax+k*static_cast<size_t>(childNode),rectMin,rectMax,k)) {
            numInRange+=this->rangeCountRecur(static_cast<size_t>(childNode),rectMin,rectMax);
        }
    }
    
    return numInRange;
}


void kdTreeCPP::findmBestNN(size_t *idxRange,double *distSquared,const double *point,const size_t numPoints, const size_t m) const {
    priority_queue<pair<double,size_t> > mBestQueue;
    size_t i, curFound;

    for(i=0;i<numPoints;i++) {
        const size_t offset=m*i;
        
        //Start the recursion to fill the queue with the k-best values.
        this->mBestRecur(0,mBestQueue,point+i*k,m);
                
        curFound=m;
        do {
            pair<double,size_t>topPair=mBestQueue.top();
            
            mBestQueue.pop();
            curFound--;
                        
            *(distSquared+curFound+offset)=topPair.first;
            *(idxRange+curFound+offset)=DATAIDX[topPair.second];
        } while(curFound>0);
        //When it gets here, the queue should be empty and can thus be
        //directly reused the next time around the loop.
    }
}

void kdTreeCPP::mBestRecur(const size_t curNodeIdx, priority_queue<pair<double,size_t> > &mBestQueue, const double *point,const size_t m) const {
//First, go down the path on the nearest side of the splitting
//dimension from this point.
    double cost, *splitPoint;
    size_t splitDim;
    ptrdiff_t farNode;
    
    splitPoint=data+k*DATAIDX[curNodeIdx];
    splitDim=DISC[curNodeIdx];
        
    if(point[splitDim]<splitPoint[splitDim]) {
        ptrdiff_t lIdx=LOSON[curNodeIdx];
        if(lIdx!=-1) {
            this->mBestRecur(static_cast<size_t>(lIdx),mBestQueue,point,m);
        }

        farNode=HISON[curNodeIdx];
    } else {
        ptrdiff_t hIdx=HISON[curNodeIdx];
        if(hIdx!=-1) {
            this->mBestRecur(static_cast<size_t>(hIdx),mBestQueue,point,m);
        }

        farNode=LOSON[curNodeIdx];
    }
     
    //Next, visit this node.
    cost=dist(point,splitPoint,k);
    if(mBestQueue.empty()==true) {
        //Since no nodes have been found to this point, the best node is
        //the current one found thus far. This path will only be visited
        //at a leaf of the tree, so we can just return.
        pair<double,size_t> newPair(cost,curNodeIdx);
        mBestQueue.push(newPair);
        return;
    } else {
        //The point is only added to the queue if it is lower than the
        //maximum cost point found thus far or if there are fewer than k
        //points already in the queue.
        pair<double,size_t> curPair;
        double maxDist;
        
        curPair=mBestQueue.top();
        maxDist=curPair.first;
        
        if(maxDist>cost || mBestQueue.size() <m) {
            pair<double,size_t> newPair(cost,curNodeIdx);
            if(mBestQueue.size()==m) {
                mBestQueue.pop();
            }
            
            mBestQueue.push(newPair);
            maxDist=cost;
        }
        
        //Now, see if it is necessary to visit the other branch of the
        //tree. That is only the case if the bounding box intersects
        //with a ball centered at the point to find whose squared radius
        //is equal to maxDist or if there are fewer than k things in the
        //queue.
        if(farNode!=-1) {
            if(mBestQueue.size()<m||boundsIntersectBall(point,maxDist,BMin+k*static_cast<size_t>(farNode),BMax+k*static_cast<size_t>(farNode),k)) {
                this->mBestRecur(static_cast<size_t>(farNode),mBestQueue,point,m);
            }
        }
    }
}

void kdTreeCPP::getSubtreeIdx(const size_t nodeIdx, size_t *idxRange, size_t &numFound) const {
    //Add the current node.
    idxRange[numFound]=DATAIDX[nodeIdx];
    numFound++;

    //Add the child nodes.
    if(HISON[nodeIdx]!=-1) {//If there are children to the right
        this->getSubtreeIdx(static_cast<size_t>(HISON[nodeIdx]), idxRange, numFound);
    }

    if(LOSON[nodeIdx]!=-1) {//If there are children to the left.
        this->getSubtreeIdx(static_cast<size_t>(LOSON[nodeIdx]), idxRange, numFound);
    }
}

kdTreeCPP::~kdTreeCPP() {
    if(buffer !=NULL) {
        delete[] buffer;
    }
}

/*LICENSE:
%
%The source code is in the public domain and not licensed or under
%copyright. The information and software may be used freely by the public.
%As required by 17 U.S.C. 403, third parties producing copyrighted works
%consisting predominantly of the material produced by U.S. government
%agencies must provide notice with such work(s) identifying the U.S.
%Government material incorporated and stating that such material is not
%subject to copyright protection.
%
%Derived works shall not identify themselves in a manner that implies an
%endorsement by or an affiliation with the Naval Research Laboratory.
%
%RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
%SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY THE NAVAL
%RESEARCH LABORATORY FOR ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS
%OF RECIPIENT IN THE USE OF THE SOFTWARE.*/
