/*
 *December 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
 *(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#include "metricTreeCPP.hpp"
#include <limits>
//Needed for sqrt
#include <cmath>
#include <vector>
#include "mathFuncs.hpp"
//For memcpy
#include <cstring>

using namespace std;

//Prototypes for functions that are not in a header.
double distEuclid(const double *a, const double *b,const size_t numEl);

/*This structure is used with the sort function to sort one array according
 *to the values in another*/
struct CompVal {
    const double *compList;
    
    CompVal(const double *compList1): compList(compList1)
    {}

    bool operator()(const size_t Lhs, const size_t Rhs)const
    {
        return compList[Lhs] < compList[Rhs];
    }
};

double distEuclid(const double *a, const double *b,const size_t numEl) {
//DISTEUCLID Compute the Euclidean distance between two vectors.
    
    double temp, distVal=0;
    size_t i;
    
    for(i=0;i<numEl;i++) {
        temp=a[i]-b[i];
        distVal+=temp*temp;
    }
    
    return sqrt(distVal);
}

metricTreeCPP::metricTreeCPP() {
    buffer=NULL;
    N=0;
    k=0;
}

metricTreeCPP::metricTreeCPP(const size_t kDes, const size_t NDes) {
    char *basePtr;
    N=NDes;
    k=kDes;
    
/*To minimize the number of calls to memory allocation and deallocation
 * routines, a big chunk of memory is allocated at once and pointers
 * to parts of it for the different variables are saved.*/
    buffer=new char[sizeof(size_t)*N+sizeof(ptrdiff_t)*2*N+sizeof(double)*2*N+sizeof(double)*k*N];
    basePtr=buffer;

    DATAIDX=reinterpret_cast<size_t*>(basePtr);
    basePtr+=sizeof(size_t)*N;
    innerChild=reinterpret_cast<ptrdiff_t*>(basePtr);
    basePtr+=sizeof(ptrdiff_t)*N;
    outerChild=reinterpret_cast<ptrdiff_t*>(basePtr);
    basePtr+=sizeof(ptrdiff_t)*N;
    innerRadii=reinterpret_cast<double*>(basePtr);
    basePtr+=sizeof(double)*N;
    outerRadii=reinterpret_cast<double*>(basePtr);
    basePtr+=sizeof(double)*N;
    data=reinterpret_cast<double*>(basePtr);
}


void metricTreeCPP::buildTreeFromBatch(const double *dataBatch) {
    size_t *idx,i;
    char *buffLoc,*curPtr;
    double *adjMat;
    size_t *tempSortIdx;
    double *tempSortRow;
    
    //Copy the batch of data into the  memory for the class.
    memcpy(data,dataBatch,k*N*sizeof(double));
    
    //Allocate memory for index vectors and temporary buffers that will be
    //used for sorting during the recursion.
    buffLoc= new char[sizeof(double)*N+sizeof(size_t)*N+sizeof(size_t)*N+max(sizeof(double),sizeof(size_t))*N];
    curPtr=buffLoc;
    
    adjMat=reinterpret_cast<double*>(curPtr);
    curPtr+=sizeof(double)*N;
    idx=reinterpret_cast<size_t*>(curPtr);
    curPtr+=sizeof(size_t)*N;
    tempSortIdx=reinterpret_cast<size_t*>(curPtr);
    curPtr+=sizeof(size_t)*N;
    tempSortRow=reinterpret_cast<double*>(curPtr);
    
//Initialize the indices.
    for(i=0;i<N;i++) {
        idx[i]=i;
    }
    
    this->treeGrow(adjMat,0,idx,tempSortRow,tempSortIdx,0, N);

    delete[] buffLoc;    
}

size_t metricTreeCPP::treeGrow(double *adjMat,const size_t colOffset, size_t *idx, double *tempSortRow,size_t *tempSortIdx,const size_t curNode, size_t NSubTree) {
//TREEGROW A recursion function for creating a new metric tree.
    size_t midIdx, nextFreeNode, i;
    
    //The first index in idx is added to the tree.
    DATAIDX[curNode]=idx[colOffset];

    //If this is a leaf node.
    if(NSubTree==1) {
        innerChild[curNode]=-1;
        outerChild[curNode]=-1;
        innerRadii[curNode]=-numeric_limits<double>::infinity();
        outerRadii[curNode]=numeric_limits<double>::infinity();
        return curNode+1;
    }

    //Fill the NCurX1 subvector of adjMat with all pairwise
    //distances from the current point to the other points in idx.
    {size_t curPoint;
        NSubTree--;
        for(curPoint=0;curPoint<NSubTree;curPoint++) {
            adjMat[curPoint]=distEuclid(data+k*idx[colOffset],data+k*idx[colOffset+curPoint+1],k);
        }
    }
    
    //Next, sort the points according to their distances from the current
    //point.
    for(i=0;i<NSubTree;i++){
        tempSortIdx[i]=i;
    }
    sort(tempSortIdx,tempSortIdx+NSubTree,CompVal(adjMat));
    
    //The ordering of tempSortIdx corresponds to how the indicies from
    //idx[colOffset+1] to idx[colOffset+NSubTree] should be rearranged.
    memcpy(reinterpret_cast<size_t*>(tempSortRow),idx+colOffset+1,NSubTree*sizeof(size_t));
    for(i=0;i<NSubTree;i++) {
        idx[colOffset+1+i]=(reinterpret_cast<size_t*>(tempSortRow))[tempSortIdx[i]];
    }
    
    //Now, rearrange the distances in adjMat according to the same sort
    //order so that we can find the median distance.
    memcpy(tempSortRow,adjMat,NSubTree*sizeof(double));
    for(i=0;i<NSubTree;i++) {
        adjMat[i]=tempSortRow[tempSortIdx[i]];
    }

    //Find the first occurance of the median distance. The manipulation of
    //midIdx before being put into the function is the same as saying
    //ceil((NSubTree+1)/2)) if one were using floating point arithmetic.
    midIdx=(NSubTree+1);
    midIdx=findFirstMaxCPP(adjMat,midIdx/2+midIdx%2);
    
    //The inner and outer radii of how things are split in terms of
    //distances from the given point.
    outerRadii[curNode]=adjMat[midIdx];
    if(midIdx>0) {
        innerRadii[curNode]=adjMat[midIdx-1];
    } else {
        innerRadii[curNode]=-numeric_limits<double>::infinity();
    }

    //Continue the recursion.
    outerChild[curNode]=static_cast<ptrdiff_t>(curNode+1);
    nextFreeNode=this->treeGrow(adjMat,colOffset+1+midIdx,idx,tempSortRow,tempSortIdx,curNode+1,NSubTree-midIdx);

    //If a full partitioning of the nodes took place
    if(midIdx>0) {
        innerChild[curNode]=static_cast<ptrdiff_t>(nextFreeNode);
        nextFreeNode=this->treeGrow(adjMat,colOffset+1,idx,tempSortRow,tempSortIdx,nextFreeNode,midIdx);
    } else {
        innerChild[curNode]=-1;
    }
    
    return nextFreeNode;
}

void metricTreeCPP::searchRadius(ClusterSetCPP<size_t> &pointClust,ClusterSetCPP<double> &distClust,const double *point,const double *radius, const  size_t numPoints) const {
    //We do not initially know how many points will be in the search radius
    //and there does not appear to be any quick way to tell. Thus, we will
    //use vectors and just push points on as they are found.
    vector<size_t> idxRange;
    vector<double> distList;
    vector<size_t> clusterElsIdx;
    vector<double> clusterElsDist;
    size_t i, *clusterSizes;
    
    //Allocate temporary space
    clusterSizes=new size_t[N];
    
    //Reserve some space in the vectors so that it (hopefully) does not
    //have to call memory allocation routines very often.
    idxRange.reserve(64);
    distList.reserve(64);
    clusterElsIdx.reserve(64);
    clusterElsDist.reserve(64);
    
    for(i=0;i<numPoints;i++) {
        idxRange.clear();
        distList.clear();
        
        //Fill idxRange and distList with the values for this point.
        this->searchRadRecur(idxRange,distList,point+k*i,radius[i],0);
                
        //Append the vectors to the total collectioncollection of vectors found....
        clusterElsIdx.insert(clusterElsIdx.end(),idxRange.begin(),idxRange.end());
        clusterElsDist.insert(clusterElsDist.end(),distList.begin(),distList.end());
        
        clusterSizes[i]=idxRange.size();
    }

    //Now, place the values into the cluster sets to return.
    pointClust.initWithClusterSizes(clusterSizes,numPoints);
    distClust.initWithClusterSizes(clusterSizes,numPoints);
    memcpy(pointClust.clusterEls, &(clusterElsIdx[0]), clusterElsIdx.size()*sizeof(size_t));    
    memcpy(distClust.clusterEls, &(clusterElsDist[0]), clusterElsDist.size()*sizeof(double));

    //Free temporary space.
    delete[] clusterSizes;
}

void metricTreeCPP::searchRadRecur(vector<size_t> &idxRange,vector<double> &distList,const double *point,const double radius,const size_t curNode) const {
    double distCur;
        
    distCur=distEuclid(point,data+k*DATAIDX[curNode],k);
    if(distCur<=radius) {        
        idxRange.push_back(DATAIDX[curNode]);
        distList.push_back(distCur);
    }
    
    if(distCur+radius>=outerRadii[curNode]&&outerChild[curNode]!=-1) {
        this->searchRadRecur(idxRange,distList,point,radius,static_cast<size_t>(outerChild[curNode]));
    }

    if(distCur-radius<=innerRadii[curNode]&&innerChild[curNode]!=-1) {
        this->searchRadRecur(idxRange,distList,point,radius,static_cast<size_t>(innerChild[curNode]));
    }
}


metricTreeCPP::~metricTreeCPP() {
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
