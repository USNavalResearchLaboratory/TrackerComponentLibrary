/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#ifndef KDTREECPP
#define KDTREECPP

#include <queue>
#include "ClusterSetCPP.hpp"

class kdTreeCPP {
public:
    size_t N;//The number of data points.
    size_t k;//The number of dimensions per data point.
    
    //LOSON lists the index of the next low node. The type ptrdiff_t is a
    //signed version of size_t and allows a LOSON of -1 to be used to
    //indicate that there is no LOSON child node.
    ptrdiff_t *LOSON;
    ptrdiff_t *HISON;//Lists the index of the next high node in the tree.
    size_t *DATAIDX;//The index of the data at a node.
    size_t *DISC;//The discriminating dimension index at the node.
    size_t *subtreeSizes;//Holds the size of all of the nodes in the subtree from a given node (including the given node).
    double *BMin;//An array of bounds of the children of each node for a given level.
    double *BMax;
    double *data;//A matrix of the data points. (A matrix ordered row-first).

    kdTreeCPP();
    kdTreeCPP(const size_t kDes, const size_t NDes);    
    void buildTreeFromBatch(const double *dataBatch);
    size_t *rangeCount(const double *rectMin,const double *rectMax,const size_t numRanges) const;
    void rangeQuery(ClusterSetCPP<size_t> &rangeClust,const double *rectMin,const double *rectMax,const  size_t numRanges) const;
    void findmBestNN(size_t *idxRange, double  *distSquared,const double *point,const size_t numPoints, const size_t m) const;
    ~kdTreeCPP();
    
private:
    char * buffer;
    size_t treeGrow(double *sortBatch, size_t *idx, double *tempSortRow,const size_t level,const size_t curNode,const size_t NSubTree);
    size_t rangeCountRecur(const size_t curNode,const double *rectMin,const double *rectMax) const;
    void rangeQueryRecur(const size_t curNode, const double *rectMin, const double *rectMax,size_t *idxRange, size_t &numFound, const size_t numInRange) const;
    void getSubtreeIdx(const size_t nodeIdx, size_t *idxRange, size_t &numFound) const;
    //The returned value is the number actually found.
    void mBestRecur(const size_t curNodeIdx, std::priority_queue<std::pair<double,size_t> > &mBestQueue, const double *point,const size_t m) const;
};
#endif

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
