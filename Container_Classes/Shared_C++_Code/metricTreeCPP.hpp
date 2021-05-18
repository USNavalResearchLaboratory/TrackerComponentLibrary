/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#ifndef METRICTREECPP
#define METRICTREECPP

#include <queue>
#include <vector>
#include "ClusterSetCPP.hpp"

class metricTreeCPP {
public:
    size_t N;//The number of data points.
    size_t k;//The number of dimensions per data point.
    
    size_t *DATAIDX;//The index of the data at a node.
    //The type ptrdiff_t is a signed version of size_t and allows a value 
    //of -1 to be used to denote no children.
    ptrdiff_t *innerChild;
    ptrdiff_t *outerChild;
    double *innerRadii;
    double *outerRadii;
    double *data;//A matrix of the data points.
    
    metricTreeCPP();
    metricTreeCPP(const size_t NDes, const size_t kDes);    
    void buildTreeFromBatch(const double *dataBatch);
    void searchRadius(ClusterSetCPP<size_t> &pointClust,ClusterSetCPP<double> &distClust,const double *point,const double *radius, const  size_t numPoints) const;

    ~metricTreeCPP();
    
private:
    char * buffer;
    size_t treeGrow(double *adjMat, const size_t colOffset, size_t *idx, double *tempSortRow, size_t *tempSortIdx,const size_t curNode, size_t NSubTree);
    void searchRadRecur(std::vector<size_t> &idxRange,std::vector<double> &distList,const double *point,const double radius,const size_t curNode) const;
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
