/**CLUSTERSETCPP This class overloads the [] operator such that an instance
 *               A can access item c in cluster r using the notation
 *               A[c][r].
 * 
 * The elements of the clusters are held in thearray clusterEls. The offset
 * of the beginning of cluster c is given by offsetArray[c]. The size of
 * cluster c is given by clusterSizes[c]. the number of clusters and the
 * total number of elements are given by numClust and totalNumEl.
 *
 * If the initWithClusterSizes method is used to allocate space for the
 * class, then the allocated memory is freed when the destructor is called.
 * One can also explicitly set the elements of the class to one's own
 * buffers, and never call initWithClusterSizes, in which case the
 * destructor does not free the memory.
 *
 * The entire class is implemented in the header file.
 *
 *December 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
 **/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#ifndef CLUSTERSETCPP
#define CLUSTERSETCPP

//A class to store and access collections of indices.
//The class is initialized with an array of clusterSizes, which the class
//becomes responsible to free.

//For the accumulate function
#include <numeric>
//For memcpy
#include <cstring>
#include <cstddef>
//For sort
#include <algorithm>

template<typename T>
class ClusterSetCPP {
private:
    char *buffer;
public:
    size_t numClust;
    size_t totalNumEl;
    T *clusterEls;
    size_t *offsetArray;
    size_t *clusterSizes;
    
    ClusterSetCPP<T>() {
        numClust=0;
        buffer=NULL;
    }

    void initWithClusterSizes(const size_t *clustSizes,const size_t numClusters) {
    //Allocate space for the cluster elements and fill in the offsetArray
    //and clusterSizes arrays.
        char *curIdx;
        size_t i;
        
        numClust=numClusters;
        totalNumEl=std::accumulate(clustSizes,clustSizes+numClusters,static_cast<size_t>(0));
        //Allocate space for a ClusterSet of the specified size and then
        //separate out the pointers for the elements and offset array.
        buffer = new char[sizeof(T)*totalNumEl+sizeof(size_t)*numClusters*2];
        curIdx=buffer;
        clusterEls=reinterpret_cast<T*>(curIdx);
        curIdx+=sizeof(T)*totalNumEl;
        offsetArray=reinterpret_cast<size_t*>(curIdx);
        curIdx+=sizeof(size_t)*numClusters;
        clusterSizes=reinterpret_cast<size_t*>(curIdx);
        
        memcpy(clusterSizes,clustSizes,numClusters*sizeof(size_t));

        offsetArray[0]=0;
        for(i=1;i<numClusters;i++) {
            offsetArray[i]=offsetArray[i-1]+clustSizes[i-1];
        }
    }

    T*operator[] (const size_t idx) const {
    //Overload the [] operator to return the start of a particular cluster.
    //Additional indexation can access elements in the cluster.
    
       return clusterEls+offsetArray[idx]; 
    }
    
    ~ClusterSetCPP() {
        if(buffer!=NULL){
            delete[] buffer;
        }
    }
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
