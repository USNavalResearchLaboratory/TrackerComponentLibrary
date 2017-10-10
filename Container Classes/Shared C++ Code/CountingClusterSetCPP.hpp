/**COUNTINGCLUSTERSETCPP This header file implements two template classes:
 *              CountingClusterSetCPP and CountingClusterSetVecCPP
 * 
 **THE COUNTINGCLUSTERSETCPP CLASS
 *The CountingClusterSetCPP class overloads the [] operator such that one
 *can access arrays of size 1,2,3,4,5 using the notation A[c][r], where c
 *selects the array and r selects the element in the array. The elements of
 *the array are ordered by row and then by column. more general than this
 *class is the ClusterSetCPP class, as it allows non-sequential sizes for
 *the rows.
 * 
 *The elements of the clusters are held in the array clusterEls. 
 *
 *If the initWithNumClust method is used to allocate space for the class,
 *then the allocated memory is freed when the destructor is called.
 *One can also explicitly set the elements of the class to one's own
 *buffers, and never call initWithNumCLust, in which case the
 *destructor does not free the memory.
 *
 **THE COUNTINGCLUSTERSETVECCPP CLASS
 *This class overloads the () operator such that one can read multiple sets
 *of arrays of size 1,2,3,4,5 using the notation A(c,r,set), where c
 *selects the array, r selects the element in the array, and set selects
 *which set of elements to choose. The elements of the array are
 *ordered by row and then by column. Unlike the CountingClusterSetCPP
 *class, the access method only lets one read a value; it cannot be used to
 *write to a value.
 *
 *July 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
 **/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#ifndef COUNTINGCLUSTERSETCPP
#define COUNTINGCLUSTERSETCPP

template<typename T>
class CountingClusterSetCPP {
private:
    char *buffer;
public:
    size_t numClust;
    size_t totalNumEl;
    T *clusterEls;
    
    CountingClusterSetCPP<T>() {
        numClust=0;
        buffer=NULL;
    }

    void initWithNumClust(const size_t numClusters) {
        //Allocate space for the cluster elements.
        numClust=numClusters;
        totalNumEl=(numClust*(numClust+1))/2;
        
        //Allocate space for a CountingClusterSet of the specified size.
        buffer = new char[sizeof(T)*totalNumEl];
        clusterEls=reinterpret_cast<T*>(buffer);
    }

    T*operator[] (const size_t idx) const {
       //Overload the [] operator to return the start of a particular cluster.
       //Additional indexation can access elements in the cluster.
        const size_t clustStartIdx=(idx*(idx+1))/2;
        
        return clusterEls+clustStartIdx; 
    }

    ~CountingClusterSetCPP() {
        if(buffer!=NULL){
            delete[] buffer;
        }
    }
};

template<typename T>
class CountingClusterSetVecCPP {
private:
    char *buffer;
public:
    size_t numClust;
    size_t totalNumEl;
    size_t numSets;
    T *clusterEls;

    CountingClusterSetVecCPP<T>() {
        numClust=0;
        numSets=0;
        buffer=NULL;
    }

    void initWithNumClust(const size_t numClusters, const size_t numOfSets) {
        //Allocate space for the cluster elements.
        numClust=numClusters;
        numSets=numOfSets;
        totalNumEl=(numClust*(numClust+1))/2;
        
        //Allocate space for a CountingClusterSet of the specified size.
        buffer = new char[sizeof(T)*totalNumEl*numSets];
        clusterEls=reinterpret_cast<T*>(buffer);
    }

    T operator() (const size_t clustIdx, const size_t idxInClust, const size_t setIdx) const {
        //Overload the () operator to index a particular cluster, the
        //element in the cluster, and the set.
        const size_t clustStartOffset=(clustIdx*(clustIdx+1))/2;
        const size_t setStartOffset=setIdx*totalNumEl;
        
        return *(clusterEls+clustStartOffset+setStartOffset+idxInClust); 
    }

    ~CountingClusterSetVecCPP() {
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
