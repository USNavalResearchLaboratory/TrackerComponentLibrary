/**KDTREECPPINT An interface class between the Matlab kd-tree class and a
 *              more efficient C++ kd tree class. This function is meant to
 *              be called by the kdTree class in Matlab; not directly by
 *              the user.
 *
 *As the data of the true C++ class is stored in the CPPData input that is
 *passed to this function, passing garbage for the CPPData input can cause
 *Matlab to crash.
 *
 *The function is called as
 *newTree.CPPData=kdTreeCPPInt('kdTreeCPP',k,N);
 *or
 *N=kdTreeCPPInt('getN',CPPData);
 *or
 *k=kdTreeCPPInt('getk',CPPData);
 *or
 *kdTreeCPPInt('buildTreeFromBatch',CPPData,dataBatch);
 *or
 *retSet=kdTreeCPPInt('rangeQuery',CPPData,rectMin,rectMax);
 *or
 *numInRange=kdTreeCPPInt('rangeCount',CPPData,rectMin,rectMax);
 *or
 *[idxRange, distSquared]=kdTreeCPPInt('findmBestNN',CPPData,point,m);
 *or
 *[LOSON,HISON,DATAIDX,DISC,subtreeSizes,BMin,BMax,data]=kdTreeCPPInt('getAllData',CPPData);
 *or
 *kdTreeCPPInt('~kdTreeCPP',CPPData);
 *
 *December 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

//For strcmp
#include <cstring>
#include "ClusterSetCPP.hpp"
#include "MexValidation.h"
#include "kdTreeCPP.hpp"
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    char cmd[64];
    kdTreeCPP *theTree;
 
    if(nrhs>5) {
        mexErrMsgTxt("Too many inputs.");
    }
    
    //Get the command string that is passed.
    mxGetString(prhs[0], cmd, sizeof(cmd));
    
    //prhs[0] is assumed to be the string telling
    if(!strcmp("kdTreeCPP", cmd)){
        size_t k, N;
        mxArray *retPtr;
        
        if(nrhs!=3) {
            mexErrMsgTxt("Wrong Number of inputs.");
        }
        
        k=getSizeTFromMatlab(prhs[1]);
        N=getSizeTFromMatlab(prhs[2]);
        
        theTree = new kdTreeCPP(k,N);
        
        //Convert the pointer to a Matlab matrix to return.
        retPtr=ptr2Matlab<kdTreeCPP*>(theTree);
        
        //Lock this mex file so that it can not be cleared until the object
        //has been deleted (This avoids a memory leak).
        mexLock();
        //Return the pointer to the tree
        plhs[0]=retPtr;
    } else if(!strcmp("buildTreeFromBatch",cmd)) {
        double *dataBatch;
        
        if(nrhs!=3) {
            mexErrMsgTxt("Wrong Number of inputs.");
        }
        
        //Get the pointer back from Matlab.
        theTree=Matlab2Ptr<kdTreeCPP*>(prhs[1]);   
        
        checkRealDoubleArray(prhs[2]);
        dataBatch=mxGetDoubles(prhs[2]);
        
        theTree->buildTreeFromBatch(dataBatch);
    } else if(!strcmp("rangeQuery",cmd)) {
        size_t numRects;
        ClusterSetCPP<size_t> rangeClust;
        double *rectMin, *rectMax;
        mxArray *clustParams[3];

        if(nrhs!=4) {
            mexErrMsgTxt("Wrong Number of inputs.");
        }
        
        //Get the inputs
        theTree=Matlab2Ptr<kdTreeCPP*>(prhs[1]);
        checkRealDoubleArray(prhs[2]);
        checkRealDoubleArray(prhs[3]);
        rectMin=mxGetDoubles(prhs[2]);
        rectMax=mxGetDoubles(prhs[3]);
        numRects=mxGetN(prhs[2]);

        //Run the search; rangeCluster now contains the results.
        theTree->rangeQuery(rangeClust,rectMin,rectMax,numRects);
        
        //Put the results into an instance of the ClusterSet container 
        //class in Matlab.
        clustParams[0]=unsignedSizeMat2Matlab(rangeClust.clusterEls,rangeClust.totalNumEl,1);
        clustParams[1]=unsignedSizeMat2Matlab(rangeClust.clusterSizes,rangeClust.numClust,1);
        clustParams[2]=unsignedSizeMat2Matlab(rangeClust.offsetArray,rangeClust.numClust,1);
        
        //Return a ClusterSet containing the appropriate data.
        mexCallMATLAB(1, plhs, 3,  clustParams, "ClusterSet");        
    } else if(!strcmp("rangeCount",cmd)) {
        size_t numRects;
        double *rectMin, *rectMax;
        size_t *rangeCounts;
        
        if(nrhs!=4) {
            mexErrMsgTxt("Wrong Number of inputs.");
        }
        
        //Get the inputs
        theTree=Matlab2Ptr<kdTreeCPP*>(prhs[1]);
        checkRealDoubleArray(prhs[2]);
        checkRealDoubleArray(prhs[3]);
        rectMin=mxGetDoubles(prhs[2]);
        rectMax=mxGetDoubles(prhs[3]);
        numRects=mxGetN(prhs[2]);
        
        //Run the search
        rangeCounts=theTree->rangeCount(rectMin,rectMax,numRects);
        
        //Process the output
        plhs[0]=unsignedSizeMat2Matlab(rangeCounts,numRects, 1);
    } else if(!strcmp("findmBestNN",cmd)){
        double *point;
        size_t m, numPoints;
        mxArray *idxRangeMATLAB,*distSquaredMATLAB;
        size_t *idxRange;
        double *distSquared;
        
        if(nrhs!=4) {
            mexErrMsgTxt("Wrong Number of inputs.");
        }
        
        //Get the inputs
        theTree=Matlab2Ptr<kdTreeCPP*>(prhs[1]);
        checkRealDoubleArray(prhs[2]);        
        point=mxGetDoubles(prhs[2]);
        numPoints=mxGetN(prhs[2]);
        m=getSizeTFromMatlab(prhs[3]);

        //Allocate space for the return variables.
        idxRangeMATLAB=allocUnsignedSizeMatInMatlab(m, numPoints);
        
        distSquaredMATLAB=mxCreateNumericMatrix(m,numPoints,mxDOUBLE_CLASS,mxREAL);
        if(sizeof(size_t)==4) {//32 bit
            idxRange=(size_t*)mxGetUint32s(idxRangeMATLAB);
        } else {//64 bit
            idxRange=(size_t*)mxGetUint64s(idxRangeMATLAB);
        }

        distSquared=mxGetDoubles(distSquaredMATLAB);
        
        theTree->findmBestNN(idxRange,distSquared, point, numPoints, m);

        plhs[0]=idxRangeMATLAB;
        if(nlhs>1){
            plhs[1]=distSquaredMATLAB;
        }
    } else if(!strcmp("~kdTreeCPP", cmd)){
        theTree=Matlab2Ptr<kdTreeCPP*>(prhs[1]);

        delete theTree;
        //Unlock the mex file allowing it to be cleared.
        mexUnlock();
    } else if(!strcmp("getAllData", cmd)){
        size_t N;
        size_t k;
        
        if(nrhs!=2) {
            mexErrMsgTxt("Wrong Number of inputs.");
        }
        
        theTree=Matlab2Ptr<kdTreeCPP*>(prhs[1]);
        N=theTree->N;
        k=theTree->k;
        
        switch(nlhs) {
            case 8:
                plhs[7]=doubleMat2Matlab(theTree->data,k, N);
            case 7:
                plhs[6]=doubleMat2Matlab(theTree->BMax,k, N);
            case 6:
                plhs[5]=doubleMat2Matlab(theTree->BMin,k, N);
            case 5:
                plhs[4]=unsignedSizeMat2Matlab(theTree->subtreeSizes,N, 1);
            case 4:
                plhs[3]=unsignedSizeMat2Matlab(theTree->DISC,N, 1);
            case 3:
                plhs[2]=unsignedSizeMat2Matlab(theTree->DATAIDX,N, 1);
            case 2:
                plhs[1]=signedSizeMat2Matlab(theTree->HISON,N, 1);
            default:
                plhs[0]=signedSizeMat2Matlab(theTree->LOSON,N, 1);
        }
    }else if(!strcmp("getk", cmd)) {
        if(nrhs!=2) {
            mexErrMsgTxt("Wrong Number of inputs.");
        }
        
        theTree=Matlab2Ptr<kdTreeCPP*>(prhs[1]);
        plhs[0]=unsignedSizeMat2Matlab(&(theTree->k),1,1);
    }else if(!strcmp("getN", cmd)) {
        if(nrhs!=2) {
            mexErrMsgTxt("Wrong Number of inputs.");
        }

        theTree=Matlab2Ptr<kdTreeCPP*>(prhs[1]);
        plhs[0]=unsignedSizeMat2Matlab(&(theTree->N),1,1);
    }else {
        mexErrMsgTxt("Invalid string passed to kdTreeCPPInt.");
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
