/**METRICTREECPPINT An interface class between the Matlab metric tree class
 *                  and a more efficient C++ metric tree class. This
 *                  function is meant to be called by the metric class in
 *                  Matlab; not directly by the user.
 *
 * *As the data of the true C++ class is stored in the CPPData input that is
 *passed to this function, passing garbage for the CPPData input can cause
 *Matlab to crash.
 *
 *The calling convention is
 *newTree.CPPData=metricTreeCPPInt('metricTreeCPP',k,N);
 *or
 *metricTreeCPPInt('buildTreeFromBatch',CPPData,dataBatch);
 *or
 *[retSet,distSet]=metricTreeCPPInt('searchRadius',CPPData,point,radius);
 *or
 *[DATAIDX,innerChild,outerChild,innerRadii,outerRadii,data]=metricTreeCPPInt('getAllData',CPPData);
 *or
 *N=metricTreeCPPInt('getN',CPPData);
 *or
 *k=metricTreeCPPInt('getk',CPPData);
 *or
 *metricTreeCPPInt('~metricTreeCPP',CPPData);
 *
 *December 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

//For strcmp
#include <cstring>
#include "ClusterSetCPP.hpp"
#include "MexValidation.h"
#include "metricTreeCPP.hpp"
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    char cmd[64];
    metricTreeCPP *theTree;
    
    if(nrhs>4) {
        mexErrMsgTxt("Too many inputs.");   
    }
    
    //Get the command string that is passed.
    mxGetString(prhs[0], cmd, sizeof(cmd));
    
    //prhs[0] is assumed to be the string telling
    if(!strcmp("metricTreeCPP", cmd)){
        size_t k, N;
        mxArray *retPtr;
        
        k=getSizeTFromMatlab(prhs[1]);
        N=getSizeTFromMatlab(prhs[2]);
        
        theTree =  new metricTreeCPP(k,N);
        
        //Convert the pointer to a Matlab matrix to return.
        retPtr=ptr2Matlab<metricTreeCPP*>(theTree);
        
        //Lock this mex file so that it can not be cleared until the object
        //has been deleted (This avoids a memory leak).
        mexLock();
        //Return the pointer to the tree
        plhs[0]=retPtr;
    } else if(!strcmp("buildTreeFromBatch",cmd)) {
        double *dataBatch;
        
        //Get the pointer back from Matlab.
        theTree=Matlab2Ptr<metricTreeCPP*>(prhs[1]);   
        
        checkRealDoubleArray(prhs[2]);
        dataBatch=mxGetDoubles(prhs[2]);
        
        theTree->buildTreeFromBatch(dataBatch);
    } else if(!strcmp("searchRadius",cmd)) {
        size_t numPoints;
        ClusterSetCPP<size_t> pointClust;
        ClusterSetCPP<double> distClust;
        double *point;
        double *radius;
        mxArray *clustParams[3];
        
        //Get the inputs
        theTree=Matlab2Ptr<metricTreeCPP*>(prhs[1]);
        checkRealDoubleArray(prhs[2]);
        checkRealDoubleArray(prhs[3]);
        point=mxGetDoubles(prhs[2]);
        radius=mxGetDoubles(prhs[3]);
        
        numPoints=mxGetN(prhs[2]);
        if(mxGetM(prhs[2])!=theTree->k){
            mexErrMsgTxt("Invalid point size passed.");
        }
        
        //Run the search; rangeCluster now contains the results.
        theTree->searchRadius(pointClust,distClust,point,radius, numPoints);
        
        //Put the results into an instance of the ClusterSet container 
        //class in Matlab.
        clustParams[0]=unsignedSizeMat2Matlab(pointClust.clusterEls,pointClust.totalNumEl,1);
        clustParams[1]=unsignedSizeMat2Matlab(pointClust.clusterSizes,pointClust.numClust,1);
        clustParams[2]=unsignedSizeMat2Matlab(pointClust.offsetArray,pointClust.numClust,1);
        
        //Return a ClusterSet containing the appropriate data.
        mexCallMATLAB(1, &(plhs[0]), 3,  clustParams, "ClusterSet");
        
        //If the distances should also be returned.
        if(nlhs>1) {
            clustParams[0]=doubleMat2Matlab(distClust.clusterEls,distClust.totalNumEl,1);
            clustParams[1]=unsignedSizeMat2Matlab(distClust.clusterSizes,distClust.numClust,1);
            clustParams[2]=unsignedSizeMat2Matlab(distClust.offsetArray,distClust.numClust,1);
        
            //Return a ClusterSet containing the appropriate data.
            mexCallMATLAB(1, &(plhs[1]), 3,  clustParams, "ClusterSet");
        }
    } else if(!strcmp("~metricTreeCPP", cmd)){
        theTree=Matlab2Ptr<metricTreeCPP*>(prhs[1]);

        delete theTree;
        //Unlock the mex file allowing it to be cleared.
        mexUnlock();
    } else if(!strcmp("getAllData", cmd)){
        size_t N;
        size_t k;
        
        theTree=Matlab2Ptr<metricTreeCPP*>(prhs[1]);
        N=theTree->N;
        k=theTree->k;
        
        plhs[0]=unsignedSizeMat2Matlab(theTree->DATAIDX,N, 1);
        if(nlhs>1) {
            plhs[1]=signedSizeMat2Matlab(theTree->innerChild,N, 1);
            if(nlhs>2) {
                plhs[2]=signedSizeMat2Matlab(theTree->outerChild,N, 1);
                if(nlhs>3) {
                    plhs[3]=doubleMat2Matlab(theTree->innerRadii,N, 1);
                    if(nlhs>4) {
                        plhs[4]=doubleMat2Matlab(theTree->outerRadii,N, 1);
                        if(nlhs>5) {
                            plhs[5]=doubleMat2Matlab(theTree->data,k, N);
                        }
                    }
                }
            }   
        }
    } else if(!strcmp("getN", cmd)) {
        theTree=Matlab2Ptr<metricTreeCPP*>(prhs[1]);
        plhs[0]=unsignedSizeMat2Matlab(&(theTree->N),1,1);
    } else if(!strcmp("getk", cmd)) {
        theTree=Matlab2Ptr<metricTreeCPP*>(prhs[1]);
        plhs[0]=unsignedSizeMat2Matlab(&(theTree->k),1,1);
    }else {
        mexErrMsgTxt("Invalid string passed to metricTreeCPPInt.");
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
