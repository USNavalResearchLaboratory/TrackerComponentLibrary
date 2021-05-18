/**RELVECADD Special relativistic addition of velocity vectors. In the
 *           inertial reference coordinate system, an observer moves with
 *           constant velocity vObsFrame. In the coordinate system of the
 *           observer, an object moves with constant velocity vObjInFrame.
 *           This computes the velocity of the observed object in the
 *           inertial reference frame. Under special relativity, it is not
 *           just vObsFrame+vObjInFrame as it is in Newtonian mechanics.
 *
 *INPUTS:  vObsFrame The 3XN set of N velocity vectors in meters per second
 *                  of the observer with respect to the inertial reference
 *                  coordinate system.
 *      vObjInFrame The 3XN set of N velocity vectors in meters per second
 *                  of the object with respect to the observer's coordinate
 *                  system.
 *
 *OUTPUTS: vObjRef  The 3XN set of velocity vectors in vObjInFrame
 *                 transformed into the reference coordinate system.
 *
 *The formulae for special relativistic velocity addition is derived in
 *Chapter 1.4 of [1]. The magnitudes of vObsFrame and vObjInFrame must both
 *be less than the speed of light (Constants.speedOfLight).
 *
 *The algorithm can be compiled for use in Matlab  using the 
 *CompileCLibraries function.
 *
 *The algorithm is run in Matlab using the command format
 *vObjRef=relVecAdd(vObsFrame,vObjInFrame)
 *
 *REFERENCES:
 *[1] G. Ludyk, Einstein in Matrix Form: Exact Derivation of the Theory of
 *    Special and General Relativity without Tensors. Heidelberg: Springer,
 *    2013.
 *
 *March 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
 **/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

/*This header is required by Matlab.*/
#include "mex.h"
#include "MexValidation.h"
#include "relFuncs.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    double c, *vObsFrame, *vObjInFrame, *vObjRef;
    mxArray *retMATLAB;
    size_t i, numVec;
    
    if(nrhs!=2) {
        mexErrMsgTxt("Incorrect number of inputs.");
        return;
    }
    
    if(nlhs>1) {
        mexErrMsgTxt("Wrong number of outputs.");
        return;
    }

    numVec=mxGetN(prhs[0]);
    if(numVec!=mxGetN(prhs[1])||mxGetM(prhs[0])!=3||mxGetM(prhs[1])!=3) {
        mexErrMsgTxt("The input vectors have the wrong dimensionality.");
        return;
    }
    
    checkRealDoubleArray(prhs[0]);
    vObsFrame=mxGetDoubles(prhs[0]);
    checkRealDoubleArray(prhs[1]);
    vObjInFrame=mxGetDoubles(prhs[1]);
    
    c=getScalarMatlabClassConst("Constants","speedOfLight");
    
    //Allocate space for the return values.
    retMATLAB=mxCreateDoubleMatrix(3,numVec,mxREAL);
    vObjRef=mxGetDoubles(retMATLAB);

    for(i=0;i<numVec;i++) {
        relVecAddC(c,vObsFrame+3*i,vObjInFrame+3*i,vObjRef+3*i);
    }
    
    plhs[0]=retMATLAB;
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
