/**NDIM2INDEX A C++ implementation of a function such that given the number
 *           of dimensions of a multidimensional array/matrix and a list of
 *           indices for each dimension, the linear index of the point
 *           addressed by the indices is found. This function is like the
 *           sub2ind function that comes with Matlab, except the sub2ind
 *           function cannot take a vector of indices. Rather, to address n
 *           dimensions, it takes n inputs. That is not useful if one
 *           wishes to write a function where the number of dimensions
 *           needed in a hypermatrix is not known a priori, such as when
 *           handling multivariate polynomials.
 *
 *INPUTS: dims A 1XnumDim or numDimX1 vector holding the size of each of
 *             the nDim dimensions. The values are >=1. Alternatively, a
 *             scalar value can be passed if all of the dimensions have the
 *             same size.
 *     indices A numDim1XnumIdx matrix of indices where each of the
 *             corresponding numIdx linear indices is desired. The
 *             indexation is Matlab-style, starting from 1, not 0.
 *             Normally, numDim1=numDim. However, if numDim1<numDim, then
 *             it is assumed that the omitted indices are all ones. This
 *             helps deal with singleton dimensions in some instances.
 *
 *OUTPUTS: idx A numIdxX1 vector of linear indices corresponding to each of
 *             the multidimensional index sets given in indices.
 *
 *The function index2NDim is the inverse of this function.
 *
 *The algorithm can be compiled for use in Matlab using the 
 *CompileCLibraries function.
 *
 *The algorithm is run in Matlab using the command format
 *indices=nDim2Index(dims,idx);
 *
 *June 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#include "mex.h"
/* This header validates inputs and includes a header needed to handle
 * Matlab matrices.*/
#include "MexValidation.h"
//For std::copy
#include <algorithm>
//For floor
#include <cmath>

void mexFunction(const int nlhs, mxArray *plhs[], const int nrhs, const mxArray *prhs[]) {
    size_t numDim,numDim1;
    double *dims, *multIdx, *indices,*indicesIn,*indicesCopy=NULL;
    bool dimIsScalar;
    size_t i,j, numIdx;
    mxArray *retMat;
    double *idx;
    
    if(nrhs!=2) {
        mexErrMsgTxt("Incorrect number of inputs.");
        return;
    }
    
    if(nlhs>1) {
        mexErrMsgTxt("Incorrect number of outputs.");
        return;  
    }
    
    //If an empty matrix is passed, return an empty matrix.
    if(mxIsEmpty(prhs[1])) {
        plhs[0]=mxCreateDoubleMatrix(0,0,mxREAL);
        return;
    }
    
    numDim=mxGetNumberOfElements(prhs[0]);
    numDim1=mxGetM(prhs[1]);
    numIdx=mxGetN(prhs[1]);
    dimIsScalar=mxIsScalar(prhs[0]);

    if(numDim1>numDim&&!dimIsScalar) {
        mexErrMsgTxt("Inconsistent input dimensionalities given.");
        return;    
    }

    checkRealDoubleArray(prhs[0]);
    checkRealDoubleArray(prhs[1]);
    
    dims=mxGetDoubles(prhs[0]);
    indicesIn=mxGetDoubles(prhs[1]);

    multIdx=reinterpret_cast<double*>(mxMalloc(numDim1*sizeof(double)));
    multIdx[0]=1;
    if(dimIsScalar==true) {
        for(i=1;i<numDim1;i++) {
            multIdx[i]=multIdx[i-1]*dims[0];
        }
    } else {
        for(i=1;i<numDim1;i++) {
            multIdx[i]=multIdx[i-1]*dims[i-1];
        }
    }

    if(numDim1<numDim) {
        indicesCopy=reinterpret_cast<double*>(mxMalloc(numDim*numIdx*sizeof(double)));

        for(j=0;j<numIdx;j++) {
            std::copy(indicesIn+j*numDim1,indicesIn+(j+1)*numDim1,indicesCopy+numDim*j);
            for(i=numDim1;i<numDim;i++) {
                indicesCopy[i+numDim*j]=1;
            }
        }

        numDim1=numDim;
        indices=indicesCopy;
    } else{
        indices=indicesIn;     
    }
    
    retMat=mxCreateDoubleMatrix(numIdx,1,mxREAL);
    idx=mxGetDoubles(retMat);
    
    for(j=0;j<numIdx;j++) {
        idx[j]=1;
        
        for(i=0;i<numDim1;i++) {
            idx[j]+=(indices[i+numDim1*j]-1)*multIdx[i];
        }
    }

    plhs[0]=retMat;
    mxFree(multIdx);
    if(indicesCopy!=NULL) {
        mxFree(indicesCopy);
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
