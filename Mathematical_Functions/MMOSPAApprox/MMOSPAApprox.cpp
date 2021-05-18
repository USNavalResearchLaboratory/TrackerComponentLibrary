/**MMOSPAAPPROX  Find the approximate minimum mean optimal sub-pattern
*                assignment (MMOSPA) estimate from a set of weighted 
*                discrete sets of target estimates using 2D assignment in a
*                forward-backward algorithm.
*
* INPUTS:   x        An xDim X numTar X numHyp hypermatrix that holds
*                    numHyp hypotheses each consisting or numTar targets
*                    (or generic vectors) with xDim dimensions per target
*                    (per generic vector). One should not pass x values
*                    that are a combination of position and velocity
*                    components, because the OSPA error in that instance
*                    has no meaning. In general, the components of each
*                    target vector for each hypothesis will be position
*                    only.
*           w        A numHyp X 1 vector of the probabilities of each of
*                    the numHyp hypotheses in x. The elements must all be
*                    positive and sum to one.
*           numScans An optional parameter >=1 specifying how many forward
*                    scans of the approximate algorithm to perform. Even 
*                    though more scans can improve the estimate, global
*                    convergence is not guaranteed. The default is 1 if
*                    this parameter is not provided.
* 
* OUTPUTS:  MMOSPAEst The approximate xDim X numTar MMOSPA estimate.
*           orderList A numTarXnumHyp matrix specifying the ordering of the
*                     targets in each hypothesis that went into the
*                     approximate MMOSPA estimate.
* 
* Given a set of numHyp hypotheses, the standard expected value minimizes
* the mean squared error. The standard expected value is just
* the weighted sum of the hypothese for all of the targets. In the MMOSPA
* estimate, a weighted sum is used, but the ordering of the target states
* is modified to minimize the expected value of a specific version of the
* OSPA metric. orderList holds indices of the reordered states going into
* the MMOSPA estimate such that the MMOSPA estimate uses  
* x(:,orderList(:,curHyp),curHyp) instead of x(:,:,curHyp) for the ordered
* matrix of target states in hypothesis curHyp. Finding the MMOSPA estimate
* is generally NP-hard for more than two hypotheses. Thus, this algorithm
* approximates the MMOSPA estimate.
*
* The basic algorithm uses sequential 2D assignment going forward to
* approximate the MMOSPA estimate. If desired, assignments can be
* reevaluated in additional backward-forward passes to try to obtain an
* approximation close to the true MMOSPA estimate.
*
* The algorithm as well as the concept of MOSPA error are described in
* detail in [1].
*
* The algorithm can be compiled for use in Matlab  using the 
* CompileCLibraries function.
*
* The algorithm is run in Matlab using the command format
* [MMOSPAEst,orderList]=MMOSPAApprox(x,w,numScans);
*
*REFERENCES:
*[1] D. F. Crouse, "Advances in displaying uncertain estimates of multiple
*    targets," in Proceedings of SPIE: Signal Processing, Sensor Fusion,
*    and Target Recognition XXII, vol. 8745, Baltimore, MD, Apr. 2013.
*
* November 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.*/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

/*This header is required by Matlab*/
#include "mex.h"
#include "MexValidation.h"
#include "MMOSPAApproxCPP.hpp"
#include <algorithm>

using namespace std;

void mexFunction(const int nlhs, mxArray *plhs[], const int nrhs, const mxArray *prhs[]) {
    size_t xDim,numTar,numHyp,numScans;
    mxArray *MMOSPAEstMATLAB,*orderListMATLAB;//These will hold the values to be returned.
    const mwSize *xDims;
    double *MMOSPAEst,*x,*w;
    size_t *orderList;
    
    if(nrhs<2){
        mexErrMsgTxt("Not enough inputs.");
        return;
    }
    
    if(nrhs<3){
        numScans=1;
    }else {
        numScans=getSizeTFromMatlab(prhs[2]);
        
        if(numScans<1) {
            mexErrMsgTxt("Invalid number of scans specified.");
            return;
        }
    }
    
    if(nrhs>3) {
        mexErrMsgTxt("Too many inputs.");
        return;
    }
    
    /*Verify the validity of the x and w parameters.*/
    checkRealDoubleArray(prhs[0]);
    checkRealDoubleArray(prhs[1]);
    
    xDims=mxGetDimensions(prhs[0]);
    xDim=xDims[0];
    numTar=xDims[1];//Matlab should always provide at least a 2D dimension vector.
    
    if(mxGetNumberOfDimensions(prhs[0])<3){
        numHyp=1;
    } else if(mxGetNumberOfDimensions(prhs[0])>3) {
        mexErrMsgTxt("The first parameter has too many dimensions.");
        return;
    } else {
        numHyp=xDims[2];
    }
    
    //Check the dimensionality of the second input
    if(mxGetNumberOfDimensions(prhs[1])>2){
        mexErrMsgTxt("The second parameter has too may dimensions.");
        return;
    }
    
    xDims = mxGetDimensions(prhs[1]);
    if(xDims[0]!=numHyp||xDims[1]!=1) {
        mexErrMsgTxt("The dimensionality of the second parameter is inconsistent.");
        return;
    }
    
    //Allocate space for the return variables.
    MMOSPAEstMATLAB = mxCreateNumericMatrix(xDim,numTar,mxDOUBLE_CLASS,mxREAL);
    orderListMATLAB=allocUnsignedSizeMatInMatlab(numTar,numHyp);

    MMOSPAEst=mxGetDoubles(MMOSPAEstMATLAB);
    if(sizeof(size_t)==4) {//32 bit
        orderList=reinterpret_cast<size_t*>(mxGetUint32s(orderListMATLAB));
    } else {//64 bit
        orderList=reinterpret_cast<size_t*>(mxGetUint64s(orderListMATLAB));
    }

/*Get the matrices*/
    x=mxGetDoubles(prhs[0]); 
    w=mxGetDoubles(prhs[1]);
    
    //Run the algorithm
    MMOSPAApproxCPP(MMOSPAEst,
                    orderList,
                    x,
                    w,
                    xDim,
                    numTar,
                    numHyp,
                    numScans);
    
/*Set the outputs*/
    plhs[0]=MMOSPAEstMATLAB;
    if(nlhs>1) {
        /*Convert C++ indices to Matlab indices*/
        for_each(orderList, orderList+numTar*numHyp, increment<size_t>);
        plhs[1]=orderListMATLAB;
    } else {
        mxDestroyArray(orderListMATLAB);
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
