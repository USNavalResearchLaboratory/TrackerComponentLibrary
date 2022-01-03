/**INVCONDNUMBERAPPROX Find the approximate inverse condition number of a
*      real, square matrix using the l1 norm. The output is very similar
*       to the rcond function built into Matlab, when given a real matrix.
*       It approximates the inverse condition number without directly
*       inverting the matrix X. See the comments to the Matlam
*       implementation for more information.
*
* The algorithm can be compiled for use in Matlab  using the 
* CompileCLibraries function.
*
* The algorithm is run in Matlab using the command format
* rCond=invCondNumberApprox(X)
*
*June 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
*/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#include "mex.h"

/* This header validates inputs and includes a header needed to handle
 * Matlab matrices.*/
#include "MexValidation.h"
#include "basicMatOpsEigen.hpp"

double approxInvMat1Norm(const Eigen::MatrixXd &X, Eigen::MatrixXd &v);
double invCondNumberApprox(const Eigen::MatrixXd &X);

void mexFunction(const int nlhs, mxArray *plhs[], const int nrhs, const mxArray *prhs[]) {

    if(nrhs<1){
        mexErrMsgTxt("Incorrect number of inputs.");
        return;
    }
    
    if(nlhs>1) {
        mexErrMsgTxt("Invalid number of outputs.");
        return;
    }
    
    if(mxIsEmpty(prhs[0])) {//Special case.
        const double infVal=std::numeric_limits<double>::infinity();
        plhs[0]=doubleMat2Matlab(&infVal,1,1);
        return;
    }
    
    checkDoubleArray(prhs[0]);
    const size_t N=mxGetM(prhs[0]);
    
    if(mxGetN(prhs[0])!=N) {
        mexErrMsgTxt("The input matrix must be square.");
        return;
    }

    double *X=mxGetDoubles(prhs[0]);
    const Eigen::Map<Eigen::MatrixXd> XEigen(X,N,N);
    
    const double condNum=invCondNumberApprox(XEigen);
    
    plhs[0]=doubleMat2Matlab(&condNum,1,1);
}

/*LICENSE:
*
*The source code is in the public domain and not licensed or under
*copyright. The information and software may be used freely by the public.
*As required by 17 U.S.C. 403, third parties producing copyrighted works
*consisting predominantly of the material produced by U.S. government
*agencies must provide notice with such work(s) identifying the U.S.
*Government material incorporated and stating that such material is not
*subject to copyright protection.
*
*Derived works shall not identify themselves in a manner that implies an
*endorsement by or an affiliation with the Naval Research Laboratory.
*
*RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
*SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY THE NAVAL
*RESEARCH LABORATORY FOR ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS
*OF RECIPIENT IN THE USE OF THE SOFTWARE.*/
 