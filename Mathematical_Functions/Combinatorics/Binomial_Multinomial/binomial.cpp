/**BINOMIAL A C++ implementation of a function to evaluate the binomial
*          coefficient n choose k, which is the total number of ways of
*          choosing k items from a total selection of n where the ordering
*          of the selection does not matter. See the Matlab implementation
*          for more comments.
*
*The algorithm can be compiled for use in Matlab  using the 
*CompileCLibraries function.
*
*The algorithm is run in Matlab using the command format
*val=binomial(n,k);
*
*December 2024 David F. Crouse, Naval Research Laboratory, Washington D.C.*/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release. */

#ifdef _MSC_VER
#pragma warning( push )
#pragma warning( disable : 4514 )
#endif

#include "mex.h"

#ifdef _MSC_VER
#pragma warning( pop )
//Get rid of useless Spectre mitigation warnings.
#pragma warning( disable : 5045 )
#endif

/* This header validates inputs and includes a header needed to handle
 * Matlab matrices.*/
#include "MexValidation.h"
#include "combinatorialFuns.hpp"

void mexFunction(const int nlhs, mxArray *plhs[], const int nrhs, const mxArray *prhs[]) {

    if(nrhs<1||nrhs>2){
        mexErrMsgTxt("Incorrect number of inputs.");
        return;
    }
    
    if(nlhs>1) {
        mexErrMsgTxt("Wrong number of outputs.");
        return;
    }

    const size_t n=getSizeTFromMatlab(prhs[0]);
    const size_t k=getSizeTFromMatlab(prhs[1]);

    if(k>n) {
        plhs[0]=mxCreateDoubleScalar(0);
        return;
    }

    double val;
    if(n<=67) {
        val=static_cast<double>(binomialFromTableCPP(n,k));
    } else {
        val=binomialCPP(n,k);
    }

    plhs[0]=mxCreateDoubleScalar(val);
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