/**EXACTPAIRSUM A C++ implementation of a function to evaluates a+b such
*       that s is the floating point value given by a+b with fixed
*       precision and e is a value such that, with infinite precision,
*       (s+e)=(a+b). See the comments in the Matlab implementation of this
*       file for more information on the inputs to the function and how it
*       works.
*
* The algorithm can be compiled for use in Matlab  using the 
* CompileCLibraries function.
*
* The algorithm is run in Matlab using the command format
* [s,e]=exactPairSum(a,b)
*
*February 2025 David F. Crouse, Naval Research Laboratory, Washington D.C.
*/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#ifdef _MSC_VER
#pragma warning( push )
#pragma warning( disable : 4514 )
#endif

#include "mex.h"

#ifdef _MSC_VER
#pragma warning( pop )
#endif

/* This header validates inputs and includes a header needed to handle
 * Matlab matrices.*/
#include "MexValidation.h"
#include <cmath>

template<class T>//T is a floating point data type.
void exactPairSum(const T a,const T b, T &s, T &e) {
    s=a+b;
    if(std::abs(a)>std::abs(b)) {
        const T temp=s-a;
        e=b-temp;
    } else {
        const T temp=s-b;
        e=a-temp;
    }
}

void mexFunction(const int nlhs, mxArray *plhs[], const int nrhs, const mxArray *prhs[]) {
    mxArray *sMat;
    mxArray *eMat;

    if(nrhs!=2) {
        mexErrMsgTxt("Incorrect number of inputs.");
        return;
    }

    if(nlhs>2) {
        mexErrMsgTxt("Incorrect number of outputs.");
        return;  
    }

    checkRealDoubleVector(prhs[0]);
    checkRealDoubleVector(prhs[1]);

    const double *a=mxGetDoubles(prhs[0]);
    const double *b=mxGetDoubles(prhs[1]);
    const size_t aLen=mxGetNumberOfElements(prhs[0]);
    const size_t bLen=mxGetNumberOfElements(prhs[1]);

    const bool aScal=(aLen==1);
    const bool bScal=(bLen==1);
    
    if(aScal&&bScal) {
        sMat=mxCreateDoubleMatrix(1, 1, mxREAL);
        eMat=mxCreateDoubleMatrix(1, 1, mxREAL);
        double *s=mxGetDoubles(sMat);
        double *e=mxGetDoubles(eMat);

        exactPairSum<double>(*a,*b, *s, *e);
    } else if(!aScal&&!bScal) {
        if(aLen!=bLen) {
            mexErrMsgTxt("a and b have an inconsistent number of elements.");
        }
        sMat=mxCreateDoubleMatrix(aLen, 1, mxREAL);
        eMat=mxCreateDoubleMatrix(aLen, 1, mxREAL);
        double *s=mxGetDoubles(sMat);
        double *e=mxGetDoubles(eMat);

        for(size_t k=0;k<aLen;k++) {
            exactPairSum<double>(a[k],b[k],s[k],e[k]);
        }
    } else {
        const size_t theLen=std::max(aLen,bLen);

        sMat=mxCreateDoubleMatrix(theLen, 1, mxREAL);
        eMat=mxCreateDoubleMatrix(theLen, 1, mxREAL);
        double *s=mxGetDoubles(sMat);
        double *e=mxGetDoubles(eMat);

        if(!bScal) {
            for(size_t k=0;k<theLen;k++) {
                 exactPairSum<double>(*a,b[k],s[k],e[k]);
            }
        } else {
            for(size_t k=0;k<theLen;k++) {
                 exactPairSum<double>(a[k],*b,s[k],e[k]);
            }
        }
    }

    plhs[0]=sMat;
    if(nlhs>1) {
        plhs[1]=eMat;
    } else {
        mxFree(eMat);
    }
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
