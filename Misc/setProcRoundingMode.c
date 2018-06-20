/**SETPROCROUNDINGMODE Set the rounding mode used in the processor. This
 *                     specifies how results should be truncated to a
 *                     finite number of bits when performing floating point
 *                     arithmetic. Control of the rounding mode is useful
 *                     when implementing routines utilizing interval
 *                     algebra. Note that this does not appear to change
 *                     the rounding mode used in other mex files called
 *                     from Matlab. When running things in parallel, the 
 *                     rounding mode will have to be set for each parallel
 *                     process. For example, set it at the start of a
 *                     parfor loop, not outside the loop.
 *
 *INPUTS: roundMode An integer specifying the rounding mode to use.
 *                  Possible values are
 *                  0 Rounding is done towards negative infinity.
 *                  1 Rounding is done towards zero.
 *                  2 Rounding is done to the nearest value.
 *                  3 Rounding is done towards positive infinity.
 *
 *OUTPUTS: retVal This is zero if the rounding direction was successfully
 *                set and a nonzero value (corresponding to the outout of
 *                fesetround in C) otherwise.
 *
 *The rounding modes are standardized in [1]. The native Matlab version of
 *this function calls the undocumented system_dependent('setround',dir)
 *function in Matlab. The version that can be compiled in C just calls the
 *standard C function fesetround.
 *
 *The algorithm can be compiled for use in Matlab  using the 
 *CompileCLibraries function.
 *
 *The algorithm is run in Matlab using the command format
 *retVal=setProcRoundingMode(roundMode);
 *
 *EXAMPLE:
 *  setProcRoundingMode(0);
 *  1+eps(0)>1
 *  %The result should be false (0).
 *  setProcRoundingMode(3);
 *  1+eps(0)>1
 *  %The result should be true (1).
 *
 *REFERENCES:
 *[1] IEEE Standard 754-1985 for Binary Floating-point Arithmetic, IEEE,
 *    (1985).
 *
 *September 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

/*This header is required by Matlab.*/
#include "mex.h"
#if _MSC_VER
//This header is needed for _controlfp_s under Windows in Visual Studio
#include <float.h>
#else
/*This header is needed for fesetround*/
#include <fenv.h>
#endif
#include "MexValidation.h"
#pragma STDC FENV_ACCESS ON
#pragma fenv_access (on)

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
//Depending on the compiler, one of these pragmas will probably be
//necessary. The other should generate a warning about it not being
//supported. Without these pragmas, some compilers will use optimization
//routines that make this function not work.

    int roundMode;
    int retVal;
    if(nrhs>1||nrhs>2) {
        mexErrMsgTxt("Wrong number of inputs");
    }
    
    if(nlhs>1) {
        mexErrMsgTxt("Wrong number of outputs");   
    }
    
    roundMode=getIntFromMatlab(prhs[0]);

    //Most versions of Windows have their own hard-to-find method of
    //controlling the rounding direction of the processor rather than
    //actually supporting the C99 standard and providing fesetround.
    #if _MSC_VER
    {
        int err;
        unsigned int control_word;

        switch(roundMode) {
            case 0:
                err = _controlfp_s(&control_word, _RC_DOWN, _MCW_RC);
                break;
            case 1:
                err = _controlfp_s(&control_word, _RC_CHOP, _MCW_RC);
                break;
            case 2:
                err = _controlfp_s(&control_word, _RC_NEAR, _MCW_RC);
                break;
            case 3:
                err = _controlfp_s(&control_word, _RC_UP, _MCW_RC);
                break;
            default:
            mexErrMsgTxt("Unknown rounding mode specified.");
        }
    }
    #else
    switch(roundMode) {
        case 0:
            retVal=fesetround(FE_DOWNWARD);
            break;
        case 1:
            retVal=fesetround(FE_TOWARDZERO);
            break;
        case 2:
            retVal=fesetround(FE_TONEAREST);
            break;
        case 3:
            retVal=fesetround(FE_UPWARD);
            break;
        default:
            mexErrMsgTxt("Unknown rounding mode specified.");
    }
    #endif
    plhs[0]=intMat2MatlabDoubles(&retVal,1,1);
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
