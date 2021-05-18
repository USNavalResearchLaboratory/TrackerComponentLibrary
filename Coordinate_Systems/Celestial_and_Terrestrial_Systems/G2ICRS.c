/**G2ICRS Rotate vectors from the orientation of the International
 *        Astronomical Union's (IAU's) 1958 system of galactic coordinates
 *        to that of the International Celestial Reference System (ICRS).
 *        This function works with Cartesian vectors and with two spherical
 *        angles (azimuth and elevation) without range. The ICRS is aligned
 *        with the GCRS and BCRS.
 *
 *INPUTS: x The NXnumVec collection of vectors to convert. N can be 2, or
 *          3. If the vectors are 2D, then they are assumed to be azimuth
 *          and elevation in radians. 3D vectors are assumed to be
 *          Cartesian position.
 *
 *OUTPUTS: xICRS The vectors rotated into the galactic coordinates. If the
 *               input was 2D azimuth and elevation, the output will be the
 *               same. If the input was Cartesian, then the output will be
 *               Cartesian.
 *
 *This is mostly a wrapper for the function iauG2icrs in the International
 *Astronomical Union's (IAU) Standard's of Fundamental Astronomy library.
 *
 *The original coordinate system is defined in [1]. However, as documented
 *in the IAU's function iauG2icrs, the Galactic coordinate system has been
 *redefined in terms of more accurate modern data.
 *
 *The algorithm can be compiled for use in Matlab  using the 
 *CompileCLibraries function.
 *
 *The algorithm is run in Matlab using the command format
 *xICRS=G2ICRS(x,Jul1,Jul2);
 *
 *REFERENCES:
 *[1] A. Blaauw, C. S. Gum, J. L. Pawsey, and G. Westerhout, "The new
 *    I.A.U. system of galactic coordinates (1958 revision)," Monthly Notes
 *    of the Royal Astronomical Society, vol. 121, no. 2, pp. 123-131,
 *    1960.
 *
 *April 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

/*This header is required by Matlab.*/
#include "mex.h"
/*This header is for the SOFA library.*/
#include "sofa.h"
#include "MexValidation.h"
//For sqrt
#include <math.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    size_t curVec, numRow, numVec;
    double *xVec, *retData;
    mxArray *retMat;
            
    if(nrhs!=1){
        mexErrMsgTxt("Wrong number of inputs");
    }
    
    if(nlhs>1) {
        mexErrMsgTxt("Wrong number of outputs.");
    }
    
    checkRealDoubleArray(prhs[0]);
    
    numRow = mxGetM(prhs[0]);
    numVec = mxGetN(prhs[0]);
    
    if(!(numRow==2||numRow==3)) {
        mexErrMsgTxt("The input vector has a bad dimensionality.");
    }

    xVec=mxGetDoubles(prhs[0]);
    
    //Allocate space for the return vectors.
    retMat=mxCreateDoubleMatrix(numRow,numVec,mxREAL);
    retData=mxGetDoubles(retMat);
    
    if(numRow==3) {
        for(curVec=0;curVec<numVec;curVec++) {
            double curR, curAz, curEl;
            double *curInput=xVec+numRow*curVec;
            double *curOutput=retData+numRow*curVec;
       
            //If a 3D vector was given, convert it to spherical
            //coordinates.
            iauC2s(curInput,&curAz,&curEl);
            curR=sqrt(curInput[0]*curInput[0]+curInput[1]*curInput[1]+curInput[2]*curInput[2]);

            //Perform the rotation
            iauG2icrs (curAz, curEl, &curAz, &curEl);
        
            //Convert back to Cartesian coordinates
            iauS2c(curAz, curEl, curOutput);
            
            //Adjust for the magnitude of the vector.
            curOutput[0]*=curR;
            curOutput[1]*=curR;
            curOutput[2]*=curR;
        }
    } else {
        for(curVec=0;curVec<numVec;curVec++) {
            double curAz, curEl;
            double *curInput=xVec+numRow*curVec;
            double *curOutput=retData+numRow*curVec;
        
           //If a 2D vector was given, just get the azimuth and elevation.
            curAz=curInput[0];
            curEl=curInput[1];
            
            //Perform the rotation, saving the result in the return vector.
            iauG2icrs(curAz, curEl, curOutput, curOutput+1);
        }
    }
    //Set the return value.
    plhs[0]=retMat;
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
