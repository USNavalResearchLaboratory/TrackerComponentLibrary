/**ICRS2ECLIPTIC Convert a location vector from the International
 *               Celestial Reference System (ICRS) to eliptic coordinates
 *               either using the IAU 2006 precession model or the Vondrak
 *               400 millennia precession model.
 *
 *INPUTS: x The NXnumVec collection of vectors in the ICRS to convert. N
 *          can be 2, or 3. If the vectors are 2D, then they are assumed to
 *          be azimuth and elevation in radians. 3D vectors are assumed to
 *          be Cartesian position.
 * Jul1, Jul2 Two parts of a Julian date given in terrestrial time (TT).
 *          The units of the date are days. The full date is the sum of
 *          both terms. The date is broken into two parts to provide more
 *          bits of precision. It does not matter how the date is split.
 *   method An optional parameter specifying which algorithm is to be used.
 *          Possible values are
 *          0 (The default if omitted or an empty matrix is passed) Use the
 *            IAU 2006 precession model.
 *          1 Use the long-term (Vondrak) precession model.
 *
 *OUTPUTS: xG The vectors rotated into the ecliptic coordinate system. If
 *            the input was 2D azimuth and elevation, the output will be
 *            the same. If the input was Cartesian, then the output will be
 *            Cartesian.
 *
 *This function is a Matlab interface for the relevant functions in the
 *International Astronomical Union's (IAU) Standard's of Fundamental
 *Astronomy library.
 *
 *The ecliptic is defined in the IERS Conventions [1] to be the "the
 *plane perpendicular to the mean heliocentric orbital angular momentum
 *vector of the Earth-Moon barycentre in the BCRS".
 *
 *The algorithm can be compiled for use in Matlab  using the
 *CompileCLibraries function.
 *
 *The algorithm is run in Matlab using the command format
 *eclipVec=ICRS2Eliptic(vec,TT1,TT2,method);
 *
 *REFERENCES:
 *[1] G. Petit and B. Luzum, IERS Conventions (2010), International Earth
 *    Rotation and Reference Systems Service Std. 36, 2010.
 *
 *July 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

 /*This header is required by Matlab.*/
#include "mex.h"
/*This header is for the SOFA library.*/
#include "sofa.h"
#include "MexValidation.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    size_t curVec, numRow, numVec;
    double *xVec, *retData;
    double TT1, TT2, julEpoch;
    mxArray *retMATLAB;
    size_t method=0;
    
    if(nrhs<3||nrhs>4){
        mexErrMsgTxt("Wrong number of inputs");
        return;
    }
    
    if(nlhs>1) {
        mexErrMsgTxt("Wrong number of outputs");
        return;
    }
    
    numVec=mxGetN(prhs[0]);
    numRow=mxGetM(prhs[0]);
    if(!(numRow==2||numRow==3)) {
       mexErrMsgTxt("vec has a bad dimensionality.");
       return;
    }
    
    checkRealDoubleArray(prhs[0]);
    xVec=mxGetDoubles(prhs[0]);
    
    TT1=getDoubleFromMatlab(prhs[1]);
    TT2=getDoubleFromMatlab(prhs[2]);
    
    if(nrhs>3) {
       method=getSizeTFromMatlab(prhs[3]);
    }
    
    //Allocate the return vector
    retMATLAB=mxCreateDoubleMatrix(numRow,numVec,mxREAL);
    retData=mxGetDoubles(retMATLAB);
    
    switch(method) {
        case 0://Use the IAU 2006 model
            if(numRow==3){
            	for(curVec=0;curVec<numVec;curVec++) {
                    double curR, curAz, curEl;
                    double *curInput=xVec+numRow*curVec;
                    double *curOutput=retData+numRow*curVec;

                    //If a 3D vector was given, convert it to spherical
                    //coordinates.
                    iauC2s(curInput,&curAz,&curEl);
                    curR=sqrt(curInput[0]*curInput[0]+curInput[1]*curInput[1]+curInput[2]*curInput[2]);
                    
                    //Perform the conversion
                    iauEqec06(TT1,TT2,curAz,curEl,&curAz,&curEl);

                    //Convert back to Cartesian coordinates
                    iauS2c(curAz, curEl, curOutput);
                    
                    //Adjust for the magnitude of the vector.
                    curOutput[0]*=curR;
                    curOutput[1]*=curR;
                    curOutput[2]*=curR;
                }
            }else {
                for(curVec=0;curVec<numVec;curVec++) {
                    double curAz, curEl;
                    double *curInput=xVec+numRow*curVec;
                    double *curOutput=retData+numRow*curVec;

                   //If a 2D vector was given, just get the azimuth and elevation.
                    curAz=curInput[0];
                    curEl=curInput[1];

                    //Perform the rotation, saving the result in the return vector.
                    iauEqec06(TT1,TT2,curAz, curEl, curOutput, curOutput+1);
                }
            }
            break;
        case 1:
            julEpoch=iauEpj(TT1, TT2);
            if(numRow==3){
            	for(curVec=0;curVec<numVec;curVec++) {
                    double curR, curAz, curEl;
                    double *curInput=xVec+numRow*curVec;
                    double *curOutput=retData+numRow*curVec;

                    //If a 3D vector was given, convert it to spherical
                    //coordinates.
                    iauC2s(curInput,&curAz,&curEl);
                    curR=sqrt(curInput[0]*curInput[0]+curInput[1]*curInput[1]+curInput[2]*curInput[2]);
                    
                    //Perform the conversion
                    iauLteqec(julEpoch,curAz,curEl,&curAz,&curEl);

                    //Convert back to Cartesian coordinates
                    iauS2c(curAz, curEl, curOutput);
                    
                    //Adjust for the magnitude of the vector.
                    curOutput[0]*=curR;
                    curOutput[1]*=curR;
                    curOutput[2]*=curR;
                }
            }else {
                for(curVec=0;curVec<numVec;curVec++) {
                    double curAz, curEl;
                    double *curInput=xVec+numRow*curVec;
                    double *curOutput=retData+numRow*curVec;

                   //If a 2D vector was given, just get the azimuth and elevation.
                    curAz=curInput[0];
                    curEl=curInput[1];

                    //Perform the rotation, saving the result in the return vector.
                    iauLteqec(julEpoch,curAz, curEl, curOutput, curOutput+1);
                }
            }
            break;
        default:
            mexErrMsgTxt("Unknown method specified.");
    }
    //Set the return value.
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
