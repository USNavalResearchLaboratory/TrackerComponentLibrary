/**GCRS2MOD Rotate a vector from the geocentric celestial reference system
 *          (GCRS) to the mean of date (MOD) coordinate system, which is
 *          the coordinate system using the mean equinox and
 *          ecliptic of date, IAU 2006/2000A model. The transformation is
 *          performed by adding the precession and frame bias.
 *
 *INPUTS: xVec The 3XN matrix of N 3X1 Cartesian vectors that are to be
 *             rotated from the GCRS coordinate system into the MOD
 *             coordinate system.
 * TT1, TT2 Jul1,Jul2 Two parts of a Julian date given in TT. The units
 *             of the date are days. The full date is the sum of both
 *             terms. The date is broken into two parts to provide more
 *             bits of precision. It does not matter how the date is
 *             split.
 *
 *OUTPUTS: xRot The 3XN matrix of the N 3X1 input vector rotated into the
 *              MOD coordinate system.
 *       rotMat The 3X3 rotation matrix such that
 *              xRot(:,i)=rotMat*xVec(:,i).
 *
 *This uses functions in the International Astronomical Union's (IAU)
 *Standard's of Fundamental Astronomy (SOFA) library to obtain the product
 *of the precession rotation matrix and the frame rotation bias matrix. One
 *goes from GCRS to mean of date by applying a frame bias and then
 *precession. Thus this function removes those rotations. The rotations are
 *discussed in the documentation for the SOFA library as well as in [1]
 *among other sources.
 *
 *The algorithm can be compiled for use in Matlab  using the 
 *CompileCLibraries function.
 *
 *The algorithm is run in Matlab using the command format
 *[xRot,rotMat]=GCRS2MOD(xVec,TT1,TT2);
 *
 *Different celestial coordinate systems are compared in [2].
 *
 *REFERENCES:
 *[1] G. Petit and B. Luzum, IERS Conventions (2010), International Earth
 *    Rotation and Reference Systems Service Std. 36, 2010.
 *[2] D. F. Crouse, "An Overview of Major Terrestrial, Celestial, and
 *    Temporal Coordinate Systems for Target Tracking," Formal Report,
 *    Naval Research Laboratory, no. NRL/FR/5344--16-10,279, 10 Aug. 2016,
 *    173 pages.
 *
 *March 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

/*This header is required by Matlab.*/
#include "mex.h"
/*This header is for the SOFA library.*/
#include "sofa.h"
/*This is for input validation*/
#include "MexValidation.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    double *origVec;
    double TT1,TT2;
    double rotMat[3][3];//To hold a rotation matrix product.
    size_t numItems,i;
    mxArray *retMATLAB;
    double *retVec;

    if(nrhs!=3) {
        mexErrMsgTxt("Incorrect number of inputs.");
        return;
    }
    
    if(nlhs>2) {
        mexErrMsgTxt("Wrong number of outputs.");
        return;
    }
    
    numItems=mxGetN(prhs[0]);
    
    if(mxGetM(prhs[0])!=3) {
       mexErrMsgTxt("xVec has the wrong dimensionality. It must be an 3XN matrix.");
       return;
    }
    checkRealDoubleArray(prhs[0]);
    origVec=mxGetDoubles(prhs[0]);
    
    TT1=getDoubleFromMatlab(prhs[1]);
    TT2=getDoubleFromMatlab(prhs[2]);

    {
    double dpsi,deps,epsa;
    double rb[3][3];
    double rp[3][3];
    double rn[3][3];
    double rbpn[3][3];
    
    iauPn06a(TT1, TT2,
              &dpsi, &deps, &epsa,
              rb,//frame bias matrix
              rp,//precession matrix
              rotMat,//bias-precession matrix
              rn,//nutation matrix
              rbpn);//GCRS-to-true matrix

    }
    
    retMATLAB=mxCreateDoubleMatrix(3,numItems,mxREAL);
    retVec=mxGetDoubles(retMATLAB);
    
    for(i=0;i<numItems;i++) {
        //Multiply the original vectors by the matrix to put it into the GCRS.
        iauRxp(rotMat, origVec+3*i, retVec+3*i);
    }
    
    plhs[0]=retMATLAB;

    if(nlhs>1) {
        double *retData;
        size_t j;
        plhs[1]=mxCreateDoubleMatrix(3,3,mxREAL);
        retData=mxGetDoubles(plhs[1]);
        for(i=0;i<3;i++) {
            for(j=0;j<3;j++) {
                retData[3*i+j]=rotMat[j][i];
            }
        }
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
