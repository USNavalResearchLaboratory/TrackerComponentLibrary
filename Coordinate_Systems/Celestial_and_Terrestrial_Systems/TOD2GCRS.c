/**TOD2GCRS Rotate a vector from the true equator and equinox of date
 *          coordinate system (TOD) using the IAU 2006/2000A model, where a
 *          location is known as an "apparent place", to the geocentric
 *          celestial reference system (GCRS). The transformation is
 *          performed by removing the precession, nutation, and frame bias.
 *
 *INPUTS:  xVec The 3XN matrix of N 3X1 Cartesian vectors that are to be
 *              rotated from the TOD into the GCRS coordinate system.
 * TT1, TT2 Jul1,Jul2 Two parts of a Julian date given in TT. The units of
 *              the date are days. The full date is the sum of both terms.
 *              The date is broken into two parts to provide more bits of
 *              precision. It does not matter how the date is split.
 *         dXdY dXdY=[dX;dY] are the celestial pole offsets with respect to
 *              the IAU 2006/2000A precession/nutation model in radians If
 *              this parameter is omitted or an empty matrix is passed, the
 *              value from the function getEOP will be used.
 *
 *OUTPUTS: xRot The 3XN matrix of the N 3X1 input vector rotated into the
 *              GCRS coordinate system.
 *       rotMat The 3X3 rotation matrix such that
 *              xRot(:,i)=rotMat*xVec(:,i).
 *
 *This uses functions in the International Astronomical Union's (IAU)
 *Standard's of Fundamental Astronomy (SOFA) library to obtain the product
 *of the nutation and precession rotation matrices  and the frame rotation
 *bias matrix. One goes from GCRS to TOD by applying a frame bias and then
 *precession and a nutation. Thus this function removes those rotations.
 *The rotations are discussed in the documentation for the SOFA library as
 *well as in among other sources.
 *
 *The correction for using dXdY is the most accurate formula in [2].
 *
 *The algorithm can be compiled for use in Matlab  using the 
 *CompileCLibraries function.
 *
 *The algorithm is run in Matlab using the command format
 *[xRot,rotMat]=TOD2GCRS(xVec,TT1,TT2);
 *or
 *[xRot,rotMat]=TOD2GCRS(xVec,TT1,TT2,dXdY);
 *
 *REFERENCES:
 *[1] G. Petit and B. Luzum, IERS Conventions (2010), International Earth
 *    Rotation and Reference Systems Service Std. 36, 2010.
 *[2] G. H. Kaplan, "Celestial pole offsets: Conversion from (dx,dy) to
 *    (dpsi,deps)," U.S. Naval Observatory, Tech. Rep., May 2005. [Online].
 *    Available: http://aa.usno.navy.mil/publications/reports/dXdY_to_dpsideps.pdf
 *
 *March 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

/*This header is required by Matlab.*/
#include "mex.h"
/*This header is for the SOFA library.*/
#include "sofa.h"
#include "MexValidation.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {    double *origVec;
    double TT1,TT2;
    double dX, dY;
    double rotMat[3][3];//To hold a rotation matrix product.
    size_t numItems,i;
    mxArray *retMATLAB;
    double *retVec;

    if(nrhs<3||nrhs>4) {
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

    //Get the celestial pole offsets
    if(nrhs>3&&!mxIsEmpty(prhs[3])) {//If they are provided
        size_t elInVec=mxGetM(prhs[3])*mxGetN(prhs[3]);
        double *dXdY;
        if(elInVec!=2) {
            mexErrMsgTxt("dXdY has the wrong dimensionality. It must contain two elements.");
            return;
        }
        checkRealDoubleArray(prhs[3]);
        dXdY=mxGetDoubles(prhs[3]);
        dX=dXdY[0];
        dY=dXdY[1];
    }else {//get the from the function getEOP, if theya re not provided.
        mxArray *retVals[2];
        double *dXdY;
        mxArray *JulUTCMATLAB[2];
        double JulUTC[2];
        int retVal;
        
        //Get the time in UTC to look up the parameters by going to TAI and
        //then UTC.
        retVal=iauTttai(TT1, TT2, &JulUTC[0], &JulUTC[1]);
        if(retVal!=0) {
            mexErrMsgTxt("An error occurred computing TAI.");
        }
        retVal=iauTaiutc(JulUTC[0], JulUTC[1], &JulUTC[0], &JulUTC[1]);
        switch(retVal){
            case 1:
                mexWarnMsgTxt("Dubious Date entered.");
                break;
            case -1:
                mexErrMsgTxt("Unacceptable date entered");
                break;
            default:
                break;
        }
        
        JulUTCMATLAB[0]=doubleMat2Matlab(&JulUTC[0],1,1);
        JulUTCMATLAB[1]=doubleMat2Matlab(&JulUTC[1],1,1);

        //Get the Earth orientation parameters for the given date.
        mexCallMATLAB(2,retVals,2,JulUTCMATLAB,"getEOP");
        mxDestroyArray(JulUTCMATLAB[0]);
        mxDestroyArray(JulUTCMATLAB[1]);
        
        checkRealDoubleArray(retVals[1]);
        if(mxGetM(retVals[1])!=2||mxGetN(retVals[1])!=1) {
            mxDestroyArray(retVals[0]);
            mxDestroyArray(retVals[1]);

            mexErrMsgTxt("Error using the getEOP function.");
            return;
        }
        
        dXdY=mxGetDoubles(retVals[1]);//The celestial pole offsets are not used.
        dX=dXdY[0];
        dY=dXdY[1];
        
        //Free the returned arrays.
        mxDestroyArray(retVals[0]);
        mxDestroyArray(retVals[1]);
    }
    
    {
        double dpsi,deps,epsa;
        double rb[3][3];
        double rp[3][3];
        double rn[3][3];
        double rbp[3][3];
        double rbpn[3][3];
        double XYZVec[3];
        double dZ;

        iauPn06a(TT1, TT2,
                  &dpsi, &deps, &epsa,
                  rb,//frame bias matrix
                  rp,//precession matrix
                  rbp,//bias-precession matrix
                  rn,//nutation matrix without dXdY correction.
                  rbpn);//GCRS-to-true matrix without dXdY correction

        //Now, we have to put the corrections for dXdY into the nutation matrix.
        //First, invert BPN by taking the Transpose. The result is B'P'N'.
        iauTr(rbpn, rbpn);
        //Next, get the pole coordinates by multiplying the inverted rbpn by
        //[0;0;1].
        XYZVec[0]=0;
        XYZVec[1]=0;
        XYZVec[2]=1;
        iauRxp(rbpn, XYZVec, XYZVec);
        //XYZVec now holds the pole coordinates X, Y, Z.
        dZ=-(XYZVec[0]/XYZVec[2])*dX-(XYZVec[1]/XYZVec[2])*dY;
        //Now multiply  P*B*[dX;dY;dZ] to get dX',dY'dZ'.
        XYZVec[0]=dX;
        XYZVec[1]=dY;
        XYZVec[2]=dZ;
        iauRxp(rbp, XYZVec, XYZVec);
        //Add in the correction terms
        dpsi+=XYZVec[0]/sin(epsa);
        deps+=XYZVec[1];
        //Use the corrected terms to get the full, corrected nutation matrix.
        iauPn06(TT1, TT2,
                dpsi, deps, &epsa,
                rb,//frame bias matrix
                rp,//precession matrix
                rbp,//bias-precession matrix
                rn,//nutation matrix with dXdY correction.
                rotMat);//GCRS-to-true matrix with dXdY correction
        
        iauTr(rotMat, rotMat);
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
