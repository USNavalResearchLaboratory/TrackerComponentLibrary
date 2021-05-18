/**ABERCORR Rotate vectors for stellar aberration, approximately
 *           including the effects of the gravitational potential of the
 *           Sun if the distances between the observers and the Sun are
 *           provided.
 *
 *INPUTS: vOrig A 3XN set of N vectors pointing from the observers to the
 *              objects being observed without corruption due to
 *              aberration. The units of the vector do not matter.
 *              Normally, one would expect vOrig to be composed of unit
 *              vectors as aberration affects the pointing direction of the
 *              vectors.
 *       obsVel The 3XN matrix of the velocity of each observer with
 *              respect to the stellar coordinate system origin in meters
 *              per second. For example, if using the barycentric celestial
 *              reference system (BCRS), they are the set of velocity
 *              vectors of the observers with respect to the barycenter.
 *      sunDist An optional 3XN vector of the scalar distances in meters
 *              from the Sun to the observer. If this parameter is
 *              included, an approximate correction due to the
 *              gravitational potential of the Sun is included.
 *
 *OUTPUTS: vAberr The 3XN set of N vectors rotated to deal with aberration
 *                effects.
 *
 *If all parameters are provided, the function is essentially a wrapper for
 *the function iauAb in the International Astronomical Union's Standards of
 *Fundamental Astronomy library with an adjustment so that the vector in
 *question need not have unit magnitude. If the distance to the Sun is
 *omitted, then the standard special relativistic aberration correction
 *without any gravitational effects is applied. The standard special
 *relativistic aberration correction is described in Chapter 7 of [1]
 *
 *Note that if the vectors provided are meant to be apparent distances,
 *rather than just unit vectors representing directions, the transformation
 *does not adjust the magnitudes of the vectors to account for special
 *relativistic contraction of space due to the motion of the observer with
 *respect to the origin. The magnitudes of the output vectors equal the
 *magnitudes of the input vectors.
 *
 *The algorithm can be compiled for use in Matlab  using the 
 *CompileCLibraries function.
 *
 *The algorithm is run in Matlab using the command format
 *vAberr=aberCorr(vOrig,obsVel,sunDist);
 *or
 *vAberr=aberCorr(vOrig,obsVel);
 *
 *REFERENCES:
 *[1] S. E. Urban and K. P.Seidelmann, Eds.,Explanatory Supplement to the
 *    Astronomical Almanac, 3rd ed. Mill Valley, CA: University Science
 *    Books, 2013.
 *
 *March 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
 **/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

/*This header is required by Matlab.*/
#include "mex.h"
/*This header is for the SOFA library.*/
#include "sofa.h"
//This is for the velocity addition for the correction without
//gravitational effects.
#include "relFuncs.h"
#include "MexValidation.h"
//For the sqrt function.
#include <math.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    double *vOrig;
    double *obsVel;
    mxArray *retMATLAB;
    double c, *retVec;
    size_t numVec;
    
    if(nrhs<2||nrhs>3) {
        mexErrMsgTxt("Incorrect number of inputs.");
        return;
    }
    
    if(nlhs>1) {
        mexErrMsgTxt("Wrong number of outputs.");
        return;
    }

    numVec=mxGetN(prhs[0]);
    if(numVec!=mxGetN(prhs[1])||mxGetM(prhs[0])!=3||mxGetM(prhs[1])!=3) {
        mexErrMsgTxt("The input vectors have the wrong dimensionality.");
        return;
    }
    
    checkRealDoubleArray(prhs[0]);
    vOrig=mxGetDoubles(prhs[0]);
    checkRealDoubleArray(prhs[1]);
    obsVel=mxGetDoubles(prhs[1]);
    
    c=getScalarMatlabClassConst("Constants","speedOfLight");
    
    //Allocate space for the return values.
    retMATLAB=mxCreateDoubleMatrix(3,numVec,mxREAL);
    retVec=mxGetDoubles(retMATLAB);
    
    //If the third parameter is provided, then use the algorithm from the IAU.
    if(nrhs>2) {
        double AU, *sunDist;
        size_t curVec;
        
        //Needed to convert units.
        AU=getScalarMatlabClassConst("Constants","AstronomialUnit");
        
        //If the dimensionality is wrong.
        if(!((mxGetM(prhs[2])==1&&mxGetN(prhs[2])==numVec)||(mxGetM(prhs[2])==numVec&&mxGetN(prhs[2])==1))) {
            mxDestroyArray(retMATLAB);
            mexErrMsgTxt("The input vectors have the wrong dimensionality.");
            return;
        }
        checkRealDoubleArray(prhs[0]);
        sunDist=mxGetDoubles(prhs[2]);

        for(curVec=0;curVec<numVec;curVec++) {
            double vecMag,unitVec[3];
            double v[3];
            double s;
            double bm1;
            
            //Get a unit direction vector and magnitude.
            iauPn(vOrig+3*curVec, &vecMag, unitVec);
            
            //Convert the velocity to units of the speed of light.
            v[0]=obsVel[3*curVec]/c;
            v[1]=obsVel[3*curVec+1]/c;
            v[2]=obsVel[3*curVec+2]/c;
            
            //The distance to the sun in AU.
            s=sunDist[curVec]/AU;
            //The reciprocal of the Lorentz factor.
            bm1=sqrt(1-(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]));
            
            //Perform the correction with the IAU's function.
            iauAb(unitVec, v, s, bm1,retVec+3*curVec);
            
            //Set the magnitude back to what it was.
            retVec[3*curVec]*=vecMag;
            retVec[3*curVec+1]*=vecMag;
            retVec[3*curVec+2]*=vecMag;
        }
    } else {
    //If the distance to the sun is not given, then just perform a normal
    //special relativistic correction.
        size_t curVec;
        
        for(curVec=0;curVec<numVec;curVec++) {
            double vecMag,lightVec[3];
            
            //Get a unit direction vector and magnitude.
            iauPn(vOrig+3*curVec, &vecMag, lightVec);
            
            //The light is in direction lightVec with speed c
            lightVec[0]*=c;
            lightVec[1]*=c;
            lightVec[2]*=c;
            
            //The light travels at speed c with true direction uPosition.
            //Add the velocity vector of the observer to that of light.
            relVecAddC(c,obsVel+3*curVec,lightVec,retVec+3*curVec);
            //Because one vector had a magnitude of c, the returned vector
            //must have the same magnitude. Restore the previous magnitude.
            retVec[3*curVec]*=vecMag/c;
            retVec[3*curVec+1]*=vecMag/c;
            retVec[3*curVec+2]*=vecMag/c;
        }
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
