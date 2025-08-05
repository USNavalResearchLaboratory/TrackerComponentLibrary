/**GETCOMPRESSEDIDXMAPPINGS  A C++ implementation of a function that maps a
*      set of indices with gaps to a set of indices with no gaps. Given a
*       matrix holding a set of hypotheses of assigned measurements at each
*       dimension of an S-dimensional assignment, obtain a mapping of the
*       indices of the measurements at each dimension to one that has no
*       gaps. For example, if in the first dimension, hypotheses contain
*       measurements 1, 3, 4, and 5, these would be remapped to indices 1,
*       2, 3, and 4 (the gap between 1 and 3 is eliminated). At the same
*       time, create a matrix to perform the inverse mapping from the
*       compressed indices to the spread out indices. The zero index is
*       assumed to be the missed detection hypothesis and thus is not
*       considered for remapping (it implicitly does not change). Also,
*       this function can handle the case of multiple groups of hypotheses
*       with regions between them that should be ignored. See the comments
*       in the Matlab implementation of this file for more information on
*       the inputs to the function and how it works.
*
* The algorithm can be compiled for use in Matlab  using the 
* CompileCLibraries function.
*
* The algorithm is run in Matlab using the command format
* [numMeasPerScan,measMap,invMeasMap]=getCompressedIdxMappings(theHypIdxHist,numSec,numPartsPerSec,secLength)
*
*January 2025 David F. Crouse, Naval Research Laboratory, Washington D.C.
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
#include <algorithm>
#include <numeric>
/*
 *INPUTS: numDims The size of the first dimension of theHypIdxHist.
 *     numHyps The size of the second dimension of theHypIdxHist.
 * theHypIdxHist A numDimsXnumHyps matrix of indices. Negative indices
 *             are skipped by the algorithm.
 * numSec, numPartsPerSec, secLength If all of the hypotheses are
 *            continugous, then set these to 0, NULL, 0. Otherwise, the
 *            data is given in numSec sections of length secLength, where
 *            the ith section contains numPartsPerSec(i) hypotheses.
 * maxHypsPerScan The largest value in theHypIdxHist.
 *
 *OUTPUTS: maxMeasPerScan The maximum value in the numMeasPerScan output.
 *         numMeasPerScan A length numDims array where numMeasPerScan[i] is
 *            the number of unique indices in theHypIdxHist(i,:). Negative
 *            indices present in theHypIdxHist do not count.
 *    measMap One must pass enough memory to hold maxHypsPerScanXnumDims
 *            matrix. However, what is filled in is actually a
 *            maxMeasPerScanXnumDims matrix of how compressed
 *            indices map to the global indices in theHypIdxHist.
 *            measMap[i,j] is how the ith compressed index in the jth
 *            dimension maps to a measurement index in jth dimension of
 *            theHypIdxHist.
 * invMeasMap This is the inverse map of measMap. This is a
 *            maxMeasPerScanXnumDims matrix such that if one takes an
 *            index idx from the jth dimension of theHypIdxHist, it
 *            maps to the compressed indices as invMeasMap(idx,j)
 *            (think [row,column]), assuming that idx>=0.
 *
 *January 2025 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
void getCompressedIdxMappings(const size_t numDims,const size_t maxHyps,const ptrdiff_t *theHypIdxHist,const size_t numSec,const size_t *numPartsPerSec,const size_t secLength,const size_t maxHypsPerScan, size_t &maxMeasPerScan, size_t *numMeasPerScan, size_t *measMap,size_t *invMeasMap) {
    //Zero the elements of invMeasMap. It will initially be treated as a
    //boolean array and then the actual indices will get filled in.
    std::fill_n(invMeasMap,maxHypsPerScan*numDims,static_cast<size_t>(0));

    if(numPartsPerSec==NULL) {
        for(size_t curHyp=0;curHyp<maxHyps;curHyp++) {
            const size_t offset=curHyp*numDims;
            for(size_t curDim=0;curDim<numDims;curDim++) {
                const ptrdiff_t idx=theHypIdxHist[offset+curDim];
    
                if(idx>=0) {
                    invMeasMap[maxHypsPerScan*curDim+idx]=1;
                }
            }
        }
    } else {
        for(size_t curSec=0;curSec<numSec;curSec++) {
            size_t startIdx=curSec*secLength;
            for(size_t hypIdx=startIdx;hypIdx<startIdx+numPartsPerSec[curSec];hypIdx++) {
                const size_t hypOffset=hypIdx*numDims;
                for(size_t curDim=0;curDim<numDims;curDim++) {
                    const ptrdiff_t idx=theHypIdxHist[hypOffset+curDim];
                    if(idx>=0) {
                        invMeasMap[maxHypsPerScan*curDim+idx]=1;
                    }
                }
            }
        }
    }

    //invMeasMap is now essentially a boolean matrix marking which
    //measurement indices are present at each scan. Sum over the first
    //dimension to fill in numMeasPerScan.
    maxMeasPerScan=0;
    for(size_t curDim=0;curDim<numDims;curDim++) {
        size_t offset=curDim*maxHypsPerScan;
        numMeasPerScan[curDim]=std::accumulate(invMeasMap+offset, invMeasMap+offset+maxHypsPerScan,static_cast<size_t>(0));
        maxMeasPerScan=std::max(maxMeasPerScan,numMeasPerScan[curDim]);
    }

    //Fill in a maxMeasPerScanXnumDims matrix measMat that will give the
    //index of the ith measurement from the specified scan. At the same
    //time, replace the nonnegative indices in invMeasMat with the
    //corresponding indices in measMat. This lets one map from the
    //larger set of measurement values at each scan in theHypIdxHist to
    //a compressed set of values (no gaps).
    for(size_t curDim=0;curDim<numDims;curDim++) {
        const size_t numMeasCur=numMeasPerScan[curDim];
        const size_t dimOffsetI=curDim*maxHypsPerScan;
        const size_t dimOffsetM=curDim*maxMeasPerScan;
        if(numMeasCur>0) {
            //Fill in the first numMeasCur items in the curDimth column of
            //measMap with the indices of the nonzero values in invMeasMap.
            size_t numFound=0;
            for(size_t k=0;k<maxHypsPerScan;k++) {
                if(invMeasMap[k+dimOffsetI]) {
                    measMap[numFound+dimOffsetM]=k;
                    numFound++;
                    if(numFound==numMeasCur) {
                        break;
                    }
                }
            }

            //Replace the boolean values with actual indices.
            for(size_t curMeas=0;curMeas<numMeasCur;curMeas++) {
                const size_t idx=measMap[curMeas+dimOffsetM];
                invMeasMap[dimOffsetI+idx]=curMeas;
            }
        }
    }
}
 
void mexFunction(const int nlhs, mxArray *plhs[], const int nrhs, const mxArray *prhs[]) { 
    size_t numSec;
    size_t *numPartsPerSec=NULL;
    size_t secLength;
    size_t *invMeasMap, *numMeasPerScan, *measMap;
    size_t maxMeasPerScan=0;
    size_t maxHypsPerScan=0;

    if(nrhs<1||(nrhs!=1&&nrhs!=4)) {
        mexErrMsgTxt("Incorrect number of inputs.");
        return;
    }

    if(nlhs>3) {
        mexErrMsgTxt("Incorrect number of outputs.");
        return;  
    }

    checkRealArray(prhs[0]);

    size_t numDims;
    size_t maxHyps;
    //Copying it as size_t makes it check for being <0. Once copied, we
    //have to subtract 1 from all of the elements to put the indexation
    //into C++ format.
    ptrdiff_t *theHypIdxHist=reinterpret_cast<ptrdiff_t*>(copySizeTMatrixFromMatlab(prhs[0],&numDims,&maxHyps));
    const size_t totalNumElements=numDims*maxHyps;

    //totalNumElements should equal numDims*maxHyps.
    for(size_t k=0;k<totalNumElements;k++) {
        theHypIdxHist[k]--;//Convert to C++ indices.
    }

    const bool hasSkipSections=nrhs>1;

    if(hasSkipSections) {
        numSec=getSizeTFromMatlab(prhs[1]);
        if(numSec>maxHyps) {
            mexErrMsgTxt("numSec cannot be larger than maxHyps.");
        }

        if(mxGetNumberOfElements(prhs[2])<numSec) {
            mexErrMsgTxt("numPartsPerSec must be at least numSec in length (extra elements are ignored).");
        }

        secLength=getSizeTFromMatlab(prhs[3]);

        if(numSec*secLength>maxHyps) {
            mexErrMsgTxt("numSec, secLength and the size of theHypIdxHist are inconsistent.");
        }
        
        size_t numElsCopied;
        numPartsPerSec=copySizeTArrayFromMatlab(prhs[2],&numElsCopied);
        

        for(size_t k=0;k<numSec;k++) {
            if(numPartsPerSec[k]>secLength) {
                 mexErrMsgTxt("An entry in numPartsPerSec is too long.");
            }
        }

        maxHypsPerScan=0;
        for(size_t curSec=0;curSec<numSec;curSec++) {
            const size_t startIdx=curSec*secLength;
            
            for(size_t hypIdx=startIdx;hypIdx<(startIdx+numPartsPerSec[curSec]);hypIdx++){
                const ptrdiff_t numHypsCur=1+(*std::max_element(theHypIdxHist+hypIdx*numDims,theHypIdxHist+(hypIdx+1)*numDims));
                if(numHypsCur>0) {
                    maxHypsPerScan=std::max(maxHypsPerScan,static_cast<size_t>(numHypsCur));
                }
            }
        }

        //Allocate space for the return variables.
        invMeasMap=new size_t[maxHypsPerScan*numDims];
        numMeasPerScan=new size_t[numDims];
        measMap=new size_t[maxHypsPerScan*numDims];

        getCompressedIdxMappings(numDims, maxHyps, theHypIdxHist, numSec, numPartsPerSec, secLength, maxHypsPerScan, maxMeasPerScan, numMeasPerScan, measMap,invMeasMap);
        
        mxFree(numPartsPerSec);
    } else {
        ptrdiff_t maxElVal=1+(*std::max_element(theHypIdxHist,theHypIdxHist+totalNumElements));
        if(maxElVal>0) {
            maxHypsPerScan=static_cast<size_t>(maxElVal);
        }

        //Allocate space for the return variables.
        invMeasMap=new size_t[maxHypsPerScan*numDims];
        numMeasPerScan=new size_t[numDims];
        measMap=new size_t[maxHypsPerScan*numDims];

        getCompressedIdxMappings(numDims,maxHyps,theHypIdxHist,0,NULL,0,maxHypsPerScan, maxMeasPerScan, numMeasPerScan, measMap, invMeasMap);
    }
    mxFree(theHypIdxHist);
    //Copy the values to Matlab to be returned.
    plhs[0]=sizeTMat2MatlabDoubles(numMeasPerScan,1,numDims);
    delete[] numMeasPerScan;
    if(nlhs>1) {
        //The indices inside of measMap must all be incremented by 1 to
        //work in Matlab.
        for(size_t k=0;k<maxMeasPerScan*numDims;k++) {
            measMap[k]++;
        }

        plhs[1]=sizeTMat2MatlabDoubles(measMap,maxMeasPerScan,numDims);
        delete[] measMap;
        if(nlhs>2) {
            //The indices inside of invMeasMapMat must all be incremented
            //by 1 to work in Matlab.
            for(size_t k=0;k<maxHypsPerScan*numDims;k++) {
                invMeasMap[k]++;
            }

            plhs[2]=sizeTMat2MatlabDoubles(invMeasMap,maxHypsPerScan,numDims);
        }
        delete[] invMeasMap;
    } else{
        delete[] measMap;
        delete[] invMeasMap;
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
