/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#ifndef MMOSPAALGS
#define MMOSPAALGS
#include "ShortestPathCPP.hpp"

void MMOSPAApproxCPP(double *MMOSPAEst,
                     size_t *orderList,
                     const double *x,
                     const double *w,
                     const size_t xDim,
                     const size_t numTar,
                     const size_t numHyp,
                     size_t numScan);
/*MMOSPAAPPROXCPP A C++ implementation of aan approximate minimum MOSPA
 *                optimization algorithm.
 *
 *INPUTS: MMOSPAEst  An array with xDim*numtar elements to hold the
 *                   approxiamte MMOSPA estimate.
 *        orderList  An array with numTar*numHyp element to hold the
 *                   ordering of the targets going into the approximate
 *                   MMOSPA estimate.
 *        x         An xDim * numTar * numHyp array that holds
 *                  numHyp hypotheses each consisting or numTar targets
 *                  (or generic vectors) with xDim dimensions per target
 *                  (per generic vector).
 *        w         A numHyp X 1 vector of the probabilities of each of
 *                  the numHyp hypotheses in x. The elements must all be
 *                  positive and sum to one.
 *        xDim      The size of the first dimension that x splits into.
 *        numTar    The size of the second dimensions that x splits into
 *        numHyp    The size of the third dimensions that x splits into.
 *        numScan   The number of forward scans of the approximation
 *                  algoirthm to perform.
 *
 *OUTPUT: The output is placed in MMOSPAEst and orderList and is
 *        The MMOSPA estimate and the ordering idices of the hypotheses.
 *
 * The algorithm as well as the concept of MOSPA error are described in
 * detail in 
 * D. F. Crouse, "Advances in displaying uncertain estimates of multiple
 * targets," in Proceedings of SPIE: Signal Processing, Sensor Fusion, and
 * Target Recognition XXII, vol. 8745, Baltimore, MD, Apr. 2013.
 *
 * November 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */ 

#endif

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
