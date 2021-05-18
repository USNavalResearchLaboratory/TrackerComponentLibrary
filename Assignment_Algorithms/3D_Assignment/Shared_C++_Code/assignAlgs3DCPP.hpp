/**ASSIGNALGS3DCPP This is a header for C++ language functions to solve the
 *              axial operations research 3D assignment problem.
 * 
 *Better understanding of the algorithms can usually be obtained from
 *looking at the Matlab implementations.
 *
 *November 2020 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#ifndef ASSIGNALGS3DCPP
#define ASSIGNALGS3DCPP

//Defines the size_t and ptrdiff_t types
#include <cstddef>

/**ASSIGN3DBBBUFFERSIZE Obtain the size (in bytes) of the temporary buffer
 *                      needed by the assign3DBBCPP function.
 *
 *INPUTS: nDims A length-3 vector holding the size of each dimension of the
 *              3D cost matrix.
 *    boundType This selects the bound to use for the branch and bound
 *              method. These correspond to the method input of the
 *              assign3DLB function.
 *   initMethod A parameter indicating how the branch-and-bound algorithm
 *              should be initialized. Possible values are:
 *            0 Do not use an initial estimate.
 *            1 Use an initial solution (and lower bound) via algorithm 0
 *              of the assign3D function.
 *
 *October 2020 David F. Crouse, Naval Research Laboratory, Washington D.C.
 **/
size_t assign3DBBBufferSize(const size_t *nDims,const int boundType,const int initMethod);
/**ASSIGN3DBBCPP Solve the axial operations research 3D assign problem
 *               using a branch-and-bound algorithm.
 *
 *INPUTS: minCostTuples A length 3*nDims[0] array holding space for each of
 *               the nDims[0] 3D tuples. The assigned tuplesa re put in
 *               here.
 *         COrig The 3D cost matrix. The size of each dimensions is in
 *               nDims. The is a nDIms[0]*nDims[1]*nDims[2] matrix. Values
 *               are stored, by row, then column, then third dimension.
 *         nDims A length 3 vector where nDims[i] is the size of the ith
 *               dimension of COrig.
 *      maximize If true, the minimization problem is transformed into a
 *               maximization problem.
 *     boundType This selects the bound to use for the branch and bound
 *               method. These correspond to the method input of the
 *               assign3DLB function.
 *    initMethod A parameter indicating how the branch-and-bound algorithm
 *               should be initialized. Possible values are:
 *               0 Do not use an initial estimate.
 *               1 Use an initial solution (and lower bound) via algorithm
 *                 0 of the assign3D function.
 *       maxIter If initMethod=1, then this is the maximum number of
 *               iterations to allow the initialization routine in assign3D
 *               to perform. If initMethod~=1, then this parameter is
 *               ignored.
 *        epsVal If initMethod=1, then this is both the AbsTold and RelTol
 *               inputs to assign3D to determine whether the function has
 *               converged in terms of the relative duality gap. If
 *               initMethod~=1, then this parameter is ignored.
 *    tempBuffer A pointer to a buffer with sufficient space to hold at
 *               least the number of bytes given by the function
 *               assign3DBBBufferSize.
 *
 *OUTPUTS: Tuples are placed in minCostTuples. The return value is the cost
 *         value of the optimal assignment found.
 *
 *October 2020 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
double assign3DBBCPP(ptrdiff_t *minCostTuples, const double *const COrig,const size_t *const nDims,const int maximize,const int boundType,const int initMethod,const size_t maxIter,const double epsVal, void *tempBuffer);

#endif

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
*OF RECIPIENT IN THE USE OF THE SOFTWARE.
 */
