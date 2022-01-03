/**ASSIGN3DC This file contains C language functions implementing a
 *          Lagrangian relaxation-based approximate solution to the axial
 *          operations research 3D assignment problem. See the comments to
 *          the Matlab implementation of assign3D for more details on the
 *          algorithm.
 *
 *March 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

//The header for the implementations here.
#include "assignAlgs3D.h"

//The header for the 2D assignment algorithms in C.
#include "assignAlgs2D.h"

//For some basic operations on matrices that are simple to perform in
//Matlab, but can be a bit tedious in C.
#include "basicMatOps.h"

//For sqrt and pow and INFINITY.
#include <math.h>

//For uint8_t and other types.
#include <stdint.h>

//For memset.
#include <string.h>

/*If a compiler does not support INFINITY in C99, then it must be
 * explicitly defined.*/
#ifndef INFINITY
static const uint64_t infVal=0x7ff0000000000000;
#define INFINITY (*(double*)(&infVal))
#endif

//BLAS routines for basic algebra.
#include "blas.h"

//For isfinite
#include <math.h>

static ptrdiff_t assign3DCBasic(ptrdiff_t * tuples,double * fStar, double *qStar,double *u, void *tempSpace,const size_t *nDims,const double *C,const int subgradMethod,const size_t maxIter,const double AbsTol,const double RelTol,const double param1,const double param2,const size_t param3);
static size_t updateBestFeasSol3DCBufferSize(const size_t *nDims);
static ptrdiff_t updateBestFeasSol3DC(ptrdiff_t *tuples, double * fStar,void *tempBuffer2DAssign,const size_t *nDims, const double *C, const ptrdiff_t *gamma1, const double qStar, const double AbsTol, const double RelTol);

static size_t updateBestFeasSol3DCBufferSize(const size_t *nDims) {
 /**UPDATEBESTFEASSOL3DCBUFFERSIZE Given the dimensions of the assignment
 *     matrix, return the minimum size of the input tempBuffer needed (in
 *     bytes) for the updateBestFeasSol3DC function.
 *
 *February 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */   
    const size_t n1=nDims[0];
    const size_t n3=nDims[2];
    
    return assign2DCBufferSize(n3,n1)+(n1+n3)*sizeof(ptrdiff_t)+(n1+n3+n1*n3)*sizeof(double);
}

size_t assign3DCBufferSize(const size_t *nDims,const int subgradMethod) {
/**ASSIGN3DCBUFFERSIZE Given the dimensions of the assignment matrix,
 *      return the minimum size of the input tempBuffer needed (in bytes)
 *      for the assign3DC and assign3DCBasic algorithms.
 *
 *February 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
    const size_t n1=nDims[0];
    const size_t n2=nDims[1];
    const size_t n3=nDims[2];
    size_t buffSize=sizeof(ptrdiff_t)*(n1*n2+n1)+sizeof(double)*(n3+n1*n2)+updateBestFeasSol3DCBufferSize(nDims);

    if(subgradMethod>=3) {//If it is using space dilation.
        buffSize+=sizeof(double)*(3*n3+3*n3*n3);
    }

    return buffSize;
}

ptrdiff_t assign3DC(ptrdiff_t * tuples,double * fStar, double *qStar,double *u,void *tempSpace,const size_t *nDims,double *C,bool maximize,const int subgradMethod,const size_t maxIter,const double AbsTol,const double RelTol,const double param1,const double param2,const size_t param3) {
/**ASSIGN3DC Approximate the solution to the operations research axial 3D
 *         assignment problem using a dual-primal Lagrangian relaxation
 *         algorithm. This function is the same as assign3DCBasic except
 *         it has a maximize input and the C matrix is modified if
 *         maximize=true.
 *
 *See assign3DCBasic for more comments as to the inputs of this function.
 *
 *February 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
const size_t n1=nDims[0];
const size_t n2=nDims[1];
const size_t n3=nDims[2];
const size_t numEls=n1*n2*n3;
 
//The scalar special case.
if(numEls==1) {
    tuples[0]=0;
    tuples[1]=0;
    tuples[2]=0;
    
    *fStar=C[0];
    *qStar=C[0];
    u[0]=0;
    return 0;
}

if(maximize==true) {
    size_t i;
    ptrdiff_t retVal;
    
    for(i=0;i<numEls;i++) {
        C[i]=-C[i];
    }
    
    retVal=assign3DCBasic(tuples,fStar,qStar,u,tempSpace,nDims,C,subgradMethod,maxIter,AbsTol,RelTol,param1,param2,param3);
    
    //Account for the difference between the minimization and maximization
    //problems.
    *fStar=-*fStar;
    *qStar=-*qStar;
    
    for(i=0;i<n3;i++) {
        u[i]=-u[i];
    }
    
    return retVal;
} else {
    return assign3DCBasic(tuples,fStar,qStar,u,tempSpace,nDims,C,subgradMethod,maxIter,AbsTol,RelTol,param1,param2,param3);
}
}

static ptrdiff_t assign3DCBasic(ptrdiff_t *tuples,double *fStar, double *qStar,double *u, void *tempSpace,const size_t *nDims,const double *C,const int subgradMethod,const size_t maxIter,const double AbsTol,const double RelTol,const double param1,const double param2,const size_t param3) {
/**ASSIGN3DBASIC Approximate the solution to the minimization-only
 *        operations research axial 3D assignment problem using a dual-
 *        primal Lagrangian relaxation  technique. The optimization problem
 *        being solved is
 *        minimize
 *        sum_{i=1}^{n1}sum_{j=1}^{n2}sum_{k=1}^{n3}C_{i,j,k}*rho_{i,j,k}
 *        subject to
 *        sum_{i=1}^{n1}sum_{j=1}^{n2}rho_{i,j,k}<=1 for all k
 *        sum_{i=1}^{n1}sum_{k=1}^{n3}rho_{i,j,k}<=1 for all j
 *        sum_{j=1}^{n2}sum_{k=1}^{n3}rho_{i,j,k} =1 for all i
 *        rho_{i,j,k} = 0 or 1
 *        assuming that n1<=n2<=n3, and C is an n1Xn2Xn3 cost matrix.
 *
 *INPUTS: tuples A pointer to an array of ptrdiff_t values to hold the 3Xn1
 *         matrix of assigned tuples.
 *   fStar A pointer to a double that will hold the best (lowest) primal
 *         solution value found by this function.
 *   qStar A pointer to a double that will hold the best (highest)
 *         dual solution value found by this function.
 *       u A pointer to an array of n3 doubles to hold the dual values.
 * tempSpace A pointer to a buffer that is used as temporary space by
 *         this function. The buffer must be at least
 *         assign3DCBufferSize(nDims,subgradMethod) bytes in size.
 *   nDims A length 3 vector of ptrdiff_t values holding the dimensions of
 *         C. All values should be >=1.
 *       C The original n1Xn2Xn3 3D cost matrix of doubles with
 *         n1>=n2>=n3. Elements are stored by column as is standard in
 *         Matlab and Fortran.
 * subgradMethod A parameter indicating the subgradient optimization
 *         algorithm to use. Possible values are:
 *         0 Polyak's Method from [2]. (The default if omitted or an empty
 *           matrix is passed) In this method, for minimization, the dual
 *           variables are updated as uNew=u+gamma*(fStar-q)/norm(g)^2*g,
 *           where g is the subgradient, q is the current dual cost, and g
 *           is the subgradient vector. gamma is a design parameter.
 *         1 Bragin's Method from [3]. The dual update for minimization is
 *           uNew=u+alpha*g
 *           where during the first step, alpha=(fStar-q)/norm(g)^2 and for
 *           subsequent steps it is
 *           alpha=(1-1/(M*k^(1-1/k^r)))*alphaPrev*norm(gPrev)/norm(g)
 *           where alphaPrev is alpha for the previous step, k is the step
 *           number (starting at 0), gPrev is the subgradient before the
 *           current step and M and r are design parameters.
 *         2 Bertsekas' Heuristic Method from [4]. The update is
 *           uNew=u+alpha*g
 *           where alpha=((1+a/beta^b)*qMax-q)/norm(g)^2 where qMax is the
 *           highest dual cost value encountered thus far (when
 *           minimizing), a and b are design parameters, and beta is a
 *           value that increases or decreases in the algorithm depending
 *           on whether or not there was an improvement in the best dual
 *           cost found after the last step.
 *         3 Shor's Space Dilation Algorithm from [5].
 *           The update is
 *           uNew=u+alpha*H*g;
 *           where alpha=(2*M/(M+1))*(fStar-q)/norm(d), d=g and H is a
 *           matrix that also depends on M and d, which is a design
 *           parameter. Also, after NR iterations, the recursive turms that
 *           go into H reset.
 *         4 The r Space Dilation Algorithm from [5]. This is the same as 3
 *           except d=g-gPrev except for the first iteration and if
 *           g=gPrev, in which case d=g.
 * maxIter The maximum number of iterations to perform. This should be
 *         >=1.
 *  AbsTol The absolute duality gap to use for convergence determiniation.
 *         Convergence is declared if (fStar-qStar)<=AbsTol.
 *  RelTol The relative duality gap to use for convergence determiniation.
 *         Convergence is declared if (fStar-qStar)<=RelTol*abs(qStar).
 * param1, param2, param3 Parameters for the selected subgradient algorithm.
 *         There are given in order for each method:
 *         Method 0: 'gamma'.
 *         Method 1: 'M' and 'r' with defaults 2.8 and 0.06.
 *         Method 2: 'a and 'b' with defaults 0.3 and 1.5.
 *         Method 3,4: 'M', 'normBound', and 'NR' with defaults 2, n3, and
 *                    eps(). normBound is such that if norm(B.'*d) in the
 *                    computation of the matrix H is <= normBound, then B
 *                    is reset to the identity matrix.
 *
 *OUTPUTS: The outputs are put in tuples, fStar, qStar and u. The return
 *         value is an exitCode, which takes the following values:
 *               -3 A non-finite number arose in the dual variables, so the
 *                  algorithm stopped.
 *               -2 The algorithm did not converge within the alotted
 *                  number of iterations. However, a valid feasible
 *                  assignment was found.
 *               -1 A subproblem in computing the dual or in obtaining a
 *                  feasible primal solution was infeasible. A feasible
 *                  valid assignment might not have been found.
 *               >=0 Values that are zero or positive indicate the
 *                  convergence was obtained and the returned value is the
 *                  number of iterations.
 *
 *This function implements the algorithm of [1], but modified so that it
 *does not have any unconstrained indices and offering different
 *subgradient methods. We are solving the "operations research" 3D
 *assignment problem rather than the "data fusion" 3D assignment problem of
 *[1].
 *
 *REFERENCES:
 *[1] K. Pattipati, S. Deb, Y. Bar-Shalom, and R. B. Washburn Jr., "A
 *   new relaxation algorithm and passive sensor data association," IEEE
 *   Transactions on Automatic Control, vol. 37, no. 2, pp. 198-213, Feb.
 *   1992.
 *[2] B. T. Polyak, "Minimization of unsmooth functionals," USSR
 *   Computational Mathematics and Mathematical Physics, vol. 9, no. 3, pp.
 *   14-29, 1969.
 *[3] M. A. Bragin, P. B. Luh, J. H. Yan, N. Yu, and G. A. Stern,
 *   "Convergence of the surrogate Lagrangian relaxation method,? Journal
 *   of Optimization Theory and Applications, vol. 164, no. 1, pp. 173-201,
 *   Jan. 2015.
 *[4] D. P. Bertsekas, Nonlinear Programming, 3rd ed. Belmont, MA: Athena
 *   Scientific, 2016, Chapter 7.5.
 *[5] N. Z. Shor, "Utilization of the operation of space dilation in the
 *   minimization of convex functions," Cybernetics, vol. 6, no. 1, pp. 7-
 *   15, Dec. 1972.
 *
 *February 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
    const ptrdiff_t n1=(ptrdiff_t)nDims[0];
    const ptrdiff_t n2=(ptrdiff_t)nDims[1];
    const ptrdiff_t n3=(ptrdiff_t)nDims[2];
    const ptrdiff_t n1n2=(ptrdiff_t)n1*n2;
    //These typed constant values are needed for the input to some BLAS
    //functions, because the functions take everything as pointers.
    const ptrdiff_t onePtrDiff=1;
    const char TChar='T';
    const char NChar='N';
    const double oneDouble=1;
    const double zeroDouble=0;
    const bool unconstU=(n1==n2)&&(n2==n3);
    //beta is used if Bertsekas' Heuristic Method is selected as the
    //subgradient method.
    double beta=1;
    size_t k;
    double q;
    //Parameters for the 2D assignment function.
    double *u2D, *v2D;
    ptrdiff_t *row4col;
    void *tempBuffer2DAssignFeas, *tempBuffer2DAssign;
    //These are needed if Bragin's method is used.
    double alphaPrev=0;
    double gNormPrev=0;
    //These are needed if a space dilation algorithm is selected.
    double *gPrev=NULL;
    double *B=NULL;
    double *H=NULL;
    double *d=NULL;
    double *dPrev=NULL;
    double *R=NULL;
    double rho=0;
    void *curBuffPos=tempSpace;
    //Divide the buffer tempSpace among the variables used in this
    //function.

    ptrdiff_t * gamma2=(ptrdiff_t*)curBuffPos;//Size n1*n2
    curBuffPos=(void*)((uint8_t*)curBuffPos+nDims[0]*nDims[1]*sizeof(ptrdiff_t));
    
    ptrdiff_t * gamma1=(ptrdiff_t*)(curBuffPos);//Size n1
    curBuffPos=(void*)((uint8_t*)curBuffPos+nDims[0]*sizeof(ptrdiff_t));
    
    double * g=(double*)(ptrdiff_t*)curBuffPos;//Size n3
    curBuffPos=(void*)((uint8_t*)curBuffPos+nDims[2]*sizeof(double));
    
    double * d2=(double*)curBuffPos;//Size n1*n2;
    curBuffPos=(void*)((uint8_t*)curBuffPos+nDims[0]*nDims[1]*sizeof(double));

    //The size of tempBuffer2DAssignFeas is determined by the size needed by
    //the updateBestFeasSol3DC function, because it needs more than
    //assign2D and its inputs. updateBestFeasSol3DC needs
    //least assign2DCBufferSize(n3,n1)+(n1+n3)*sizeof(ptrdiff_t)+
    //(n1+n3+n1*n3)*sizeof(double) bytes.
    if(subgradMethod>2) {//If a space dilation algorithm is used.
        gPrev=(double*)curBuffPos;//gPrev is length n3.
        B=gPrev+n3;//B is an n3Xn3 matrix.
        H=B+n3*n3;//H is an n3Xn3 matrix.
        d=H+n3*n3;//d is a length n3X1 vector.
        dPrev=d+n3;//dPrev is a length n3X1 vector.
        R=dPrev+n3;//R is an n3*n3 matrix.
        tempBuffer2DAssignFeas=(void*)(R+n3*n3);
        
        //Initialize.
        rho=(param1-1)/(param1+1);//param1 is M with space dilation.
    } else {
        tempBuffer2DAssignFeas=curBuffPos;
    }
    //We will not be calling assign2DC at the same time as
    //updateBestFeasSol3DC in this function, so the buffer needed by the
    //updateBestFeasSol3DC can overlap with that used by the assign2DC
    //function and its inputs u2D, v2D, and row4col. Also, the buffer
    //needed by updateBestFeasSol3DC is LARGER than that needed by
    //assign2DC and its inputs, so it is the determining factor is how much
    //memory this function needs in its tempSpace input.
    u2D=(double*)tempBuffer2DAssignFeas;//u2D is length n2.
    v2D=u2D+n2;//v2D is length n1.
    row4col=(ptrdiff_t*)(v2D+n1);//row4col is length n2.
    //assign2DC requires this buffer to be at least
    //assign2DCBufferSize(n2, n1) bytes in size.
    tempBuffer2DAssign=(void*)(row4col+n2);

    *qStar=-(double)INFINITY;
    *fStar=(double)INFINITY;

    //Initialize the dual variables to zero.
    memset(u,0,(size_t)n3*sizeof(double));
    
    for(k=0;k<maxIter;k++) {
        ptrdiff_t i1,i2,i3;
        double gNorm2;
//////
//DUAL COST AND SUBGRADIENT UPDATE.
//////
        //These loops essentially do:
        //[d2,gamma2]=min(bsxfun(@plus,C,reshape(u,[1,1,n3])),[],3);
        //d2=d2';
        //The transpose on d2 is necessary, because assign2DC requires that
        //the number of rows of the input be >= the number of columns and
        //we know that n2>=n1.
        {
            ptrdiff_t *curGamma2=gamma2;

            for(i2=0;i2<n2;i2++) {
                for(i1=0;i1<n1;i1++) {
                    double minVal=C[i1+n1*i2]+u[0];
                    ptrdiff_t minIdx=0;

                    for(i3=1;i3<n3;i3++) {
                        double curVal=C[i1+n1*i2+n1n2*i3]+u[i3];

                        if(curVal<minVal) {
                            minVal=curVal;
                            minIdx=i3;
                        }
                    }

                    //Store in a transposed order.
                    d2[i2+n2*i1]=minVal;

                    *curGamma2=minIdx;
                    curGamma2++;
                }
            }
        }

        {
            double minVal;
            //This function modifies (the transposed) d2, but it does not
            //matter, because it isn't used again in this loop.
            if(assign2DC(false, d2, &minVal, row4col,gamma1, tempBuffer2DAssign, v2D, u2D, (size_t)n2, (size_t)n1)) {
                //If the 2D assignment problem is infeasible.
                return -1;
            }
            
            q=minVal-sumVectorD(u,(size_t)n3);//The dual cost.
        }

        //Keep track of the maximum q value. This is used for testing
        //convergence.
        if(q>*qStar) {
            double costGap;
            *qStar=q;
            
            costGap=*fStar-*qStar;//The duality gap.
            
            //The correctness of this on the first iterations requires
            //proper handling of NaNs.
            if(costGap<=AbsTol||(costGap<fabs(*qStar)*RelTol)) {
                return (ptrdiff_t)(k+1);//The algorithm converged.
            }
        }
        
        //Compute the subgradient.
        //g=-1*ones(n3,1);
        for (i3=0;i3<n3;i3++) {
            g[i3]=-1;
        }
        for(i1=0;i1<n1;i1++) {
            i2=gamma1[i1];
            i3=gamma2[i1+n1*i2];
            
            g[i3]+=1;
        }
//////
//TEST THE GRADIENT AND FINISH THE LOOP.
//////

        //BLAS function for the dot product to get g'.*g.
        gNorm2=ddot(&n3,g,&onePtrDiff,g,&onePtrDiff);
        
        //If no constraints are violated, then we have a feasible solution
        //and the algorithm has converged to a local or global minimum
        //point.
        if(gNorm2==0) {
            //If the dual cost is less than the best primal solution found
            //with a heuristic, then we use the tuples for that solution.
            //Otherwise, the best primal solution to this point is
            //returned.
            
            if(q<*fStar) {
                *fStar=q;
                
                //Record the tuples
                for(i1=0;i1<n1;i1++) {
                    i2=gamma1[i1];
                    i3=gamma2[i1+n1*i2];

                    tuples[3*i1]=i1;
                    tuples[3*i1+1]=i2;
                    tuples[3*i1+2]=i3;
                }
            }
            
            return (ptrdiff_t)(k+1);//The algorithm converged.
        }

        {
            const ptrdiff_t retVal=updateBestFeasSol3DC(tuples, fStar, tempBuffer2DAssignFeas, nDims, C, gamma1, *qStar, AbsTol, RelTol);
            //If the primal converged.
            if(retVal==0) {
                return (ptrdiff_t)(k+1);
            }
            
            if(retVal<0) {
                //If the feasible assignment problem is infeasible.
                return -1;
            }
        }

        //Perform the subgradient update.
        switch(subgradMethod) {
            case 0://Polyak's method
            {//param1 is gammaVal
                const double alpha=param1*((*fStar-q)/gNorm2);
                
                //u=u+alpha*g
                daxpy(&n3,&alpha,g,&onePtrDiff,u,&onePtrDiff);
            }
            break;
            case 1://Bragin's Method
            {
                const double gNorm=sqrt(gNorm2);
                double alpha;
                
                if(k==0) {
                    alpha=((*fStar-q)/gNorm2);
                } else {
                    double kD=(double)k;
                    //param1 is M.
                    //param2 is r.
                    
                    alpha=(1-1/(param1*pow(kD,1-pow(kD,-param2))))*alphaPrev*(gNormPrev/gNorm);
                }
                //u=u+alpha*g
                daxpy(&n3,&alpha,g,&onePtrDiff,u,&onePtrDiff);
                
                alphaPrev=alpha;
                gNormPrev=gNorm;
            }
            break;
            case 2://Bertsekas' Heuristic Method
            {
                double alpha;
                if(q<*qStar) {
                    beta++;
                } else if(beta-1<1) {
                    beta=1;
                } else {
                    beta--;
                }
                //param1 is a and param2 is b
                alpha=fabs(((1+param1)/(pow(beta,param2)))*(*qStar)-q)/gNorm2;
                
                //u=u+alpha*g
                daxpy(&n3,&alpha,g,&onePtrDiff,u,&onePtrDiff);
            }
            break;
            default://Shor's Space Dilation Algorithm or the r Space
                    //Dilation Algorithm.
            {
                double normVal;
                double alpha;

                if(subgradMethod==3||k==0) {//Shor's algorithm
                    //d=g;
                    dcopy(&n3,g,&onePtrDiff,d,&onePtrDiff);
                } else {
                    bool allEqual=true;
                    //Check if g==gPrev
                    for(i3=0;i3<n3;i3++) {
                        if(g[i3]!=gPrev[i3]) {
                            allEqual=false;
                            break;
                        }
                    }
                    
                    //d=g;
                    dcopy(&n3,g,&onePtrDiff,d,&onePtrDiff);
                    
                    if(allEqual==false) {
                        const double negOne=-1;
                        //d=-gPrev+d;
                        
                        daxpy(&n3,&negOne,gPrev,&onePtrDiff,d,&onePtrDiff);
                    }
                }
                
                //param3 is NR
                if(k%param3==0) {
                    //B=eye(n3,n3);
                    identMatD((size_t)n3,B);

                    normVal=dnrm2(&n3,d,&onePtrDiff);
                    //We set H to zero. otherwise, if the unitialized H
                    //happened to contain NaNs or Inf terms, then the dgemm
                    //functon setting it would propagate NaNs (because
                    //there is a zero time H).
                    memset(H,0,(size_t)(n3*n3)*sizeof(double));
                } else {

                    //H(:,1)=B'*dPrev; -- we are just using part of the H matrix
                    //to store this temporary result. We can't store it
                    //back into dPrev -- the dgemv function does not
                    //support that.
                    dgemv(&TChar,&n3,&n3,&oneDouble,B,&n3,dPrev,&onePtrDiff,&zeroDouble,H,&onePtrDiff);

                    //normVal=norm(H(:,1)); --H(:,1) holds B'*dPrev
                    normVal=dnrm2(&n3,H,&onePtrDiff);

                    //We overwrite H(:,1) with zeta. zeta=B'*dPrev/normVal;
                    {
                        const double temp=1/normVal;
                        dscal(&n3,&temp,H,&onePtrDiff);
                    }

                    //Store zeta in dPrev.
                    dcopy(&n3,H,&onePtrDiff,dPrev,&onePtrDiff);

                    //H=eye(n3,n3);
                    identMatD((size_t)n3,H);
                    //H=H+(rho-1)*(zeta*zeta.');
                    {
                        const double temp=rho-1;
                        dger(&n3,&n3,&temp,dPrev,&onePtrDiff,dPrev,&onePtrDiff,H,&n3);
                    }

                    //R=B*H;
                    dgemm(&NChar,&NChar,&n3,&n3,&n3,&oneDouble,B,&n3,H,&n3,&zeroDouble,R,&n3);
                    //B=R;
                    {
                        const ptrdiff_t n3n3=n3*n3;
                        dcopy(&n3n3,R,&onePtrDiff,B,&onePtrDiff);
                    }

                    //dPrev=B'*d
                    dgemv(&TChar,&n3,&n3,&oneDouble,B,&n3,d,&onePtrDiff,&zeroDouble,dPrev,&onePtrDiff);

                    //normVal=norm(dPrev)
                    normVal=dnrm2(&n3,dPrev,&onePtrDiff);
                    if(normVal<=param2) {//param2 is normBound
                        //B=eye(n3,n3);
                        identMatD((size_t)n3,B);
                        normVal=dnrm2(&n3,d,&onePtrDiff);
                    }
                }
                //normVal now holds norm(B.'*d).

                //H=(B*B.')/norm(B.'*d);
                {
                    const double temp=1/normVal;
                    
                    dgemm(&NChar,&TChar,&n3,&n3,&n3,&temp,B,&n3,B,&n3,&zeroDouble,H,&n3);
                }

                //normVal=norm(d);
                normVal=dnrm2(&n3,d,&onePtrDiff);

                //alpha=(2*M/(M+1))*((fStar-q)/norm(d));
                alpha=(2*param1/(param1+1))*((*fStar-q)/normVal);

                //u=u+alpha*H*d;
                dgemv(&NChar,&n3,&n3,&alpha,H,&n3,d,&onePtrDiff,&oneDouble,u,&onePtrDiff);

                //dPrev=d;
                dcopy(&n3,d,&onePtrDiff,dPrev,&onePtrDiff);
                
                if(subgradMethod!=3) {
                    //gPrev=g
                    dcopy(&n3,g,&onePtrDiff,gPrev,&onePtrDiff);
                }
            }
        }
        
        //if(any(~isfinite(u)))
        for(i3=0;i3<n3;i3++) {
            if(!isfinite(u[i3])) {
                //This can sometimes occur with big problems and poor stepsizes.
                return -3;
            }
        }
        //Unless the constraints are equality constraints, the dual
        //variables have to be clipped.
        if(unconstU==false) {
            for(i3=0;i3<n3;i3++) {
                if(u[i3]<0) {
                    u[i3]=0;
                }
            }
        }
    }
    //The maximum number of iterations was hit.  
    return -2;
}

static ptrdiff_t updateBestFeasSol3DC(ptrdiff_t *tuples, double * fStar,void *tempBuffer2DAssign,const size_t *nDims, const double *C, const ptrdiff_t *gamma1, const double qStar, const double AbsTol, const double RelTol) {
/**UPDATEBESTFEASSOL3DC This function implements a subroutine to
 *        heuristically obtain a feasible solution given a partial set of
 *        assignments. Here, it is just for the 3D assignment problem.
 *        Additionally, this function updates the best feasible solution
 *        cost fStar and its associated set of tuples.
 *
 *INPUTS: tuples The set of tuples associated with the current best primal
 *               cost. This is a pointer to space enough to hold a 3Xn1
 *               matrix of tuples as ptrdiff_t values. This gets updated.
 *         fStar A pointer to the best (lowest) primal cost value
 *               encountered thus far as a double. This gets updated.
 * tempBuffer2DAssign A buffer for intermediate values. This must be at
 *               least updateBestFeasSol3DCBufferSize(nDims) bytes in
 *               length.
 *         nDims A length-3 array of size_t values that are the
 *               dimensions of the matrix C.
 *             C The original n1Xn2Xn3 3D cost matrix of doubles with
 *               n1>=n2>=n3. Elements are stored by column.
 *        gamma1 The length-n1 set of assignment of elements in the first
 *               index of C to those in the second index of C. The indices
 *               are stored in a ptrdiff_t type.
 *         qStar The best dual cost value encountered this far as a double.
 *        AbsTol The absolute threshold for determining convergence of the
 *               algorithm according to the duality gap.
 *        RelTol The threshold for determining convergence of the algorithm
 *               according to the relative duality gap.
 *
 *OUTPUTS: Outputs are placed in tuples and fStar. The return value is
 *                 an exitCode, which takes the following values:
 *                 -1 The subproblem posted is not feasible. 
 *                  0 Convergence achieved.
 *                  1 No errors occurred, but convegrence was not achieved.
 *
 *February 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
    const size_t n1=nDims[0];
    const size_t n2=nDims[1];
    const size_t n3=nDims[2];
    const size_t n1n2=n1*n2;
    double minVal;
    size_t i1, i2, i3;
    void *curBuffPos=tempBuffer2DAssign;
    
    //Allocate the memory in tempBuffer2DAssign. Once divided, each
    //variable is disjoint, so they are declared using the restrct keyword. 
    double * CFeas=(double *)curBuffPos;//CFeas is n1Xn3
    curBuffPos=(void*)((uint8_t*)curBuffPos+n1*n3*sizeof(double));
    
    double * u2D=(double*)curBuffPos;//u2D is n3X1
    curBuffPos=(void*)((uint8_t*)curBuffPos+n3*sizeof(double));
    
    double * v2D=(double*)curBuffPos;//v2D is n1X1
    curBuffPos=(void*)((uint8_t*)curBuffPos+n1*sizeof(double));
    
    ptrdiff_t * row4col=(ptrdiff_t*)(curBuffPos);//Length n3
    curBuffPos=(void*)((uint8_t*)curBuffPos+n3*sizeof(ptrdiff_t));

    ptrdiff_t * gammaTilde3=(ptrdiff_t*)curBuffPos;//Length n1
    curBuffPos=(void*)((uint8_t*)curBuffPos+n1*sizeof(ptrdiff_t));
    
    //tempBuffer2DAssign is at least size assign2DCBufferSize(n1,n3).
    tempBuffer2DAssign=curBuffPos;

    //We store the elements in CFeas in a transposed manner, because
    //assign2DC requires that the number of rows is >= the number of
    //columns and we know that n1<=n3.
    for(i1=0;i1<n1;i1++) {
        i2=(size_t)gamma1[i1];

        for(i3=0;i3<n3;i3++) {
            CFeas[i3+n3*i1]=C[i1+n1*i2+n1n2*i3];
        }
    }

    if(assign2DC(false, CFeas, &minVal, row4col, gammaTilde3, tempBuffer2DAssign, v2D, u2D, (size_t)n3, (size_t)n1)) {
        //If the 2D assignment problem is infeasible.
        return -1;
    }

    if(minVal<*fStar) {
        double costGap;
        
        *fStar=minVal;
        
        for(i1=0;i1<n1;i1++) {
            tuples[3*i1]=(ptrdiff_t)i1;
            tuples[3*i1+1]=gamma1[i1];
            tuples[3*i1+2]=gammaTilde3[i1];
        }
        
        //The duality gap.
        costGap=*fStar-qStar;
        
        if(costGap<=AbsTol||costGap<fabs(qStar)*RelTol) {
            return 0;//The algorithm converged.
        } else {
            return 1;
        }
    } else {
        return 1;
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
