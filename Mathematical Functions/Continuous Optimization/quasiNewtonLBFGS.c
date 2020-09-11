/**QUASINEWTONLBFGS Perform unconstrained nonlinear optimization using the
 *                  limited-memory Broyden-Fletcher-Goldfarb-Shanno (BFGS)
 *                  algorithm. The algorithm performs unconstrained
 *                  minimization of a nonlinear function without one having
 *                  to provide a Hessian matrix. It is a type of
 *                  Quasi-Newton algorithm. It minimizes a function f(x),
 *                  where x is a vector. If the C option in the line-search
 *                  routine is not zero, then the function being minimized
 *                  is f(x)+C*norm(x,1), where C>0 (or the norm can be
 *                  taken using a subset of the elements of x), in which
 *                  case this is the orthant-wise limited-memory
 *                  Quasi-Newton method. Unlike the traditional BFGS
 *                  algorithm, this algorithm does not store the entire
 *                  inverse Hessian matrix estimate, thus making it
 *                  significantly more efficient for very large problems.
 *                  If one wishes to zero a vector (not a scalar), then
 *                  NewtonsMethod is more appropriate as this function
 *                  assumes the Hessian matrix is symmetric.
 *
 *INPUTS: f A handle to the function (and its gradient) over which the
 *          minimization is to be performed. The function [fVal,gVal]=f(x)
 *          takes the NX1 x vector returns the real scalar function value
 *          fVal and gradient gVal at the point x. 
 *       x0 The NX1-dimensional point from which the minimization starts.
 *  numCorr The number of corrections to approximate the inverse Hessian
 *          matrix. The default if omitted or an empty matrix is passed is
 *          6. The L-BFGS documentation recommends not using fewer than 3.
 *  epsilon The parameter determining the accuracy of the desired solution.
 *          The function terminates when norm(g) < epsilon*max([1, norm(x)])
 *          where g is the gradient. The default if omitted or an empty
 *          matrix is passed is 1e-6.
 * deltaTestDist The number of iterations back to use to compute the
 *          decrease of the objective function if a delta-based convergence
 *          test is performed. If zero, then no delta-based convergence
 *          testing is done. The default if omitted or an empty matrix is
 *          passed is zero.
 *    delta The delta for the delta convergence test. This determines the
 *          minimum rate of decrease of the objective function. Convergence
 *          is determined if (f'-f)/f<delta, where f' is the value of the
 *          objective function f deltaTestDist iterations ago, and f is the
 *          current objective function value. The default if this parameter
 *          is omitted or an empty matrix is passed is 0.
 * lineSearchParams An optional structure whose members specify tolerances
 *          for the line search. The parameters are described as in the
 *          lineSearch function. Possible members are algorithm, C, fTol,
 *          wolfeTol, xTol, minStep, maxStep, maxIter, l1NormRange. If this
 *          parameter is omitted or an empty matrix is passed, then the
 *          default values as described in the lineSearch function are
 *          used. If any member of this structure is omitted or is assigned
 *          an empty matrix, then the default value will be used.
 * maxIterations The maximum number of iterations to use for the overall
 *          L-BFGS algorithm. The default if this parameter is omitted or
 *          an empty matrix is past is 1000. If this parameter is set to
 *          zero, the algorithm will continue either until convergence is
 *          obtained or until an error occurs.
 * progressCallback An optional Matlab function handle that is called
 *          during each step of the algorithm with information on the
 *          progress. This can be used to analyze how well the optimization
 *          is working. The function takes 8 arguments and returns one
 *          argument. The format is
 *          progressCallback(x,...%The state
 *                           g,...%The gradient
 *                           fx,...%The function value
 *                           xnorm,...%The l2 norm of the state
 *                           gnorm,...%The l2 norm of the gradient
 *                                 ...%The line-search step used for this
 *                                 ...%iteration
 *                           step,...
 *                           k,...%The iteration number.
 *                             ...%The number of function evaluations
 *                             ...%called for this iteration.
 *                           ls);
 *          If the return value is 0, then the optimization process will
 *          continue. If the return value is nonzero, then the optimization
 *          process will not continue. The return value should be a real
 *          scalar value. If this parameter is omitted or an empty
 *          matrix is passed, then no callback is performed.
 *
 *OUTPUTS: xMin The value of x at the minimum point found. If exitCode is
 *              negative, then this value might be invalid.
 *         fMin The cost function value at the minimum point found. If
 *              exitCode is negative, then this value might be invalid.
 *     exitCode A value indicating the termination condition of the
 *              algorithm. Negative values indicate errors Possible
 *              values are:
 *                  0 The algorithm terminated successfully.
 *                  1 Termination according to the delta stopping
 *                    criterion occurred.
 *                  2 The initial value already minimizes the objective
 *                    function.
 *              -1023 A logical error in the code occurred.
 *              -1022 Insufficient memory.
 *              -1021 The optimization was cancelled by the user callback
 *                    progress function returning a nonzero value.
 *              -1001 A finite precision error occurred or no line-search
 *                    step satisfies the sufficient decrease and curvature
 *                    conditions.
 *              -1000 The line-search step size became less than minStep.
 *               -999 The line-search step size became larger than maxStep.
 *               -998 The maximum number of line-search iterations was
 *                    reached.
 *               -997 The maximum number of overall iterations was reached.
 *               -996 The relative width of the interval of uncertainty is
 *                    at most xTol
 *               -995 A negative line-search step occurred.
 *               -994 The current search direction increases the objective
 *                    function.
 *
 *This is a Matlab interface for the C implementation of the L-BFGS
 *algorithm of
 *https://github.com/chokkan/liblbfgs
 *which is based on the Fortran L-BFGS algorithm of
 *http://www.ece.northwestern.edu/~nocedal/lbfgs.html
 *but extends the work. The basic algorithm is described in [1] and [2].
 *However, the C implementation contains additional controls over the type
 *of line search perfromed and the convergence conditions. The above BFGS
 *library has been slightly modified to use Matlab's memory allocation and
 *deallocation routines.
 *
 *The algorithm can be compiled for use in Matlab  using the 
 *CompileCLibraries function.
 *
 *The algorithm is run in Matlab using the command format
 *[xMin,fMin,exitCode]=quasiNewtonLBFGS(f,x0);
 *or if more options are used
 *[xMin,fMin,exitCode]=quasiNewtonLBFGS(f,x0,numCorr,epsilon,deltaTestDist,delta,lineSearchParams,maxIterations,progressCallback);
 *
 *Examples:
 *The first example is that used in the lineSearch file. 
 * f=@(x)deal((x(1)+x(2)-3)*(x(1)+x(2))^3*(x(1)+x(2)-6)^4,... %The function
 *            [(-6+x(1)+x(2))^3*(x(1)+x(2))^2*(54+8*x(1)^2+x(2)*(-45+8*x(2))+x(1)*(-45+16*x(2)));
 *            (-6+x(1)+x(2))^3*(x(1)+x(2))^2*(54+8*x(1)^2+x(2)*(-45+8*x(2))+x(1)*(-45+16*x(2)))]);%And the gradient as the second return.
 * %Note that the deal function is used to make an anonymous function have
 * %two outputs.
 * x0=[0.5;0.25];
 * [xMin,fMin,exitCode]=quasiNewtonLBFGS(f,x0)
 *The optimum point found is such that sum(xMin) is approximately
 *1.73539450 with a minimum function value of approximately -2.1860756.
 *
 *The second example requires that the cost function be placed in a
 *separate file. The cost function is
 * function [fx,g]=objFun(x)
 *    n=length(x);
 *    fx=0;
 *   
 *    g=zeros(n,1);
 *    for i=0:(n-1)
 *        if(mod(i,2)==1)
 *            continue;
 *        end
 *        t1=1-x(i+1);
 *        t2=10*(x(i+1+1)-x(i+1)^2);
 *        g(i+1+1)=20*t2;
 *        g(i+1)=-2*(x(i+1)*g(i+1+1)+t1);
 *        fx =fx+t1^2+t2^2;
 *    end
 * end
 *It is used as
 * n=100;
 * x0=zeros(n,1);
 * for i=0:(n-1)
 *    if(mod(i,2)==1)
 *        continue;
 *    end
 *
 *    x0(i+1)=-1.2;
 *    x0(i+1+1)=1;
 * end
 * [xMin,fMin,exitCode]=quasiNewtonLBFGS(f,x0)
 *whereby the optimal solution is all ones with a minimum function value of
 *zero. This second example is the same as that provided with the L-BFGS
 *library in C.
 *
 *REFERENCES:
 *[1] Liu, D. C.; Nocedal, J. (1989). "On the Limited Memory Method for
 *    Large Scale Optimization". Mathematical Programming B 45 (3):
 *    503-528. doi:10.1007/BF01589116.
 *[2] Byrd, Richard H.; Lu, Peihuang; Nocedal, Jorge; Zhu, Ciyou (1995).
 *    "A Limited Memory Algorithm for Bound Constrained Optimization".
 *    SIAM Journal on Scientific and Statistical Computing 16 (5):
 *    1190-1208. doi:10.1137/0916069.
 *[3] J. J. Mor� and D. J. Thuente, "Line search algorithms with
 *    guaranteed sufficient decrease," ACM Transactions on Mathematical
 *    Software, vol. 20, no. 3, pp. 286-307, Sep. 1994.
 *
 *August 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

/*This header is required by Matlab.*/
#include "mex.h"

/*The header for the L-BFGS library with the definition making sure that
 *double precision floating point values are used.*/
#define  LBFGS_FLOAT 64
#include "lbfgs.h"
#include "MexValidation.h"
#include <limits.h>

//Prototype for the callback function wrappers.
static double MatlabCallback(void *MatlabFunctionHandles,
         //The initial state; this must have been allocated using mxCalloc.
            const double *x,
            double *gVal,//The gradient vector is filled in here.
            const int n,//The dimensionality of the state/ gradient.
            const double step//Unused input
            );

static int progressMatlabCallback(
    void *MatlabFunctionHandles,
    const double *x,//The state
    const double *g,//The gradient
    const double fx,//The function value
    const double xnorm,//The l2 norm of the state
    const double gnorm,//The l2 norm of the gradient
    const double step,//The line-search step used for this iteration
    int n,//The dimensionality of the state
    int k,//The iteration number.
    int ls//The number of function evaluations called for this iteration.
    );

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    size_t xDim;
    double *x;
    lbfgs_parameter_t param;
    double fVal;
    int exitCode;
    const mxArray *MatlabFunctionHandles[2];//To hold the callback functions.

    if(nrhs<2||nrhs>9){
        mexErrMsgTxt("Wrong number of inputs");
    }
    
    if(nlhs>3) {
        mexErrMsgTxt("Wrong number of outputs.");
    }

    //Check that a function handle was passed.
    if(!mxIsClass(prhs[0],"function_handle")) {
        mexErrMsgTxt("The first input must be a function handle.");
    }
    MatlabFunctionHandles[0]=prhs[0];

    //Check that a valid initial estimate x was passed. 
    checkRealDoubleArray(prhs[1]);
    xDim=mxGetM(prhs[1]);
    if(xDim<1||mxGetN(prhs[1])!=1) {
        mexErrMsgTxt("The point x has the wrong dimensionality."); 
    }
    
    //The function in the library only uses integers for dimensions.
    if(xDim>INT_MAX) {
        mexErrMsgTxt("The problem has too many dimensions to be solved using this function.");
    }
    
    //Allocate space for the state.
    x=mxCalloc(xDim, sizeof(double));
    //Copy the passed initial estimate.
    memcpy(x, mxGetDoubles(prhs[1]), sizeof(double)*xDim);

    //Set all of the default values for the param structure. If other
    //inputs are provided, then these will be changed.
    param.max_iterations=1000;
    param.past=0;
    param.delta=0;
    param.epsilon=1e-6;
    param.m=6;
    //Default parameters for the line search
    param.linesearch=LBFGS_LINESEARCH_MORETHUENTE;
    param.orthantwise_c=0;
    param.ftol=1e-6;
    param.wolfe=0.9;
    param.gtol=param.wolfe;
    param.xtol=1e-16;
    param.min_step=1e-20;
    param.max_step=1e20;
    param.max_linesearch=20;
    param.orthantwise_start=0;
    param.orthantwise_end=(int)xDim-1;
    
    //Now load optional parameters.
    if(nrhs>2&&!mxIsEmpty(prhs[2])) {
        param.m=getIntFromMatlab(prhs[2]);
        if(param.m<0) {
            mexErrMsgTxt("param.m must be a positive integer.");
        }
    }
    
    if(nrhs>3&&!mxIsEmpty(prhs[3])) {
        param.epsilon=getDoubleFromMatlab(prhs[3]);
    }
    
    if(nrhs>4&&!mxIsEmpty(prhs[4])) {
        param.past=getIntFromMatlab(prhs[4]);
        if(param.past<0) {
            mexErrMsgTxt("param.past must be a positive integer.");
        }
    }
    
    if(nrhs>5&&!mxIsEmpty(prhs[5])) {
        param.delta=getDoubleFromMatlab(prhs[5]);
    }
    
    if(nrhs>6&&!mxIsEmpty(prhs[6])) {
        //We have to check for every possible structure member.
        mxArray *theField;
        
        if(!mxIsStruct(prhs[6])) {
            mexErrMsgTxt("The line search parameters must be given in a structure.");
        }
        
        theField=mxGetField(prhs[6],0,"algorithm");
        if(theField!=NULL) {//If the field is present.
            int algorithm=getIntFromMatlab(theField);

            switch(algorithm){
                case 0:
                    param.linesearch=LBFGS_LINESEARCH_MORETHUENTE;
                    break;
                case 1:
                    param.linesearch=LBFGS_LINESEARCH_BACKTRACKING_ARMIJO;
                    break;
                case 2:
                    param.linesearch=LBFGS_LINESEARCH_BACKTRACKING_WOLFE;
                    break;
                case 3:
                    param.linesearch=LBFGS_LINESEARCH_BACKTRACKING_STRONG_WOLFE;
                    break;
                default:
                    mexErrMsgTxt("Unknown line search algorithm specified.");      
            }  
        }
        
        theField=mxGetField(prhs[6],0,"C");
        if(theField!=NULL) {//If the field is present.
            param.orthantwise_c=getDoubleFromMatlab(theField);
        }
                
        theField=mxGetField(prhs[6],0,"fTol");
        if(theField!=NULL) {//If the field is present.
            param.ftol=getDoubleFromMatlab(theField);
        }
                        
        theField=mxGetField(prhs[6],0,"wolfeTol");
        if(theField!=NULL) {//If the field is present.
            param.wolfe=getDoubleFromMatlab(theField);
            param.gtol=param.wolfe;
        }
                                
        theField=mxGetField(prhs[6],0,"xTol");
        if(theField!=NULL) {//If the field is present.
            param.xtol=getDoubleFromMatlab(theField);
        }
                                        
        theField=mxGetField(prhs[6],0,"minStep");
        if(theField!=NULL) {//If the field is present.
            param.min_step=getDoubleFromMatlab(theField);
        }
        
        theField=mxGetField(prhs[6],0,"maxStep");
        if(theField!=NULL) {//If the field is present.
            param.max_step=getDoubleFromMatlab(theField);
        }
        
        theField=mxGetField(prhs[6],0,"maxIter");
        if(theField!=NULL) {//If the field is present.
            param.max_linesearch=getIntFromMatlab(theField);
        
            if(param.max_linesearch<0) {
               mexErrMsgTxt("param.max_linesearch must be positive."); 
            }
        }
        
        theField=mxGetField(prhs[6],0,"l1NormRange");
        if(theField!=NULL) {//If the field is present.
            double *indices, indexMin, indexMax;
        
            checkRealDoubleArray(theField);
            if(mxGetM(theField)*mxGetN(theField)!=2) {
                mexErrMsgTxt("The size of l1NormRange is incorrect."); 
            }

            indices=mxGetDoubles(theField);

            if(indices[0]!=floor(indices[0])||indices[1]!=floor(indices[1])) {
                mexErrMsgTxt("The indices in l1NormRange must be integers."); 
            }

            if(indices[0]<indices[1]) {
                indexMin=indices[0];
                indexMax=indices[1];
            } else {
                indexMin=indices[1];
                indexMax=indices[0];
            }

            if(indexMin<1||indexMax>xDim) {
                mexErrMsgTxt("Invalid range given in l1NormRange.");
            }

            param.orthantwise_start=(int)(indexMin)-1;
            param.orthantwise_end=(int)(indexMax)-1;
        }
    }
    
    if(nrhs>7&&!mxIsEmpty(prhs[7])) {
        param.max_iterations=getIntFromMatlab(prhs[7]);
        if(param.max_iterations<0) {
            mexErrMsgTxt("param.max_iterations must be a positive integer.");
        }
    }
    
    if(nrhs>8&&!mxIsEmpty(prhs[8])) {
        if(!mxIsClass(prhs[8],"function_handle")) {
            mexErrMsgTxt("The callback function must be a function handle.");
        }
        MatlabFunctionHandles[1]=prhs[8];
    } else {
        MatlabFunctionHandles[1]=NULL;
    }
    
    //Check the inputs for validity.
    if(param.epsilon < 0) {
        mexErrMsgTxt("The epsilon parameter must be nonnegative.");    
    }
    
    if(param.delta < 0.) {
        mexErrMsgTxt("The delta parameter must be nonnegative."); 
    }
    
    //The function tolerance must be nonegative and if a Wolfe parameter is
    //used, the Wolfe parameter must have the proper relation to the
    //function tolerance.
    if(param.ftol<0||((param.linesearch!=LBFGS_LINESEARCH_BACKTRACKING_ARMIJO)&&(param.wolfe <= param.ftol || 1 <= param.wolfe|| param.wolfe<0))) {
        mexErrMsgTxt("Invalid Wolfe parameter given.");
    }
    
    if(param.xtol < 0) {
        mexErrMsgTxt("Invalid XTol specified");
    }
            
    if(param.max_linesearch <= 0) {
        mexErrMsgTxt("Invalid maxIter specified");
    }
    
    if(param.orthantwise_c < 0) {
       mexErrMsgTxt("Invalid C value specified"); 
    }
    
    if(param.min_step < 0) {
       mexErrMsgTxt("Invalid minimum step size specified");  
    }
    
    if(param.max_step < param.min_step) {
        mexErrMsgTxt("The maximum step size is less than the minimum step size.");  
    }
    
    if(param.orthantwise_c!=0&&param.linesearch==LBFGS_LINESEARCH_MORETHUENTE) {
       mexErrMsgTxt("The More and Thuente algorithm cannot be used with C~=0.");  
    }

    //Now, the algorithm can be run.
    if(MatlabFunctionHandles[1]==NULL) {
        //If there is no callback function.
        exitCode=lbfgs((int)xDim,
               x,
               &fVal,
               &MatlabCallback,
               NULL,
               (void *)MatlabFunctionHandles,
               &param);
    }
    else {
        //If there is a callback function.
        exitCode=lbfgs((int)xDim,
                       x,
                       &fVal,
                       &MatlabCallback,
                       &progressMatlabCallback,
                       (void *)MatlabFunctionHandles,
                       &param);
    }
    
    //Set the outputs
    //First set x.
    plhs[0]=mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
    mxFree(mxGetDoubles(plhs[0]));
    mxSetDoubles(plhs[0], x);
    mxSetM(plhs[0], xDim);
    mxSetN(plhs[0], 1);
    
    if(nlhs>1) {
        //Return the minimum function value.
        plhs[1]=doubleMat2Matlab(&fVal,1,1);
        
        if(nlhs>2) {
            //Return the error code
            plhs[2]=intMat2MatlabDoubles(&exitCode,1,1);
        }
    }
}

static double MatlabCallback(void *MatlabFunctionHandles,
        //The initial state; this must have been allocated using mxCalloc.
            const double *x,
            double *gVal,//The gradient vector is filled in here.
            const int n,//The dimensionality of the state/ gradient.
            const double step//Unused input
            ) {
//The function returns the fval=f(x), the function value, and the gradient
//using the given Matlab function handle.
    mxArray *rhs[2];
    mxArray *lhs[2];
    double *oldPtr;
    double fVal;
    
    //feval in Matlab will take the function handle and the state as
    //inputs.
    //The first function handle is f. 
    rhs[0]=((mxArray**)MatlabFunctionHandles)[0];
    rhs[1]=mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);

    //Set the matrix data to x 
    oldPtr=mxGetDoubles(rhs[1]);
    //x will not be modified, but the const must be typecast away to use
    //the mxSetDoubles function.
    mxSetDoubles(rhs[1], (double*)x);
    mxSetM(rhs[1], (size_t)n);
    mxSetN(rhs[1], 1);

    //Get the function value and gradient.
    mexCallMATLAB(2,lhs,2,rhs,"feval");
    
    //Get the function value.
    fVal=getDoubleFromMatlab(lhs[0]);

    //Copy the gradient into gVal, checking for errors.
    verifySizeReal((size_t)n,1,lhs[1]);
    memcpy(gVal, mxGetDoubles(lhs[1]), sizeof(double)*(size_t)n);

    //Get rid of the returned Matlab matrices.
    mxDestroyArray(lhs[0]);
    mxDestroyArray(lhs[1]);
    
    //Set the data pointer back to what it was during allocation that
    //mxDestroyArray does not have a problem. 
    mxSetDoubles(rhs[1],oldPtr);
    mxSetM(rhs[1], 0);
    mxSetN(rhs[1], 0);

    //Get rid of the temporary natrix.
    mxDestroyArray(rhs[1]);

    return fVal;
}


static int progressMatlabCallback(
    void *MatlabFunctionHandles,
    const double *x,//The state
    const double *g,//The gradient
    const double fx,//The function value
    const double xnorm,//The l2 norm of the state
    const double gnorm,//The l2 norm of the gradient
    const double step,//The line-search step used for this iteration
    int n,//The dimensionality of the state
    int k,//The iteration number.
    int ls//The number of function evaluations called for this iteration.
    )
{
    mxArray *rhs[9];
    mxArray *lhs[1];
    size_t i;
    int retVal;
    
    //Allocate temporary variables to pass to Matlab.
    //The second function handle is the callback function.
    rhs[0]=((mxArray**)MatlabFunctionHandles)[1];
    rhs[1]=doubleMat2Matlab(x,(size_t)n,1);
    rhs[2]=doubleMat2Matlab(g,(size_t)n,1);
    rhs[3]=doubleMat2Matlab(&fx,1,1);
    rhs[4]=doubleMat2Matlab(&xnorm,1,1);
    rhs[5]=doubleMat2Matlab(&gnorm,1,1);
    rhs[6]=doubleMat2Matlab(&step,1,1);
    rhs[7]=intMat2MatlabDoubles(&k,1,1);
    rhs[8]=intMat2MatlabDoubles(&ls,1,1);
    
    mexCallMATLAB(1,lhs,9,rhs,"feval");
    
    retVal=getIntFromMatlab(lhs[0]);
    
    //Get rid of the returned Matlab matrix
    mxDestroyArray(lhs[0]);
    
    //Free the temporary variables.
    for(i=1;i<9;i++) {
        mxDestroyArray(rhs[i]);
    }
    
    //Returning a nonzero value will terminate the optimization
    return retVal;
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
