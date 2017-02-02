/**LINESEARCH Given a function f taking an NX1 dimensional vector x and a
 *            descent direction d, find the positive value alpha such that
 *            f(x+alpha*d) is minimized to an extent necessary for
 *            satisfying certain convergence criterion in quasi-Newton
 *            methods. Some of the algorithms can also minimize 
 *            f(x+alpha*d)+C*norm(x+alpha*d,1), where C is some parameter,
 *            which is useful in l1-norm optimization problems. The
 *            algorithms assume that the function f has only one finite
 *            minimum along the descent direction.
 *
 *INPUTS: f        A handle to the function (and its gradient) over which
 *                 the line search minimization is to be performed. The
 *                 function [fVal,gVal]=f(x) takes the NX1 x vector returns
 *                 the real scalar function value fVal and gradient gVal at
 *                 the point x. 
 *        x        The NX1-dimensional point from which the line search
 *                 starts.
 *        d        The NX1 vector defining the initial descent direction.
 *                 The minimization is performed such that f(x+alpha*d) is
 *                 minimized where alpha is a scalar parameter.
 *       algorithm A parameter specifying the line search algorithm to
 *                 use. The methods differ in the method by which they
 *                 determine that a sufficient reduction in the objective
 *                 function has occurred. Possible values are:
 *                 0) (The default if omitted or an empty matrix is
 *                    passed). Use the algorithm proposed by Moré and
 *                    Thuente in [1] to determine a sufficient decrease in
 *                    the function value from f(x). This method brackets
 *                    the uncertainty interval. A sufficient decrease in
 *                    the objective function is determined if the relative
 *                    width of the interval of uncertainty is at most XTol
 *                    or if the sufficient decrease condition of Armijo's
 *                    rule (see algorithm 2) that
 *                    f(x+alpha*d)<=fVal + fTol*alpha*gVal'*d
 *                    holds and the absolute curvature condition that
 *                    abs(g(x+alpha *d)'*d) <= wolfeTol*abs(g(x)'*d)
 *                    (See algorithm 3) holds. This algorithm will not work
 *                    with an orphant-wise parameter C being nonzero.
 *                 1) Use a backtracking-based algorithm with Armijo's
 *                    rule for determing a sufficient decrease in the
 *                    function value from f(x). Armijo's rule is orignally
 *                    from [2] and states that the result of a line search
 *                    is sufficient for use in a gradient-descent-based
 *                    optimization algorithm if
 *                    f(x+alpha*d)<=fVal + fTol*alpha*gVal'*d
 *                    where fVal and gVal are defined as in the first input
 *                    to this function and fTol is a tolerance value
 *                    (another input to this function). If the input
 *                    parameter C is not zero, then the optimization 
 *                    includes an L1-norm term as in the orphant-wise
 *                    method of [6].
 *                 2) Use a backtracking-based algorithm with the regular
 *                    Wolfe condition for determining a sufficient decrease
 *                    in the objective function. This condition is Armijo's
 *                    rule plus the curvature condition that
 *                    g(x+alpha*d)'*d >= wolfeTol*gVal'*d,
 *                    where gVal is defined as in the first input to this
 *                    function and wolfeTol is a tolerance parameter
 *                    (another input to this function). This new constraint
 *                    is based on a condition in [3,4]. Specifically,
 *                    condition iii of [3]. Note that 0<fTol<wolfeTol<1 to
 *                    satisfy the Wolfe conditions for a sufficient step.
 *                    If the input parameter C is not zero, then the
 *                    optimization includes an L1-norm term as in the
 *                    orphant-wise method of [6].
 *                 3) Use a backtracking-based algorithm with the strong
 *                    Wolfe condition for determining a sufficient decrease
 *                    in the objective function. This condition is Armijo's
 *                    rule plus the absolute curvature condition that
 *                    abs(g(x+alpha *d)'*d) <= wolfeTol*abs(g(x)'*d)
 *                    This modification of the Wolfe condition is described
 *                    in many places including in Chapter 3 of [5].
 *                    Note that 0<fTol<wolfeTol<1 to satisfy the Wolfe
 *                    conditions for s sufficient step. If the input
 *                    parameter C is not zero, then the optimization
 *                    includes an L1-norm term as in the orphant-wise
 *                    method of [6].
 *             C  If this parameter is given and is not zero or an empty
 *                matrix, then the cost function to minimize over alpha
 *                is f(x+alpha*d)+C*norm(x+alpha*d,1); where C must be
 *                positive. That is, C is the weight of the l1 norm term
 *                added to the cost function.
 *   stepSizeInit The initial stepsize for the search. If this is omitted
 *                or an empty matrix is passed, then the initial stepsize
 *                estimate is stepSizeInit = 1; This is equivalent to
 *                saying that one first tries to take a step to x+alpha*d
 *                with alpha=1 as alpha is the stepsize.
 *           fTol The tolerance condition for Armijo's rule. If this
 *                parameter is omitted or an empty matrix is apssed, the
 *                default value of 1e-6; is used.
 *       wolfeTol The wole tolerance condition for the gradient. If this
 *                parameter is omitted or an empty matrix is passed, the
 *                default value of 0.9 is used.
 *           xTol The tolerance condition for the bounding region of the
 *                line search. This is only used by algorithm 0. If this
 *                parameter is omitted or an empty matrix is passed, the
 *                default value of 1e-16 is used.
 *        minStep The minimum allowable step size. If omitted or an empty
 *                matrix is passed, the default value of 1e-20 is used.
 *        maxStep The maximum allowable step size. If omitted or an empty
 *                matrix is passed, the default value of 1e20 is used.
 *        maxIter The maximum number of iterations. If omitted or an empty
 *                matrix is passed, the default value of 20 is used.
 *    l1NormRange If C!=0, then this is used, if provided. This is a 2X1 or
 *                1X2 vector where the l1 norm in a cost function of
 *                f(x+alpha*d)+C*norm(x+alpha*d,1) is nor taken over all
 *                elements as written but rather as
 *                y=x+alpha*d;
 *                f(y)+C*norm(y(l1NormRange),1)
 *                If this parameter is omitted or an empty matrix is
 *                passed, then all elements of y are used in the l1 norm.
 *                The indices should be given as floating point doubles
 *                (the default format in Matlab). The index range can only
 *                go up to 2^32-1.
 *
 *OUTPUTS: xMin     The value of x at the minimum point found.
 *         fMin     The cost function value at the minimum point found.
 *         gMin     The value of the gradient of the cost function at the
 *                  minimum point found.
 *         alpha    The step size found such that xMin=x+alpha*d
 *         exitCode A value indicating the termination condition of the
 *                  algorithm. Negative values can indicate errors, but the
 *                  results are often still usable. Possible values are
 *                      0 or a positive number: The algorithm terminated
 *                        successfully.
 *                  -1023 A logical error in the code occurred.
 *                  -1001 A finite precision error occurred or no
 *                        line-search step satisfies the sufficient
 *                        decrease and curvature conditions.
 *                  -1000 The line-search step size became less than
 *                        minStep.
 *                   -999 The line-search step size became larger than
 *                        maxStep.
 *                   -998 The maximum number of iterations was reached.
 *                   -996 The relative width of the interval of uncertainty
 *                        is at most xTol
 *                   -995 A negative line-search step occurred.
 *                   -994 The current search direction increases the
 *                        objective function.
 *
 *Line searches play a pivotal role in multivariate optimization routines.
 *They are related to single variable optimization methods, such as in the
 *function goldenSectionSearch or in Matlab fminbnd function, in that they
 *try to minimize a function over one parameter (in this case a
 *multivariate function being minimized over a line). However, line search
 *methods generally try less to bound the minimum than to determine a point
 *that sufficiently decreases the objective function to assure convergence
 *of the overall multivariate optimization problem over time.
 *
 *Thus, outside of the context of larger multivariate optimization
 *problems, line search algorithms are not very good at performing
 *univariate optimization and functions like goldenSectionSearch would be
 *better for finding the minimum in a particular direction. On the other
 *hand, univariate optimization algorithms can usually be used for
 *performing a line search, because with many objective function, one can
 *often use an initial step size of 1 as the maximum step size and then
 *find the minimum over the bounded line.
 *
 *The line search algorithms here are from the libLBFGS library. The Matlab
 *function handle is wrapped into a C function and passed as a function
 *handle to the library. The library can be obtained from
 *http://www.chokkan.org/software/liblbfgs/index.html
 *and is under the MIT license.
 *
 *Algorithm 3 is generally the best backtracking algorithm to use, and thus
 *is the best algorithm to use when an L1-norm parameter is involved.
 *Algorithm 0 is often the best algorithm to use when no L1-norm parameter
 *is involved.
 *
 *The algorithm can be compiled for use in Matlab using the 
 *CompileCLibraries function.
 *
 *The algorithm is run in Matlab using the command format
 *[xMin,fMin,gMin,alpha,exitCode]=lineSearch(f,x,d);
 *or with all parameters as
 *[xMin,fMin,gMin,alpha,exitCode]=lineSearch(f,x,d,algorithm,C,stepSizeInit,fTol,wolfeTol,xTol,minStep,maxStep,maxIter,l1NormRange);
 *
 *As an example, we will take the sample problem used in
 *goldenSectionSearch and turn it into a bivariate minimization problem
 *where a particular line over which we want to minimize is the same as in
 *the goldenSectionSearch example. This will let us see how a line search
 *algorithms helps minimize a function, but does not necessarily find the
 *global minimum.
 *
 * f=@(x)deal((x(1)+x(2)-3)*(x(1)+x(2))^3*(x(1)+x(2)-6)^4,... %The function
 *            [(-6+x(1)+x(2))^3*(x(1)+x(2))^2*(54+8*x(1)^2+x(2)*(-45+8*x(2))+x(1)*(-45+16*x(2)));
 *            (-6+x(1)+x(2))^3*(x(1)+x(2))^2*(54+8*x(1)^2+x(2)*(-45+8*x(2))+x(1)*(-45+16*x(2)))]);%And the gradient as the second return.
 * %Note that the deal function is used to make an anonymous function have
 * %two outputs.
 * x=[4;0];
 * d=[1;0];
 * [xMin,fMin,gMin,alpha,exitCode]=lineSearch(f,x,d)
 *
 *In the above example, the starting point is x=[4;0] and only the x
 *coordinate changes due to the direction of d, so the y coordinate can be
 *effectively ignored. The local minimizing point found is
 *x=[6.091159463293529;0] with a minimum function value of
 *fMin=0.048242338941319. However, this is not the absolute local minimum.
 *The absolute local minimum is at x=[6;0], where fMin=0. Thus, unlike
 *bounding a region and using goldenSectionSearch, this function decreases
 *the cost sufficient according to certain conditions, but it does not find
 *a local minimum.
 *
 *REFERENCES:
 *[1] J. J. Moré and D. J. Thuente, "Line search algorithms with
 *    guaranteed sufficient decrease," ACM Transactions on Mathematical
 *    Software, vol. 20, no. 3, pp. 286-307, Sep. 1994.
 *[2] L. Armijo, "Minimization of functions having Lipschitz continuous
 *    first partial derivatives," Pacific Journal of Mathematics, vol. 16,
 *    no. 1, pp. 1-3, Nov. 1966.
 *[3] P. Wolfe, "Convergence conditions for ascent methods," SIAM Review,
 *    vol. 11, no. 2, pp. 226-235, Apr. 1969.
 *[4] P. Wolfe, "Convergence conditions for ascent methods. II: Some
 *    corrections," SIAM Review, vol. 13, no. 2, pp. 185-188, Apr. 1971.
 *[5] J. Nocedal and S. Wright, Numerical Optimization, 2nd ed. New
 *    York: Springer, 2006.
 *[6] A. Galen and G. Jianfeng, "Scalable training of l1-regularized 
 *    loglinear models," in Proceedings of the 24th International
 *    Conference on Machine Learning, Corvallis, OR, 20-24 Jun. 2007, pp.
 *    33-40.
 *
 *July 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
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
 
//Prototype for the callback function wrapper.
static double MatlabCallback(void *MatlabFunctionHandle,
            const double *x,
            double *gVal,
            const int n,
            const double step
            );

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    size_t xDim;
    //This will contain the information needed for the callback to the
    //Matlab function for the function value and gradient. 
    callback_data_t callbackData;
    double *x, *gVal;//The estimate and the gradient.
    double *xp, *gp;//Working space
    double *d;//A vector in the descent (search) direction
    double fVal;
    lbfgs_parameter_t param;
    double stepSize;
    int exitCode;
    
    if(nrhs<3||nrhs>13){
        mexErrMsgTxt("Wrong number of inputs");
    }
    
    if(nlhs>5) {
        mexErrMsgTxt("Wrong number of outputs.");
    }
    
    //Check that a function handle was passed.
    if(!mxIsClass(prhs[0],"function_handle")) {
        mexErrMsgTxt("The first input must be a function handle.");
    }
    
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
    
    //Check that a valid descent direction was passed.
    checkRealDoubleArray(prhs[2]);
    if(mxGetM(prhs[2])!=xDim||mxGetN(prhs[2])!=1) {
        mexErrMsgTxt("The descent direction d has the wrong dimensionality."); 
    }

    //Make a function wrapper for the Matlab callback to get the function
    //value and gradient.
    callbackData.n=(int)xDim;
    callbackData.instance=(void *)prhs[0];//The Matlab callback function handle.
    callbackData.proc_evaluate=&MatlabCallback;
    callbackData.proc_progress=NULL;
    
    //Allocate space for the state, the gradient, and the descent
    //direction. These should be allocated using mxCalloc so that the
    //function MatlabCallback does not have to needlessly copy data.
    x=mxCalloc(xDim, sizeof(double));
    gVal=mxCalloc(xDim, sizeof(double));
    d=mxCalloc(xDim,sizeof(double));
    //Allocate working space as well. This will have to be initialized with
    //the initial value of x and the gradient.
    xp=mxCalloc(xDim, sizeof(double));
    gp=mxCalloc(xDim, sizeof(double));
    
    //Copy the passed initial estimate and descent direction
    memcpy(x, mxGetData(prhs[1]), sizeof(double)*xDim);
    memcpy(d, mxGetData(prhs[2]), sizeof(double)*xDim);
    
    //Determine the algorithm to use.
    if(nrhs>3&&!mxIsEmpty(prhs[3])) {
        int algorithm=getIntFromMatlab(prhs[3]);

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
                mexErrMsgTxt("Unknown algorithm specified.");      
        }  
    } else {
        //Use the More-Thuente algorithm by default.
        param.linesearch=LBFGS_LINESEARCH_MORETHUENTE;
    }
    
    if(nrhs>4&&!mxIsEmpty(prhs[4])) {
        param.orthantwise_c=getDoubleFromMatlab(prhs[4]);
    } else {
        param.orthantwise_c=0;
    }
    
    //Get the initial stepsize
    if(nrhs>5&&!mxIsEmpty(prhs[5])) {
        stepSize=getDoubleFromMatlab(prhs[5]);
    } else {//Use the default stepsize.
        stepSize=1;
    }

    //Zero unused parameters in the param structure.
    param.max_iterations=0;
    param.past=0;
    param.delta=0;
    param.epsilon=0;
    param.m=0;
    
    //Fill in the parameters in the param structures related to line
    //searches.
    if(nrhs>6&&!mxIsEmpty(prhs[6])) {
        param.ftol=getDoubleFromMatlab(prhs[6]);
    } else {
        param.ftol=1e-6;
    }
    
    if(nrhs>7&&!mxIsEmpty(prhs[7])) {
        param.wolfe=getDoubleFromMatlab(prhs[7]);
    } else {
        param.wolfe=0.9;
    }
    param.gtol=param.wolfe;
    
    if(nrhs>8&&!mxIsEmpty(prhs[8])) {
        param.xtol=getDoubleFromMatlab(prhs[8]);
    } else {
        param.xtol=1e-16;
    }
    
    if(nrhs>9&&!mxIsEmpty(prhs[9])) {
        param.min_step=getDoubleFromMatlab(prhs[9]);
    } else {
        param.min_step=1e-20;
    }
    
    if(nrhs>10&&!mxIsEmpty(prhs[10])) {
        param.max_step=getDoubleFromMatlab(prhs[10]);
    } else {
        param.max_step=1e20;
    }
    
    if(nrhs>11&&!mxIsEmpty(prhs[11])) {
        param.max_linesearch=getIntFromMatlab(prhs[11]);
        
        if(param.max_linesearch<0) {
           mexErrMsgTxt("param.max_linesearch must be positive."); 
        }
    } else {
        param.max_linesearch=20;
    }
    
    if(nrhs>12&&!mxIsEmpty(prhs[12])) {
        double *indices, indexMin, indexMax;
        
        checkRealDoubleArray(prhs[12]);
        if(mxGetM(prhs[12])*mxGetN(prhs[12])!=2) {
            mexErrMsgTxt("The size of l1NormRange is incorrect."); 
        }
        
        indices=mxGetData(prhs[12]);
        
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
    } else {
        param.orthantwise_start=0;
        param.orthantwise_end=(int)xDim-1;
    }
    
    //Check for the validity of the parameters used.
    
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
    
    //Fill in the initial value of the function and the gradient.
    fVal=MatlabCallback(callbackData.instance,x,gVal,(int)xDim,0);
    //Put the current position and gradient vectors into the temporary
    //storage space.
    memcpy(xp, x, sizeof(double)*xDim);
    memcpy(gp, gVal, sizeof(double)*xDim);
        
    switch(param.linesearch){
        case LBFGS_LINESEARCH_MORETHUENTE:
            exitCode=line_search_morethuente((int)xDim,//The size of x.
                     x,//The current (returned final) estimate.
                     &fVal,//The returned final function value.
                     gVal,//The gradient of the function at x.
                     d,//A vector in the search (descent) direction. 
                     &stepSize,//The initial value of the stepsize
                     xp,//Working space the size of x (for x).
                     gp,//working space the size of g (for the gradient).
                     NULL,//Unused. Just pass NULL.
                     &callbackData,//Callback to get the function/ gradient.
                     &param//structure of parameter values.
                     );
            break;
        case LBFGS_LINESEARCH_BACKTRACKING_ARMIJO:
        case LBFGS_LINESEARCH_BACKTRACKING_WOLFE:
        case LBFGS_LINESEARCH_BACKTRACKING_STRONG_WOLFE:
            //If there is no l1 norm optimization term.
            if(param.orthantwise_c==0) {
                exitCode=line_search_backtracking((int)xDim,//The size of x
                         x,//The current (returned final) estimate.
                         //The returned final function value.
                         &fVal,
                         gVal,//The gradient of the function at x.
                         d,//A vector in the search (descent) direction. 
                         &stepSize,//The initial value of the stepsize
                         xp,//Working space the size of x (for x).
                         //Working space the size of g (for the gradient).
                         gp,
                         NULL,//Unused. Just pass NULL.
                         //Callback to get the function/ gradient.
                         &callbackData,
                         &param//structure of parameter values.
                         );
            } else {
                //Allocate space for the orthant.
                double *wp=mxCalloc(xDim, sizeof(double));

                exitCode=line_search_backtracking_owlqn((int)xDim,//The size of x
                         x,//The current (returned final) estimate.
                         //The returned final function value.
                         &fVal,
                         gVal,//The gradient of the function at x.
                         d,//A vector in the search (descent) direction. 
                         &stepSize,//The initial value of the stepsize
                         xp,//Working space the size of x (for x).
                         //Working space the size of g (for the gradient).
                         gp,
                         //Working space for the orthant.
                         wp,
                         //Callback to get the function/ gradient.
                         &callbackData,
                         &param//structure of parameter values.
                         );
                
                //Free scratch space.
                mxFree(wp);
            }
            break;
        default:
            mexErrMsgTxt("Error in choosing the proper algorithm.");      
    }
    //Free scratch space
    mxFree(xp);
    mxFree(gp);
    mxFree(d);
    
    //Set the return values.
    //First set x.
    plhs[0]=mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
    mxFree(mxGetPr(plhs[0]));
    mxSetPr(plhs[0], x);
    mxSetM(plhs[0], xDim);
    mxSetN(plhs[0], 1);
    
    if(nlhs>1) {
        //Return the minimum function value.
        plhs[1]=doubleMat2Matlab(&fVal,1,1);
        
        if(nlhs>2) {
            //Return the gradient vector.
            plhs[2]=mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
            mxFree(mxGetPr(plhs[2]));
            mxSetPr(plhs[2], gVal);
            mxSetM(plhs[2], xDim);
            mxSetN(plhs[2], 1);
            
            if(nlhs>3) {
                //Return the size of the step.
                plhs[3]=doubleMat2Matlab(&stepSize,1,1);

                if(nlhs>4) {
                    //Return the error code
                    plhs[4]=intMat2MatlabDoubles(&exitCode,1,1);
                }
            }
        } else {
            mxFree(gVal);
        }
    } else {
       mxFree(gVal); 
    }
}

static double MatlabCallback(void *MatlabFunctionHandle,
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
    rhs[0]=(mxArray*)MatlabFunctionHandle;
    rhs[1]=mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);

    //Set the matrix data to x 
    oldPtr=mxGetPr(rhs[1]);
    //x will not be modified, but the const must be typecast away to use
    //the mxSetPr function.
    mxSetPr(rhs[1], (double*)x);
    mxSetM(rhs[1], (size_t)n);
    mxSetN(rhs[1], 1);

    //Get the function value and gradient.
    mexCallMATLAB(2,lhs,2,rhs,"feval");

    //Get the function value.
    fVal=getDoubleFromMatlab(lhs[0]);

    //Copy the gradient into gVal, checking for errors.
    verifySizeReal((size_t)n,1,lhs[1]);
    memcpy(gVal, mxGetData(lhs[1]), sizeof(double)*(size_t)n);

    //Set the data pointer back to what it was during allocation that
    //mxDestroyArray does not have a problem. 
    mxSetPr(rhs[1],oldPtr);
    mxSetM(rhs[1], 0);
    mxSetN(rhs[1], 0);

    //Get rid of the temporary natrix.
    mxDestroyArray(rhs[1]);

    return fVal;
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
