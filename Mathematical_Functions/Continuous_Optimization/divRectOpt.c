/**DIVRECTOPT The dividing rectangles optimization algorithm. This is a
 *           derivative-free global optimization algorithm where the
 *           multidimensional search space can be bounded. This algorithm
 *           tries to minimize a given function of a vector parameter.
 *
 *INPUTS: f A handle to the function over which the minimization is to be
 *          performed. The function [fVal,violatesConstraints]=f(x) takes
 *          the NX1 x vector and returns the real scalar function value
 *          fVal. The parameter violatesConstraints is used to mark
 *          forbidden areas, for example, where the function might not be
 *          defined. violatesConstraints=false if there are not constraint
 *          violations and it is true if there is a constraint violation.
 * lowerBounds,upperBounds NX1 vectors of lower bounds defining lower and
 *          upper bounds on the components of x at the optimal points. This
 *          defines the search space. The elements must be finite and
 *          cannot be NaNs.
 *  options A structure where all elements are optional. The elements set
 *          options that affect the algorithm and its convergence. Default
 *          values are used if a field is missing. Possible fields in the
 *          structure are:
 *          algorithm (default 0) This selects which of the two direct
 *                    optimization algorithms should be posed. Possible
 *                    values are
 *                    0: Use the algorithm of [1], which is probably
 *                       better for functions with many local minima.
 *                    1: Use the algorithm of [2], which favors local
 *                       search and thus might be better for functions
 *                       with few local minima.
 *             maxDiv (default 5000) The maximum number of divisions of
 *                    the hyperrectangles allowed. This must be >0.
 *            maxIter (default 100) The maximum number of iterations.
 *           maxFEval (default 500) The maximum number of function
 *                     evaluations.
 *            epsilon (default 1e-4) The Jones' factor. It is the epsilon
 *                    term in Definitions 3.1 and 4.1 in [1]. It affects
 *                    the convergence rate of the algorithm by avoiding
 *                    oversampling near points with low function values
 *                    and biasing the sampling towards a global search.
 *                    If a positive value of epsilon is specified, then
 *                    the same value is used for all iterations. In [1],
 *                    it is suggested that epsilon be set to the desired
 *                    solution accuracy (in f) for good convergence,
 *                    though other convergence conditions are discussed in
 *                    [2]. If a negative value is specified, then epsilon
 *                    is adaptively changed using the formula
 *                    epsilon = max(1.D-4*abs(fMin),abs(epsilon)).
 *                    The magnitude of epislon must be between 0 and 1,
 *                    not inclusive of 1, though 0 tends to work. 
 *         epsilonAbs (default 0) Definitions 3.1 and 4.1 in [1] involve
 *                    comparisons to fMin-epsilon*abs(fMin), where fMin is
 *                    the current minimum function value found. If
 *                    epsilonAbs is nonzero, then the comparison is
 *                    actually taken with respect to
 *                    min(fMin-epsilon*abs(fMin),epsilonAbs)
 *                    epsilonAbs must be positive.
 *       volumeRelTol (default 1e-10) Convergence is declared if the
 *                    volume of the search region is this fraction of the
 *                    original volume (must be between 0 and 1)
 *        sigmaRelTol (default 0) Convergence is declared if the hypercube
 *                    measure is smaller than this value. The notion of
 *                    the hypercube measure sigma(S) is discussed in [2].
 *            fGlobal If the actual minimum value of f is known, even
 *                    though the value of x providing that minimum is
 *                    unknown, then it can be specified. Otherwise, this
 *                    parameter is not used.
 *      fGlobalRelTol (default 1e-12 if fGlobal is provided) This is only
 *                    used if fGlobal is provided. Convergence is declared
 *                    if
 *                    (fMin - fGlobal)/max(1,abs(fGlobal)) < fGlobalRelTol
 *
 *OUTPUTS: x The NX1 optimal point found.
 *      fVal The value of the function at the optimal point found.
 *  exitCode A return value indicating the status of the algorithm on
 *          termination. Possible values are:
 *          1 Function terminated, because  the number of function
 *            evaluations done is larger than maxFEval.
 *          2 Function terminated, because the number of iterations equals
 *            maxIter.
 *          3 The function value found is within fGlobalRelTol of the
 *            global optimum given in fGlobal.
 *          4 The volume of the hyperrectangle with fVal at its center is
 *            less than volumeRelTol times the volume of the original
 *            hyperrectangle.
 *          5 The measure of the hyperrectangle with fVal at its center is
 *            less than sigmaRelTol.
 *
 *Global optimization algorithms such as this are very good at handling
 *functions with many local minima. Often, one might want to use a global
 *optimization algorithm to get an initial estimate close to the global
 *solution and then do a few iterations of, for example, Newton's method,
 *to find the global optimum.
 *
 *This function is an interface to the DIRECT library that is part of the
 *NLOpt set of nonlinear optimization routines available at [4], which
 *implements [1] and [2].  Though the aforementioned code online provides
 *a Matlab interface, it was decided to include the dividing rectangles
 *algorithm without the rest of the library, in part because some of the
 *other functions include weak copyleft provisions in the form of the LGPL.
 *Additionally, the implementation above did not allow one to specify the
 *maximum number of hyperrectangle divisions (it was statically fixed) nor
 *did it use the memory allocation routines in Matlab that mex files are
 *supposed to (but do not have to) use. Thus, the code was slightly
 *modified (but is still included with the rest of the third party code).
 *
 *Though [1] and [2] use a number of exmaple, functions, they formulae are
 *not explicitely given. Rather, one must consult a few references, which
 *can be tedious to obtain. On the other hand, numrous exmaples functions
 *used in [1] are given explicitely in [3]. Examples of the algorithm are
 *thus:
 * %Branin's RCOS function
 * %Note that the deal function is used to make an anonymous function have
 * %two outputs.
 * f=@(x)deal((x(2)-(5/(4*pi^2))*x(1)^2+(5/pi)*x(1)-6)^2+10*(1-(1/(8*pi)))*cos(x(1))+10,0);
 * lowerBounds=[-5;10];
 * upperBounds=[0;15];
 * [x,fVal,exitCode]=divRectOpt(f,lowerBounds,upperBounds)
 *The globally optimal value is given in [3] as 0.397887357729739. The
 *returned fVal is close.
 *
 * %Six-Hump Camel
 * f=@(x)deal((x(1)^2*(4-2.1*x(1)^2+x(1)^4/3)+x(1)*x(2)+x(2)^2*(-4+4*x(2)^2)),0);
 * lowerBounds=[-3;-2];
 * upperBounds=[3;2];
 * [x,fVal,exitCode]=divRectOpt(f,lowerBounds,upperBounds)
 *The globally optimal value is given in [3] as -1.0316284535. The returned
 *fVal is close.
 *
 * %Two-Dimensional Schubert
 * f=@(x)deal(sum((1:5).*cos(((1:5)+1)*x(1)+(1:5)))*sum((1:5).*cos(((1:5)+1)*x(2)+(1:5))),0)
 * lowerBounds=[-10;-10];
 * upperBounds=[10;10];
 * options.maxFEval=5000;
 * options.maxIter=200;
 * [x,fVal,exitCode]=divRectOpt(f,lowerBounds,upperBounds,options)
 *Here, increasing the number of function evaluations and iterations is
 *important to getting a solution that is close to the optimal value as
 *this function has 760 local minima (and 18 global minima). The global
 *minimum is -186.730908831024 as per [3].
 *
 *REFERENCES:
 *[1] D. R. Jones, C. D. Peritunen, and B. E. Stuckman, "Lipschitzian
 *    optimization without the Lipschitz constant," Journal of Optimization
 *    Theory and Application, vol. 79, no. 1, pp. 157-181, Oct. 1993.
 *[2] J. M. Gablonsky and C. T. Kelley, "A locally-biased form of the
 *    DIRECT algorithm," Journal of Global Optimization, vol. 21, no. 1,
 *    pp. 27-37, Sep. 2001.
 *[3] M. Bj�kman and K. Holmstr�m, "Global optimization using the DIRECT
 *    algorithm in Matlab," The Electronic International Journal Advanced
 *    Modeling and Optimization, vol. 1, no. 2, pp. 17-37, 1999.
 *    [Online]. Available: http://camo.ici.ro/journal/v1n2.htm
 *[4] S. G. Johnson. (2014, 20 May) The NLopt nonlinear-optimization
 *    package. [Online]. Available: http://ab-initio.mit.edu/nlopt
 *
 *October 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

/*This header is required by Matlab.*/
#include "mex.h"

//For the dividing rectangle optimization routines
#include "direct.h"
#include "MexValidation.h"
#include <limits.h>

static double MatlabCallback(int n,//The dimensionality of the point.
                      const double *x,//The array of the point to evaluate.
                      int *undefined_flag,//Set to 1 on return if x violates constraints; otherwise unused.
                      void *data);//Any data that the user passed as f_data for the function.

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    const mxArray *MatlabFunctionHandle;
    direct_algorithm theAlgorithm=DIRECT_ORIGINAL;
    const direct_objective_func f=&MatlabCallback;
    size_t numDim;
    double *lowerBounds;
    double *upperBounds;
    int maxFEval=500;
    int maxIter=100;
    double epsilon=1e-4;
    double epsilonAbs=0;
    double volumeRelTol=1e-10;
    double sigmaRelTol=0;
    double fGlobal=DIRECT_UNKNOWN_FGLOBAL;
    double fGlobalRelTol=DIRECT_UNKNOWN_FGLOBAL;
    int maxDiv=5000;
    //To hold return values
    direct_return_code retVal;
    double *x;
    double fVal;
    
    if(nrhs<3||nrhs>4){
        mexErrMsgTxt("Wrong number of inputs");
    }

    if(nlhs>3) {
        mexErrMsgTxt("Wrong number of outputs.");
    }
    
    //Check that a function handle was passed.
    if(!mxIsClass(prhs[0],"function_handle")) {
        mexErrMsgTxt("The first input must be a function handle.");
    }
    MatlabFunctionHandle=prhs[0];
    
    //Check the lower and upper bounds
    checkRealDoubleArray(prhs[1]);
    checkRealDoubleArray(prhs[2]);
    numDim=mxGetM(prhs[1]);
    if(numDim<1||mxGetN(prhs[1])!=1) {
        mexErrMsgTxt("The lower bounds have the wrong dimensionality.");
    }
    if(mxGetM(prhs[2])!=numDim||mxGetN(prhs[1])!=1) {
        mexErrMsgTxt("The upper bounds have the wrong dimensionality.");
    }
    
    //The function in the library only uses integers for dimensions.
    if(numDim>INT_MAX) {
        mexErrMsgTxt("The problem has too many dimensions to be solved using this function.");
    }
    
    lowerBounds=mxGetDoubles(prhs[1]);
    upperBounds=mxGetDoubles(prhs[2]);
    
    //If a structure of options was given
    if(nrhs>3&&!mxIsEmpty(prhs[3])) {
        //We have to check for every possible structure member.
        mxArray *theField;
        
        if(!mxIsStruct(prhs[3])) {
            mexErrMsgTxt("The options must be given in a structure.");
        }
    
        theField=mxGetField(prhs[3],0,"algorithm");
        
        if(theField!=NULL) {//If the field is present.
            int algType=getIntFromMatlab(theField);
            
            //Algorithm type
            switch(algType) {
              case 0:
                  theAlgorithm=DIRECT_ORIGINAL;
                  break;
              case 1:
                  theAlgorithm=DIRECT_GABLONSKY;
                  break;
              default:
                  mexErrMsgTxt("Invalid algorithm specified");
            }
        }
        
        theField=mxGetField(prhs[3],0,"maxDiv");
        if(theField!=NULL) {//If the field is present.
            maxDiv=getIntFromMatlab(theField);
            if(maxDiv<1) {
                mexErrMsgTxt("Invalid maxDiv specified");
            }
        }
        
        theField=mxGetField(prhs[3],0,"maxIter");
        if(theField!=NULL) {//If the field is present.
            maxIter=getIntFromMatlab(theField);
            if(maxIter<1) {
                mexErrMsgTxt("Invalid maxIter specified");
            }
        }
        
        theField=mxGetField(prhs[3],0,"maxFEval");
        if(theField!=NULL) {//If the field is present.
            maxFEval=getIntFromMatlab(theField);
            if(maxIter<1) {
                mexErrMsgTxt("Invalid maxFEval specified");
            }
        }
        
        theField=mxGetField(prhs[3],0,"epsilon");
        if(theField!=NULL) {//If the field is present.
            epsilon=getDoubleFromMatlab(theField);
            if(fabs(epsilon)>=1) {
                mexErrMsgTxt("Invalid epsilon specified");
            }
        }
        
        theField=mxGetField(prhs[3],0,"epsilonAbs");
        if(theField!=NULL) {//If the field is present.
            epsilonAbs=getDoubleFromMatlab(theField);
            if(epsilonAbs>=1||epsilonAbs<0) {
                mexErrMsgTxt("Invalid epsilonAbs specified");
            }
        }
        
        theField=mxGetField(prhs[3],0,"volumeRelTol");
        if(theField!=NULL) {//If the field is present.
            volumeRelTol=getDoubleFromMatlab(theField);
            if(volumeRelTol>=1||volumeRelTol<0) {
                mexErrMsgTxt("Invalid volumeRelTol specified");
            }
        }
        
        theField=mxGetField(prhs[3],0,"sigmaRelTol");
        if(theField!=NULL) {//If the field is present.
            sigmaRelTol=getDoubleFromMatlab(theField);
            if(sigmaRelTol>=1||sigmaRelTol<0) {
                mexErrMsgTxt("Invalid sigmaRelTol specified");
            }
        }
        
        theField=mxGetField(prhs[3],0,"fGlobal");
        if(theField!=NULL) {//If the field is present.
            fGlobal=getDoubleFromMatlab(theField);
            fGlobalRelTol=1e-12;//Default value if fGlobal is given.
        }
        
        theField=mxGetField(prhs[3],0,"fGlobalRelTol");
        if(theField!=NULL) {//If the field is present.
            fGlobalRelTol=getDoubleFromMatlab(theField);
            if(fGlobalRelTol<0) {
                mexErrMsgTxt("Invalid fGlobalRelTol specified");
            }
        }      
    }
    
    //Allocate space for the return values
    x=mxCalloc(numDim, sizeof(double));

    retVal=direct_optimize(f,//Objective function pointer
            (void*)MatlabFunctionHandle,//Data that passed to the objective function. Here, the function handle passed by the user.
            (int)numDim,//The dimensionality of the problem
            lowerBounds,//Array of the lower bound of each dimension.
            upperBounds,//Array of the upper bound of each dimension.
            x,//An array that will be set to the optimum value on return.
            &fVal,//On return, set to the minimum value.
            maxFEval,//Maximum number of function evaluations.
            maxIter,//Maximum number of iterations
            epsilon,//Jones' epsilon parameter (1e-4 is recommended)
            epsilonAbs,//An absolute version of magic_eps
            //Relative tolerance on the hypercube volume (use 0 if none). This is
            //a percentage (0-1) of the volume of the original hypercube
            volumeRelTol,
            //Relative tolerance on the hypercube measure (0 if none)
            sigmaRelTol,
            //A pointer to a variable that can terminate the optimization. This is
            //in the library's code so that an external thread can write to this
            //and force the optimization to stop. Here, we are not using it, so we
            //can pass NULL.
            NULL,
            fGlobal,//Function value of the global optimum, if known. Use DIRECT_UNKNOWN_FGLOBAL if unknown.
            fGlobalRelTol,//Relative tolerance for convergence if fglobal is known. If unknown, use DIRECT_UNKNOWN_FGLOBAL.
            NULL,//There is no output to a file.
            theAlgorithm,//The algorithm to use
            maxDiv);//The maximum number of hyperrectangle divisions

    //If an error occurred
    if(retVal<0) {
        mxFree(x);
        switch(retVal){
            case DIRECT_INVALID_BOUNDS:
                mexErrMsgTxt("Upper bounds are <= lower bounds for one or more dimensions.");
            case DIRECT_MAXFEVAL_TOOBIG:
                mexErrMsgTxt("maxFEval is too large.");
            case DIRECT_INIT_FAILED:
                mexErrMsgTxt("An initialization step failed.");
            case DIRECT_SAMPLEPOINTS_FAILED:
                 mexErrMsgTxt("An error occurred creating sampling points.");
            case DIRECT_SAMPLE_FAILED:
                mexErrMsgTxt("An error occurred while sampling the function.");
            case DIRECT_OUT_OF_MEMORY:
                mexErrMsgTxt("Out of memory");
            case DIRECT_INVALID_ARGS:
                mexErrMsgTxt("Invalid arguments provided");
            //These errors should not occur either because they should
            //never be called or because retVal>=0 so these are not
            //actually errors.
            case DIRECT_FORCED_STOP:
            case DIRECT_MAXFEVAL_EXCEEDED:
            case DIRECT_MAXITER_EXCEEDED:
            case DIRECT_GLOBAL_FOUND:
            case DIRECT_VOLTOL:
            case DIRECT_SIGMATOL:
            default:
                mexErrMsgTxt("An unknown error occurred.");
            
        }
    }
        
    //Set the return values
    plhs[0]=mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
    mxFree(mxGetDoubles(plhs[0]));
    mxSetDoubles(plhs[0], x);
    mxSetM(plhs[0], numDim);
    mxSetN(plhs[0], 1);

    if(nlhs>1) {
        plhs[1]=doubleMat2Matlab(&fVal,1,1);
        if(nlhs>2) {
            plhs[2]=intMat2MatlabDoubles(&retVal,1,1);
        }
    }
}

static double MatlabCallback(int n, const double *x, int *undefined_flag, void *data) {
    mxArray *rhs[2];
    mxArray *lhs[2];
    double *oldPtr;
    double fVal;
    bool violatesConstraints;

    //feval in Matlab will take the function handle and the state as
    //inputs.
    //The first function handle is f. 
    rhs[0]=(mxArray*)data;
    rhs[1]=mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);

    //Set the matrix data to x 
    oldPtr=mxGetDoubles(rhs[1]);
    
    //x will not be modified, but the const must be typecast away to use
    //the mxSetDoubles function.    
    mxSetDoubles(rhs[1],(double*)x);
    mxSetM(rhs[1], (size_t)n);
    mxSetN(rhs[1], 1);
    
    //Get the function value and gradient.
    mexCallMATLAB(2,lhs,2,rhs,"feval");
    
    //Get the function value.
    fVal=getDoubleFromMatlab(lhs[0]);
      
    violatesConstraints=getBoolFromMatlab(lhs[1]);
    
    if(violatesConstraints) {
        *undefined_flag=1;
    }
    
    //Get rid of the returned Matlab Matrices.
    mxDestroyArray(lhs[0]);
    mxDestroyArray(lhs[1]);
    
    //Set the data pointer back to what it was during allocation so that
    //mxDestroyArray does not have a problem. 
    mxSetDoubles(rhs[1],oldPtr);
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
