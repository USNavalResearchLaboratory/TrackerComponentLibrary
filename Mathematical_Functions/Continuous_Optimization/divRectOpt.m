function [x,fVal,exitCode]=divRectOpt(f,lowerBounds,upperBounds,options)
%%DIVRECTOPT The dividing rectangles optimization algorithm. This is a
%           derivative-free global optimization algorithm where the
%           multidimensional search space can be bounded. This algorithm
%           tries to minimize a given function of a vector parameter.
%
%INPUTS: f A handle to the function over which the minimization is to be
%          performed. The function [fVal,violatesConstraints]=f(x) takes
%          the NX1 xvector and returns the real scalar function value
%          fVal. The parameter violatesConstraints is used to mark
%          forbidden areas, for example, where the function might not be
%          defined. violatesConstraints=false if there are not constraint
%          violations and it is true if there is a constraint violation.
% lowerBounds,upperBounds NX1 vectors of lower bounds defining lower and
%          upper bounds on the components of x at the optimal points. This
%          defines the search space. The elements must be finite and
%          cannot be NaNs.
%  options A structure where all elements are optional. The elements set
%          options that affect the algorithm and its convergence. Default
%          values are used if a field is missing. Possible fields in the
%          structure are:
%          algorithm (default 0) This selects which of the two direct
%                    optimization algorithms should be posed. Possible
%                    values are
%                    0: Use the algorithm of [1], which is probably
%                       better for functions with many local minima.
%                    1: Use the algorithm of [2], which favors local
%                       search and thus might be better for functions
%                       with few local minima.
%             maxDiv (default 5000) The maximum number of divisions of
%                    the hyperrectangles allowed. This must be >0.
%            maxIter (default 100) The maximum number of iterations.
%           maxFEval (default 500) The maximum number of function
%                     evaluations.
%            epsilon (default 1e-4) The Jones' factor. It is the epsilon
%                    term in Definitions 3.1 and 4.1 in [1]. It affects
%                    the convergence rate of the algorithm by avoiding
%                    oversampling near points with low function values
%                    and biasing the sampling towards a global search.
%                    If a positive value of epsilon is specified, then
%                    the same value is used for all iterations. In [1],
%                    it is suggested that epsilon be set to the desired
%                    solution accuracy (in f) for good convergence,
%                    though other convergence conditions are discussed in
%                    [2]. If a negative value is specified, then epsilon
%                    is adaptively changed using the formula
%                    epsilon = max(1.D-4*abs(fMin),abs(epsilon)).
%                    The magnitude of epislon must be between 0 and 1,
%                    not inclusive of 1, though 0 tends to work. 
%         epsilonAbs (default 0) Definitions 3.1 and 4.1 in [1] involve
%                    comparisons to fMin-epsilon*abs(fMin), where fMin is
%                    the current minimum function value found. If
%                    epsilonAbs is nonzero, then the comparison is
%                    actually taken with respect to
%                    min(fMin-epsilon*abs(fMin),epsilonAbs)
%                    epsilonAbs must be positive.
%       volumeRelTol (default 1e-10) Convergence is declared if the
%                    volume of the search region is this fraction of the
%                    original volume (must be between 0 and 1)
%        sigmaRelTol (default 0) Convergence is declared if the hypercube
%                    measure is smaller than this value. The notion of
%                    the hypercube measure sigma(S) is discussed in [2].
%            fGlobal If the actual minimum value of f is known, even
%                    though the value of x providing that minimum is
%                    unknown, then it can be specified. Otherwise, this
%                    parameter is not used.
%      fGlobalRelTol (default 1e-12 if fGlobal is provided) This is only
%                    used if fGlobal is provided. Convergence is declared
%                    if
%                    (fMin - fGlobal)/max(1,abs(fGlobal)) < fGlobalRelTol
%
%OUTPUTS: x The NX1 optimal point found.
%      fVal The value of the function at the optimal point found.
%  exitCode A return value indicating the status of the algorithm on
%          termination. Possible values are:
%          1 Function terminated, because  the number of function
%            evaluations done is larger than maxFEval.
%          2 Function terminated, because the number of iterations equals
%            maxIter.
%          3 The function value found is within fGlobalRelTol of the
%            global optimum given in fGlobal.
%          4 The volume of the hyperrectangle with fVal at its center is
%            less than volumeRelTol times the volume of the original
%            hyperrectangle.
%          5 The measure of the hyperrectangle with fVal at its center is
%            less than sigmaRelTol.
%
%Global optimization algorithms such as this are very good at handling
%functions with many local minima. Often, one might want to use a global
%optimization algorithm to get an initial estimate close to the global
%solution and then do a few iterations of, for example, Newton's method,
%to find the global optimum.
%
%This function is an interface to the DIRECT library that is part of the
%NLOpt set of nonlinear optimization routines available at [4], which
%implements [1] and [2].  Though the aforementioned code online provides
%a Matlab interface, it was decided to include the dividing rectangles
%algorithm without the rest of the library, in part because some of the
%other functions include weak copyleft provisions in the form of the LGPL.
%Additionally, the implementation above did not allow one to specify the
%maximum number of hyperrectangle divisions (it was statically fixed) nor
%did it use the memory allocation routines in Matlab that mex files are
%supposed to (but do not have to) use. Thus, the code was slightly
%modified (but is still included with the rest of the third party code).
%
%Though [1] and [2] use a number of exmaple, functions, they formulae are
%not explicitely given. Rather, one must consult a few references, which
%can be tedious to obtain. On the other hand, numrous exmaples functions
%used in [1] are given explicitely in [3]. Examples of the algorithm are
%thus:
% %Branin's RCOS function
% %Note that the deal function is used to make an anonymous function have
% %two outputs.
% f=@(x)deal((x(2)-(5/(4*pi^2))*x(1)^2+(5/pi)*x(1)-6)^2+10*(1-(1/(8*pi)))*cos(x(1))+10,0);
% lowerBounds=[-5;10];
% upperBounds=[0;15];
% [x,fVal,exitCode]=divRectOpt(f,lowerBounds,upperBounds)
%The globally optimal value is given in [3] as 0.397887357729739. The
%returned fVal is close.
%
% %Six-Hump Camel
% f=@(x)deal((x(1)^2*(4-2.1*x(1)^2+x(1)^4/3)+x(1)*x(2)+x(2)^2*(-4+4*x(2)^2)),0);
% lowerBounds=[-3;-2];
% upperBounds=[3;2];
% [x,fVal,exitCode]=divRectOpt(f,lowerBounds,upperBounds)
%The globally optimal value is given in [3] as -1.0316284535. The returned
%fVal is close.
%
% %Two-Dimensional Schubert
% f=@(x)deal(sum((1:5).*cos(((1:5)+1)*x(1)+(1:5)))*sum((1:5).*cos(((1:5)+1)*x(2)+(1:5))),0)
% lowerBounds=[-10;-10];
% upperBounds=[10;10];
% options.maxFEval=5000;
% options.maxIter=200;
% [x,fVal,exitCode]=divRectOpt(f,lowerBounds,upperBounds,options)
%Here, increasing the number of function evaluations and iterations is
%important to getting a solution that is close to the optimal value as
%this function has 760 local minima (and 18 global minima). The global
%minimum is -186.730908831024 as per [3].
%
%REFERENCES:
%[1] D. R. Jones, C. D. Peritunen, and B. E. Stuckman, "Lipschitzian
%    optimization without the Lipschitz constant," Journal of Optimization
%    Theory and Application, vol. 79, no. 1, pp. 157-181, Oct. 1993.
%[2] J. M. Gablonsky and C. T. Kelley, "A locally-biased form of the
%    DIRECT algorithm," Journal of Global Optimization, vol. 21, no. 1,
%    pp. 27-37, Sep. 2001.
%[3] M. Bj�kman and K. Holmstr�m, "Global optimization using the DIRECT
%    algorithm in Matlab," The Electronic International Journal Advanced
%    Modeling and Optimization, vol. 1, no. 2, pp. 17-37, 1999.
%    [Online]. Available: http://camo.ici.ro/journal/v1n2.htm
%[4] S. G. Johnson. (2014, 20 May) The NLopt nonlinear-optimization
%    package. [Online]. Available: http://ab-initio.mit.edu/nlopt
%
%October 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

error('This function is only implemented as a mexed C or C++ function. Please run CompileCLibraries.m to compile the function for use.')

end

%LICENSE:
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
%OF RECIPIENT IN THE USE OF THE SOFTWARE.
