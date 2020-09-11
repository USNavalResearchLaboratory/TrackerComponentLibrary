function [xMin,fMin,exitCode]=quasiNewtonLBFGS(f,x0,numCorr,epsilon,deltaTestDist,delta,lineSearchParams,maxIterations,progressCallback)
%%QUASINEWTONLBFGS Perform unconstrained nonlinear optimization using the
%                  limited-memory Broyden-Fletcher-Goldfarb-Shanno (BFGS)
%                  algorithm. The algorithm performs unconstrained
%                  minimization of a nonlinear function without one having
%                  to provide a Hessian matrix. It is a type of
%                  Quasi-Newton algorithm. It minimizes a function f(x),
%                  where x is a vector. If the C option in the line-search
%                  routine is not zero, then the function being minimized
%                  is f(x)+C*norm(x,1), where C>0 (or the norm can be
%                  taken using a subset of the elements of x), in which
%                  case this is the orthant-wise limited-memory
%                  Quasi-Newton method. Unlike the traditional BFGS
%                  algorithm, this algorithm does not store the entire
%                  inverse Hessian matrix estimate, thus making it
%                  significantly more efficient for very large problems.
%                  If one wishes to zero a vector (not a scalar), then
%                  NewtonsMethod is more appropriate as this function
%                  assumes the Hessian matrix is symmetric.
%
%INPUTS: f A handle to the function (and its gradient) over which the
%          minimization is to be performed. The function [fVal,gVal]=f(x)
%          takes the NX1 x vector returns the real scalar function value
%          fVal and gradient gVal at the point x. 
%       x0 The NX1-dimensional point from which the minimization starts.
%  numCorr The number of corrections to approximate the inverse Hessian
%          matrix. The default if omitted or an empty matrix is passed is
%          6. The L-BFGS documentation recommends not using fewer than 3.
%  epsilon The parameter determining the accuracy of the desired solution.
%          The function terminates when norm(g) < epsilon*max([1, norm(x)])
%          where g is the gradient. The default if omitted or an empty
%          matrix is passed is 1e-6.
% deltaTestDist The number of iterations back to use to compute the
%          decrease of the objective function if a delta-based convergence
%          test is performed. If zero, then no delta-based convergence
%          testing is done. The default if omitted or an empty matrix is
%          passed is zero.
%    delta The delta for the delta convergence test. This determines the
%          minimum rate of decrease of the objective function. Convergence
%          is determined if (f'-f)/f<delta, where f' is the value of the
%          objective function f deltaTestDist iterations ago, and f is the
%          current objective function value. The default if this parameter
%          is omitted or an empty matrix is passed is 0.
% lineSearchParams An optional structure whose members specify tolerances
%          for the line search. The parameters are described as in the
%          lineSearch function. Possible members are algorithm, C, fTol,
%          wolfeTol, xTol, minStep, maxStep, maxIter, l1NormRange. If this
%          parameter is omitted or an empty matrix is passed, then the
%          default values as described in the lineSearch function are
%          used. If any member of this structure is omitted or is assigned
%          an empty matrix, then the default value will be used.
% maxIterations The maximum number of iterations to use for the overall
%          L-BFGS algorithm. The default if this parameter is omitted or
%          an empty matrix is past is 1000. If this parameter is set to
%          zero, the algorithm will continue either until convergence is
%          obtained or until an error occurs.
% progressCallback An optional Matlab function handle that is called
%          during each step of the algorithm with information on the
%          progress. This can be used to analyze how well the optimization
%          is working. The function takes 8 arguments and returns one
%          argument. The format is
%          progressCallback(x,...%The state
%                           g,...%The gradient
%                           fx,...%The function value
%                           xnorm,...%The l2 norm of the state
%                           gnorm,...%The l2 norm of the gradient
%                                 ...%The line-search step used for this
%                                 ...%iteration
%                           step,...
%                           k,...%The iteration number.
%                             ...%The number of function evaluations
%                             ...%called for this iteration.
%                           ls);
%          If the return value is 0, then the optimization process will
%          continue. If the return value is nonzero, then the optimization
%          process will not continue. The return value should be a real
%          scalar value. If this parameter is omitted or an empty
%          matrix is passed, then no callback is performed.
%
%OUTPUTS: xMin The value of x at the minimum point found. If exitCode is
%              negative, then this value might be invalid.
%         fMin The cost function value at the minimum point found. If
%              exitCode is negative, then this value might be invalid.
%     exitCode A value indicating the termination condition of the
%              algorithm. Negative values indicate errors Possible
%              values are:
%                  0 The algorithm terminated successfully.
%                  1 Termination according to the delta stopping
%                    criterion occurred.
%                  2 The initial value already minimizes the objective
%                    function.
%              -1023 A logical error in the code occurred.
%              -1022 Insufficient memory.
%              -1021 The optimization was cancelled by the user callback
%                    progress function returning a nonzero value.
%              -1001 A finite precision error occurred or no line-search
%                    step satisfies the sufficient decrease and curvature
%                    conditions.
%              -1000 The line-search step size became less than minStep.
%               -999 The line-search step size became larger than maxStep.
%               -998 The maximum number of line-search iterations was
%                    reached.
%               -997 The maximum number of overall iterations was reached.
%               -996 The relative width of the interval of uncertainty is
%                    at most xTol
%               -995 A negative line-search step occurred.
%               -994 The current search direction increases the objective
%                    function.
%
%This is a Matlab interface for the C implementation of the L-BFGS
%algorithm of
%https://github.com/chokkan/liblbfgs
%which is based on the Fortran L-BFGS algorithm of
%http://www.ece.northwestern.edu/~nocedal/lbfgs.html
%but extends the work. The basic algorithm is described in [1] and [2].
%However, the C implementation contains additional controls over the type
%of line search perfromed and the convergence conditions. The above BFGS
%library has been slightly modified to use Matlab's memory allocation and
%deallocation routines.
%
%The algorithm can be compiled for use in Matlab  using the 
%CompileCLibraries function.
%
%The algorithm is run in Matlab using the command format
%[xMin,fMin,exitCode]=quasiNewtonLBFGS(f,x0);
%or if more options are used
%[xMin,fMin,exitCode]=quasiNewtonLBFGS(f,x0,numCorr,epsilon,deltaTestDist,delta,lineSearchParams,maxIterations,progressCallback);
%
%Examples:
%The first example is that used in the lineSearch file. 
% f=@(x)deal((x(1)+x(2)-3)*(x(1)+x(2))^3*(x(1)+x(2)-6)^4,... %The function
%            [(-6+x(1)+x(2))^3*(x(1)+x(2))^2*(54+8*x(1)^2+x(2)*(-45+8*x(2))+x(1)*(-45+16*x(2)));
%            (-6+x(1)+x(2))^3*(x(1)+x(2))^2*(54+8*x(1)^2+x(2)*(-45+8*x(2))+x(1)*(-45+16*x(2)))]);%And the gradient as the second return.
% %Note that the deal function is used to make an anonymous function have
% %two outputs.
% x0=[0.5;0.25];
% [xMin,fMin,exitCode]=quasiNewtonLBFGS(f,x0)
%The optimum point found is such that sum(xMin) is approximately
%1.73539450 with a minimum function value of approximately -2.1860756.
%
%The second example requires that the cost function be placed in a
%separate file. The cost function is
% function [fx,g]=objFun(x)
%    n=length(x);
%    fx=0;
%   
%    g=zeros(n,1);
%    for i=0:(n-1)
%        if(mod(i,2)==1)
%            continue;
%        end
%        t1=1-x(i+1);
%        t2=10*(x(i+1+1)-x(i+1)^2);
%        g(i+1+1)=20*t2;
%        g(i+1)=-2*(x(i+1)*g(i+1+1)+t1);
%        fx =fx+t1^2+t2^2;
%    end
% end
%It is used as
% n=100;
% x0=zeros(n,1);
% for i=0:(n-1)
%    if(mod(i,2)==1)
%        continue;
%    end
%
%    x0(i+1)=-1.2;
%    x0(i+1+1)=1;
% end
% [xMin,fMin,exitCode]=quasiNewtonLBFGS(f,x0)
%whereby the optimal solution is all ones with a minimum function value of
%zero. This second example is the same as that provided with the L-BFGS
%library in C.
%
%REFERENCES:
%[1] Liu, D. C.; Nocedal, J. (1989). "On the Limited Memory Method for
%    Large Scale Optimization". Mathematical Programming B 45 (3):
%    503-528. doi:10.1007/BF01589116.
%[2] Byrd, Richard H.; Lu, Peihuang; Nocedal, Jorge; Zhu, Ciyou (1995).
%    "A Limited Memory Algorithm for Bound Constrained Optimization".
%    SIAM Journal on Scientific and Statistical Computing 16 (5):
%    1190-1208. doi:10.1137/0916069.
%[3] J. J. Mor� and D. J. Thuente, "Line search algorithms with
%    guaranteed sufficient decrease," ACM Transactions on Mathematical
%    Software, vol. 20, no. 3, pp. 286-307, Sep. 1994.
%
%August 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
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
