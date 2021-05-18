function [xk,objfunvals,xkk] = solveLASSOProblem(A,b,alpha,varargin)
%%SOLVELASSOPROBLEM This function solves the Least Absolute Shrinkage and
%               Selection Operator (LASSO) problem. This optimization
%               problem is
%               minimize (1/2)*norm(A*x-b,2)^2+alpha*norm(x,1)
%               Thus, the LASSO problem is the same as the ridge regression
%               problem except the second term uses an l_1 norm instead of
%               an l_2 norm. This type of optimization plays a notable role
%               in compressive sensing algorithms. The optimization problem
%               is solved using the Fast-Iterative Shrinkage and
%               Thresholding Algorithm (FISTA) for sparse vectors.
%
%INPUTS:      A An mXn matrix. This can be complex.
%             b An mX1 vector. This can be complex.
%         alpha A real scalar value.
%      varargin Arguments that are specified as fieldname, value pairs.
%               Possible pairs are:
%               FIELD NAME         VALUE TYPE     DESCRIPTION
%               -----------------  -------------- -----------------------
%               'maxiter'          numeric scalar The maximum number of
%                                                 iterations (default 2000)
%               'tol'              numeric scalar Tolerence for termination
%                                                 (default 1e-6).
%               'tau'              numeric scalar Descent parameter used to
%                                                 to control the rate of
%                                                 minimization (default:
%                                                 svds(A,1).^2).
%               'backtracking'     bool           Performs backtracking if
%                                                 turned on.
%               'verbose'          bool           Enables verbose mode.
%               'stoppingcriteria' char           'cauchy' or 'gradient'.
%               'init'             vector         n x 1 vector for warm
%                                                 start.
%
%OUTPUTS: xk The nX1 vector solution to the optimization problem.
% objfunvals If verbose mode is activated, then this will hold the value of
%            the objective function at every iteration. in a numIterX1
%            vector. Otherwise, this is an empty matrix.
%        xkk An nXnumIter matrix holding the value of the estimate at each
%            step.
%
%The algorithm implemented is that of [1]. The algorithm typically needs a
%very large number of iterations to converge. However, the number of
%iterations is not particularly dependent on the scale of the problem,
%making the solver good for large problems.
%
%EXAMPLE:
%A trivial example is
% A=[10,  4,   -3, 6;
%     0,  -18, 24, -10;
%     5,  7,   3,  12];
% xTrue=[0;2;0;12];
% alpha=2;
% b=A*xTrue;
% xk=solveLASSOProblem(A,b,alpha)
%The xk found should have the required sparsity, but the l1 penalty term
%means that the xk found will not equal xTrue; Rather, it will be
%xk=[0;2.006;0;11.9922]. Increasing alpha forces the sparsity up. For
%example, using alpha=3200 changes the result to xk=[0;2.6581;0;0] and
%using alpha=4000 changes xk to all zeros. Making alpha close to 0 reduces
%the result to approach pinv(A)*b.
%
%REFERENCES:
%[1] A. Beck and M. Teboulle, "A fast iterative shrinkage-thresholding
%    algorithm for linear inverse problems," SIAM Journal on Imaging
%    Sciences, vol. 2, no. 1, pp. 183-202, 2009.
%
%February 2016 Hatim F. Alqadah, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    %% PROCESS INPUTS HERE
    if(nargin < 3)
        alpha = [];
    end
    if(nargin < 2)
        error('invalid number of arguments');
    end
    
    % check argument type and sizes
    assert(ismatrix(A),'A should be a matrix');
    assert(isvector(b),'b should be a vector');
    [m,n] = size(A);
    assert(size(b,1) == m,'invalid vector length on right hand side');
    
    % if matrix or data is complex we need to seperate out
    % the real and imaginary components
    is_complex = ~(isreal(A) && isreal(b));
    if(is_complex)
        % reformulate into a purely real problem
        A = [real(A) -imag(A);imag(A) real(A)];
        b = [real(b);imag(b)];
    end
    
    % set default options and process user defined ones
    opts.maxiter = 2000;
    opts.tol = 1e-6;
    opts.stopcrit = 'cauchy';
    opts.backtracking = false;
    opts.verbose = false;
    opts.restart = false;
    opts.return_as_dense = true;
    opts.init = [];
    opts.tau = [];
    
    [opts] = addCelListToStruct(opts,varargin);
    clear varargin
    
    if(~isempty(opts.init))
        assert(length(opts.x0) == n,'given warm start has incorrect length');
        xk = x0;
        if(is_complex)
            xk = [real(xk);imag(xk)];
        end
    else
        if(is_complex)
            xk = sparse(2*n,1);
        else
            xk = sparse(n,1);
        end
    end
    
    %% PREPROCESSING
    
    % set initial quantities
    xp = xk;
    nn = 0;
    tk = 1;
    tp = 1;
    At = A';
    Atb = At*b;
    has_converged = false;
    keep_history = nargout > 2;
    objfunvals = [];
    xkk = [];
    
    if(isempty(alpha))
        alpha = .5*norm(Atb,'inf');
    end
    
    if(isempty(opts.tau))
        %Use iterative method to find the Lipschitz constant
        opts.tau = svds(A,1).^2;
    end  
    
    %Setup reporting functions if quiet mode is disabled
    if(opts.verbose)
        fprintf('\n*******************************************\n');
        fprintf('* FISTA: SPARSE VECTOR RECONSTRUCTION \n');
        if(is_complex)
        	fprintf('* CMPLX     : %s\n','true');
        else
            fprintf('* CMPLX     : %s\n','false');
        end
        fprintf('* MAXITER   : %d\n',opts.maxiter);
        fprintf('* TOL       : %e\n',opts.tol);
        if(opts.restart)
            fprintf('* RESTART   : %s\n','true');
        else
            fprintf('* RESTART   : %s\n','false');
        end
        fprintf('* STOPCRIT  : %s\n',opts.stopcrit);
        fprintf('*******************************************\n');
        msgstr = ['Current Precision: ',num2str(inf)];
        h = waitbar(0,msgstr);   
        precision = zeros(opts.maxiter,1);
        objfun = @(x) ObjFunc(x,A,b,alpha);
        objfunvals = zeros(opts.maxiter,1);
    end
    if(keep_history)
        xkk = zeros(length(xk),opts.maxiter);
    end
    
    %% MAINLOOP
    
    while(nn < opts.maxiter)
        if opts.verbose
            objfunvals(nn+1) = objfun(xk); %#ok<*AGROW>
        end
        if(keep_history)
            xkk(:,nn+1) = xk;
        end
        
        yk = xk + ((tp-1)/tk)*(xk - xp);
        grd = At*(A*yk) - Atb;
        xp = xk;
        if(opts.backtracking)
            % solve the local regularization problem with an appropriate
            % step size
            [xk,opts.tau] = LineSearch(A,b,grd,alpha,yk,opts.tau);
        else
            xk = sparse(shrink(yk - (1/opts.tau)*grd,alpha/opts.tau));
        end

        % check stopping criterion
        switch(lower(opts.stopcrit))
            case 'cauchy'
                prec = norm(xk(:) - yk(:))/max(norm(xk),1);
            case 'subgradient'
                % we need to estimate the subgradient at current pt
                % optimality occurs at zero
                prec = ComputeOptimalityResidue(A,At,xk,Atb,alpha);
        end
        if(prec < opts.tol)
            has_converged = true;
            break;
        end
        
        % keep going
        nn = nn + 1;
        tp = tk;
        tk = 0.5*(1+sqrt(1+4*tk*tk));
        
        % check if a restart is necessary
        if(opts.restart && ((yk - xk)'*(xk - xp) > 0))
            % forget all previous iterations
            % and reset momentum back to zero
            tp = 1;
            tk = 1;
        end    
        if(opts.verbose)
            precision(nn) = prec;
            waitbar(nn/opts.maxiter,h,['Current precision: ',num2str(prec)]);
        end
    end
    if(~has_converged)
        warning('MATLAB:no_convergence','method failed to converge within given time steps');
    end
    
    %% Post-processing
    if(opts.verbose)
        waitbar(1,h);
        objfunvals = objfunvals(1:nn);
        precision = precision(1:nn);
        iter = 0:nn-1;
        figure;
        semilogy(iter,objfunvals,'--r','linewidth',1.2);
        grid on;
        title(['Objective Function for \alpha = ',num2str(alpha)]);
        xlabel('ITERATION');
        ylabel('FOBJ VAL');
        
        figure;
        semilogy(iter,precision,'--g','linewidth',1.2);
        grid on;
        title(['Precision for \alpha = ',num2str(alpha)]);
        xlabel('ITERATION');
        ylabel('PREC LVL');
        close(h);
    end
    if(opts.return_as_dense)
        xk = full(xk);
        if keep_history
            xkk = full(xkk(:,1:nn));
        end
    end
    if(is_complex)
        xk = complex(xk(1:n),xk(n+1:end));
        if(keep_history)
            xkk = complex(xkk(1:n,:),xkk(n+1:end,:));
        end
    end

end
%% UTILITY
% defines the standard shrinkage operator. 
function z = shrink(x,mu)
    z = sign(x).*max(abs(x)-mu,0);
end

% Experimental stopping criteria particlarly for proximal gradient methods.
% It is not clear if this can help issues do non-strict convexivity of the
% data fidelity term
function prec = ComputeOptimalityResidue(A,At,xk,Atb,alpha)
    % first compute the gradient at the current point
    gradF = At*(A*xk) - Atb;
    
    % compute the unique portion of the sub-gradient
    subGrad = alpha*sign(xk);
    
    % now for the points that are zero minimize the
    % optimality residue
    inds = find(subGrad == 0);
    
    % compute solution to the constrained
    % sub optimization problem
    g = gradF(inds);
    y = zeros(size(g));
    
    for i=1:50
        yp = y;
        lambda = abs(y+g)/alpha;
        y = (y - .5*(y+g))./(1 - .5*lambda);
        if(norm(y-yp)/max([norm(y),1]) < 1e-5)
            break;
        end
    end
    subGrad(inds) = y;
    
    % evaluate the maximum sub-gradient component
    prec = norm(gradF + subGrad,'inf');
end

% returns the objective function value of the unconstrained l2-l1 problem
function val = ObjFunc(x,A,b,alpha)
   val = .5*norm(A*x - b).^2 + alpha*norm(x,1);
end

% A straight-forward line search method for first order methods
function [xk,fxk,L] = LineSearch(A,b,grd,alpha,yk,L)
    beta = 1.2;
    stopBacktrack = false;
    fyk = .5*norm(A*yk - b).^2 + alpha*norm(yk,1); 
    while(~stopBacktrack)
        % predict the next value
        xk = sparse(shrink(yk - (1/L)*grd,alpha/L));
        fxk = .5*norm(A*xk - b).^2 + alpha*norm(xk,1);               

        % find the Lipschitz constant that ensures the Taylor
        % expansion is valid
        q =  fyk + (xk-yk)'*grd + 0.5*L*norm(xk-yk).^2 + alpha*norm(xk,1);
        if(fxk <= q)
            stopBacktrack = true;
        else
            % increase the Lipschitz constant estimate 
            L = L*beta;
        end
    end
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
