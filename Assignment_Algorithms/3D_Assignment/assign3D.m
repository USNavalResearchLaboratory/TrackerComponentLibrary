function [tuples,fStar,qStar,u,exitCode]=assign3D(C,maximize,subgradMethod,subgradParams,maxIter,AbsTol,RelTol)
%%ASSIGN3D Approximate the solution to the operations research axial 3D
%         assignment problem using a dual-primal Lagrangian relaxation 
%         technique. Such problems are NP-hard. The optimization problem
%         being solved is
%         minimize (or maximize)
%         \sum_{i=1}^{n1}\sum_{j=1}^{n2}\sum_{k=1}^{n3}C_{i,j,k}*\rho_{i,j,k}
%         subject to
%         \sum_{i=1}^{n1}\sum_{j=1}^{n2}\rho_{i,j,k}<=1 for all k
%         \sum_{i=1}^{n1}\sum_{k=1}^{n3}\rho_{i,j,k}<=1 for all j
%         \sum_{j=1}^{n2}\sum_{k=1}^{n3}\rho_{i,j,k} =1 for all i
%         \rho_{i,j,k} = 0 or 1
%         assuming that n1<=n2<=n3, and C is an n1Xn2Xn3 cost matrix.
%         This is equivalent to the optimization problem
%         min (or max) sum_{i} C(i,phi_2(i),phi_3(i))
%         where phi_2 and phi_3 are length n1 arrangements of n2 and n3
%         items over which the minimization is performed.
%
%INPUTS: C An n1Xn2Xn3 cost hypermatrix with n1<=n2<=n3. C cannot contain
%          any NaNs and the largest finite element minus the smallest
%          element is a finite quantity (does not overflow) when performing
%          minimization and where the smallest finite element minus the
%          largest element is finite when performing maximization.
%          Forbidden assignments can be given costs of +Inf for
%          minimization and -Inf for maximization. During minimization,
%          there should be no elements with -Inf cost and no elements with
%          +Inf cost during maximization.
% maximize If true, the minimization problem is transformed into a
%          maximization problem. The default if this parameter is omitted
%          or an empty matrix is passed is false.
% subgradMethod A parameter indicating the subgradient optimization
%          algorithm to use. Possible values are:
%          0 Polyak's Method from [2]. (The default if omitted or an empty
%            matrix is passed) In this method, for minimization, the dual
%            variables are updated as uNew=u+gamma*(fStar-q)/norm(g)^2*g,
%            where g is the subgradient, q is the current dual cost, and g
%            is the subgradient vector. gamma is a design parameter.
%          1 Bragin's Method from [3]. The dual update for minimization is
%            uNew=u+alpha*g
%            where during the first step, alpha=(fStar-q)/norm(g)^2 and for
%            subsequent steps it is
%            alpha=(1-1/(M*k^(1-1/k^r)))*alphaPrev*norm(gPrev)/norm(g)
%            where alphaPrev is alpha for the previous step, k is the step
%            number (starting at 0), gPrev is the subgradient before the
%            current step and M and r are design parameters.
%          2 Bertsekas' Heuristic Method from [4]. The update is
%            uNew=u+alpha*g
%            where alpha=((1+a/beta^b)*qMax-q)/norm(g)^2 where qMax is the
%            highest dual cost value encountered thus far (when
%            minimizing), a and b are design parameters, and beta is a
%            value that increases or decreases in the algorithm depending
%            on whether or not there was an improvement in the best dual
%            cost found after the last step.
%          3 Shor's Space Dilation Algorithm from [5].
%            The update is
%            uNew=u+alpha*H*g;
%            where alpha=(2*M/(M+1))*(fStar-q)/norm(d), d=g and H is a
%            matrix that also depends on M and d, which is a design
%            parameter. Also, after NR iterations, the recursive turms that
%            go into H reset.
%          4 The r Space Dilation Algorithm from [5]. This is the same as 3
%            except d=g-gPrev except for the first iteration and if
%            g=gPrev, in which case d=g.
% subgradParams This is a structure that takes the design parameters for
%          the selected cubgradient algorithm. Possible fields of the
%          structure as well as default values depend on the selected
%          subgradient method and are:
%          Method 0: 'gamma' with default 1.
%          Method 1: 'M' and 'r' with defaults 2.8 and 0.06.
%          Method 2: 'a and 'b' with defaults 0.3 and 1.5.
%          Method 3,4: 'M', 'NR', and 'normBound' with defaults 2, n3, and
%                     eps(). normBound is such that if norm(B.'*d) in the
%                     computation of the matrix H is <= normBound, then B
%                     is reset to the identity matrix.
%  maxIter The maximum number of iterations to perform. The default value
%          if this parameter is omitted or an empty matrix is passed is 20.
%   AbsTol The absolute duality gap to use for convergence determiniation.
%          Convergence is declared if (fStar-qStar)<=AbsTol, where if
%          minimizing, fStar is the smallest value of the primal function
%          found and qStar is the highest value of the dual cost function
%          found. The default if omitted or an empty matrix is passed is
%          1e-10.
%   RelTol The relative duality gap to use for convergence determiniation.
%          Convergence is declared if (fStar-qStar)<=RelTol*abs(qStar). The
%          default if omitted or an empty matrix is passed is 0.05.
%
%OUTPUTS: tuples A 3Xn1 matrix where tuples(:,i) is the ith assigned tuple
%                in C as described above. An empty matrix is returned if
%                the problem is infeasible or the heuristic to obtain a
%                feasible solution failed.
%          fStar The cost of the assignment in tuples. This is the sum of
%                the values of the assigned elements in C.
%          qStar The value of the best dual solution found. The absolute
%                duality gap is fStar-qStar and the relative duality gap is
%                abs(fStar-qStar)/abs(qStar)
%              u The n3X1 vector of dual variables.
%       exitCode A parameter indicating how the algorithm terminated.
%                Possible values are:
%                -3 A non-finite number arose in the dual variables, so the
%                   algorithm stopped.
%                -2 The algorithm did not converge within the alotted
%                   number of iterations. However, a valid feasible
%                   assignment was found.
%                -1 A subproblem in computing the dual or in obtaining a
%                   feasible primal solution was infeasible.The output
%                   tuples will be an empty matrix.
%                >=0 Values that are zero or positive indicate the
%                   convergence was obtained and the returned value is the
%                   number of iterations.
%
%This function implements the algorithm of [1], but modified so that it
%does not have any unconstrained indices and offering different subgradient
%methods. We are solving the "operations research" 3D assignment problem
%rather than the "data fusion" 3D assignment problem of [1].
%
%EXAMPLE:
%Here is an example using a random problem:
% n1=10;
% n2=12;
% n3=15;
% C=randn(n1,n2,n3);
% maximize=false;
% [tuples,fStar,qStar,u,exitCode]=assign3D(C,maximize)
%
%REFERENCES:
%[1] K. Pattipati, S. Deb, Y. Bar-Shalom, and R. B. Washburn Jr., "A
%    new relaxation algorithm and passive sensor data association," IEEE
%    Transactions on Automatic Control, vol. 37, no. 2, pp. 198-213, Feb.
%    1992.
%[2] B. T. Polyak, "Minimization of unsmooth functionals," USSR
%    Computational Mathematics and Mathematical Physics, vol. 9, no. 3, pp.
%    14-29, 1969.
%[3] M. A. Bragin, P. B. Luh, J. H. Yan, N. Yu, and G. A. Stern,
%    "Convergence of the surrogate Lagrangian relaxation method," Journal
%    of Optimization Theory and Applications, vol. 164, no. 1, pp. 173-201,
%    Jan. 2015.
%[4] D. P. Bertsekas, Nonlinear Programming, 3rd ed. Belmont, MA: Athena
%    Scientific, 2016, Chapter 7.5.
%[5] N. Z. Shor, "Utilization of the operation of space dilation in the
%    minimization of convex functions," Cybernetics, vol. 6, no. 1, pp. 7-
%    15, Dec. 1972.
%
%February 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(maximize))
   maximize=false; 
end

if(nargin<3||isempty(subgradMethod))
    subgradMethod=0; 
end

if(nargin<4)
    subgradParams=[];
end

if(nargin<5||isempty(maxIter))
    maxIter=20;
end

if(nargin<6||isempty(AbsTol))
    AbsTol=1e-10;
end

if(nargin<7||isempty(RelTol))
    RelTol=0.05;
end

nVals=size(C);

%Deal with the special case of C being scalar.
if(isscalar(C))
    tuples=[1;1;1];
    fStar=C(1);
    qStar=fStar;
    u=0;
    exitCode=0;
    return;
end

if(isempty(C))%The empty matrix special case.
    tuples=[];
    fStar=0;
    qStar=0;
    u=[];
    exitCode=0;
    return;
end

if(maximize)
    C=-C;
end

if(length(nVals)<3||~(nVals(1)<=nVals(2)&&nVals(2)<=nVals(3)))
    error('It is required that size(C,1)<=size(C,2)<=size(C,3)')
end

n1=nVals(1);
n2=nVals(2);
n3=nVals(3);

if((n1==n2)&&(n2==n3))
    unconstU=true;
else
    unconstU=false;
end

%Get the parameters and initial values for the chosen subgradient method.
switch(subgradMethod)
    case 0%Polyak's method
        gammaParam=1;
        
        if(nargin>3&&~isempty(subgradParams))
            if(isfield(subgradParams,'gamma'))
                gammaParam=subgradParams.gamma;
            end
        end
    case 1%Bragin's Method
        M=2.8;
        r=0.06;
        gNormPrev=[];
        alphaPrev=[];
        
        if(nargin>3&&~isempty(subgradParams))
            if(isfield(subgradParams,'M'))
                M=subgradParams.M;
            end
            
            if(isfield(subgradParams,'r'))
                r=subgradParams.r;
            end
        end
    case 2%Bertsekas' Heuristic Method
        a=0.3;
        b=1.5;
        beta=1;%Initialize
        
        if(nargin>3&&~isempty(subgradParams))
            if(isfield(subgradParams,'a'))
                a=subgradParams.a;
            end
            
            if(isfield(subgradParams,'b'))
                b=subgradParams.b;
            end
        end
    case 3%Shor's Space Dilation Algorithm.
        NR=n3;
        M=2;
        rho=(M-1)/(M+1);
        normBound=eps();

        B=eye(n3,n3);%Allocate and initialize
        dPrev=[];

        if(nargin>3&&~isempty(subgradParams))
            if(isfield(subgradParams,'NR'))
                NR=subgradParams.NR;
            end

            if(isfield(subgradParams,'M'))
                M=subgradParams.M;
            end

            if(isfield(subgradParams,'normBound'))
                normBound=subgradParams.normBound;
            end
        end
    case 4%The r Space Dilation Algorithm.
        NR=n3;
        M=2;
        rho=(M-1)/(M+1);
        normBound=eps();
        
        B=eye(n3,n3);%Allocate and initialize
        gPrev=[];
        dPrev=[];
        
        if(nargin>3&&~isempty(subgradParams))
            if(isfield(subgradParams,'NR'))
                NR=subgradParams.NR;
            end
            
            if(isfield(subgradParams,'M'))
                M=subgradParams.M;
            end
            
            if(isfield(subgradParams,'normBound'))
                normBound=subgradParams.normBound;
            end
        end
    otherwise
        error('Invalid subgradient method specified.')
end

u=zeros(n3,1);%Allocate and initialize the dual variables.
%Allocate space for the return tuples.
tuples=zeros(3,n1);

%qStar is the best (maximum) dual cost encountered thus far.
qStar=-Inf;

%fStar is the best (minimum) feasible primal cost encountered
%thus far.
fStar=Inf;
%If the loop exits without convergence, then exitCode is not modified.
exitCode=-2;
for k=0:(maxIter-1)
%%%%%%
%DUAL COST AND SUBGRADIENT UPDATE. SEE FIG. 3 IN [1].
%%%%%%
    %Note that d3=C;
    [d2,gamma2]=min(bsxfun(@plus,C,reshape(u,[1,1,n3])),[],3);
   
    %gamma1 is the columns for rows (nonzero elements in omega).
    [gamma1, ~, minVal]=assign2D(d2,false);
    
    %If the subproblem is infeasible.
    if(isempty(gamma1))
        tuples=[];
        fStar=[];
        qStar=[];
        u=[];
        exitCode=-1;
        return
    end

    q=minVal-sum(u);%The dual cost.
    
    %Keep track of the maximum q value. This is used for testing
    %convergence. 
    if(q>qStar)
        qStar=q;

        %The duality gap.
        costGap=fStar-qStar;
        
        %The correctness of this on the first iterations requires proper
        %handling of NaNs.
        if(costGap<=AbsTol||(costGap<abs(qStar)*RelTol))
            exitCode=k+1;%The algorithm converged.
            break; 
        end
    end

    %Compute the subgradient
    g=-1*ones(n3,1);
    for i1=1:n1
        i2=gamma1(i1);
        i3=gamma2(i1,i2);

        g(i3)=g(i3)+1;
    end
    
%%%%%%
%TEST THE DUAL VARIABLES AND FINISH THE LOOP. See Fig. 2 in [1].
%%%%%%
    gNorm2=dot(g,g);
    
    %If no constraints are violated, then we have a feasible solution and
    %the algorithm has converged to a local or global minimum point.
    if(gNorm2==0)
        %If the dual cost is less than the best primal solution found with
        %a heuristic, then we use the tuples for that solution. Otherwise,
        %the best primal solution to this point is returned.
        if(q<fStar)
            fStar=q;
            %Record the tuples.
            for i1=1:n1
                i2=gamma1(i1);
                i3=gamma2(i1,i2);
                tuples(:,i1)=[i1;i2;i3];
            end
        end

        exitCode=k+1;%The algorithm converged.
        break;
    end

    [tuples,fStar,retVal]=updateBestFeasSol3D(C,gamma1,tuples,fStar,qStar,AbsTol,RelTol);

    %If the primal converged.
    if(retVal==0)
        exitCode=k+1;
        break; 
    end

    %If an error occurred obtaining a feasible solution.
    if(retVal<0)
        tuples=[];
        fStar=[];
        qStar=[];
        u=[];
        exitCode=retVal;
        return;
    end

    %Perform the subgradient update.
    switch(subgradMethod)
        case 0%Polyak's method
            alpha=gammaParam*((fStar-q)/gNorm2);
            u=u+alpha*g;
        case 1%Bragin's Method
            gNorm=sqrt(gNorm2);
            if(k==0)
                alpha=((fStar-q)/gNorm2);
            else
                alpha=(1-1/(M*k^(1-k^(-r))))*alphaPrev*(gNormPrev/gNorm);
            end
            u=u+alpha*g;
            
            alphaPrev=alpha;
            gNormPrev=gNorm;
        case 2%Bertsekas' Heuristic Method
            %Note that we do not have to keep track of qStar value from the
            %previous iteration. Though in [1], the comparison is with the
            %previous qStar, if q>qStarPrev, then qStar will be set to q,
            %so the comparison q<qStarPrev has the same result as q<qStar
            %and if q<=qStarPrev, then qStarPrev=qStar, so the comparison
            %is still correct.
            if(q<qStar)
                beta=beta+1;
            else
                beta=max(beta-1,1);
            end
            
            alpha=abs(((1+a)/(beta^b))*qStar-q)/gNorm2;
            
            u=u+alpha*g;
        case 3%Shor's Space Dilation Algorithm.
            d=g;
            
            %The H matrix is initialized on the first iteration and reset
            %to the identity matrix every NR iterations.
            if(mod(k,NR)==0)
                B=eye(n3,n3);
                normVal=norm(d);%=norm(B'*d)
            else
                zeta=B.'*dPrev/norm(B.'*dPrev);
                R=eye(n3,n3)+(rho-1)*(zeta*zeta.');
                B=B*R;
                normVal=norm(B.'*d);
                if(normVal<=normBound)
                    B=eye(n3,n3);
                    normVal=norm(d);
                end
            end
            H=(B*B.')/normVal;
            
            alpha=(2*M/(M+1))*((fStar-q)/norm(d));
            
            u=u+alpha*H*d;
            
            dPrev=d;
        case 4%The r Space Dilation Algorithm.
            if(k==0||all(g==gPrev))
                d=g;
            else
                d=g-gPrev;
            end

            %The H matrix is initialized on the first iteration and reset
            %to the identity matrix every NR iterations.
            if(mod(k,NR)==0)
                B=eye(n3,n3);
                normVal=norm(d);%=norm(B'*d)
            else
                zeta=B.'*dPrev/norm(B.'*dPrev);
                R=eye(n3,n3)+(rho-1)*(zeta*zeta.');
                B=B*R;
                normVal=norm(B.'*d);
                if(norm(B.'*d)<=normBound)
                    B=eye(n3,n3);
                    normVal=norm(d);
                end
            end
            H=(B*B.')/normVal;

            alpha=(2*M/(M+1))*((fStar-q)/norm(d));

            u=u+alpha*H*d;

            gPrev=g;
            dPrev=d;
        otherwise
            error('Invalid subgradient method specified.')
    end

    if(any(~isfinite(u)))
        %This can sometimes occur with big problems and poor stepsizes.
        exitCode=-3;
        break;
    end

    %Unless the constraints are equality constraints, the dual variables
    %have to be clipped.
    if(unconstU==false)
        u(u<0)=0;
    end
end

%Adjust the gains for the case where the initial cost matrix is transformed
%so that maximization can be performed.
if(maximize)
    fStar=-fStar;
    qStar=-qStar;
    u=-u;
end
end

function [optTuples,fStar,exitCode]=updateBestFeasSol3D(C,gamma1,optTuples,fStar,qStar,AbsTol,RelTol)
%%UPDATEBESTFEASSOL3D This function implements a subroutine to
%         heuristically obtain a feasible solution given a partial set of
%         assignments. Here, it is just for the 3D assignment problem.
%         Additionally, this function updates the best feasible solution
%         cost fStar and its associated set of tuples optTuples.
%
%INPUTS: C The original n1Xn2Xn3 3D cost matrix.
%   gamma1 gamma1(i1) gives the value of n2 for each value of i1.
% optTuples The set of tuples associated with the current best primal cost.
%    fStar The current best (lowest) primal cost function value found.
%    qStar The current best (highest) dual cost function value found.
%   AbsTol The absolute threshold for determining convergence of the
%          algorithm according to the duality gap.
%   RelTol The threshold for determining convergence of the algorithm
%          according to the relative duality gap.
%
%OUTPUTS: optTuples The updated set of tuples associated with fStar.
%             fStar The updated best cost function found thus far.
%          exitCode A value indicating the terminating state of the
%                   function. Possible values are:
%                  -1 The subproblem posted is not feasible. 
%                   0 Convergence achieved.
%                   1 No errors occurred, but convegrence was not achieved.
%
%February 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.

    n1=size(C,1);
    n3=size(C,3);

    CFeas=zeros(n1,n3);
    for i1=1:n1
        CFeas(i1,:)=reshape(C(i1,gamma1(i1),:),[1,n3]);
    end

    [gammaTilde3,~,f]=assign2D(CFeas,false);
    
    %If the subproblem is not feasible.
    if(isempty(gammaTilde3))
        exitCode=-1;
        return;
    end

    if(f<fStar)
        fStar=f;
        
        for i1=1:n1
           optTuples(:,i1)=[i1;gamma1(i1);gammaTilde3(i1)];
        end
        
        %The duality gap
        costGap=fStar-qStar;

        if(costGap<=AbsTol||(costGap<abs(qStar)*RelTol))
            exitCode=0;%The algorithm converged.
        else
            exitCode=1;
        end
    else
        exitCode=1;
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
