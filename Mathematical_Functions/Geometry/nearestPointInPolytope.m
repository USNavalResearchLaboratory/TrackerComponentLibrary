function [x,exitCode]=nearestPointInPolytope(P,z,algorithm,algParams)
%%NEARESTPOINTINPOLYTOPE Given a polytope defined as the convex hull of a
%          set of points, determine the point on the hull that is closest
%          to z. That is, solve the optimization problem:
%          minimize norm(x-z)^2
%          such that x=sum_{i=1}^m w(m)*P(:,m)   and w(m)>=0
%          This is useful in implementing conjugate subgradient algorithms.
%
%INPUTS: P A nXm set of m vectors that define the polytope. It
%          is okay if extra vectors within the polytope are present.
%        z An optional nX1 point used as the point of reference to which
%          the closest point on the polytope is desired. If omitted, the
%          origin is used.
% algorithm A parameter selecting the algorithm to use. Possible values
%          are:
%          0 (The default if omitted or an empty matrix is passed) Use
%            Wolfe's algorithm of [1].
%          1 Reformulate the problem as a convex semidefinite quadratic
%            programming problem (see below) and solve it using the
%            convexQuadProgSOCP function. 
% algParams Optional parameters for the algorithm. For algorithm=0,
%          tolerances values Z1, Z2 and Z3 in [1] can be specified. Also,
%          max_iters can be used to specify the maximum number of
%          iterations of the algorithm. the default is 100*m. The defaults
%          for Z1, Z2, and Z3 are those given in [1]: 1e-12, 1e-10, and
%          1e-10. For algorithm=1, algParams is the params input to the
%          convexQuadProgSOCP function.
%
%OUTPUTS: x The nX1 vector on the polytope that is closest to z or z if z
%           is in the polytope.
%  exitCode A value indicating how the algorithm terminated. For
%           algorithm=1, this is the value of the info output of the
%           function convexQuadProgSOCP. For algorithm=0, the values are
%           as follows:
%           0 The algorithm converged according to the optimality criterion
%             using Z1 in Step 1c in [1]. 
%           1 The algorithm stopped in step 1d in [1], which usually
%             indicates that Z1 was chosen to be too small.
%           2 Finite precision errors made it such that there were no zero
%             elements in w in step 3e in [1].
%           -1 The maximum number of iterations passed.
%
%The algorithm of [1] is generally very fast, but appears to have
%exponential theoretical worst-case complexity.
%
%The equivalent quadratic programming problem that is used if algorithm=1
%is minimize norm(P*w)^2
%   subject to sum(w)=1 and w>=0
%The vector x is recovered via x=P*w.
%
%EXAMPLE:
%Here, we demonstrate that both algorithms produce essentially the same
%result, but Wolfe's algorithm is more accurate and faster.
% P=2*rand(5,100)-1;
% z=[10;0;0;0;0];
% algorithm=0;
% x0=nearestPointInPolytope(P,z,algorithm);
% algorithm=1;
% algParams.eps=1e-12;
% x1=nearestPointInPolytope(P,z,algorithm,algParams);
% %Because of how both solutions are obtained, they should be both on the
% %polytope. However, one will see that Wolfe's solution is generally
% %better than the quadratic programming solution:
% norm(x0-z)
% norm(x1-z)%Larger error.
%
%REFERENCES:
%[1] P. Wolfe, "Finding the nearest point in a polytope," Mathematical
%    Programming, vol. 11, no. 1, pp. 128-149, Dec. 1976.
%[2] J. De Loera, J. Haddock, and R. Luis, "The minimum Euclidean-norm
%    point in a convex polytope: Wolfe's combinatorial algorithm is
%    exponential," ArXiv, 3 Nov. 2017. [Online]. Available:
%    https://arxiv.org/abs/1710.02608
%
%January 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin>1&&~isempty(z))
    P=bsxfun(@minus,P,z);
end

if(nargin<3||isempty(algorithm))
   algorithm=0; 
end

switch(algorithm)
    case 0
        %Default tolerance values.
        Z1=1e-12;
        Z2=1e-10;
        Z3=1e-10;
        
        %A default number of iterations.
        m=size(P,2);
        maxIter=100*m;
        
        if(nargin>3&&~isempty(algParams))
            if(~isstruct(algParams))
               error('algParams must be a structure.') 
            end
            
            if(isfield(algParams,'Z1'))
               Z1=algParams.Z1;
            end
            
            if(isfield(algParams,'Z2'))
               Z2=algParams.Z1;
            end
            
            if(isfield(algParams,'Z3'))
               Z3=algParams.Z3;
            end
            
            if(isfield(algParams,'max_iters'))
                maxIter=algParams.max_iters;
            end
        end
        
        [x,exitCode]=WolfeAlg(P,Z1,Z2,Z3,maxIter);
    case 1%The convex quadratic programming formulation.
        m=size(P,2);

        if(nargin<4)
            algParams=[];
        end
        
        G=P'*P;
        
        C=zeros(m+1,m);
        %The equality constraint.
        numEqConst=1;
        C(1,:)=ones(1,m);
        
        %The inequality constraints.
        C(2:(m+1),1:m)=eye(m);
        b=zeros(m+1,1);
        b(1)=1;
        
        a=zeros(m,1);
        
        [w,~,exitCode]=convexQuadProgSOCP(G,a,C',b,numEqConst,algParams);
        
        %Force it to satisfy the constraints.
        w(w<0)=0;
        w=w/sum(w);
        
        x=P*w;
    otherwise
        error('Unknown Algorithm Specified.')
end
%Undo any shift due to the point in question being moved.
if(nargin>1&&~isempty(z))
    x=bsxfun(@plus,x,z);
end

end

function [x,exitCode]=WolfeAlg(P,Z1,Z2,Z3,maxIter)
%%WOLFEALG This function implement's Wolfe's algorithm in [1] for finding
%          the point in a polytope that is closest to the origin.
%
%REFERENCES:
%[1] P. Wolfe, "Finding the nearest point in a polytope," Mathematical
%    Programming, vol. 11, no. 1, pp. 128-149, Dec. 1976.
%
%January 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.

n=size(P,1);
m=size(P,2);

%Indices of the columns in P that are in the columns constituting Q.
SIdx=zeros(m,1);
w=zeros(m,1);

PMags2=sum(P.*P,1);

bVec=zeros(n+1,1);
bVec(1)=1;

%Step 0
[~,J]=min(PMags2);
SIdx(1)=J;
numInSet=1;
w(1)=1;

exitCode=-1;
skip1=false;
for curIter=1:maxIter
    if(skip1==false)
        %Step 1a
        x=P(:,SIdx(1:numInSet))*w(1:numInSet);

        %Step 1b
        [~,J]=min(sum(bsxfun(@times,x,P),1));

        %Step 1c
        if(x'*P(:,J)>x'*x-Z1*max(PMags2(J),PMags2(SIdx(1:numInSet))))
            exitCode=0;
            break;
        end

        %Step 1d
        if(any(J==SIdx(1:numInSet)))
            exitCode=1;
            break;
        end

        %Step 1e
        numInSet=numInSet+1;
        SIdx(numInSet)=J;
        w(numInSet)=0;
    end
    
    %Step 2a. We are solving Equation 4.1 using Algorithm C.
    u=[ones(1,numInSet);
       P(:,SIdx(1:numInSet))]\bVec;
    v=u/sum(u);

    %Step 2b
    if(any(v>Z2))
        w=v;
        skip1=false;
        continue;
    end
    
    %Step 3a
    POS=find(w(1:numInSet)-v>Z3);
    
    %Step3b
    theta=min(1,min(w(POS)./(w(POS)-v(POS))));
    
    %Step 3c
    w(1:numInSet)=theta*w(1:numInSet)+(1-theta)*v;
    
    %Step 3d This only needs to work on w(1:numInSet), but is is simpler in
    %Matlab to do it over all of the w vector that was preallocated.
    w(w<=Z2)=0;
    
    %Step 3e We choose to delete the first zero component.
    idx=find(w(1:numInSet)==0,1);
    if(isempty(idx))
        exitCode=2;
        break;
    end
    w(idx)=w(numInSet);
    SIdx(idx)=SIdx(numInSet);
    numInSet=numInSet-1;
    skip1=true;
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
