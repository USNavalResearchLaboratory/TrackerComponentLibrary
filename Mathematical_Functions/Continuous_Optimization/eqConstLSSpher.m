function x=eqConstLSSpher(A,b,alpha,epsRed)
%%EQCONSTLSSPHER Find x to minimize norm(A*x-b,2) under the constraint
%                that norm(x,2)=alpha. This is essentialy contraining x to
%                the surface of a sphere of radius alpha. This only solves
%                real systems.
%
%INPUTS: A A real mXn matrix with m>=n.
%        b A real mX1 vector.
%    alpha The equality constraint value. If this parameter is omitted or
%          an empty matrix is passed, the default of 1 is used.
%   epsRed Let xU be the unconstrained solution to the problem. A possible
%          constrained solutions is considered valid if
%          abs(norm(x)^2-alpha^2)<=epsRed*abs(norm(xU)^2-alpha^2)
%          If no such solutions are found, then xU is returned. The default
%          for this parameter if omitted or an empty matrix is passed is
%          1e-9. This parameter should be between 0 and 1.
%
%OUTPUTS: x The optimal value of x subject to the spherical constraint.
%    lambda The Lagrangian multiplier used in the optimization. lambda>=0.
%           If lambda=0, then the constraint on x did not have to be
%           enforced.
%
%This implements a modified version of the algorithm of Chapter 6.2.1 of
%[1]. The algorithm of Chapter 6.2.1 of [1] solves the optimization
%constrained such that norm(x,2)<=alpha. We wish to solve the equality
%constrained problem. To do so, we always enforce the constraint. An
%equation has to be solved for scalar zeros over lambda. To do so, we
%multiply both sids by the denominators resulting in a polynomial system.
%The system is then solved. However, some solutions might not satisfy the
%constraint (due to cancellation in denominators). Thus, candidate
%solutions that do not improve the error in the magnitude by a sufficient
%amount compared to the unconstrained solution are discarded. If no
%solutions are left, the the unconstrained solution is used.
%
%EXAMPLE:
%This is a simple example where the x returned by the
%inequality-constrained algorithm is too small.
% A=magic(8)+20*eye(8);
% b=(1:8).';
% x=inv(A)*b;
% norm(x)%Norm <1.
% xConst=eqConstLSSpher(A,b);
% norm(xConst)%Norm=1
% %It is not just normalizing the vector; elements are not scaled by a
% %constant.
% x./xConst
%
%REFERENCES:
%[1] G. H. Golub and C. F. Van Loan, Matrix Computations, 4th ed.
%    Baltimore: Johns Hopkins University Press, 2013.
%
%December 2020 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(alpha))
    alpha=1;
end

if(nargin<4||isempty(epsRed))
    epsRed=1e-9;
end

r=rank(A);

[U,Sigma,V]=svd(A,0);
sigma=diag(Sigma);

%Sums are only up to r, so get rid of the extra elements.
U=U(:,1:r);
V=V(:,1:r);
sigma=sigma(1:r);

bTilde=U'*b;
xUnconst=V*(bTilde./sigma);

alpha2Unconst=abs(norm(xUnconst)^2-alpha^2);

numVal=(sigma.*bTilde).^2;
denomVal=sigma.^2;

numDim=length(sigma);
polynoms=[ones(numDim,1),2*denomVal,denomVal.^2];

%Construct the polynomial to solve.
polyNom=0;
for k1=1:numDim
    curPoly=numVal(k1);
    for k2=1:numDim
        if(k2==k1)
            continue;
        end
        curPoly=conv(curPoly,polynoms(k2,:));
    end
    polyNom=polySum(curPoly,polyNom);
end

%The final term
curPoly=-alpha^2;
for k2=1:numDim
    curPoly=conv(curPoly,polynoms(k2,:));
end
polyNom=polySum(curPoly,polyNom);
lambdaVals=roots(polyNom).';
%Get rid of imaginary solutions.
lambdaVals=lambdaVals(imag(lambdaVals)==0);
numSol=length(lambdaVals);

%Due to certain values corresponding to zero denominators, there can be
%more solutions in lambda than are valid. We must eliminate all solutions
%that do not produce x values with the correct magnitude. First, get all
%possibly valid solutions regardless of the magnitude.
x=zeros(numDim,numSol);
numKept=0;
for k=1:numSol
    xCur=V*((sigma.*bTilde)./(sigma.^2+lambdaVals(k)));
    normErr=abs(norm(xCur)^2-alpha^2);
    if(all(isfinite(xCur))&&normErr<=epsRed*alpha2Unconst)
        numKept=numKept+1;
        x(:,numKept)=xCur;
    end
end

if(numKept==0)
    %If nothing was kept, then just use the unconstrained solution. 
    x=xUnconst;
    return 
end
x=x(:,1:numKept);

%Of all of the solutions kept, take the one that minimizes the original
%optimization problem.
if(numKept>1)
    minCost=norm(A*x(:,1)-b,2);
    minIdx=1;
    for k=2:numKept
       curCost=norm(A*x(:,k)-b,2);
       
       if(curCost<minCost)
           minIdx=k;
           minCost=curCost;
       end
    end
    x=x(:,minIdx);
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
