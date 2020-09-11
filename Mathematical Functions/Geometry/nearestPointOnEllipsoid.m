function zp=nearestPointOnEllipsoid(zIn,A,pIn,gammaVal,epsVal)
%%NEARESTPOINTONELLIPSOID Given an ellipsoid in an arbitrary number of
%        dimensions such that a point zp on the surface of the ellipsoid
%        satisfies the equation (zIn-zp)'*A*(zIn-zp)=gammaVal find the
%        point on the ellipsoid that is closest to another point p. All
%        inputs are real.
%
%INPUTS: zIn The numDimX1 center of the ellipsoid.
%          A A numDimXnumDim symmetric, positive definite matrix that
%            specifies the size and shape of the ellipse or ellipsoid,
%            where a point zp is on the ellipse/ellipsoid if
%            (zIn-zp)'*A*(zIn-zp)=gammaVal.
%        pIn A numDimX1 point.
%   gammaVal The threshold for declaring a point to be on the ellipsoid. If
%            this parameter is omitted or an empty matrix is passed, the
%            default value of 1 is used.
%     epsVal This parameter can usually be omitted as the default often
%            works. When going through potential solutions, this is the
%            maximum value of abs(diff'*A*diff-gammaVal) for which a
%            potential solution is considered correct. The default if this
%            parameter is omitted or an empty matrix is passed is
%            2^10*eps(gammaVal).
%
%OUTPUTS: zp The numDimX1 point on the ellipse that is closest to p. When
%            multiple solutions exist, only one is chosen. If an empty
%            matrix is returned, then finite precision errors caused the
%            roots function to return no valid solutions.
%
%The solution in 3D is outlined in Chapter 10.5.2 of [1]. However, Equation
%10.8 is missing a lambda times the gradient term (See Chapter 10.5.1 for a
%similar scenario where such a term is necessary and present). The solution
%in this file is to directly solve the polynomial in the unnumbered
%equation after Equation 10.9 in [1] for lambda. Given lambda, by
%substituting Equation 10.9 into 10.10, one can obtain the x y and z
%components.
%
%However, this function is not limited to 3D solutions. For other
%dimensionalities, the same steps are taken to form the polynomial after
%Equation 10.9 in [1]. Thus, this function implements the general
%polynomial and solves it. As only one solution to the polynomial can be
%correct, all real solutions to lambda are tried and the one having the
%lowest cost that is within epsVal of gammaVal is used. When multiple
%finite solutions exist (e.g. in the center of the ellipsoid), then only
%one is chosen.
%
%The algorithm in Chapter 10.5.2 of [1] only handles axis-aligned
%ellipsoids. This function deals with general ellipsoids by shifting
%everything to the origin, doing an eigenvalue/ eigenvector decomposition
%of A to rotate things to be axis-aligned and then undoing that in the end.
%The value gammaVal is handled by replacing A with A/gammaVal.
%(zp-z)*A*(zp-z)=gammaVal
%
%EXAMPLE 1:
%Here, we find the nearest point on an ellipsoid in 3D.
% A=[27,  4,  10;
%     4, 21, 16;
%    10, 16, 15];
% z=[12;24;36];
% p=[1000;-1000;2000];
% gammaVal=2.5;
% zp=nearestPointOnEllipsoid(z,A,p,gammaVal)
%One gets zp=[11.606225228425036;22.829459247176374;37.594179506992610]
%One can verify that (z-zp)'*A*(z-zp)=gammaVal within finite precision
%limits.
%
%EXAMPLE 2:
%This example is used in the function nearestPointInEllipsoid. There, since
%the point chosen in inside of the 2D ellipse, the function returns that
%point. However, this function returns the point on the outside of the
%ellipse that is closest.
% A=[1,0;
%    0,5];
% z=[0;10];
% p=[0.5;10];
% gammaVal=1;
% zp=nearestPointOnEllipsoid(z,A,p,gammaVal);
% figure(1)
% clf
% hold on
% axis([-6,6,-6+10,6+10])
% axis square
% drawEllipse(z,A,gammaVal,'g','linewidth',2)
% plot([p(1),zp(1)],[p(2),zp(2)],'--c')
% scatter(p(1),p(2),'ok','linewidth',2)
% scatter(zp(1),zp(2),'xr','linewidth',2)
%
%REFERENCES:
%[1] P. J. Schneider and D. H. Eberly, Geometric Tools for Computer
%    Graphics. Amsterdam: Morgan Kaufmann Publishers, 2003.
%
%April 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4||isempty(gammaVal))
    gammaVal=1;
end

if(nargin<5||isempty(epsVal))
    epsVal=2^10*eps(gammaVal);
end

B=A/gammaVal;
[V,D]=eig(B);

%The point in the zero-centered, diagonalized coordinate system.
p=V'*(pIn-zIn);
a2=1./diag(D);

numDims=length(a2);

%Do the first term.
prodVal=-1;
for i=1:numDims
    prodVal=conv(prodVal,conv([2,a2(i)],[2;a2(i)]));
end
polyEq=prodVal;

%Construct and add the subsequent terms.
for curDim=1:numDims
    prodVal=a2(curDim)*p(curDim)^2;
    for i=[1:(curDim-1),(curDim+1):numDims]
        prodVal=conv(prodVal,conv([2,a2(i)],[2;a2(i)]));
    end
    polyEq=polySum(polyEq,prodVal);
end

lambda=roots(polyEq).';
lambda=lambda(imag(lambda)==0);
if(isempty(lambda))
    zp=[];
    return;
end

numSol=length(lambda);
zp=zeros(numDims,numSol);
for curDim=1:numDims
    zp(curDim,:)=a2(curDim)*p(curDim)./(a2(curDim)+2*lambda);
end

zp=bsxfun(@plus,V*zp,zIn);

minDiffVal=Inf;
minSol=1;
for curSol=1:numSol
    diffMag=norm(pIn-zp(:,curSol));
    if(diffMag<minDiffVal)
        diff=zIn-zp(:,curSol);
        %Verify that it is actually a solution.
        if(abs(diff'*A*diff-gammaVal)<epsVal)
            minSol=curSol;
            minDiffVal=diffMag;
        end
    end
end

zp=zp(:,minSol);

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
