function zp=nearestPointOnQuadricSurface(y,A,bIn,c,epsVal)
%%NEARESTPOINTONQUADRICSURFACE Given a point y, find the nearest point to
%       it on a quadric surface. The surface is defined by
%       x'*A*x+bIn'*x+c=0. All values are real.
%
%INPUTS: y A numDimX1 point.
%        A A numDimXnumDim symmetric matrix. It does not need to be
%          invertible.
%      bIn A numDimX1 vector.
%        c A scalar value.
%
%OUTPUTS: zp The closes point on the surface to y, or an empty matrix if
%            finite precision errors caused the roots function to return no
%            valid solutions.
%
%The solution in 3D is given in Chapter 10.5.1 of [1]. However, this
%function is not limited to 3D solutions. For other dimensionalities, the
%polynomial in terms of t is found and explicitly solved.
%
%EXAMPLE 1:
% A=[82,  7, 12, 17;
%     7, 72, 17, 22;
%    12, 17, 62, 27;
%    17, 22, 27, 52];
% bIn=[4;-3;1;0];
% c=-6;
% y=[1;2;13;-8];
% zp=nearestPointOnQuadricSurface(y,A,bIn,c)
%One will get
%zp=[0.003766587429133;0.073180186270540;0.301170048336470;-0.324551318196678];
%and can verify that it satisfies the equation zp'*A*zp+bIn'*zp+c=0 within
%finite precision bounds.
%
%EXAMPLE 2:
%In this example, we verify that we get the same solution as
%nearestPointOnEllipsoid when both are passed equivalent problems.
% A=[82,  7, 12, 17;
%     7, 72, 17, 22;
%    12, 17, 62, 27;
%    17, 22, 27, 52];
% bIn=[0;0;0;0];
% c=-6;
% y=[1;2;13;-8];
% z1=nearestPointOnEllipsoid([0;0;0;0],A,y,-c)
% z2=nearestPointOnQuadricSurface(y,A,bIn,c)
%One will see the both solutions are essentially equal with 
%z2=[0.028873225664696;0.048315933230922;0.310279553948698;-0.325946108741623]
%and can verify that it satisfies the equation zp'*A*zp+bIn'*zp+c=0 within
%finite precision bounds.
%
%REFERENCES:
%[1] P. J. Schneider and D. H. Eberly, Geometric Tools for Computer
%    Graphics. Amsterdam: Morgan Kaufmann Publishers, 2003.
%
%April 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<5||isempty(epsVal))
    epsVal=2^12*eps(1);
end

[R,D]=eig(A);
a=R'*y;
b=R'*bIn;
d=diag(D);

numDims=length(d);

%Do the term times c first
prodVal=c;
for i=1:numDims
    prodVal=conv(prodVal,conv([2*d(i);1],[2*d(i);1]));
end
polyEq=prodVal;

%Next, add in the terms from the equations with two inv(I+2*t*D) entries,
for curDim=1:numDims
    prodVal=d(curDim)*conv([-b(curDim);a(curDim)],[-b(curDim);a(curDim)]);
    for i=[1:(curDim-1),(curDim+1):numDims]
        prodVal=conv(prodVal,conv([2*d(i);1],[2*d(i);1]));
    end
    polyEq=polySum(polyEq,prodVal);
end

%Finally, add in the terms from the equation containing only one
%inv(I+2*t*D)  entry.
for curDim=1:numDims
    prodVal=b(curDim)*conv([-b(curDim);a(curDim)],[2*d(curDim);1]);
    for i=[1:(curDim-1),(curDim+1):numDims]
        prodVal=conv(prodVal,conv([2*d(i);1],[2*d(i);1]));
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
I=eye(numDims,numDims);
for curSol=1:numSol
    zp(:,curSol)=(I+2*lambda(curSol)*A)\(y-lambda(curSol)*bIn);
    %ztest=R*inv(I+2*lambda(curSol)*D)*(a-t*b
end

minDiffVal=Inf;
minSol=1;
for curSol=1:numSol
    diffMag=norm(y-zp(:,curSol));
    if(diffMag<minDiffVal)
        %Verify that it is actually a solution.
        if(abs(zp(:,curSol)'*A*zp(:,curSol)+bIn'*zp(:,curSol)+c)<epsVal)
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
