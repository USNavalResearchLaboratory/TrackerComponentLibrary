function [x,b]=pseudoRangeLoc(rPseudo,s,RelDiff)
%%PSEUDORANGELOC Given a set of numDim=1 pseudorange measurements, such as
%                one might encounter in a GPS receiver, determine the
%                location of the receiver as well as its clock bias. The 
%                ith pseudorange measurement is assumed to be of the form
%                rPseudo(:,i)=norm(x-s(:,i))+b (possibly with noise)
%                where b is the bias, x is the unknown receiver location
%                and s is the location of the satellite. This algorithm can
%                be used to obtain initial estimates to more sophisticated
%                GPS algorithms.
%
%INPUTS: rPseudo The numMeasX1 or 1XnumMeas set of pseudoRange
%                measurements. There must be at least numDim+1
%                measurements.
%              s The numDimXnumMeas set of reference locations (satellite
%                locations in GPS applications).
%        RelDiff The algorithm will produce multiple solutions in some
%                instances, such as when all sensors are coplanar. Also,
%                given large biases, finite precision errors, and the
%                presence of noise, it can be difficult to determine which
%                solution is correct. This is a threshold. If
%                abs(residVals(1)-residVals(2))<=RelDiff*residVals(2)
%                where residVals(1) is the solution with the least
%                disagreement with rPseudo and residVals(2) is the solution
%                with the highest disagreement with rPseudo, then both
%                solutions will be returned. Otherwise, only the solution
%                with the best agreement will be returned. The default if
%                omitted or an empty matrix is passed is RelDiff=1/2. To
%                always return both solutions, set RelDiff>=1.
%
%OUTPUTS: x The numDimXnumSol target location solution or set of solutions.
%           If the discriminant term in the algorithm is negative
%           (unlikely), then this and the next output will be empty
%           matrices.
%         b The 1XnumSol set of biases.
%
%This function implements the algorithm of [1], which is also given in
%Appendix D of [2]. Often, when a plausible false solution is present, it
%will be far from the Earth and can be eliminated.
%
%EXAMPLE 1:
%This is the 1D example in [1].
% s=[-4,4];
% rPseudo=[4;2];
% [x,b]=pseudoRangeLoc(rPseudo,s);
% %One will get r=1 and xDeltaT=-1.
% rPseudoCalc=sqrt(bsxfun(@minus,x,s).^2).'+b;
% max(abs(rPseudoCalc-rPseudo))
%The error in the computed ranges will be on the order of finite precision
%limitiations.
%
%EXAMPLE 2:
%This is a simple 3D example. We find the solution and then verify that the
%computed pseudo ranges given the position and bias equal those provided.
% s=[10, 8,5,2;
%    -4,12,5,7;
%     6,-3,5,-13];
% xTrue=[1;2;3];
% biasTrue=1/2;
% rPseudo=sqrt(sum(bsxfun(@minus,xTrue,s).^2,1))+biasTrue;
% [x,b]=pseudoRangeLoc(rPseudo,s);
% rPseudoCalc=sqrt(sum(bsxfun(@minus,x,s).^2,1))+b;
% max(abs(rPseudoCalc(:)-rPseudo(:)))
%The calculated ranges agree within finite precision limits. One will also
%see that x and b agree with xTrue and biasTrue.
%
%REFERENCES:
%[1] S. Bancroft, "An algebraic solution of the GPS equations," IEEE
%    Transactions on Aerospace and Electronic Systems, vol. 21, no. 7, pp.
%    56-59, Jan. 1985.
%[2] J. Sanz Subirana, J. M. Juan Zornoza, and M. Hernandez-Pajarez, "GNSS
%    data processing, Volume I: Fundamentals and algorithms," European
%    Space Agency, Tech. Rep. ESA TM-23/1, May 2013.
%
%September 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(RelDiff))
    RelDiff=1/2;
end

numDim=size(s,1);
numSat=size(s,2);

if(numSat<numDim+1)
    error('There must be at least numDim+1 measurements for observability.')
end

rPseudo=rPseudo(:);

%Equation 5 in [1].
A=[s',rPseudo];

r=zeros(numSat,1);
for k=1:numSat
    %Equation 8 in [1].
    r(k)=minkowskiInnerProd(A(k,:),A(k,:))/2;
end

%Equation 9 in [1].
B=pinv(A);
oneVec=ones(numSat,1);

%Equation 10 in [1].
u=B*oneVec;
%Equation 11 in [1].
v=B*r;

E=minkowskiInnerProd(u,u);
F=minkowskiInnerProd(u,v)-1;
G=minkowskiInnerProd(v,v);

%Just directly incorporate the 2 from the quadratic in Equation 15 in [1]
%into F.
F=2*F;

discrim=F^2-4*E*G;
if(discrim<0)
    x=[];
    b=[];
    return;
end

rootTerm=sqrt(discrim);
denom=2*E;

lambda=zeros(2,1);
lambda(1)=(-F+rootTerm)/denom;
lambda(2)=(-F-rootTerm)/denom;

y=zeros(numDim+1,1);
residVals=zeros(2,1);
for k=1:2
    %Equation 16
    y(:,k)=lambda(k)*u+v;

    x=y(1:numDim,k);
    b=-y(numDim+1,k);
    
    %The calculated measurements given a target at x with bias b. This is
    %just Equation 2 in [1].
    calcMeas=sqrt(sum(bsxfun(@minus,x,s).^2,1))+b;
    
    residVals(k)=norm(rPseudo(:)-calcMeas(:));
end

if(residVals(1)<residVals(2))
    minIdx=1;
    otherIdx=2;
else
    minIdx=2;
    otherIdx=1;
end

if(abs(residVals(1)-residVals(2))<=RelDiff*residVals(otherIdx))
    %Return both solutions.
    x=y(1:numDim,[minIdx,otherIdx]);
    b=-y(numDim+1,[minIdx,otherIdx]);
else
    %Return only one solution
    x=y(1:numDim,minIdx);
    b=-y(numDim+1,minIdx);
end
end

function val=minkowskiInnerProd(a,b)
%%MINKOWSKIINNERPROD The Minkowski inner product defined in Equation 4 in
%                    [1]. This is also known as the Lorentz inner product.
%
%REFERENCES:
%[1] S. Bancroft, "An algebraic solution of the GPS equations," IEEE
%    Transactions on Aerospace and Electronic Systems, vol. 21, no. 7, pp.
%    56-59, Jan. 1985.
%
%September 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.

vec=a.*b;
vec(end)=-vec(end);
val=sum(vec);

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
