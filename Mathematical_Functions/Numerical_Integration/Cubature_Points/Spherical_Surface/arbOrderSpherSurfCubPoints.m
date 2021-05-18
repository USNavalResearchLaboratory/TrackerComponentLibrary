function [xi,w]=arbOrderSpherSurfCubPoints(numDim,n)
%%ARBORDERSPHERSURFCUBPOINTS Generate cubature points for integrating over
%               the surface of the unit (hyper-)sphere of an arbitrary
%               order.
%
%INPUTS: numDim An integer specifying the dimensionality of the points to
%               be generated.
%             n A positive integer such that 2n-1 is the highest degree to
%               which the cubature points are accurate.
%
%OUTPUTS: xi A numDim X numCubaturePoints matrix containing the cubature
%            points (Each "point" is a vector).
%          w A numCubaturePoints X 1 vector of the weights associated with
%            the cubature points.
%
%The algorithm of Theorem 2.7-3 in Chapter 2.7 of [1] is used to generate
%the points.
%
%REFERENCES:
%[1] A.H. Stroud, Approximate Calculation of Multiple Integrals. Cliffs,
%    NJ: Prentice-Hall, Inc., 1971.
%
%October 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release

r=[-1;1];

%Allocate space for the points.
numGegenbauerPoints=n;
y=zeros(numGegenbauerPoints,numDim-1);
Ak=zeros(numGegenbauerPoints,numDim-1);

%Find all of the yk and Ak values
for k=1:(numDim-1)
    [yCur,AkCur]=quadraturePoints1D(n,3,(k-1)/2);
    
    y(:,k)=yCur(:);
    Ak(:,k)=AkCur;
end

%Allocate space
numPoints=2*numGegenbauerPoints;
xi=zeros(numDim,numPoints);
w=zeros(numPoints,1);

dims=numGegenbauerPoints*ones(numDim-1,1);
curCubPoint=1;
%Go through all of the ranges
for curR=1:2
    curIdx=1;
    i=index2NDim(dims,curIdx);
    
    %Now, we go through all the tuples
    while(~isempty(i))
        xi(numDim,curCubPoint)=r(curR)*y(i(numDim-1),numDim-1);
        wProd=Ak(i(numDim-1),numDim-1);
        
        prodRecur=1;
        for k=(numDim-1):-1:2
            prodRecur=prodRecur*sqrt(1-y(i(k),k)^2);
            xi(k,curCubPoint)=r(curR)*prodRecur*y(i(k-1),k-1);
            wProd=wProd*Ak(i(k-1),k-1);
        end
        prodRecur=prodRecur*sqrt(1-y(i(1),1)^2);
        xi(1,curCubPoint)=r(curR)*prodRecur;
        w(curCubPoint)=wProd;
        
        curCubPoint=curCubPoint+1;
        curIdx=curIdx+1;
        i=index2NDim(dims,curIdx);
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
