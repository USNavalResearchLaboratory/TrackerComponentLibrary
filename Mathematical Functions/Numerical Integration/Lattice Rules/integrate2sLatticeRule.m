function [intVal,E]=integrate2sLatticeRule(f,numFDim,xi,w,errParams)
%%INTEGRATE2SLATTICERULE Perform numerical integration of a function given
%                   lattice points/ weights an embedded 2^s lattice rule,
%                   as are described in Chapter 10 of [1]. Usually, one
%                   will obtain the points using the lattice points and
%                   weights using the standardLattice2SRulesfunction and
%                   will pass them here.
%
%INPUTS: f A handle to the function over which numerical integration is to
%          be performed. This must take a column vector as the input.
%  numFDim The number of dimensions of the output of f.
%   xi,w, errParams The output parameters of the standardLattice2SRules
%          function that are to be used fornumerical integration. This
%          includes lattice points xi, associated weights w and a set of
%          parameters that are needed to estimate the error of the
%          integral.
%
%OUTPUTS: intVal The estimated value of the integral from numerical
%                integration. This is a numFDimX1 vector.
%              E The estimated error of the integral. If a non-periodic
%                function is used with non-periodized points, then the
%                error is likely to be too optimistic.
%
%This implements part of the algorithm of Chapter 10.4 of [1]. The error
%computation part is described in Chapter 10.3.
%
%EXAMPLE 1:
%Here, we consider integration over a region that is not the unit cube. We
%consider the function sum(abs(x-1/2)) in five dimensions with bounds as
%follows:
% f=@(x)sum(abs(x-1/2),1);
% numFDim=1; %(The output of the function is 1D.
% lowerBounds=[-2;-1;0;4;3];
% upperBounds=[2;1;1;6;8];
% %The exact solution against which we will compare is easily found:
% exactSol=915;
% %Let us just use the first formula available.
% methodIdx=1;
% numDim=5;
% %We first consider it without periodizing the points.
% [xi,w,errParams]=standardLattice2SRules(numDim,methodIdx,[],upperBounds,lowerBounds);
% %The integral estimate is
% [est1,E]=integrate2sLatticeRule(f,numFDim,xi,w,errParams);
% error1=abs(exactSol-est1)
% E1=E
% %Now we try it again periodizing the points (since the function over
% %which integration is performed is not periodic).
% periodizeAlg=[];
% periodizeAlg.method=0;
% [xi,w,errParams]=standardLattice2SRules(numDim,methodIdx,periodizeAlg,upperBounds,lowerBounds);
% %The integral estimate is
% [est2,E]=integrate2sLatticeRule(f,numFDim,xi,w,errParams);
% error2=abs(exactSol-est2)
% E2=E
%In the above examples, it can be seen that periodizing the points improves
%the accuracy and that the error estimate of the unperiodized function is
%too optimistic, whereas the error estimate of the periodized function is
%too pessimistic.
%
%EXAMPLE 2:
%In this example, we consider a function with a vector output.
% f=@(x)[x(1)^2+x(2);exp(x(1)^2*x(2))];
% numFDim=2;%(The output of the function is 1D.
% %This time, we will assume the bounds of integration are just from 0 to 1
% %in each dimension.
% exactSol=[5/6;1.2070216633553179822478097023613];
% methodIdx=1;
% numDim=2;
% %We first consider it without periodizing the points.
% [xi,w,errParams]=standardLattice2SRules(numDim,methodIdx);
% %The integral estimate is
% [est1,E]=integrate2sLatticeRule(f,numFDim,xi,w,errParams);
% error1=abs(exactSol-est1)
% E1=E
% %Now we try it again periodizing the points (since the function over
% %which integration is performed is not periodic).
% periodizeAlg=[];
% periodizeAlg.method=0;
% [xi,w,errParams]=standardLattice2SRules(numDim,methodIdx,periodizeAlg);
% %The integral estimate is
% [est2,E]=integrate2sLatticeRule(f,numFDim,xi,w,errParams);
% error2=abs(exactSol-est2)
% E2=E
%Again,  periodization improved the results and the error estimate when
%from being too optimistic to too pessimistic.
%
%REFERENCES:
%[1] I. H. Sloan and S. Joe, Lattice Methods for Multiple integration.
%    Oxford: Clarendon Press, 1994.
%
%November 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%Assume has already been periodized and transformed to the correct region.

%We assume a 2^s lattice rule.
n=2;

s=size(xi,1);
numPoints=size(xi,2);

m=numPoints/2^s;

if(fix(m)~=m)
   error('The number of points is inconsistent with an N^s lattice rule.')
end

%The computation of the Q terms is done without the normalizing by the
%number of points.
w=w*numPoints;

%Q holds the integration result for a lattice rule having up to n^s*m
%points. There are s+1 points to deal with the presence of Q_0.
Q=zeros(numFDim,s+1);
%QParenth holds the integration result for a lattice rule having up to
%n^(s-1)*m points. The points for QParenth are embedded in those for q.
QParenth=zeros(numFDim,s);
for curPoint=1:numPoints
    val=w(curPoint)*f(xi(:,curPoint));
    
    omega=errParams.omegaStart(curPoint);
    Q(:,1+(omega:s))=bsxfun(@plus,Q(:,1+(omega:s)),val);
    sel=errParams.includeQParenth(:,curPoint);
    QParenth(:,sel)=bsxfun(@plus,QParenth(:,sel),val);
end

%The extraction of the error parameter and the integral value is as in
%Chapter 10.4 of [1].
for i=0:s
    Q(:,i+1)=Q(:,i+1)/(m*n^i);
end
E=0;
for i=1:s
    QParenth(:,i)=QParenth(:,i)/(m*n.^(s-1));
    E=E+(Q(:,s+1)-QParenth(:,i)).^2;
end
E=sqrt(E./s);

intVal=Q(:,s+1);

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
