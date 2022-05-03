function [x,info]=minSumL2NormLinear(F,g,algorithm,params)
%%MINSUML2NORMLINEAR This function solves the optimization problem:
%                    minimize sum_{i=1}^p norm(F(:,:,i)*x+g(:,i)) for x.
%
%INPUTS: F An mXnXp set of p real mXn matrices.
%        g An mXp set of p real vectors.
%algorithm An optional parameter specifying which algorithm to use to solve
%          the problem. The problem is solved as a second order cone
%          problem. Possible values are:
%          0 (The default if omitted or an empty matrix is passed) Use the
%            function splittingConicSolver.
%          1 Use the function secondOrderConeInteriorPoint.
%   params Parameters that affect how the function splittingConicSolver,
%          which is used by this function, works. See the comments to
%          splittingConicSolver for more information
%
%OUTPUTS: x The nX1 solution to the problem or an empty matrix is no
%           solution is possible.
%      info The problem is formulated as a second-order cone optimization
%           problem. This is the exitCoe return parameter from
%           secondOrderConeInteriorPoint if  algorithm=0 and this is the
%           information regarding the termination state returned by the
%           splittingConicSolver function if algorithm=1.
%
%As described in [1], the problem can be formulated as a second order cone
%optimization problem.
%
%EXAMPLE 1:
% F=zeros(2,2,3);
% F(:,:,1)=[1,2;
%           3,4];
% F(:,:,2)=[6,7;
%           8,12];
% F(:,:,3)=[1,-1;
%           -1,1];
% g=zeros(2,3);
% g(:,1)=[-45;-89];
% g(:,2)=[-155;-268];
% g(:,3)=[24;-24];
% [x,info]=minSumL2NormLinear(F,g)
%One will find that x=[-1;23];
%
%EXAMPLE 2:
%In some instances, one might want to use F matrices having differing
%numbers of rows between the norms. This can be done by just inserting zero
%rows and zero entries in g so that everything is the same size. For
%example, if we take the case of example 1, but only use the first row of
%the second  matrix, we get
% F=zeros(2,2,3);
% F(:,:,1)=[1,2;
%           3,4];
% F(:,:,2)=[6,7;
%           0,0];
% F(:,:,3)=[1,-1;
%           -1,1];
% g=zeros(2,3);
% g(:,1)=[-45;-89];
% g(:,2)=[-155;0];
% g(:,3)=[24;-24];
% [x,info]=minSumL2NormLinear(F,g)
%And the same answer as before is obtained.
%
%REFERENCES:
%[1] S. Lobo, Miguel, L. Vandenberg, S. Boyd, and H. Lebret, "Applications
%    of second-order cone programming," Linear Algebra and its
%    Applications, vol. 284, no. 1, pp. 193-228, 1998.
%
%February 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(algorithm))
   algorithm=0; 
end

if(nargin<4||isempty(params))
   params=[]; 
end

m=size(F,1);
n=size(F,2);
p=size(F,3);

%p auxiliary variables are used.
c=[zeros(n,1);ones(p,1)];

%Allocate space
A=zeros(m*p+p,n+p);
b=zeros(n+p,1);%Allocate space

e=zeros(n+p,1);

curStart=1;
for curConst=1:p
    span=curStart:(curStart+m);
    e(n+curConst)=1;
    A(span,:)=-[e';
               F(:,:,curConst),zeros(m,p)];
    b(span)=[0;g(:,curConst)];
    
    e(n+curConst)=0;
    curStart=curStart+m+1;
end

cone.q=(m+1)*ones(p,1);

switch(algorithm)
    case 0
        [x,~,~,info]=splittingConicSolver(A,b,c,cone,params);
    case 1
        [~,x,~,info]=secondOrderConeInteriorPoint(A',-c,b,cone.q,params);
    otherwise
        error('Unknown algorithm specified')
end

%Get rid of the auxiliary variables
if(~isempty(x))
    x=x(1:n);
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
