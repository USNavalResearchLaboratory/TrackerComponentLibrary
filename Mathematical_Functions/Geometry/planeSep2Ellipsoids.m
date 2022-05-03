function [a,b,info]=planeSep2Ellipsoids(xHat,R,params)
%%PLANESEP2ELLIPSOIDS Given two ellipsoids, find a hyperplane that
%                 separates them or determine that no such plane exists
%                 (the ellipsoids intersect). The ith ellipsoid is defined
%                 as the points x=xHat(:,i)+R(:,:,i)*u, where u is a vector
%                 with magnitude less than or equal to 1. The comments
%                 below relate the parameters of this problem to the
%                 typical ellipsoid parameterization used for gating in
%                 tracking applications.
%
%INPUTS: xHat A numDimX2 set of the real centers of the two ellipsoids.
%           R A numDimXnumDimX2 set of real shape matrices of two
%             ellipsoids.
%      params An optional parameter. This function uses the
%             splittingConicSolver function to find a splitting (hyper)plane
%             between two ellipsoids, it if exists.
%
%OUTPUTS: a, b The numDimX1 vector a and the scalar b define a hyperplane 
%              a'*x=b such that one ellipsoid is on one side and the other
%              is on the other side. If the ellipsoids intersect, then
%              empty matrices are returned.
%         info This is the termination status returned by the subroutine
%              splittingConicSolver. See the comments to
%              splittingConicSolver for more information.
%
%The form of an ellipsoid used for gating with a Gaussian approximation of
%a PDF in tracking is (x-xHat)*PInv*(x-xHat)<=gamma, where xHat is the
%center of the ellipsoid (the mean of the distribution), PInv
%is the inverse  covariance matrix of the distribution and gamma is a
%threshold that sets the probability region. gamma should be the value of
%a chi-square distribution that gives the desired confidence
%bounds (CDF value); the number of degrees of freedom of the chi-square
%distribution equals the dimensionality of x.
%
%We would like to express the points in the ellipsoid
%(x-xHat)*RInv*(x-xHat)<=gamma in the form xHat+PInv*u where u is a vector
%such that norm(u)<=1. To do this, we note that the form
%(x-xHat)*RInv*(x-xHat)<=gamma can be rewritten as
%u0'*RInv*u0<=gamma where x=u0+xHat defines points in the ellipsoid. Now,
%we can rewrite the quadratic form as norm(Q'*u0/sqrt(gamma))<=1 where
%Q=chol(RInv,'lower'). Thus, we can define u=Q'*u0/sqrt(gamma), which means
%that the points in the ellipsoid are x=xHat+sqrt(gamma)*inv(Q')*u. That
%has the desired form. Thus, to put the elllipsoid in the form used by this
%function, one uses R=sqrt(gamma)*inv(chol(RInv,'upper')).
%
%This function essentially implements the solution to s implified version
%of problem 4.25 in Ch. 4 of [1]. The solution is actually given on 
%https://inst.eecs.berkeley.edu/~ee127a/book/login/exa_ell_sep.html
%This notes that the hyperplane a'*x=b separates two ellipsoids if and only
%if b1=a'*xHat(:,1)+norm(R(:,:,1)'*a)<=b<=a'*xHat(:,2)-norm(R(:,:,2)'*a)=b2
%and such a hyperplane can only exist if
%a'*(xHat(:,2)-xHat(:,1)>=norm(R(:,:,1)'*a)+norm(R(:,:,2)'*a)
%The optimization problem arising from this is
% p=min_a norm(R(:,:,1)'*a)+norm(R(:,:,2)'*a)
% subject to a'*(xHat(:,2)-xHat(:,1))=1
%and the ellipsoids are separable if and only if p<=1. An approrpiate value
%of b is given as the average
%(a'*xHat(:,1)+norm(R(:,:,1)'*a)+a'*xHat(:,2)+norm(R(:,:,2)'*a))/2.
%
% This can be rewritten as
% min_a t1+t2
% subject to a'*x2-a'*x1=1
%             t1>=norm(R1'*a)
%             t2>=norm(R2'*a)
%which is a second order cone program, which can be solved using the
%splittingConicSolver function.
%
%EXAMPLE 1:
%In this example, we use 2D ellipses. The ellipses are not separable. We
%shall define the initial problem in the form typically used in tracking,
%and then transform it to the form used by this function.
% xHat=zeros(2,2);%The distribution means
% xHat(:,1)=[0;0];
% xHat(:,2)=[5;6];
% PInv=zeros(2,2,2);%The distribution's inverse covariance matrices.
% PInv(:,:,1)=inv([1,  -0.5;
%                -0.5, 1]);
% PInv(:,:,2)=inv([2,  0.5;
%                 0.5,1]);
% gammaVal=ChiSquareD.invCDF(0.9997,2);%The 99.97% confidence region.
% %Plot the ellipses
% figure(1)
% clf
% drawEllipse(xHat,PInv,gammaVal)
% R=zeros(2,2,2);%The matrices for this function.
% R(:,:,1)=sqrt(gammaVal)*inv(chol(PInv(:,:,1),'upper'));
% R(:,:,2)=sqrt(gammaVal)*inv(chol(PInv(:,:,2),'upper'));
% [a,b,info]=planeSep2Ellipsoids(xHat,R)
%
%EXAMPLE 2:
%In this example, the ellipsoids are separable and the line separating them
%is plotted.
% xHat=zeros(2,2);%The distribution means
% xHat(:,1)=[-0.5;-0.5];
% xHat(:,2)=[5;6];
% PInv=zeros(2,2,2);%The distribution's inverse covariance matrices.
% PInv(:,:,1)=inv([1,  -0.5;
%                -0.5, 1]);
% PInv(:,:,2)=inv([2,  0.5;
%                 0.5,1]);
% gammaVal=ChiSquareD.invCDF(0.9997,2);%The 99.97% confidence region.
% %Plot the ellipses
% figure(1)
% clf
% hold on
% drawEllipse(xHat,PInv,gammaVal)
% R=zeros(2,2,2);%The matrices for this function.
% R(:,:,1)=sqrt(gammaVal)*inv(chol(PInv(:,:,1),'upper'));
% R(:,:,2)=sqrt(gammaVal)*inv(chol(PInv(:,:,2),'upper'));
% [a,b,info]=planeSep2Ellipsoids(xHat,R)
% %We will plot the separating line.
% xVals=linspace(-6,12,500);
% yVals=(-a(1)/a(2))*xVals+b/a(2);
% plot(xVals,yVals)
%
%REFERENCES:
%[1] S. Boyd and L. Vandenberghe, "Convex optimization," Cambridge, United
%    Kingdom, 2004.
%
%February 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3)
   params=[];%Use the default values.
   params.max_iters=50e3;
end

x1=xHat(:,1);
x2=xHat(:,2);
R1=R(:,:,1);
R2=R(:,:,2);

n=size(x1,1);%The dimensionality of the problem.

d1=[zeros(n,1);1;0];
A1=[R1',zeros(n,2)];

d2=[zeros(n,1);0;1];
A2=[R2',zeros(n,2)];

c=[zeros(n,1);1;1];
A=[[(x2-x1)',0,0];
    -d1';
    -A1;
    -d2';
    -A2];

b=[1;zeros(2*(n+1),1)];

cone.f=1;
cone.q=[n+1;n+1];

[x,~,~,info]=splittingConicSolver(A,b,c,cone,params);

a=x(1:n);

p=norm(R1'*a)+norm(R2'*a);

if(p<=1)%If the ellipsoids are separable.
    b1=a'*x1+norm(R1'*a);
    b2=a'*x2-norm(R2'*a);
    b=(b1+b2)/2;
else
    a=[];
    b=[];
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
