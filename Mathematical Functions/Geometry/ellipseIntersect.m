function [pList,intersectionCode]=ellipseIntersect(param1,param2,param3,param4,param5,param6)
%%ELLIPSEINTERSECT Determine whether two 2D ellipses overlap. If the
%                  ellipses intersect, then return the point(s) of
%                  intersection. If they do not, indicate whether one
%                  ellipse is enclosed in the other.
%
%INPUTS:  The function behaves differently depending upon whether two or
%         six parameters are passed. If two parameters are passed, then
%         these are
%         param1: A coefficient vector, denoted subsequently as a, that
%                 specifies the coefficients of an ellipse corresponding to
%                 the form
%                 a(1)*x^2+a(2)*y^2+a(3)*x*y+a(4)*x+a(5)*y+a(6)=0
%                 Note that for a to represent an ellipse a1*a2-(a3/2)^2>0.
%         param2: A coefficient vector, like param1, which specified the
%                 coefficients of another 2D ellipse.
%         If six parameters are passed, then instead of being polynomials
%         of an ellipse, the parameters specify ellipses in terms of
%         symmetric matrices, means (the center of the ellipse) and
%         thresholds. The parameters are
%         param1, param2, param3, which shall subsequently be referred to
%         as mu, P and gamma are a 2X1 vector, a symmetrix positive
%         definite matrix, and a scalar threshold. the first 3 parameters
%         specify an ellipse such that (z-mu)'*P*(z-mu)=gamma.
%         param5, param6, and param7 are respectively the mean, symmetric
%         matrix and scalar threshold specifying a second ellipse. 
%
%OUTPUTS: pList A matrix where each column is an [x;y] point specifying an
%               intersection location of the two ellipses. If the ellipses
%               do not intersect, then pList is empty.
% intersectionCode Intersection code describes how the ellipses intersect.
%               0  means that they do not intersect.
%               -1 means that one ellipse completely engulfs the other.
%               A positive integer is the number of intersecting points.
%
%The algorithm is based in part on information in [1]. The basic idea is
%that the two ellipses can be written as two bivariate quadratic
%polynomials and solutions to intersection will correspond to real
%solutions to those polynomials.
%
%REFERENCES:
%[1] D. Eberly. (2011, 26 Sep.) Intersection of ellipses. [Online].
%    Available: http://www.geometrictools.com/Documentation/IntersectionOfEllipses.pdf
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin==2)
    a=param1;
    b=param2;
    mu1=Poly2Center(a);
    mu2=Poly2Center(b);
else
    a=MatVec2Poly(param1,param2,param3);
    b=MatVec2Poly(param4,param5,param6);
    mu1=param1;
    mu2=param4;
end

%Allocate space for the maximum number of solutions.
pList=zeros(2,4);

d0=a(1)*b(3)-b(1)*a(3);
d1=a(1)*b(2)-b(1)*a(2);
d2=a(1)*b(4)-b(1)*a(4);
d3=a(1)*b(5)-b(1)*a(5);
d4=a(1)*b(6)-b(1)*a(6);
d5=a(3)*b(2)-b(3)*a(2);
d6=a(3)*b(5)-b(3)*a(5);
d7=a(3)*b(6)-b(3)*a(6);
d8=a(2)*b(4)-b(2)*a(4);
d9=a(4)*b(5)-b(4)*a(5);
d10=a(4)*b(6)-b(4)*a(6);

c0=d2*d10-d4^2;
c1=d0*d10+d2*(d7+d9)-2*d3*d4;
c2=d0*(d7+d9)+d2*(d6-d8)-d3^2-2*d1*d4;
c3=d0*(d6-d8)+d2*d5-2*d1*d3;
c4=d0*d5-d1^2;

vVals=roots([c4;c3;c2;c1;c0]);

%Discard imaginary solutions.
vVals=vVals(vVals==real(vVals));
numV=length(vVals);

curP=0;
for curV=1:numV
    v=vVals(curV);
    alpha2=a(1);
    alpha1=a(3)*v+a(4);
    alpha0=a(5)*v+a(6)+a(2)*v^2;
    
    uVals=roots([alpha2;alpha1;alpha0]);
    %Discard imaginary solutions.
    uVals=uVals(uVals==real(uVals));
    
    if(~isempty(uVals))
        %Only one of the solutions, if real, can satisfy the other
        %equation. Ideally, the other equation should be zero. However, due
        %to rounding errors, we shall just choose the smallest solution.
        [~,idx]=min(abs(g(b,uVals,v)));
        
        curP=curP+1;
        pList(:,curP)=[uVals(idx);v];
    end
end

%Only return the intersection points.
pList=pList(:,1:curP);

%In some instances, such as when there is only one intersection point,
%pList will contain duplicates. Thus, we want to remove duplicate columns.
pList=unique(pList','rows')';

if(isempty(pList))
    intersectionCode=0;

    %Determine whether one ellipse is enclosed in the other.
    val=a(1)*mu2(1)^2+a(2)*mu2(2)^2+a(3)*mu2(1)*mu2(2)+a(4)*mu2(1)+a(5)*mu2(2)+a(6);
    if(val<=0)
        %If the center of the second ellipse is in the first.
        intersectionCode=-1;
        return;
    end
    
    val=b(1)*mu1(1)^2+b(2)*mu1(2)^2+b(3)*mu1(1)*mu1(2)+b(4)*mu1(1)+b(5)*mu1(2)+b(6);
    if(val<=0)
        %If the center of the first ellipse is in the second.
        intersectionCode=-1;
        return;
    end    
else
    intersectionCode=size(pList,2);
end

end

function vals=g(b,u,v)
    beta2=b(1);
    beta1=b(3)*v+b(4);
    beta0=b(5)*v+b(6)+b(2)*v.^2;
    
    vals=beta2.*u.^2+beta1.*u+beta0;
end

function a=MatVec2Poly(mu,P,gamma)
%%MATVEC2POLY Given an ellipse expressed in terms of an uncertainty region
%             in two dimensions using a matrix (such as a covariance
%             matrix) P, a 2X1 vector mu (the center of the ellipse or the
%             mean of an uncertainty region), and a threshold gamma, 
%             determine the corresponding coeficients of the polynomial
%             associated with the ellipse.
%
%INPUTS:  mu A 2X1 vector specifying the center of the ellipse.
%          P A 2X2 symmetric, positive definite matrix.
%      gamma A scalar threshold.
%
%OUTPUTS: The coefficients of an ellipse in the form corresponding to
%         a(1)*x^2+a(2)*y^2+a(3)*x*y+a(4)*x+a(5)*y+a(6)=0
%         Note that for a to represent an ellipse a1*a2-(a3/2)^2>0.
%
%The parameters mu P and gamma specify an ellipse given by the equation
%(z-mu)'*inv(P)*(z-mu)=gamma
%where z is a 2X1 vector.
%
%Using the return format of this function, the points inside fo the ellipse
%are such that
%a(1)*x^2+a(2)*y^2+a(3)*x*y+a(4)*x+a(5)*y+a(6)<=0
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
    
%     a(1)=P(1,1);
%     a(2)=P(2,2);
%     a(3)=2*P(1,2);
%     a(4)=-2*(mu(1)*P(1,1)+mu(2)*P(1,2));
%     a(5)=-2*(mu(1)*P(1,2)+mu(2)*P(2,2));
%     a(6)=mu(1)^2*P(1,1)+2*mu(1)*mu(2)*P(1,2)+mu(2)^2*P(2,2)-gamma;

a(1)=P(1,1); 
a(2)=P(2,2);
a(3)=2*P(1,2);
a(4)=-2*(mu(1)*P(1,1)+mu(2)*P(1,2));
a(5)=-2*(mu(1)*P(1,2)+mu(2)*P(2,2));
a(6)=mu(1)^2*P(1,1)+2*mu(1)*mu(2)*P(1,2)+mu(2)^2*P(2,2)-gamma;

end

function mu=Poly2Center(a)
%%POLY2MEAN  Given the coefficients of the polynomial of an ellipse in two
%            dimensions, obtain the mean implied by the polynomial.
%
%INPUTS:  a The coefficients of an ellipse in the form corresponding to
%           a(1)*x^2+a(2)*y^2+a(3)*x*y+a(4)*x+a(5)*y+a(6)=0
%           Note that for a to represent an ellipse a1*a2-(a3/2)^2>0.
%
%OUTPUTS: mu The mu=[x,y]; coordinate center of the ellipse.
%
%The simplest way to determine the center of the ellipse is to start with
%the matrix form as in the function MatVec2Poly, expand the polynomial form
%and then solve for the mean. That is how the equations below were
%obtained.
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.

    mu(1)=(a(3)*a(5)-2*a(2)*a(4))/(4*a(1)*a(2)-a(3)^2);
    mu(2)=(a(3)*a(4)-2*a(1)*a(5))/(4*a(1)*a(2)-a(3)^2);
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
