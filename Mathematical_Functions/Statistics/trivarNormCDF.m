function probVal=trivarNormCDF(b,mu,R)
%%TRIVARNORMCDF Evaluate the cumulative density function of the trivariate
%              normal distribution with a specified mean and covariance
%              matrix. This evaluates Pr{x1<b(1), x2<b(2), x3<b(3)) where
%              the random vector is [x1;x2;x3].
%
%INPUTS: b A 2X1 vector [b1;b2] such that b1 is the upper bound of the
%           first variable and b2 is the upper bound on the second
%           variable.
%        mu The 2X1 mean of the distribution. If this is omitted or an
%           empty matrix is passed, then the default of [0;0] is used.
%         R The 2X2 covariance matrix of the distribution. If this is
%           omitted or an empty matrix is passed, then R=eye(2,2) is used.
%
%OUTPUTS: probVal The value of the trivariate normal CDF.
%
%This function implements the second Plankett method that is described in
%[1]. Specifically, Equation 14 is used.
%
%REFERENCES:
%[1] A. Genz, "Numerical computation of rectangular bivariate and
%    trivariate normal and t probabilites," Statistics and Computing, vol.
%    14, pp. 251-260, 2004.
%
%July 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public releas

if(nargin<3||isempty(R))
    R=eye(3,3);
end

if(nargin<2||isempty(mu))
    mu=[0;0;0];
end

%The exteme for one bound.
if(any(b)==-Inf)
    probVal=0;
    return;
end

%If any of the bounds is Inf, then then problem reduces to a 2D integral.
if(b(1)==Inf)
    R=R(2:3,2:3);
    mu=mu(2:3);
    b=b(2:3);
    probVal=bivarNormCDF(b,mu,R);
    return
elseif(b(2)==Inf)
     R=R([1,3],[1,3]);
     mu=mu([1;3]);
     b=b([1;3]);
     probVal=bivarNormCDF(b,mu,R);
     return;
elseif(b(3)==Inf)
    R=R(1:2,1:2);
    mu=mu(1:2);
    b=b(1:2);
    probVal=bivarNormCDF(b,mu,R);
    return
end

%Center the distribution.
b=b-mu;

%Scale the R matrix (and the associated b values) to make the diagonals of
%R all 1.
S=diag(1./sqrt(diag(R)));
b=S*b;
R=S*R*S';

%The correlation values.
pVec=[R(2,1);R(3,1);R(3,2)];

%We want the variables of integration to be permuted to minimize 
%max(abs(p21),abs(p31)). This means that out of the above values, p21, p31,
%and p32, we want to change the ordering so that p21 and p31 are the values
%with the smallest magnitudes.
[~,idx]=sort(abs(pVec),'ascend');
pVec=pVec(idx);
b=b(idx);

b1=b(1);
b2=b(2);
b3=b(3);
p21=pVec(1);
p31=pVec(2);
p32=pVec(3);

R2D=[1,  p32;
     p32,1];
%The first term in Equation 14 in 1.
term1=GaussianD.CDF(b1)*bivarNormCDF([b2;b3],[],R2D);
f=@(t)costFun(t,b1,b2,b3,p21,p31,p32);

RelTol=1e-18;
AbsTol=1e-18;
term2=integral1DAdaptive(f,[0;1],[],[],[],RelTol,AbsTol);
probVal=term1+term2;

%In case finite precision limitations make the value invalid.
probVal=min(max(probVal,0),1);

end

function val=costFun(t,b1,b2,b3,p21,p31,p32)

t2=t.^2;
denomTerm2=1-p31^2*t2-p21^2*t2-p32^2+2*t2*p31*p21*p32;
denomTerm3=1-p21^2*t2-p31^2*t2-p32^2+2*t2*p31*p21*p32;
u2Hat=(b2*(1-p31^2*t2)-b1*t*(p21-p31*p32)-b3*(p32-p31*p21*t2))./sqrt((1-p31^2.*t2).*denomTerm2);
u3Hat=(b3*(1-p21^2*t2)-b1*t*(p31-p21*p32)-b2*(p32-p31*p21*t2))./sqrt((1-p21^2.*t2).*denomTerm3);

r=p31*t;
f2=(b1^2+b3^2-2*r*b1*b3)./(1-r.^2);
f2(isnan(f2))=Inf;
r=p21*t;
f3=(b1^2+b2^2-2*r*b1*b2)./(1-r.^2);
f3(isnan(f3))=Inf;
%The above NaN stuff should make it work properly if, for example, b1 and
%b2 are near realmax or realmin.

Phiu3Hat=GaussianD.CDF(u3Hat);
Phiu2Hat=GaussianD.CDF(u2Hat);

val=Phiu3Hat.*p21.*exp(-f3/2)./sqrt(1-p21^2.*t2)+Phiu2Hat.*p31.*exp(-f2/2)./sqrt(1-p31.^2.*t2);

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
