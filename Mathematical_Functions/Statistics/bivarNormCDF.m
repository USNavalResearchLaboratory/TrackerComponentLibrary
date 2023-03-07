function val=bivarNormCDF(b,mu,R)
%%BIVARNORMCDF Evaluate the cumulative density function of the bivariate
%              normal distribution with a specified mean and covariance
%              matrix. This evaluates Pr{x1<b(1), x2<b(2)) where the random
%              vector is [x1;x2].
%
%INPUTS: b A 2X1 vector [b1;b2] such that b1 is the upper bound of the
%           first variable and b2 is the upper bound on the second
%           variable.
%        mu The 2X1 mean of the distribution. If this is omitted or an
%           empty matrix is passed, then the default of [0;0] is used.
%         R The 2X2 covariance matrix of the distribution. If this is
%           omitted or an empty matrix is passed, then R=eye(2,2) is used.
%
%OUTPUTS: val The value of the bivariate normal CDF.
%
%This function implements the algorithm described in Section 2.3 and 2.4 of
%[1], which should be approximately accurate within the double floating
%point finite precision limits. A change of variables is used to change the
%bounds into those for the bivariate normal distribution with zero mean and
%a covariance matrix with 1s on the diagonal.
%
%REFERENCES:
%[1] A. Genz, "Numerical computation of rectangular bivariate and
%    trivariate normal and t probabilites," Statistics and Computing, vol.
%    14, pp. 251-260, 2004.
%
%July 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public releas

if(nargin<2||isempty(mu))
    mu=[0;0];
end

if(nargin<3||isempty(R))
    R=eye(2,2);
end

b1=b(1);
b2=b(2);
sigma1=sqrt(R(1,1));
sigma2=sqrt(R(2,2));
b1=(b1-mu(1))/sigma1;
b2=(b2-mu(2))/sigma2;
rho=R(1,2)/(sigma1*sigma2);

h=-b1;
k=-b2;
s=2*(rho>=0)-1;

%Points to use for cubature integration.
[xi,w]=GaussLegendrePoints1D(20);

if(abs(rho)<0.925)
    %Integrate Equation 3 in [1] using  a 20-point Gauss-Legendre formula.
    %The integral is from 0 to asin(rho), so the points must first be
    %transformed.
    upperBound=asin(rho);
    xi=(xi+1)*(upperBound/2);%Change of variables.
    %The integral in equation 3. The outer multiplying term deals with the
    %change of variables.
    intVal=(upperBound/2)*sum(bsxfun(@times,w(:).',exp(-(h.^2+k.^2-2.*h.*k.*sin(xi))./(2*cos(xi).^2))));
    
    if(isnan(intVal))
        intVal=0;
    end

    val=GaussianD.CDF(-h)*GaussianD.CDF(-k)+(1/(2*pi))*intVal;
else
    %Use the solution of Equation 6 in [1] with the modifications of
    %Section 2.4

    %The unnumbered equation before 4 in [1].
    if(s==1)
        Lhks=GaussianD.CDF(-max(h,k));
    else
        Lhks=max(0,GaussianD.CDF(-h)-GaussianD.CDF(k));
    end

    %The value of the first integral in Equation 6 in [1] with the higher-
    %order Taylor series expansion of Section 2.4. This is just another way
    %of writing the solution given in the unnumbered equation after Eq. 6.
    a=sqrt(max(0,1-rho^2));%The max helps with finite precision issues.
    b=abs(h-s*k);
    c=(4-s*h*k)/8;
    d=(12-s*h*k)/16;
    coeff=exp(-s*h*k/2);
    intVal1=coeff*(1/30)*(2*a*(15+5*(a-b)*(a+b)*c+(3*a^4-a^2*b^2+b^4)*c*d)*exp(-(b^2/(2*a^2)))-b*(15-5*b^2*c+b^4*c*d)*sqrt(2*pi)*erfc(b/(sqrt(2)*a)));
    
    if(isnan(intVal1))
        intVal1=0;
    end
    
    %Numeric integration of the second integral in Equation 6 in [1] with
    %the higher- order Taylor series expansion of Section 2.4. The integral
    %is from 0 to asin(rho), so the points must first be transformed.
    xi=(xi+1)*sqrt(max(0,1-rho^2))/2;%Max to avoid finite precision issues.
    %The integral in equation 6 (with the extra terms). The multiplying
    %term outside the sum deals with the change of variables.
    intVal2=(a/2)*sum(bsxfun(@times,w(:).',exp(-b^2./(2*xi.^2)).*(exp(-s*h*k./(1+sqrt(1-xi.^2)))./sqrt(1-xi.^2)-coeff.*(1+c.*xi.^2.*(1+d.*xi.^2)))));    
    
    if(isnan(intVal2))
        intVal2=0;
    end
    
    %Completing Equation 6.
    val=Lhks-(s/(2*pi))*(intVal1+intVal2);
end

%In case finite-precision limitations make the value invalid:
val=min(max(val,0),1);

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
