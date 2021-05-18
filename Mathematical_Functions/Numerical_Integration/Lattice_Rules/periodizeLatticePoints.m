function [xi,w]=periodizeLatticePoints(method,xi,w,param1,param2)
%%PERIODIZELATTICEPOINTS Many lattice methods for multivariate integration
%               between 0 and 1 are designed for functions that are
%               periodic in every coordinate with the interval 0-1 being
%               the first period. Non-periodic functions can also be used
%               in certain instances depending on the boundary conditions
%               at 0 and 1. For more general non-periodic functions, the
%               lattice points can be used aftering being transformed. The
%               transformation changes moves the points and changes the
%               weighting from a uniform weighting to a non-uniform
%               weighting.
%
%INPUTS: method An integer selecting the periodizing transformation to
%               perform on the lattice points. The meaning of param1 and
%               param2 in the input changes depending on the algorithm
%               selected. Possible values are
%               0 Use the sin^m method of [2]. param2 is not used and
%                 param1=m. m is an integer>=1 and small number 1-4
%                 generally work well. Even numbers are better than odd. If
%                 param1 is omitted, m=2 is used.
%               1 Use the Korovov transform mentioned in [2], but with a
%                 corrected normalization constant as the function in the
%                 paper does nto evaluate to 1 when the input is 1. param2
%                 is not used and param1 is m, which determines the order
%                 of the polynomial. m>0 and m is an integer. Results are
%                 better for even m and small values tend to work well. The
%                 default if param1 is omitted is m=2.
%               2 Use the tanh transformation described in [2]. param2 is
%                 not used and param1 is c, which changes the shape. c>0.
%                 If param1 is omitted, then c=1 is used.
%               3 Use the double exponential transform in [2]. param1 is a
%                 of the transform and param2 b. a,b>0. b is inside of a
%                 sinh function and a is inside of a tanh function. If
%                 either is omitted, then a value of pi/2 will be used.
%                 This tends to be a fairely bad method.
%            xi The numDimXnumPoints set of lattice points for
%               approximating multivariate numerical integrals of a
%               periodic function over 0-1 in all dimensions that are to be
%               transformed to work with non-periodic functions. All of the
%               components of all of the points should be in the region
%               [0,1].
%             w Optionally, a numPointsX1 set of weights for the points in
%               xi. If all of the points have equal weighting being
%               1/numPoints, which is usually the case, then an empty
%               matrix can be passed or this parameter can be omitted.
% param1,param2 The parameters whose meaning is determined by the selected
%               method. Default values will be used if these are omitted.
%
%OUTPUTS: xi The transformed lattice points that can be used for
%            multivariate integration over a non-periodic function.
%          w A numPointsX1 set of weights associated with the points. These
%            are returned even if w was not provided on the input. These
%            are important for properly evaluating the periodized integral.
%
%The transformation method used is based on the discussion in Chapter 2.12
%of [1]. Numerical integration will be a sum of the function evaluated at
%the points xi times the weights w.
%
%Note that transforming lattice points for use with non-periodic functions
%often results in lower accuracy results than lattice or cubature methods
%that have been specially derived for multivariate integration of
%non-periodic function.
%
%REFERENCES:
%[1] I. H. Sloan and S. Joe, Lattice Methods for Multiple Integration.
%    Oxford, United Kingdom: Clarendon Press, 1994.
%[2] A. Sidi, "A new variable transformation for numerical integration," in
%    Numerical Integration IV: Proceedings of the Conference at the
%    Mathematical Research Institute, Oberwolfach, Germany, 8-14 Nov. 1992,
%    pp. 359-373.
%
%November 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numPoints=size(xi,2);
%If w is not provided, the assume uniform weight.
if(nargin<3||isempty(w))
    w=(1/numPoints)*ones(numPoints,1);
end

switch(method)
    case 0%The sim^m transformation of [2]. Results are better for even m.
        if(nargin<4||isempty(param1))
            param1=2;
        end
        mMax=param1;
        phi=@(t)sinmFunction(mMax,t);
        phiDeriv=@(t)sinmFunctionDeriv(mMax,t);
    case 1%Korovov transform mentioned in [2]. m>0 and m is an integer.
        %The integral in the paper has been evaluated and the normalization
        %constant has ben corrected. Results are better for even m.
        if(nargin<4||isempty(param1))
            param1=2;
        end
        m=param1;
        phi=@(t)(1+(1-t).^m.*(t-1).*(1+t+m*t));
        phiDeriv=@(t)((1+m)*(2+m)*(1-t).^m.*t);
    case 2%The tanh transformation, c>0 in [2].
        if(nargin<4||isempty(param1))
            param1=1;
        end
        c=param1;
        phi=@(t)((1/2)*tanh(-(c/2)*(1./t-1./(1-t)))+1/2);
        phiDeriv=@(t)tanhTransformDeriv(c,t);
    case 3%The double exponential transform in [2], a,b>0.
        if(nargin<4||isempty(param1))
            param1=pi/2;
        end
        if(nargin<5||isempty(param2))
            param2=pi/2;
        end
        a=param1;
        b=param2;
        phi=@(t)((1/2)*tanh(a*sinh(b*(1./(1-t)-1./t)))+1/2);
        phiDeriv=@(t)doubleExpTransformDeriv(a,b,t);
    otherwise
        error('Unknown lattice periodization method specified')
end

for k=1:numPoints
    %Transform the point.
    w(k)=w(k)*prod(phiDeriv(xi(:,k)));
    xi(:,k)=phi(xi(:,k)); 
end
end

function dPsi=tanhTransformDeriv(c,t)
%%TANHTRANSFORMDERIV The derivative of the tanh transform in Equation 1.11
%             of [1]. Care has to be taken for values near 0 and 1, because
%             a 0*Inf arises. The asymptotic value (of 0) is substituted in
%             such instances.
%
%INPUTS:    c A constant affecting the shape of the function; c>0.
%           t The point(s) at which the derivative of the function should
%             be evaluated.
%
%OUTPUTS: dPsi The value(s) of the derivative of the tanh
%             transformation at the point(s) in t.
%
%REFERENCES:
%[1] A. Sidi, "A new variable transformation for numerical integration," in
%    Numerical Integration IV: Proceedings of the Conference at the
%    Mathematical Research Institute, Oberwolfach, Germany, 8-14 Nov. 1992,
%    pp. 359-373.
%
%November 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

dPsi=(c/4)*(1./(1-t).^2+1./t.^2).*sech((c-2*c*t)./(2*t.*(1-t))).^2;

dPsi(~isfinite(dPsi))=0;
end

function dPsi=doubleExpTransformDeriv(a,b,t)
%%DOUBLEEXPTRANSFORMDERIV The serivative of the double exponential
%             transform in Equation 1.15 of [1]. Care has to be taken for
%             values near 0 and 1, because a 0*Inf arises. The asymptotic
%             value (of 0) is substituted in
%             such instances.
%
%INPUTS:  a,b Constants affecting the shape of the function; a,b>0.
%           t The point(s) at which the derivative of the function should
%             be evaluated.
%
%OUTPUTS: dPsi The value(s) of the derivative of the doble exponential
%             transformation at the point(s) in t.
%
%REFERENCES:
%[1] A. Sidi, "A new variable transformation for numerical integration," in
%    Numerical Integration IV: Proceedings of the Conference at the
%    Mathematical Research Institute, Oberwolfach, Germany, 8-14 Nov. 1992,
%    pp. 359-373.
%
%November 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

dPsi=a*b*(1+2*(t-1).*t).*cosh(b*(1./(t-1)+1./t)).*sech(a*sinh(b*(1./(t-1)+1./t))).^2./(2*t.^2.*(t-1).^2);

dPsi(~isfinite(dPsi))=0;
end

function psi=sinmFunction(mMax,t)
%%SINMFUNCTION Evaluate the sin^m transformation defined in Equation 2.1 of
%              [1].
%
%INPUTS: mMax The order of the function. mMax>=1.
%           t The point(s) at which the function should be evaluated.
%
%OUTPUTS:psi The value(s) of the sin^m transformation at the point(s)
%             in t.
%
%REFERENCES:
%[1] A. Sidi, "A new variable transformation for numerical integration," in
%    Numerical Integration IV: Proceedings of the Conference at the
%    Mathematical Research Institute, Oberwolfach, Germany, 8-14 Nov. 1992,
%    pp. 359-373.
%
%November 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

psi=sinMThetaFunc(mMax,t)/sinMThetaFunc(mMax,1);
end

function ThetaVal=sinMThetaFunc(mMax,t)
%%SINMTHETAFUNCTION Evaluate the Theta function for the sin^m
%              transformation defined in [1], using the recursion of
%              Equation 2.3 of [1] with the initial conditions in Equation
%              2.5.
%
%INPUTS: mMax The order of the function. mMax>=1.
%           t The point(s) at which the function should be evaluated.
%
%OUTPUTS: ThetaVal The value(s) of the Theta function going into the sin^m
%             transformation evaluated at the point(s) in t.
%
%REFERENCES:
%[1] A. Sidi, "A new variable transformation for numerical integration," in
%    Numerical Integration IV: Proceedings of the Conference at the
%    Mathematical Research Institute, Oberwolfach, Germany, 8-14 Nov. 1992,
%    pp. 359-373.
%
%November 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

cosVal=cos(pi*t);
sinVal=sin(pi*t);
sinCosProd=sinVal.*cosVal;

ThetaPrev=t;
ThetaVal=(1/pi)*(1-cosVal);

for m=2:mMax
    ThetaPrev2=ThetaPrev;
    ThetaPrev=ThetaVal;
    
    ThetaVal=-1/(pi*m)*sinCosProd+((m-1)/m)*ThetaPrev2;
    sinCosProd=sinCosProd.*sinVal;
end
end


function dPsi=sinmFunctionDeriv(mMax,t)
%%SINMFUNCTIONDERIV Evaluate the derivative of the sin^m transformation 
%              defined in [1], using the derivative of the  recursion of
%              Equation 2.3 of [1] with the derivatives of the initial
%              conditions in Equation 2.5.
%
%INPUTS: mMax The order of the function. mMax>=1.
%           t The point(s) at which the derivative of the function should
%             be evaluated.
%
%OUTPUTS: dPsi The value(s) of the derivative of the sin^m 
%             transformation at the point(s) in t.
%
%REFERENCES:
%[1] A. Sidi, "A new variable transformation for numerical integration," in
%    Numerical Integration IV: Proceedings of the Conference at the
%    Mathematical Research Institute, Oberwolfach, Germany, 8-14 Nov. 1992,
%    pp. 359-373.
%
%November 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

cosVal=cos(2*pi*t);
sinVal=sin(pi*t);
sinProd=1;

dThetaPrev=1;
dThetaCur=sinVal;

for m=2:mMax
    dThetaPrev2=dThetaPrev;
    dThetaPrev=dThetaCur;
    
    dThetaCur=(-0.5*(-2+m+m*cosVal).*sinProd+(m-1)*dThetaPrev2)./m;
    sinProd=sinProd.*sinVal;
end

%Deal with the division by a constant in Equation 2.1.
dPsi=dThetaCur/sinMThetaFunc(mMax,1);
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
