function [S,C]=FresnelInt(x)
%%FRESNELINT Evaluate both types of Fresnel integrals for a real argument
%            x. The Fresnel integrals C(x)=integral_0^x cos((pi/2)*t^2 dt
%            and S(x)=integral_0^x sin((pi/2)*t^2 dt are both evaluated.
%            Fresnel integrals. Fresnel integrals arise in computations of
%            electromagnetic field intensity. Fresnel integrals also arise
%            in computing the spectrum of a linearly frequency modulated
%            (LFM) signal.
%
%INPUTS: x A real scalar, vector, or matrix of values at which one wishes
%          to evaluate the Fresnel integrals.
%
%OUTPUTS: S A matrix having the same dimension as x holding the values of
%           the Frensel sine integral.
%         C A matrix having the same dimension as x holding the values of
%           the Frensel cosine integral.
%
%This function uses the identity given in 7.3.22 of [1] to compute the
%Fresnel integral from the complex error function.
%
%REFERENCES:
%[1] Abramowitz, M. and Stegun, I. A. (Eds.). "Error Function and Fresnel
%    Integrals." in Ch. 7 in Handbook of Mathematical Functions with
%    Formulas, Graphs, and Mathematical Tables, 9th printing. New York:
%    Dover, 1972.
%
%November 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numVals=numel(x);

C=zeros(size(x));
S=zeros(size(x));

xList=x;

for curVal=1:numVals
    x=xList(curVal);

    z=(sqrt(pi)/2)*(1-1i)*x;
    val=((1+1i)/2)*erfComplex(z);
    C(curVal)=real(val);
    S(curVal)=imag(val);
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
