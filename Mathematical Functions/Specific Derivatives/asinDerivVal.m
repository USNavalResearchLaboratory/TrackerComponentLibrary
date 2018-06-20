function [vals,coeffs]=asinDerivVal(z,n,coeffs)
%%ASINDERIVVAL Return the value of the nth derivative of asin(x) with
%              respect to x  evaluated at x=z for real or imaginary z.
%
%INPUTS: z A matrix of real or complex values at which the nth derivative
%          of the arcsine function is desired. If an empty matrix is
%          passed, then just coeffs will be returned.
%        n The number of derivatives to take. n>=0.
%   coeffs The computation of the derivatives involves determining the
%          structure of a polynomial. If this function has been run before
%          for a given n value, then the polynomial can be passed back and
%          is not determined again. Otherwise, this can be omitted or an
%          empty matrix can be passed. This only makes a difference for
%          n>=12.
%
%OUTPUTS: vals The value of the nth derivative of the arcsine function
%              taken at all of the points in z. vals has the same
%              dimensions as z.
%       coeffs A vector of polynomial coefficients that can be passed back
%              to this function when evaluating with the same n to speed it
%              up. 
%
%For n=0, vals=asin(z). Subsequent derivatives are computed essentially
%using the chain rule and combining terms to get the appropriate
%polynomial. All derivatives are of the form
%polynomial*sqrt(1-z(:).^2).^(-k). For n<12, the coefficients are tabulated
%to make the function faster.
%
%May 2018 David F.Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(n==0)
    vals=asin(z);
    coeffs=[];
    return;
end

if(nargin<3||isempty(coeffs))
    if(n<12)
        switch(n)
            case 1
                coeffs=1;
            case 2
                coeffs=[1;0];
            case 3
                coeffs=[2;0;1];
            case 4
                coeffs=[6;0;9;0];
            case 5
                coeffs=[24;0;72;0;9];
            case 6
                coeffs=[120;0;600;0;225;0];
            case 7
                coeffs=[720;0;5400;0;4050;0;225];
            case 8
                coeffs=[5040;0;52920;0;66150;0;11025;0];
            case 9
                coeffs=[40320;0;564480;0;1058400;0;352800;0;11025];
            case 10
                coeffs=[362880;0;6531840;0;17146080;0;9525600;0;893025;0];
            case 11
                coeffs=[3628800;0;81648000;0;285768000;0;238140000;0;44651250;0;893025];
            otherwise
                error('Invalid n')
        end
        k=1+2*(n-1);
    else    
        %The derivative of (1-x^2)^(-(k/2)) is
        %kx(1-x^2)^(-1-k/2)
        %That is k*x*(1-x^2)^(-(k1/2)) where k1=k+2.
        xVal=[1,0];

        %This polynomial is 1-x^2.
        polyX=[-1,0,1];

        %The first derivative 
        n=n-1;
        coeffs=1;
        k=1;
        while(n>0)
            %Take the derivative of the denominator.
            coeffsNum=conv(k*xVal,coeffs);
            k=k+2;

            %Now, differentiate the numerator and multiply to get a
            %common denominator.
            coeffs=polyder(coeffs);
            coeffs=conv(coeffs,polyX);
            coeffs=polySum(coeffs,coeffsNum);
            n=n-1;
        end
    end
else
    k=1+2*(n-1);
end
    
if(~isempty(z))
    vals=zeros(size(z));
    %Evaluate the result.
    vals(:)=polyval(coeffs,z(:)).*sqrt(1-z(:).^2).^(-k);
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
