function C=changeSpherHarmonicNorm(C,origNormType,destNormType,arePolyVals)
%%CHANGESPHERHARMONNORM Spherical harmonic expansions are expressed in the
%               literature using a number of different normalizations
%               of the associated Legendre polynomials. This function lets
%               one convert the (real or complex) coefficients normalized
%               for one type of normalization to that of a different type
%               of normalization. Alternatively, it can also be used to
%               switch the normalization of a set of associated Legendre
%               polynomial values.
%
%INPUTS: C This is a length (M+2)*(M+1)/2 array holding real or complex
%          spherical harmonic coefficients or real Legendre function values
%          from degrees 0 to M. This can also be a CountingClusterSet
%          object containing the same information. As a CountingClusterSet
%          object, the value in C(n+1,m+1) is the value of degree n and
%          order m. 
% origNormType,destNormType These specify the original normalization type
%          and the desired normalization type. The normalization is
%          expressed as N(n,m) time P_{nm}(x), where P_{nm}(x) represents
%          the non-normalized associated Legendre functions without the
%          Condon-Shortley phase, as discussed below. Possible values are:
%          0 (The default if omitted or an empty matrix is passed).
%            Select fully normalized values using the definition of
%            full normalization that is used in the EGM2008 model as well
%            as other gravitational models. 
%            N(n,m)=sqrt((2*n+1)*(2-KDelta(m))*factorial(n-m)/factorial(n+m))
%            where KDelta is the Kronecker delta function.
%          1 Select fully normalized values using the definition of full
%            normalization that is used in Matlab's legendre function. Here
%            N(n,m)=sqrt((n+1/2)*factorial(n-m)/factorial(n+m))
%          2 Select fully normalized values using a definition of full
%            normalization that arises in geodesy.
%            N(n,m)=sqrt((2*n+1)*factorial(n-m)/factorial(n+m))
%          3 Select Schmidt seminormalized values using a definition of
%            Schmidt seminormalization that includes a Kronecker delta.
%            Here
%            N(n,m)=sqrt((2-KDelta(m))*factorial(n-m)/factorial(n+m))
%          4 Select Schmidt seminormalized values using a definition of
%            Schmidt seminormalization that omits a Kronecker delta.
%            N(n,m)=sqrt(factorial(n-m)/factorial(n+m))
%          5 Select values with a type of orthonormalized scaling that does
%            not use a Kronecker delta term.
%            N(n,m)=sqrt(((2*n+1)/(4*pi))*factorial(n-m)/factorial(n+m))
%          6 Select values with a type of orthonormalized scaling that used
%            a Kronecker delta term.
%            N(n,m)=sqrt(((2*n+1)/(4*pi))*(2-KDelta(m))*factorial(n-m)/factorial(n+m))
%          7 Select values with no normalization.
%            N(n,m)=1
%          8 Select values with no normalization for a definition of
%            associated Legendre polynomials that includes the Condon-
%            Shortley phase
%            N(n,m)=(-1)^m
% arePolyVals A boolean value indicating whether the values in C are
%           associated Legendre polynomial values (true) or whether they
%           are values of coefficients that go into a spherical harmonic
%           series expansion. The default if omitted or an empty matrix is
%           passed is false.
%
%OUTPUTS: C This is the same as C except the normalization has been
%           changed.
%
%The definition of associated Legendre polynomial underlying this function
%is
% P_(nm)(x)=(1-x^2)^(m/2)*Dm{P(n,x)}
%where Dm{} represents the mth-order derivative and P(n,x) is the Legendre
%polynomial having the form
%P(n,x)=(1/2^n)*sum_{k=0}^n binomial(n,n)*(x-1)^(n-k)*(x+k)^k
%This definition lacks the Condon-Shortley phase. Adding the phase is the
%same as multiplying P_(nm)(x) by (-1)^m.
%
%This function can be useful for transforming spherical Harmonic
%coefficients into a format that can be used by the spherHarmonicEval
%function.
%
%June 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(origNormType))
   origNormType=0; 
end

if(nargin<3||isempty(destNormType))
   destNormType=0; 
end

if(nargin<4||isempty(arePolyVals))
    arePolyVals=false; 
end

if(arePolyVals==false)
   temp=origNormType;
   origNormType=destNormType;
   destNormType=temp;
end

if(isfloat(C))
    M=(1/2)*(sqrt(1+8*length(C))-1)-1;
    C=CountingClusterSet(C);
    convertToFloat=true;
else%Assume it is a ClusterSet class.
    M=C.numClust-1;
    convertToFloat=false;
end

%Convert from particular associated Legendre polynomials values to
%NGA fully normalized associated Legendre values.
switch(origNormType)
    case 0%NGA fully normalized values.
        %There is no need to do anything. This is how they are given
        %above.
    case 1%Matlab fully normalized values.
        for n=0:M
            C(n+1,0+1)=C(n+1,0+1)*sqrt(2);
            C(n+1,(1:n)+1)=C(n+1,(1:n)+1)*2;
        end
    case 2%Geodesy fully normalized values.
        for n=0:M
            C(n+1,(1:n)+1)=C(n+1,(1:n)+1)*sqrt(2);
        end
    case 3%Schmidt seminormalized values type 1
        for n=0:M
            C(n+1,:)=C(n+1,:)*sqrt(1+2*n);
        end
    case 4%Schmidt seminormalized values type 2
        for n=0:M
            C(n+1,0+1)=C(n+1,0+1)*(sqrt(1+2*n));
            for m=1:n
                C(n+1,m+1)=C(n+1,m+1)*(2*sqrt(1+2*n));
            end
        end
    case 5%Orthonormalized type 1
        for n=0:M
            C(n+1,0+1)=C(n+1,0+1)*(2*sqrt(pi));
            for m=1:n
                C(n+1,m+1)=C(n+1,m+1)*(2*sqrt(2*pi));
            end
        end
    case 6%Orthonormalized type 2
        C(:)=C(:)*(2*sqrt(pi));
    case 7%Non-normalized values without the Condon-Shortley phase
        for n=0:M
            coeff=1;
            nVal=1/sqrt(1+2*n);
            
            C(n+1,0+1)=C(n+1,0+1)/(nVal*coeff);
            nVal=nVal/sqrt(2);
            for m=1:n
                coeff=coeff*sqrt(m*(1-m)+n*(1+n));
                
                C(n+1,m+1)=C(n+1,m+1)/(nVal*coeff);
            end
        end
    case 8%Non-normalized values with the Condon-Shortley phase
        for n=0:M
            coeff=1;
            signVal=1;
            nVal=1/sqrt(1+2*n);
            
            C(n+1,0+1)=signVal*C(n+1,0+1)/(nVal*coeff);
            nVal=nVal/sqrt(2);
            for m=1:n
                coeff=coeff*sqrt(m*(1-m)+n*(1+n));
                signVal=-signVal;
                
                C(n+1,m+1)=signVal*C(n+1,m+1)/(nVal*coeff);
            end
        end
    otherwise
        error('Unknown source normalization specified.')
end

%Convert from NGA fully normalized associated Legendre values to the new values.
switch(destNormType)
    case 0%NGA fully normalized values.
        %There is no need to do anything. This is how they are given
        %above.
    case 1%Matlab fully normalized values.
        for n=0:M
            C(n+1,0+1)=C(n+1,0+1)/sqrt(2);
            C(n+1,(1:n)+1)=C(n+1,(1:n)+1)/2;
        end
    case 2%Geodesy fully normalized values.
        for n=0:M
            C(n+1,(1:n)+1)=C(n+1,(1:n)+1)/sqrt(2);
        end
    case 3%Schmidt seminormalized values type 1
        for n=0:M
            C(n+1,:)=C(n+1,:)/sqrt(1+2*n);
        end
    case 4%Schmidt seminormalized values type 2
        for n=0:M
            C(n+1,0+1)=C(n+1,0+1)/(sqrt(1+2*n));
            for m=1:n
                C(n+1,m+1)=C(n+1,m+1)/(2*sqrt(1+2*n));
            end
        end
    case 5%Orthonormalized type 1
        for n=0:M
            C(n+1,0+1)=C(n+1,0+1)/(2*sqrt(pi));
            for m=1:n
                C(n+1,m+1)=C(n+1,m+1)/(2*sqrt(2*pi));
            end
        end
    case 6%Orthonormalized type 2
        C(:)=C(:)/(2*sqrt(pi));
    case 7%Non-normalized values without the Condon-Shortley phase
        for n=0:M
            coeff=1;
            nVal=1/sqrt(1+2*n);
            
            C(n+1,0+1)=nVal*coeff*C(n+1,0+1);
            nVal=nVal/sqrt(2);
            for m=1:n
                coeff=coeff*sqrt(m*(1-m)+n*(1+n));
                
                C(n+1,m+1)=nVal*coeff*C(n+1,m+1);
            end
        end
    case 8%Non-normalized values with the Condon-Shortley phase
        for n=0:M
            coeff=1;
            signVal=1;
            nVal=1/sqrt(1+2*n);
            
            C(n+1,0+1)=signVal*nVal*coeff*C(n+1,0+1);
            nVal=nVal/sqrt(2);
            for m=1:n
                coeff=coeff*sqrt(m*(1-m)+n*(1+n));
                signVal=-signVal;
                
                C(n+1,m+1)=signVal*nVal*coeff*C(n+1,m+1);
            end
        end
    otherwise
        error('Unknown source normalization specified.')
end

if(convertToFloat)
    C=C.clusterEls;
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
