function [HBar,dHBardu,d2HBardu2]=normHelmholtz(u,M,scalFactor)
%%NORMHELMHOLTZ This computes the fully normalized derived Legendre
%               functions (also known as fully normalized Helmholtz 
%               polynomials) of degree n=0 to M and for each n from order
%               m=0 to n evaluated at the point u. That is, \bar{H}^m_n(u).
%               This also finds the first and second derivatives of the
%               normalized Helmholtz polynomials with respect to the
%               parameter u, D{\bar{H}^m_n(u)} and D2{\bar{H}^m_n(u)},
%               where D{} and D2{} are respectively the first and second
%               derivative operators. All of the values can be scaled by a
%               factor of scalFac, if desired, to help prevent overflows
%               with high degrees and orders. The Helmholtz polynomials are
%               used in Pine's method for spherical harmonic evalulation.
%
%INPUTS:    u   The value at which the fully normalized Helmholtz
%               polynomials should be evaluated.
%       maxDeg  The maximum degree and order of the output. This should be
%               >=3.
%    scalFactor A scale factor to help prevent overflow of the results. A
%               value of 10^(-280) (for example) might be useful at high
%               degrees.
%
%OUTPUTS: HBar   An instance of the ClusterSet class such that
%                HBar(n+1,m+1)=scalFac*\bar{H}^m_n(u).
%       dHBardu  An instance of the ClusterSet class such that
%                dHBardu(n+1,m+1)=scalFac*D{\bar{H}^m_n(u)}.
%     d2HBardu2  An instance of the ClusterSet class such that
%                dHBardu(n+1,m+1)=scalFac*D{\bar{H}^m_n(u)}.
%
%The fully normalized derived Legendre functions (Helmholtz polynomials)
%are described in [1] and the algorithm of that paper is implemented here
%with the addition of a scale factor, if desired.
%
%The Helmholtz polynomial of degree n and order m is defined to be
%H^m_n(u)=1/(n!2^n)D_{n+m}{(u^2-1)^n}
%where the notation D_y{x} represents the yth derivative of x. A fully
%normalized Helmholtz polynomial is defined as
%\bar{H}^m_n(u)=H^m_n(u)*sqrt(k*(2*n+1)*factorial(n-m)/factorial(n+m))
%where k=1 if m=0 and k=2 otherwise. This normalization is consistent with
%the normalization that the National Geospatial Intelligence Agency (NGA)
%uses with the coefficients for the EGM96 and EGM2008 gravitational models.
%
%REFERENCES:
%[1] E. Fantino and S. Casotto, "Methods of harmonic synthesis for global
%    geopotential models and their first-, second- and third-order
%    gradients," Journal of Geodesy, vol. 83, no. 7, pp. 595-619, Jul.
%    2009.
%
%January 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    %Allocate space for the return variables. n ranges from 0 to M and m
    %ranges from 0 to n. Thus, the number of elements for each n ranges
    %from 1 to M+1. This means that a total of
    %(M+1)*(M+2)/2 elements are necessary.
    numPBarU=(M+1)*(M+2)/2;
    totalP=zeros(numPBarU,1);
    clustSizes=1:(M+1);
    HBar=ClusterSet(totalP,clustSizes);
    
    %Set the first few terms explicitly.
    HBar(0+1,0+1)=1*scalFactor;
    HBar(1+1,1+1)=sqrt(3)*scalFactor;
    HBar(1+1,0+1)=sqrt(3)*u*scalFactor;

    %Compute all terms of the form (n,n) and (n,n-1).
    for n=2:M
        %Get the (n,n) term using Equation 55.
        f=sqrt((2*n+1)/(2*n));
        HBar(n+1,n+1)=f*HBar(n-1+1,n-1+1);
        %Get the (n,n-1) term using Equation 56.
        HBar(n+1,n-1+1)=u*sqrt(2*n)*HBar(n+1,n+1);
    end
    
    %Now, compute all of the other terms using Equation 58.
    for m=0:M
        for n=(m+2):M
            g=sqrt((2*n+1)*(2*n-1)/((n+m)*(n-m)));
            h=sqrt((2*n+1)*(n-m-1)*(n+m-1)/((2*n-3)*(n+m)*(n-m)));
            
            HBar(n+1,m+1)=u*g*HBar(n-1+1,m+1)-h*HBar(n-2+1,m+1);
        end
    end
    
    %If the first derivatives are desired. The first derivatives are easily
    %computed using Equation 61.
    if(nargout>1)
        dHBardu=ClusterSet(totalP,clustSizes);
        
        %For m=0
        dHBardu(0+1,0+1)=0;
        for n=1:M
            k=sqrt(n*(n+1)/2);
            dHBardu(n+1,0+1)=k*HBar(n+1,1+1);
        end
        
        %For other n and m pairs
        for n=0:M
            for m=1:(n-1)
                k=sqrt((n-m)*(n+m+1));
                
                dHBardu(n+1,m+1)=k*HBar(n+1,m+1+1);
            end
            dHBardu(n+1,n+1)=0;
        end
    end
    
    %If the second derivatives are desired. The second derivatives are
    %easily computed using Equation 62.
    if(nargout>2)
        d2HBardu2=ClusterSet(totalP,clustSizes);
        
        %For m=0 and m=1
        d2HBardu2(0+1,0+1)=0;
        d2HBardu2(1+1,0+1)=0;
        d2HBardu2(1+1,1+1)=0;
        
        %For n=2, m=0.
        m=0;
        n=2;
        k=sqrt((n-m)*(n+m+1));
        kp=sqrt((n-(m+1))*(n+(m+1)+1));
        d2HBardu2(n+1,m+1)=k*kp*HBar(n+1,m+2+1);
        
        %For n=2, m=1;
        m=1;
        d2HBardu2(n+1,m+1)=0;
        
        %For m=0 and m=1 and all other n.
        for n=3:M            
            m=0;
            k=sqrt(n*(n+1)/2);
            kp=sqrt((n-(m+1))*(n+(m+1)+1));
            d2HBardu2(n+1,m+1)=k*kp*HBar(n+1,m+2+1);
            
            m=1;
            k=sqrt((n-m)*(n+m+1));
            kp=sqrt((n-(m+1))*(n+(m+1)+1));
            d2HBardu2(n+1,m+1)=k*kp*HBar(n+1,m+2+1);
        end
        
        %For other n and m pairs.
        for n=3:M
            for m=2:(n-2)
                k=sqrt((n-m)*(n+m+1));
                kp=sqrt((n-(m+1))*(n+(m+1)+1));
                
                d2HBardu2(n+1,m+1)=k*kp*HBar(n+1,m+2+1);
            end
            
            d2HBardu2(n-1+1,n-1+1)=0;
            d2HBardu2(n+1,n-1+1)=0;
            d2HBardu2(n+1,n+1)=0;
        end
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
