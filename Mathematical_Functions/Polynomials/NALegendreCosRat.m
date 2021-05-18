function [PBarUVals,dPBarUValsdTheta,d2PBarUValsdTheta2]=NALegendreCosRat(theta,M,scalFactor)
%NALEGENDRECOSRAT Evaluate \bar{P}_{nm}(cos(theta))/u^m for all n from
%                 0 to M and for each n for all m from 0 to n, where
%                 u=abs(sin(theta)). \bar{P}_{nm}(x) is the fully
%                 normalized associated Legendre function of x of degree n
%                 and order m. Also evaluate
%                 D{\bar{P}_{nm}(cos(theta))}/u^m and
%                 D2{\bar{P}_{nm}(cos(theta))}/u^m where D{} is the first
%                 derivative operator with respect to theta and D2{} is the
%                 second derivative operator with respect to theta. All of
%                 the values can be scaled by a factor of scalFac, if
%                 desired, to help prevent overflows with high degrees and
%                 orders.
%
%INPUTS: theta An angle in radians.
%       maxDeg The maximum degree and order of the output. This must be
%              >=0.
%   scalFactor A scale factor to help prevent overflow of the results. In
%              [1], discussed below, a value of 10^(-280) is used.
%
%OUTPUTS: PBarUVals An instance of the CountingClusterSet class such that
%                   PBarUVals(n+1,m+1)=scalFac*\bar{P}_{nm}(cos(theta))/u^m.
%  dPBarUValsdTheta An instance of the CountingClusterSet class such that
%                   dPBarUValsdTheta(n+1,m+1)=scalFac*D{\bar{P}_{nm}(cos(theta))}/u^m
% d2PBarUValsdTheta2 An instance of the CountingClusterSet class such that
%                   d2PBarUValsdTheta2(n+1,m+1)=scalFac*D2{\bar{P}_{nm}(cos(theta))}/u^m
%
%The modified forward row (MFR) algorithm of [1] is used to compute
%PBarUVals and dPBarUValsdTheta. For d2PBarUValsdTheta2, the algorithm of
%[2] is used. However, the paper contains a typo. The first term of the
%first unnumbered equation should be multiplied by an additional (1/u).
%This function uses the correct formulation, which is also described in
%Appendix F of [3].
%
%With the notation P^m_n(x) for an associated Legendre function of degree n
%and order m evaluated at the point x, one generally means 
%P^m_n(x)=(-1)^m*(1-x^)^(m/2)*D_m{P_n(x)}
%where P_n(x) is a Legendre polynomial of degree n evaluated at x and
%D_m{P_n(x)} is the mth derivative of the polynomial. These associated
%Legendre functions are the same as those used in the International
%Geomagnetic Reference Field. On the other hand, with the notation
%P_{nm}(x), one generally means (-1)^mP^m_n(x). This notation is also
%called an associated Legendre function. Matlab's built-in function
%legendre(n,x) returns all P^m_n(x) for m=0 to m=n. On the other hand, the
%notation \bar{P}_{nm}(x), refers to a fully normalized associated Legendre
%function. Multiple definitions of what is normalized exist. In this
%function, the normalization is such that
%\bar{P}_{nm}(x)=P_{nm}(x)*sqrt(k*(2*n+1)*factorial(n-m)/factorial(n+m))
%where k=1 if m=0 and k=2 otherwise. This differs from the normalization
%constant that Matlab uses in the normalized version of its built-in
%associated Legendre function. It is the same as the normalization 
%constant the NGA uses in its coefficients for the EGM96 and the
%EGM2008 gravitation models.
%
%The fully normalized associated Legendre function ratios that this
%function computes can be used in the synthesis of spherical harmonic
%coefficients.
%
%REFERENCES:
%[1] S. A. Holmes and W. E. Featherstone, "A unified approach to the
%    Clenshaw summation and the recursive computation of very high degree
%    and order normalised associated Legendre functions," Journal of
%    Geodesy, vol. 76, no. 5, pp. 279-299, May 2002.
%[2] S. A. Holmes and W. E. Featherstone, "Short note: Extending simplified
%    high-degree synthesis methods to second latitude derivatives of
%    geopotential," Journal of Geodesy, vol. 76, no. 8, pp. 447-450, Nov.
%    2002.
%[3] D. F. Crouse, "An overview of major terrestrial, celestial, and
%    temporal coordinate systems for target tracking," Naval Research
%    Laboratory, Washington, DC, Tech. Rep. NRL/FR/5344-16-10,279, 10 Aug.
%    2016.
%
%December 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    if(M<0)
       error('M must be >=0.') 
    end

    u=sin(theta);
    t=cos(theta);

    %Allocate space for the return variables. n ranges from 0 to M and m
    %ranges from 0 to n. Thus, the number of elements for each n ranges
    %from 1 to M+1. This means that a total of
    %(M+1)*(M+2)/2 elements are necessary.
    numPBarU=(M+1)*(M+2)/2;
    totalP=zeros(numPBarU,1);
    PBarUVals=CountingClusterSet(totalP);
    
    %The value of PBar_{0,0}(cos(theta)) is independent of theta and is
    %one. 
    PBarUVals(0+1,0+1)=1*scalFactor;

    if(M>0)
        %Set the seed value for PBar_{1,1}(cos(theta))/u from which the
        %other values will be computed.
        PBarUVals(1+1,1+1)=sqrt(3)*scalFactor;
        
        jTerm=1/sqrt(2);
        %First, deal with the case where n=1, m=0;
        n=1;
        m=0;
        %g is given in Equation 18 of the first Holmes and Featherstone paper.
        g=2*(m+1)/sqrt((n-m)*(n+m+1));
        PBarUVals(n+1,m+1)=jTerm*g*t*PBarUVals(n+1,m+1+1);

        %Compute the values along the main diagonal, where m=n
        %starting from m=n=2. This implements equation 28 in the first
        %Holmes and Featherstone paper for the normalized associated
        %Legendre function ratio.
        for m=2:M
            PBarUVals(m+1,m+1)=sqrt((2*m+1)/(2*m))*PBarUVals(m-1+1,m-1+1);
        end

        %Recursively compute the values using Equation 27 from the
        %first Holmes and Featherstone paper, taking into account the
        %fact that the first element of the recursion only has one term.

        %Next, evaluate the values for all other valid n and m.
        for n=2:M
            %Deal with the first element of the recursion,which is
            %where m=n-1.
            m=n-1;
            g=2*(m+1)/sqrt((n-m)*(n+m+1));
            PBarUVals(n+1,m+1)=g*t*PBarUVals(n+1,m+1+1);

            %Recursively compute the values of the rest of the m terms.
            for m=(n-2):-1:1
                g=2*(m+1)/sqrt((n-m)*(n+m+1));
                %h is given in Equation 19 of the first Holmes and
                %Featherstone paper.
                h=sqrt((n+m+2)*(n-m-1)/((n-m)*(n+m+1)));
                PBarUVals(n+1,m+1)=g*t*PBarUVals(n+1,m+1+1)-h*u^2*PBarUVals(n+1,m+2+1);
            end

            %Deal with the special m=0 case.
            m=0;
            g=2*(m+1)/sqrt((n-m)*(n+m+1));
            h=sqrt((n+m+2)*(n-m-1)/((n-m)*(n+m+1)));
            PBarUVals(n+1,m+1)=jTerm*(g*t*PBarUVals(n+1,m+1+1)-h*u^2*PBarUVals(n+1,m+2+1));
        end
    end
    
    %If the first derivative is desired.
    if(nargout>1)
        %Allocate space.
        dPBarUValsdTheta=CountingClusterSet(totalP);
        %The first derivative of PBar_{0,0}(cos(theta)) is just zero.
        dPBarUValsdTheta(1,1)=0;
        
        if(M>0)
            m=1;
            n=1;
            %From Equation 30 in the first Holmes and Featherstone paper. This
            %is the seed value from which other values will be computed.
            dPBarUValsdTheta(1+1,1+1)=m*(t/u)*PBarUVals(1+1,1+1);

            n=1;
            m=0;
            %e is given in Equation 22 of the first Holmes and Featherstone
            %paper.
            e=sqrt((n+m+1)*(n-m)/2);
            %This is Equation 30 of the first Holmes and Featherstone paper for
            %m=0.
            dPBarUValsdTheta(n+1,m+1)=-e*u*PBarUVals(n+1,m+1+1);

            %Compute the values along the main diagonal, where m=n starting 
            %from m=n=2. This implements Equation 30 in the frist Holmes and
            %Featherstone paper for the ratio of the first derivative. 
            for m=2:M
                dPBarUValsdTheta(m+1,m+1)=m*(t/u)*PBarUVals(m+1,m+1);
            end

            %Next, evaluate the values for all other valid n and m.
            for n=2:M
                %Recursively compute the values of the m terms for m>0.
                for m=(n-1):-1:1
                    e=sqrt((n+m+1)*(n-m));
                    dPBarUValsdTheta(n+1,m+1)=m*(t/u)*PBarUVals(n+1,m+1)-e*u*PBarUVals(n+1,m+1+1);
                end
                %Deal with the special m=0 case.
                m=0;
                e=sqrt((n+m+1)*(n-m)/2);
                dPBarUValsdTheta(n+1,m+1)=-e*u*PBarUVals(n+1,m+1+1);
            end
        end
    end
    
    %If the second derivative is desired
    if(nargout>2)
        %Allocate space
        d2PBarUValsdTheta2=CountingClusterSet(totalP);
        
        for n=0:M
            for m=0:n
                %From the first (un-numbered) equation in the second Holmes
                %and Featherstone paper AFTER correction.
                d2PBarUValsdTheta2(n+1,m+1)=(m^2/u^2-n*(n+1))*PBarUVals(n+1,m+1)-(t/u)*dPBarUValsdTheta(n+1,m+1);
            end
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
