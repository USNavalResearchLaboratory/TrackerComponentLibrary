function [HBar,dHBardu,d2HBardu2,d3HBardu3]=normHelmholtz(u,M,scalFactor)
%%NORMHELMHOLTZ This computes the fully normalized derived Legendre
%               functions (also known as fully normalized Helmholtz 
%               polynomials) of degree n=0 to M and for each n from order
%               m=0 to n evaluated at the point u. That is, \bar{H}^m_n(u).
%               This also finds the first, second, and third derivatives of
%               the normalized Helmholtz polynomials with respect to the
%               parameter u, D{\bar{H}^m_n(u)}, D2{\bar{H}^m_n(u)}, and
%               D3{\bar{H}^m_n(u)} where D{}, D2{}, and D3{} are
%               respectively the first, second, and third derivative
%               operators. All of the values can be scaled by a factor of
%               scalFac, if desired, to help prevent overflows with high
%               degrees and orders. The Helmholtz polynomials are used in
%               Pine's method for spherical harmonic evalulation.
%
%INPUTS: u The value at which the fully normalized Helmholtz polynomials
%          should be evaluated.
%   maxDeg The maximum degree and order of the output. This must be >=0.
% scalFactor A scale factor to help prevent overflow of the results. A
%          value of 10^(-280) (for example) might be useful at high
%          degrees. If this parameter is omitted or an empty matrix is
%          passed, then the default of 1 is used.
%
%OUTPUTS: HBar An instance of the CountingClusterSet class such that
%              HBar(n+1,m+1)=scalFac*\bar{H}^m_n(u).
%      dHBardu An instance of the CountingClusterSet class such that
%              dHBardu(n+1,m+1)=scalFac*D{\bar{H}^m_n(u)}.
%    d2HBardu2 An instance of the CountingClusterSet class such that
%              dHBardu(n+1,m+1)=scalFac*D2{\bar{H}^m_n(u)}.
%    d2HBardu3 An instance of the CountingClusterSet class such that
%              dHBardu(n+1,m+1)=scalFac*D3{\bar{H}^m_n(u)}.
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
%EXAMPLE:
%Here, we demonstrate that the deirvatives are consistent with the
%derivatives computed numerically via forward differencing.
% u=0.7047;
% M=12;
% [HBar,dHBardu,d2HBardu2,d3HBardu3]=normHelmholtz(u,M,1);
% epsVal=1e-10;
% [HBar1,dHBardu1,d2HBardu21]=normHelmholtz(u+epsVal,M,1);
% numDiff1=(HBar1(:)-HBar(:))/epsVal;
% diffErr1=max(abs((numDiff1-dHBardu(:))./dHBardu(:)))
% 
% numDiff2=(dHBardu1(:)-dHBardu(:))/epsVal;
% diffErr2=max(abs((numDiff2-d2HBardu2(:))./d2HBardu2(:)))
% 
% numDiff3=(d2HBardu21(:)-d2HBardu2(:))/epsVal;
% diffErr3=max(abs((numDiff3-d3HBardu3(:))./d3HBardu3(:)))
%One will see that the relative error between the analytic derivatives and
%the numerical derivatives is on the order of 2e-5 to 7e-6, which indicates
%a good agreement.
%
%REFERENCES:
%[1] E. Fantino and S. Casotto, "Methods of harmonic synthesis for global
%    geopotential models and their first-, second- and third-order
%    gradients," Journal of Geodesy, vol. 83, no. 7, pp. 595-619, Jul.
%    2009.
%
%January 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release. 

    if(nargin<3||isempty(scalFactor))
        scalFactor=1;
    end

    %Allocate space for the return variables. n ranges from 0 to M and m
    %ranges from 0 to n. Thus, the number of elements for each n ranges
    %from 1 to M+1. This means that a total of
    %(M+1)*(M+2)/2 elements are necessary.
    numPBarU=(M+1)*(M+2)/2;
    totalP=zeros(numPBarU,1);
    HBar=CountingClusterSet(totalP);
    
    %Set the first few terms explicitly.
    HBar(0+1,0+1)=1*scalFactor;
    if(M>0)
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
    end
    %If the first derivatives are desired. The first derivatives are easily
    %computed using Equation 61.
    if(nargout>1)
        dHBardu=CountingClusterSet(totalP);
        dHBardu(:)=NaN;
        %For m=0
        dHBardu(0+1,0+1)=0;
        if(M>0)
            for n=1:M
                k=sqrt(n*(n+1)/2);
                dHBardu(n+1,0+1)=k*HBar(n+1,1+1);
            end

            %For n=1 and m=1;
            n=1;
            m=1;
            dHBardu(n+1,m+1)=0;

            %For other n and m pairs
            for n=2:M
                for m=1:(n-1)
                    k=sqrt((n-m)*(n+m+1));

                    dHBardu(n+1,m+1)=k*HBar(n+1,m+1+1);
                end
                dHBardu(n+1,n+1)=0;
            end
        end
    end
    
    %If the second derivatives are desired. The second derivatives are
    %easily computed using Equation 62.
    if(nargout>2)
        d2HBardu2=CountingClusterSet(totalP);

        %For n=0 and n=1
        d2HBardu2(0+1,0+1)=0;
        if(M>0)
            d2HBardu2(1+1,0+1)=0;
            d2HBardu2(1+1,1+1)=0;

            if(M>1)
                %For n=2, m=0.
                n=2;
                m=0;
                k=sqrt((n-m)*(n+m+1)/2);
                kp=sqrt((n-(m+1))*(n+(m+1)+1));
                d2HBardu2(n+1,m+1)=k*kp*HBar(n+1,m+2+1);

                %For n=2, m=1,2
                d2HBardu2(2+1,1+1)=0;
                d2HBardu2(2+1,2+1)=0;

                if(M>2)
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
                    %For n=3, m=2,m=3
                    d2HBardu2(3+1,2+1)=0;
                    d2HBardu2(3+1,3+1)=0;

                    %For other n and m pairs.
                    for n=4:M
                        for m=2:(n-2)
                            k=sqrt((n-m)*(n+m+1));
                            kp=sqrt((n-(m+1))*(n+(m+1)+1));

                            d2HBardu2(n+1,m+1)=k*kp*HBar(n+1,m+2+1);
                        end

                        d2HBardu2(n+1,n-1+1)=0;
                        d2HBardu2(n+1,n+1)=0;
                    end
                end
            end
        end
    end
    
    %If the third derivatives are desired. The third derivatives are
    %easily computed using Equation 63.
    if(nargout>3)
        d3HBardu3=CountingClusterSet(totalP);
        
        %For n=0
        d3HBardu3(0+1,0+1)=0;
        if(M>0)
            %For n=1
            d3HBardu3(1+1,0+1)=0;
            d3HBardu3(1+1,1+1)=0;
            if(M>1)
                %For n=2
                d3HBardu3(2+1,0+1)=0;
                d3HBardu3(2+1,1+1)=0;
                d3HBardu3(2+1,2+1)=0;

                if(M>2)
                    %For n=3, m=0
                    n=3;
                    m=0;
                    k=sqrt(n*(n+1)/2);
                    kp1=sqrt((n-(m+1))*(n+(m+1)+1));
                    kp2=sqrt((n-(m+2))*(n+(m+2)+1));
                    d3HBardu3(n+1,m+1)=k*kp1*kp2*HBar(n+1,m+3+1);

                    %For n=3, m=1 and m=2, m=3
                    d3HBardu3(3+1,1+1)=0;
                    d3HBardu3(3+1,2+1)=0;
                    d3HBardu3(3+1,3+1)=0;

                    if(M>3)
                        %For n=4, m=0
                        n=4;
                        m=0;
                        k=sqrt(n*(n+1)/2);
                        kp1=sqrt((n-(m+1))*(n+(m+1)+1));
                        kp2=sqrt((n-(m+2))*(n+(m+2)+1));
                        d3HBardu3(n+1,m+1)=k*kp1*kp2*HBar(n+1,m+3+1);

                        %For n=4, m=1
                        n=4;
                        m=1;
                        k=sqrt((n-m)*(n+m+1));
                        kp1=sqrt((n-(m+1))*(n+(m+1)+1));
                        kp2=sqrt((n-(m+2))*(n+(m+2)+1));
                        d3HBardu3(n+1,m+1)=k*kp1*kp2*HBar(n+1,m+3+1);

                        %For n=4, m=2, 3, 4
                        d3HBardu3(4+1,2+1)=0;
                        d3HBardu3(4+1,3+1)=0;
                        d3HBardu3(4+1,4+1)=0;
                    
                        if(M>4)
                            %For m=0, m=1, and m=2 and all other n.
                            for n=5:M            
                                m=0;
                                k=sqrt(n*(n+1)/2);
                                kp1=sqrt((n-(m+1))*(n+(m+1)+1));
                                kp2=sqrt((n-(m+2))*(n+(m+2)+1));
                                d3HBardu3(n+1,m+1)=k*kp1*kp2*HBar(n+1,m+3+1);

                                for m=1:2
                                    k=sqrt((n-m)*(n+m+1));
                                    kp1=sqrt((n-(m+1))*(n+(m+1)+1));
                                    kp2=sqrt((n-(m+2))*(n+(m+2)+1));
                                    d3HBardu3(n+1,m+1)=k*kp1*kp2*HBar(n+1,m+3+1);
                                end
                            end

                            %For n=5, m=3,4,5
                            d3HBardu3(5+1,3+1)=0;
                            d3HBardu3(5+1,4+1)=0;
                            d3HBardu3(5+1,5+1)=0;

                            for n=6:M
                                %The m=0 case
                                m=0;
                                k=sqrt(n*(n+1)/2);
                                kp1=sqrt((n-(m+1))*(n+(m+1)+1));
                                kp2=sqrt((n-(m+2))*(n+(m+2)+1));
                                d3HBardu3(n+1,m+1)=k*kp1*kp2*HBar(n+1,m+3+1);

                                for m=3:(n-3)
                                    k=sqrt((n-m)*(n+m+1));
                                    kp1=sqrt((n-(m+1))*(n+(m+1)+1));
                                    kp2=sqrt((n-(m+2))*(n+(m+2)+1));
                                    d3HBardu3(n+1,m+1)=k*kp1*kp2*HBar(n+1,m+3+1);
                                end

                                d3HBardu3(n+1,n-2+1)=0;
                                d3HBardu3(n+1,n-1+1)=0;
                                d3HBardu3(n+1,n+1)=0;
                            end
                        end
                    end
                end
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
