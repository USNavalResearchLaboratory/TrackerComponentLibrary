function w=LambW(z,k)
%%LAMBDW Evaluate the primary branch and the -1 branch of the Lambert-W
%        function, which is also know as a product logarithm. The function
%        solve z=w*exp(w) for w. This implementation is meant for real
%        values of z>=-1/exp(1) and will have an error on smaller values,
%        which would normally return complex solutions.
%
%INPUTS: z A vector or matrix of parameters at which the primary branch of
%          the Lambert W function is to be evaluated. It is assumed that
%          z>=-exp(-1) for all of the elements.
%        k An optional parameter specifying the branch of the algorithm to
%          take. k=0 means the primary branch and -1 means the alternative
%          real branch. This only affects solutions for z<=0.
%
%OUTPUTS: w The solutions to z=w*exp(w) for all elements in z.
%
%Note that a loss of precision will occur for negative values on the -1
%branch around approximately abs(z)<1e-307.
%
%The origin of the algorithm is [1], where the main algorithm is Halley's
%method from Equation 5.9. For negative values on the primary branch, an
%initial value of the first few terms of Equation 4.22 is used. For
%positive values less than 1.1, an initial value of 0 is used. For positive
%values greater than or equal to 1.1, the first two perms of Equation 4.18
%are used. This same initialization (taking advantage of absolute values)
%is used for negative values on the k-1 branch. This differs from the
%suggested method in the paper in that a Padé approximation is not used for
%terms near zero. Rather, a constant initial estimate of 0 is used. 15
%iterations are performed, which should be more than enough to converge.
%
%REFERENCES:
%[1] R. M. Corless, G. H. Gonnet, D. E. G. Hare, and D. J. Jefrey, "On the
%    Lambert W Function," Advances in Computational Mathematics" vol. 5,
%    no. 1, pp. 329-359, 1996.
%
%February 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    if(nargin<2||isempty(k))
       k=0; 
    end
    
    if(k~=0&&k~=-1)
        error('The specified branch k is invalid.') 
    end

    if(any(z<-exp(-1)))
        warning('Inputs less than -exp(-1) are clipped to -exp(-1). Complex branch not implemented.')
        sel=z<-exp(-1);
        z(sel)=-exp(-1);
    end

    numIter=15;
    w=zeros(size(z));
    
    %Initialization for small positive values of z.
    zSel=z<1&z>=0;
    w(zSel)=0;
    
    %Initialize for negative values of z if the primary branch is desired.
    if(k==0)
        zSel=z<0;
        
        %Deal with finite precision errors
        p=sqrt(abs(2*(exp(1)*z(zSel)+1)));
        %Deal with finite precision errors.
        %p(imag(p)~=0)=real(p);
        w(zSel)=-1+p-p.^2/3+(11/72)*p.^3;
    end
    
    %The first two terms of equation 4.18 as an initialization.
    if(k==0)
        zSel=(z>=1.1);
    else
        zSel=(z>=1.1)|z<0;
    end
    L1=log(abs(z(zSel)));
    L2=log(abs(L1));
    w(zSel)=L1-L2;

    for curIter=1:numIter
        expVal=exp(w);

        w=w-(w.*expVal-z)./((w+1).*expVal-(w+2).*(w.*expVal-z)./(2*w+2));
    end

    %Deal with the cases that could lead to NaNs.
    w(z==-exp(-1))=-1;
    if(k==-1)
        w(z==0)=-Inf;
        w(z==-eps(0))=-751.061559539879081;
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
