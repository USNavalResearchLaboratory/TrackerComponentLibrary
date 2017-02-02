function w=LambW(z)
%%LAMBDW  Evaluate the primary branch of the Lambert-W function, which is
%         also know as a product logarithm. The function solve z=w*exp(w)
%         for w. This implementation is meant for real values of z>=-1.
%
%INPUTS:  z A vector or matrix of parameters at which the primary branch of
%           the Lambert W function is to be evaluated. It is assumed that
%           z>-exp(-1) for all of the elements.
%
%OUTPUTS: w The solutions to z=w*exp(w) for all elements in z.
%
%The origin of the algorithm is [1], where the main algorithm is Halley's
%method from Equation 5.9. For negative values, an initial value of the
%first few terms of Equation 4.22 is used. For positive values less than
%1.1, an initial value of 0 is used. For positive values greater than or
%equal to 1.1, the first two perms of Equation 4.18 are used. This differs
%from the suggested method in the paper in that a Padé approximation is not
%used for terms near zero. Rather, a constant initial estimate of 0 is
%used. 15 iterations are performed, which should be more than enough to
%converge.
%
%REFERENCES:
%[1] R. M. Corless, G. H. Gonnet, D. E. G. Hare, and D. J. Jefrey, "On the
%    Lambert W Function," Advances in Computational Mathematics" vol. 5,
%    no. 1, pp. 329-359, 1996.
%
%February 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    numIter=15;
    w=zeros(size(z));
    
    %Initialize for negative values of z.
    zSel=z<0;
    p=sqrt(2*(exp(1)*z(zSel)+1));
    w(zSel)=-1+p-p.^2/3+(11/72)*p.^3;
    
    %The first two terms of equation 4.18 as an initialization with an
    %extra ad-hoc term to make sure that L1 and L2 are always real for
    %positive values of z.
    zSel=z<1;
    w(zSel)=0;
    zSel=(z>=1.1);
    
    L1=log(z(zSel));
    L2=log(L1);
    w(zSel)=L1-L2;
    
    for curIter=1:numIter;
        expVal=exp(w);
        w=w-(w.*expVal-z)./((w+1).*expVal-(w+2).*(w.*expVal-z)./(2*w+2));
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
