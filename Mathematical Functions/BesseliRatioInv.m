function kappa= BesseliRatioInv(nu,x,maxIterInv,maxIterRat)
%%BESSELRATIOINVERSE Given the value of the modifie Bessel function ratio
%               of the first kind having the form
%               x=I_{nu}(kappa)/I_{nu-1}(kappa)
%
%INPUTS:  nu The positive integer (upper) order of the modified Bessel
%            function of the first kind in the ratio.
%          x The value of the raio I_nu(kappa)/I_{nu-1}(kappa)
% maxIterInv An optional parameter specifying the maximum number of
%            iterations to use when computing the inverse. If this
%            parameter is omitted or an empty matrix is passed, the default
%            value of 50 is used.
% maxIterRat An optional parameter specifying the maximum number of
%            iterations to use in the subroutine BesseliRatio, which is
%            used as part of the Newton's step in the inverse algoirthm. If
%            this parameter is omitted or an empty matrix is passed, then
%            the default value used in the BesseliRatio function is used.
%
%OUTPUTS: kappa The value such that x=I_{nu}(kappa)/I_{nu-1}(kappa).
%
%The algorithm is that of [1] and is based on using Newton's method after
%getting a very simple initial estimate.
%
%As an example, consider
% kappa=2;
% nu=5;
% kappaNew=BesseliRatioInv(nu,BesseliRatio(nu,kappa))
%One should find that kappaNew=kappa=2.
%
%REFERENCES:
%[1] S. Sra, "A short note on parameter approximation for von Mises-Fisher
%    distributions: and a fast implementation of Is(x)," Computational
%    Statistics, vol. 27, no. 1, pp. 177-190, Mar. 2012.
%[2] A. Banerjee, I. S. Dhillon, J. Ghosh, and S. Sra, "Clustering on the
%    unit hypersphere using von Mises-Fisher distributions," Journal of
%    Machine Learning Research, vol. 6, pp. 1345-1382, Jan. - Dec. 2005.
%
%August 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(maxIterInv))
    maxIterInv=50;
end

if(nargin<4||isempty(maxIterRat))
   maxIterRat=[]; 
end

%The approximation in Equation 4 of [1], which cites [2] as the original
%source.
p=2*nu;
kappa = x*(p-x^2)/(1-x^2);%Initial estimate

%The algorithm in [1] uses Newton's method to refine the estimate.
for curIter=1:maxIterInv
    ratVal=BesseliRatio(nu,kappa,maxIterRat);
    %Newton's method step.
    kappaNew=kappa-(ratVal-x)/(1-ratVal^2-(p-1)/kappa*ratVal);
    if(abs(kappaNew-kappa)<=eps(kappa))
        return;%If convergence has occurred.
    end
    kappa=kappaNew;
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
