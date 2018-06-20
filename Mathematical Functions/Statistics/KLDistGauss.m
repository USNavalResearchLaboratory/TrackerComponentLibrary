function KLDis=KLDistGauss(mu,P)
%%KLDISTGAUSS Obtain the Kullback-Leibler (KL) distance between two
%             multivariate Gaussian distributions.
%
%INPUTS: mu   An xDim X 2 matrix of the mean vectors of the 2 Gaussian
%             distributions whose KL distance is desired.
%        P    An xDim X xDim X2 hypermatrix of the covariance matrices of
%             the two distributions whose KL distance is desired.
%
%OUTPUTS: KLDis The KL distance between the two distributions.
%
%The KL distance between the PDF f_1(x) and the PDF f_2(x) is given by
%integral f_1(x)*log(f_1(x)/f_2(x)) dx
%The KL distance can be difference if f_1 and f_2 are switched. Here,
%mu(:,1) and P(:,:,1) correspond to f_1 and mu(:,2) and P(:,:,2) to f_2.
%
%The expression for the Kullback Leibler distance between two multivariate
%Gaussian distributions is from Theorem 1 of [1].
%
%REFERENCES:
%[1] A. R. Runnalls, "Kullback-Leibler approach to Gaussian mixture
%    reduction," IEEE Trans. Aerosp. Electron. Syst., vol. 43, no. 3, pp.
%    989-999, Jul. 2007.
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    mu1=mu(:,1);
    mu2=mu(:,2);
    P1=P(:,:,1);
    P2=P(:,:,2);
    
    KLDis=0.5*(trace(lsqminnorm(P2,(P1-P2+(mu1-mu2)*(mu1-mu2)')))+log(det(P2)/det(P1)));
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
