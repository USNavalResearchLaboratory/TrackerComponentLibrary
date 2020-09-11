function KLDis=KLDistGauss(mu1,P1,mu2,P2)
%%KLDISTGAUSS Compute the Kullback-Leibler (KL) divergence (also known as
%            the KL distance or relative entropy) between two multivariate
%            Gaussian distributions. As defined in Chapter 8.5 of [1], the
%            relative entropy  between two continuous PDFs, f and g, is
%            D(f||g)=integral_{x} f(x) log(f(x)/g(x)) dx
%            The integral is over an appropriate domain of x, in this case
%            the n-dimensional real space for an n-dimensional state. Note
%            that changing the order of f and g changes the result (the
%            operation is not commutative).
%
%INPUTS: mu1, P1 The nX1 mean and nXn covariance matrix of the first
%                Gaussian PDF.
%        mu2, P2 The nX1 mean and nXn covariance matrix of the second
%                Gaussian PDF.
%
%OUTPUTS: val The value of the Kullback-Leibler divergence between the
%             first and the second Gaussian distributions.
%
%The KL distance between the PDF f_1(x) and the PDF f_2(x) is given by
%integral f_1(x)*log(f_1(x)/f_2(x)) dx
%The KL distance can be difference if f_1 and f_2 are switched. Here,
%mu1 and P1 correspond to f_1 and mu2 and P2 to f_2.
%
%The Kullback-Leiber divergence between two Gaussian PDFs is used in
%Gaussian mixture reduction algorithms and an explicit expression is given
%in Theorem 1 in [2].
%
%REFERENCES:
%[1] T. M. Cover and J. A. Thomas, Elements of Information Theory, 2nd ed.
%    Hoboken, NJ: Wiley-Interscience, 2006.
%[2] A. R. Runnalls, "Kullback-Leibler approach to Gaussian mixture
%    reduction," IEEE Transactions on Aerospace and Electronic Systems,
%    vol. 43, no. 3, pp. 989-999, Jul. 2007.
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.
    
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
