function KLDist=KLDistDiscrete(PMF1,PMF2)
%%KLDISTDISCRETE Compute the Kullback-Leibler (KL) divergence (also known as
%            the KL distance or relative entropy) between two finite-length
%            probability mass functions (PMFs) over the same finite
%            discrete support. As defined in Chapter 17.1 of [1], the
%            relative entropy  between two PMFs is
%            D(f||g)=sum_i PMF1(i) log(PMF1(i)/PMF2(i))
%            In the case of PMF1(i)=PMF2(i)=0, we say that the logarithm
%            equals zero. Note that changing the order of PMF1 and PMF2
%            changes the result (the operation is not commutative).
%
%INPUTS: PMF1 An n1X1 or 1Xn1 set of probability mass values of the first
%             distribution. The values should be all non-negative and sum
%             to 1.
%        PMF2 An n2X1 or 1Xn2 set of probability mass values of the second
%             distribution. n2>=n1. If n2>n1, then only the first n1 values
%             are used.
%
%OUTPUTS: KLDist The value of the Kullback-Leibler divergence between the
%                first and the second PMFs.
%
%REFERENCES:
%[1] T. M. Cover and J. A. Thomas, Elements of Information Theory, 2nd ed.
%    Hoboken, NJ: Wiley-Interscience, 2006.
%
%July 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

n1=length(PMF1);

PMF1=PMF1(:);
PMF2=PMF2(:);

vals=PMF1.*(log(PMF1)-log(PMF2(1:n1)));

%Get rid of instances where PMF1 is zero or they are both zero.
vals(isnan(vals))=0;
KLDist=sum(vals);
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
