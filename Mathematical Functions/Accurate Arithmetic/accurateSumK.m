function res=accurateSumK(p,K)
%%ACCURATESUMK Sum the elements in a vector in a manner that reduces finite
%             precision errors when summing large and small values as
%             compared to sequentially summing the terms. In comparison to
%             the accurateSum function, this function allows for K-fold
%             precision (K=2 is fixed in accurateSum with algorithm 0).
%
%INPUTS: p An nX1 or 1Xn vector of real values.
%        K The multiple of the standard double precision desired for
%          intermediate results. If omitted or an empty matrix is passed,
%          the default is K=3. Values larger than n are clipped to n.
%
%OUTPUTS: res The value of the sum of the values in p.
%
%This function implements algorithm 4.8 in [1].
%
%EXAMPLE:
%The simplest example demonstrates how using higher intermediate precision
%can make sorting the input less important. We compare this function to the
%unsorted accurateSum function.
% a=[1e100;2;1;eps();eps();eps();-1e100];
% res0=accurateSumK(a)
% res1=accurateSum(a,[],false)%false makes it unsorted.
%One gets res0=3.000000000000001 and res=3. res0 is more accurate.
%
%REFERENCES:
%[1] T. Ogita, S. M. Rump, and S. Oishi, "Accurate sum and dot product,"
%    SIAM Journal on Scientific Computing, no. 6, pp. 1955-1988, 2005.
%
%November 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(K))
    K=3;
end

n=length(p);
K=min(K,n);
for k=1:(K-1)
    p=VecSum(p);
end

res=p(1);
for i=2:n
    res=res+p(i);
end
end

function p=VecSum(p)
    n=length(p);
    for k=2:n
        [p1,p2]=exactPairSum(p(k),p(k-1));
        p(k)=p1;
        p(k-1)=p2;
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
