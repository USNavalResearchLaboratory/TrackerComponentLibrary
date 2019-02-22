function val=harmonicNum(n,algorithm)
%%HARMONICNUM Determine the nth Harmonic number. The nth harmonic number is
%             \sum_{k=1}^n(1/k). Harmonic numbers arise in a number of
%             areas of statistics.
%
%INPUTS: n A scalar or matrix of Harmonic number positions. Non-integer
%          values are allowed. As are negative values. Note that negative
%          integers produce an infinite result.
% algorithm This chooses how the number is computed. Possible values are
%          0 (The default if omitted or an empty matrix is passed) Compute
%            it from the polygamma function.
%          1 Compute it by directly evaluating the \sum_{k=1}^n(1/k).
%
%OUTPUTS: val The values of the Harmonic numbers at the positions given in
%             n.
%
%The number is determined just using the identity of [1] that the nth
%harmonic number is the sum of the EulerMascheroni constant and the digamma
%function of n+1. The digamma function is psi(0,n+1) in Matlab.
%
%One example of harmonic numbers arising is in the coupon collectors
%problem. Given that there is an urn of n coupons and one draws from the
%urn with replacement (uniform probability), what is the expected number of
%coupons one has to draw to make sure that one has drawn all of the coupons
%at least once? The solution is n*harmonicNum(n).
%
%REFERENCES:
%[1] Sondow, Jonathan and Weisstein, Eric W. "Harmonic Number." From
%    MathWorld --A Wolfram Web Resource.
%    http://mathworld.wolfram.com/HarmonicNumber.html
%
%November 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(algorithm))
    algorithm=0;
end
switch(algorithm)
    case 0
        val=Constants.EulerMascheroni+polygamma(0,n+1);
	case 1
        val=0;
        for k=1:n
           val=val+1/k; 
        end
    otherwise
        error('Unknown algorithm specified.')
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
