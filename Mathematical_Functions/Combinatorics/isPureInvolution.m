function val=isPureInvolution(x)
%%ISPUREINVOLUTION Return true if x is a pure involution. That is a fixed-
%       point free involution. An involution is a permutation that is the
%       same as its inverse permutation.
%
%INPUTS: x An nX1 or 1Xn permutation of the values 1:n.
%
%OUTPUTS: val This is true if x is a pure involution and false otherwise.
%
%This function first checks whether n is odd (in which case it cannot be a
%pure involution). It then calls isInvolution and finally checks whether
%any of the values form 1 to n are in their original spots.
%
%September 2024 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

n=length(x);
if((mod(n,2)~=0)||(isInvolution(x)==false)||(any((1:n).'==x(:))))
    val=false;
else
    val=true;
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
