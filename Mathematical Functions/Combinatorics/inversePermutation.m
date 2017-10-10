function xInvPerm=inversePermutation(xPerm)
%%INVERSEPERMUTATION Given a permutation of the number 1:N, return the
%           inverse of the permutation. If something were rearranged with
%           the original permutation, the inverse permutation puts it back
%           into the original order.
%
%INPUTS: xPerm A 1XN or NX1 vector containing a permutation of the numbers
%              from 1 to N.
%
%OUTPUTS: xInvPerm The NX1 inverse of the permutation in xPerm.
%
%EXAMPLE:
% xPerm=[2;5;1;4;3];
% y=1:5;
% z=y(xPerm);
% xInvPerm=inversePermutation(xPerm);
% all(z(xInvPerm)==y)
%The result of the final line will  be 1 (true).
%
%April 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

n=length(xPerm);

xInvPerm=zeros(n,1);
for k=1:n
   xInvPerm(xPerm(k))=k; 
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
