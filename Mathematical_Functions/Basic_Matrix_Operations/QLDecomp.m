function [Q,L]=QLDecomp(A)
%%QLDECOMP Find a Q and an L such that A=Q*L, where L is a lower triangular
%          matrix.
%
%INPUTS: A An mXn matrix.
%
%OUTPUTS: Note that the return values differ depending on whether one or
%         two outputs are requested. The function is either called as
%         L=QLDecomp(A) or [Q,L]=QLDecomp(A). The Q and L are
%         Q An mXm matrix.
%         L An mXn lower-triangular matrix.
%
%This just calls the qr function, rotating the inputs and then rotates them
%back in the end to get Q and L.
%
%June 2023 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

[Q,R]=qr(rot90(A,2));
Q=rot90(Q,2);
L=rot90(R,2);

if(nargout==1)
    %If only the L matrix is desired.
    Q=L;
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
