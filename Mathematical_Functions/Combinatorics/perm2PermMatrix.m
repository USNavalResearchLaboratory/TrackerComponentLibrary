function Ip=perm2PermMatrix(p)
%%PERM2PERMMATRIX Given a permutation of the values 1:n, obtain the
%                 permutation matrix that corresponds to the permutation.
%
%INPUTS: p A 1Xn or nX1 permutation of the values 1:n.
%
%OUTPUTS: Ip A permutation matrix such that Ip*(1:n).'=p.
%
%The creation of a permutation matrix from a permutation basically involves
%rearranging the columns of an identity matrix.
%
%EXAMPLE:
% p=[3;1;4;5;2];
% Ip=perm2PermMatrix(p);
% all(Ip*(1:5).'==p)
%
%November 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

n=length(p);
Ip=zeros(n,n);
for i=1:n
    Ip(i,p(i))=1;    
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
