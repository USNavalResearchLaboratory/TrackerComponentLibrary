function a=getBalancedMPartition(n,m)
%%GETBALANCEDMPARTITION Partition the integer n into m parts that are
%           optimally balanced. A partition a(1)>=a(2)>=...>=a(m) is
%           optimally balanced when abs(a(i)-a(j))<=1 for 1<=i,j <=m.
%
%INPUTS: n The integer being partitione; n>=1.
%        m The number of parts into which n is partitioned; 1<=m. If m>n,
%          then the extra parts will be 0. 
%
%OUTPUTS: a An mX1 vector holding the partitioned values, which sum to n.
%
%This function implements Problem 3 of Section 7.2.1.4 of [1].
%
%REFERENCES:
%[1] D. E. Knuth, The Art of Computer Programming. Vol. 4A: Combinatorial
%    Algorithms, Part I, Boston: Addison-Wesley, 2011.
%
%October 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

j=(1:m).';

a=floor((n+m-j)/m);

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
