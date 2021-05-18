function num=CatalanNumber(p,q)
%%CATALANNUMBER Get the nth Catalan number C(n), if only one input is
%               provided, or find the generalized Catalan number C(p,q), if
%               two inputs are provided. The nth Catalan number C(n) is the
%               total number of possible binary trees that can be formed
%               with n nodes. The generalized Catalan numbers C(q,n) are
%               also known as ballot numbers. They represent sequences of
%               p+q ballots for which a running tabulation never favors a
%               candidate with p votes of q. q>=p. Note that C(n,n)=C(n).
%
%INPUTS: All inputs are integers >=1. If only one input is given, then C(n)
%        is the nth Catalan number. If two are given, then p<=q and C(p,q)
%        is a generalized Catalan number.
%
%OUTPUTS: num The value of the chosen Catalan number.
%
%EXAMPLE:
% num=CatalanNumber(13)
%One will get num=742900.
%
%Expressions for Catalan numbers and generalized Catalan numbers in terms
%of binomial coefficients are given in Section 7.2.1.6 of [1].
%
%REFERENCES:
%[1] D. E. Knuth, The Art of Computer Programming. Vol. 4A: Combinatorial
%    Algorithms, Part I, Boston: Addison-Wesley, 2011.
%
%October 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(q))
    q=p; 
end

if(q<p)
   error('q must be >=p.') 
end

num=((q-p+1)/(q+1))*binomial(p+q,p);

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
