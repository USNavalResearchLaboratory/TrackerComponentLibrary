function V=permuteMatrixCols(V,idx)
%%PERMUTEMATRIXCOLS This function is equivalent in Matlab to V=V(:,idx),
%   where idx is a permutation vector. That is, it rearranges the columns
%   of v according to the specified permutation. The algorithm implemented
%   here performs and in-place permutation and thus can be a template for
%   functions in programming lanagues, such as C, where one has more
%   control over the allocation and deallocation of memory. The idx vector
%   is not modified. This is just a matrix version of the function
%   permuteVector.
%
%INPUTS: V A numRowsXN matrix.
%      idx A length N permutation vector. This is a permutation of the
%          values 1:N.
%
%OUTPUTS: V The same as the input V but with the columns permuted
%           according to idx.
%
%The basic algorithm is that in [1], but corrections have been made. The
%while loop as written in [1] is usually an infinite loop, because ko is
%never changed.
%
%EXAMPLE:
%One can demonstrate the correctness of the algorithm by verifying that
% p=randperm(180);V=rand(4,180);all(all(permuteMatrixCols(V,p)==V(:,p)))
%is true.
%
%REFERENCES:
%[1] A. J. W. Duijvestin, "Correctness proof of an in-place permutation,"
%    BIT Numerical Mathematics, vol. 12, pp. 318-324, Sep. 1972.
%
%October 2020 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

G=size(V,2);
for k=1:G
    kn=idx(k);
    while(kn<k)
        kn=idx(kn);
    end
    if(kn~=k)
        wr=V(:,kn);
        V(:,kn)=V(:,k);
        V(:,k)=wr;
    end
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
