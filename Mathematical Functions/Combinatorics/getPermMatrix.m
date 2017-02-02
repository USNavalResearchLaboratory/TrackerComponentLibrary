function chi=getPermMatrix(rank,n,dim)
%%GETPERMUTATIONMATRIX Obtain the permutation matrix corresponding to a
%                      permutation of a given lexicographic order. The
%                      matrix can be used to rearrange n stacked vectors
%                      dim dimensions.
%
%INPUTS:    rank The order of the desired permutation of
%                [0 1 2 ... n-1] in lexicographic order. Note 
%                that 0<=rank<=(n!-1).
%           n    The number of elements in the desired permutation; the
%                number of that might be reordered with this matrix.
%           dim  The dimensionality of the n things being reordered.
%
%OUTPUTS:   chi  The permutation matrix (consisting of ones and zeros) that
%                can rearrange n stacked vectors of dimensionality dim.
%
%Given an (n*dim)X1 vector v of n stacked subvectors of dimension dim, 
%chi*v rearranges the subvectors in v according to the given permutation. 
%For a more computationally efficient procedure of rearranging the 
%subvectors, consider using the function getPermIndices.
%
%September 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    perm=unrankPermutation(rank,n)+1;%Get the permutation.
    chi=zeros(n*dim,n*dim);%Allocate space.
    
    %Insert little identity matrices at the correct positions in each
    %collection of dim rows.
    for k=1:n
        minIdxR=(k-1)*dim+1;
        maxIdxR=minIdxR+dim-1;
        minIdxC=(perm(k)-1)*dim+1;
        maxIdxC=minIdxC+dim-1;
        
        chi(minIdxR:maxIdxR,minIdxC:maxIdxC)=eye(dim);
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
