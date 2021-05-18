function A=permuteRowsCols(A,sigma,tau,doInversePermutation)
%%PERMUTEROWSCOLS Permute the rows and columns of a matrix in place. In
%                 Matlab, one can permute the rows and columns of a matrix
%                 using the notation A(sigma,tau), where sigma and tau are
%                 permutation indices for the rows and columns. Also, one
%                 could just write a loop over i and j that does
%                 ARet(i,j)=A(sigma(i),tau(j)); This function rearranges
%                 the rows and the columns of the matrix in-place. That is,
%                 without having to copy the matrix. While not necessarily
%                 that useful in Matlab, this can form a model for writing
%                 such a function in other languages. Also, this function
%                 can perform an inverse permutation to undo a previous
%                 reordering.
%
%INPUTS: A An mXn matrix whose rows and columns one wishes to reorder.
%    sigma A permutation of indices (from 1 to m) for how the rows should
%          be reordered (or how they were reordered if doInversePermutation
%          is true).
%      tau A permutation of indices (from 1 to n) for how the columns
%          should be reordered (or how they were reordered if
%          (doInversePermutation is true).
% doInversePermutation A parameter specifying whether the ordering given by
%          sigma and tau is desired, or whether one should use the inverse
%          permutations to, for example, undo a previous reordering. The
%          default if omitted or an empty matrix is specified is false.
%
%OUTPUTS: A The matrix A with reordered rows and columns.
%
%Reordering the rows and columns of a matrix without using a temporary
%matrix copy is a challenge. The algorithm used here is based on RENUMB in
%[1], after being modified as the book actually gave an algorithm for the
%inverse permutation, whereas one generally wants it for the actual
%permutation.
%
%REFERENCES:
%[1] A. Nijenhuis and H. S. Wilf, Combinatorial Algorithms for Computers
%    and Calculators, 2nd ed. New York: Academic press, 1978.
%
%December 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4||isempty(doInversePermutation))
    doInversePermutation=false;
end

if(doInversePermutation==false)
    [~,~,sigma]=numPermutationCycles(sigma,-1);
    [~,~,tau]=numPermutationCycles(tau,-1);
end

m=size(A,1);
n=size(A,2);

%Step A
%Determine the number of cycles NC of the permutation sigma, the sign IS of
%the permutation, and modify sigma into a tagged permutation.
[~,~,sigma]=numPermutationCycles(sigma,1);
[~,~,tau]=numPermutationCycles(tau,1);

for i=1:m
%Step B, find next cycle of sigma
    i1=-sigma(i);
    if(i1<0)
       continue;
    end
    l=0;
    
%Step C, length of cycle of sigma.
    while(i1>0)
        i1=sigma(i1);
        l=l+1;
    end
    i1=i;
    j=0;

    while(1)
        %Step D, find next cycle of tau
        while(1)
            j=j+1;
            if(j>n||tau(j)<0)
                break;
            end
        end
        
        %If we could not find a new cycle of tau.
        if(j>n)
            break;%Go back to step B.
        end
        j2=j;
        k=l;

    %Step E, start new product cycle.
        startNewProductCycle=true;
        while(startNewProductCycle)
            startNewProductCycle=false;
            j1=j2;
            t1=A(i1,j1);

        %Step F, move matrix element in one product cycle
            while(1)
                i1=abs(sigma(i1));
                j1=abs(tau(j1));
                t2=A(i1,j1);
                A(i1,j1)=t1;
                t1=t2;
                if(j1~=j2)
                    continue;
                end
                k=k-1;
                %Is it the end of aproduct cycle?
                if(i1~=i)
                    continue;
                end
                j2=abs(tau(j2));
                %Are all product cycles complete?
                if(k~=0)
                    startNewProductCycle=true;
                end
                break;
            end
        end
    end
end

%Step F, restore arrays
%Undo changes to sigma and tau. In Matlab, this does not matter, because
%the original arrays are not accessed, but if one were to code something in
%C/ Fortran based off this function, then sigma and tau would have to be
%returned to their original values.
% sigma=abs(sigma);
% tau=abs(tau);

%Undo the inversion. Again, this would only be necessary if this function
%were written in C or a similar language and used the input parameters
%sigma and tau as scratch space rather than creating copies.
% if(doInversePermutation==false)
%     [~,~,sigma]=numPermutationCycles(sigma,-1);
%     [~,~,tau]=numPermutationCycles(tau,-1);
% end

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
