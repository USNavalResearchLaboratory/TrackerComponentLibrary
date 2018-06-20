function thePerms=genAllPermutations(N)
%%GENALLPERMUTATIONS Generate all permutations of the integers 1...N. There
%                    are factorial(N) possible permutations.
%
%INPUTS: N The permutation length.
%
%OUTPUTS: thePerms An NXfactorial(N) matrix of all possible length-N
%                  permutations.
%
%This function implements a non-recursive form of heap's algorithms, which
%is presented in [1] and also given in [2].
%
%REFERENCES:
%[1] B. R. Heap, "Permutation by interchanges," The Computer Journal, vol.
%    6, no. 3, pp. 293-298, Nov. 1963.
%[2] R. Sedgewick, "Permutation generation methods," Computing Surveys,
%    vol. 9, no. 2, pp. 137-164, Jun. 1977.
%
%April 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numPerms=factorial(N);
thePerms=zeros(N,numPerms);

%The first permutation.
p=1:N;
thePerms(:,1)=p;

if(N==1)
   return; 
end

%Allocate and initialize
c=ones(N,1);

curPerm=1;
n=1;
while(1)
    if(c(n)<n)
        if(mod(n,2))%If n is odd
            k=1;
        else
            k=c(n);
        end
        
        temp=p(n);
        p(n)=p(k);
        p(k)=temp;

        curPerm=curPerm+1;
        thePerms(:,curPerm)=p;
        
        %We could alternatively check in the while loop at the top whether
        %n<=N to determine completion.
        if(curPerm==numPerms)
            break
        end
        
        c(n)=c(n)+1;
        n=1;
    else
        c(n)=1;
        n=n+1;
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
