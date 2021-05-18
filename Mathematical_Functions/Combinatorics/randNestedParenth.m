function a=randNestedParenth(n)
%%RANDNESTEDPARENTH Obtain a random character string of randomly nested
%                   parentheses. Such parantheses can represent random
%                   binary trees.
%
%INPUTS: n The number of open-closed parenthesis pairs.
%
%OUTPUTS: a A (2*n)X1 character string of n parenthesis pairs, with the
%           open parenthesis always coming before the close one.
%
%This function implements Algorithm W of Chapter 7.2.1.6 of [1].
%
%EXAMPLE:
% a=randNestedParenth(5).'
%We transpose it to make it easier to read.
%
%REFERENCES
%[1] D. E. Knuth, The Art of Computer Programming. Vol. 4A: Combinatorial
%    Algorithms, Part I, Boston: Addison-Wesley, 2011.
%
%October 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

p=n;
q=n;
m=1;
%Allocate space for the random set of nested parentheses.
a=repmat(' ',2*n,1);

while(q~=0)
    while(1)
        upperBound=(q+p)*(q-p+1);
        X=randi(upperBound-1);
        
        if(X<(q+1)*(q-p))
            q=q-1;
            a(m)=')';
            m=m+1;
            break; 
        end
        p=p-1;
        a(m)='(';
        m=m+1;
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
