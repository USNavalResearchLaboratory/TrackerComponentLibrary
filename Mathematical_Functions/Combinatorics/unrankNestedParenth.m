function pString=unrankNestedParenth(rankVal,n)
%%UNRANKNESTEDPARENTH Given a lexicographic rank from 1 to
%               CatalanNumber(n), obtain a character string of nested
%               parentheses, The ordering is such that ')' is considered
%               lexicographically smaller than ')'.
%
%INPUTS: rankVal The rank value from 1 to CatalanNumber(n).
%              n The number of parenthesis pairs that are nested. n>=1.
%
%OUTPUTS: pString The (2*n)X1 character string of the nested parentheses.
%
%This function implements ALgorithm U of Chapter 7.2.1.6 of [1].
%
%REFERENCES
%[1] D. E. Knuth, The Art of Computer Programming. Vol. 4A: Combinatorial
%    Algorithms, Part I, Boston: Addison-Wesley, 2011.
%
%October 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%Allocate space.
pString=repmat(' ',2*n,1);

q=n;
m=1;
p=1;
c=1;

while(p<n)
    p=p+1; 
    c=((4*p-2)*c)/(p+1);
end

while(q~=0)
    while(1)
        cPrime=((q+1)*(q-p)*c)/((q+p)*(q-p+1));
        if(rankVal<=cPrime)
            q=q-1;
            c=cPrime;
            pString(m)=')';
            m=m+1;
            break;
        end

        p=p-1;
        c=c-cPrime;
        rankVal=rankVal-cPrime;
        pString(m)='(';
        m=m+1;
        
        %If a rank that is greater than the maximum number of possible
        %values or that was less than 1 was given.
        if(m>2*n)
            pString=[];
            return;
        end
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
