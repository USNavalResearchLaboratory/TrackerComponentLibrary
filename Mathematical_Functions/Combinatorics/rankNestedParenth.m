function rankVal=rankNestedParenth(pString)
%%UNRANKNESTEDPARENTH Obtain the lexicographic order rank of a set of
%           nested parentheses, where in each pair, '(' always comes before
%           ')' for it to be valid, but ')' is considered lexicographically
%           smaller than ')'.
%
%INPUTS: pString A string or character array holding a valid set of nested
%                parenthese, such as pString='((()()()))'. Invalid string
%                will either cause an error or an incorrect result to be
%                produced.
%
%OUTPUTS: rankVal The lexicographic rank of the string, string from 1.
%
%This function implements the modification of Algorithm U of Section
%7.2.1.6 of [1] that is given in Problem 50 in Section 7.2.1.6.
%
%REFERENCES
%[1] D. E. Knuth, The Art of Computer Programming. Vol. 4A: Combinatorial
%    Algorithms, Part I, Boston: Addison-Wesley, 2011.
%
%October 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(isempty(pString))
    rankVal=[];
    return
end

pString=char(pString);
n=length(pString)/2;

q=n;
m=1;
p=1;
c=1;
rankVal=1;

while(p<n)
    p=p+1; 
    c=((4*p-2)*c)/(p+1);
end

while(q~=0)
    while(1)
        if(m>2*n)
            error('Invalid parenthesis character string provided.') 
        end
        
        cPrime=((q+1)*(q-p)*c)/((q+p)*(q-p+1));
        if(pString(m)==')')
            q=q-1;
            c=cPrime;
            m=m+1;
            break;
        end

        p=p-1;
        c=c-cPrime;
        rankVal=rankVal+cPrime;
        m=m+1;
    end
end
if(m~=2*n+1)
    error('Invalid parenthesis character string provided.')
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
