function [code,nCard,isLast,j]=getNextGrayCode(code,nCard)
%%GETNEXTGRAYCODE Return the next gray code value in the sequence given the
%                 current gray code value. The first code to start the
%                 sequence is all logical false (zeros), which is returned
%                 if this function is called with an empty matrix for the
%                 code and nCard=n (the length of the code). A gray code is
%                 a  binary code such that only one entry changes from 0 to
%                 1 or back each step. This function can be used to get all
%                 subsets of an n-set.
%
%INPUTS: code An nX1 or 1Xn vector consisting of zeros and ones
%             representing the current gray code value. The first code in
%             the sequence is all zeros. If getNexGrayCode is called with
%             an empty matrix for code and nCard=n, then the returned code
%             will be the first in the sequence.
%       nCard If code is empty, then this is the dimensionality of the code
%             sequence desired. Otherwise, this is the number of ones in
%             code. When getting a next code, nCard can speed things up,
%             but if omitted, it is just found as sum(code).
%
%OUTPUTS: code The next length codeLen gray code value in the sequence. If
%              the final gray code in the sequence was passed, then an
%              empty matrix is returned.
%        nCard The number of ones in the returned code, or n if the last
%              code was passed.
%       isLast True if the returned code is the last in the series and
%              passing it to getNexGrayCode would return an empty matrix.
%            j The index of the entry in code that was changed by this
%              function call. If this is the first function call and code
%              was just created, then j is an empty matrix.
%
%The algorithm is based on NEXSUB from Chapter 1 of [1]. Gray codes are
%also discussed in Chapter 7.2.1.1 of [2].
%
%REFERENCES:
%[1] A. Nijenhuis and H. S. Wilf, Combinatorial Algorithms for Computers
%    and Calculators, 2nd ed. New York: Academic press, 1978.
%[2] D. E. Knuth, The Art of Computer Programming. Vol. 4, Fascicle 2:
%    Generating all Tuples and Permutations, Upper Saddle River, NJ:
%    Addison-Wesley, 2009.
%
%October 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%If called with no code, then return the first code in the sequence.
if(isempty(code))
    n=nCard;
    code=zeros(n,1);
    isLast=false;
    nCard=0;
    j=[];
    return;
end

n=length(code);

%If no cardinality was passed, then compute it.
if(nargin<2)
   nCard=sum(code); 
end

%If the final gray code was passed, then just return empty matrices.
if(nCard==code(n)&&nCard~=0)
    code=[];
    nCard=n;
    isLast=[];
    j=[];
    return;
end

j=1;
if(mod(nCard,2)~=0)
    while(1)
        j=j+1;
        if(code(j-1)~=0)
            break;
        end
    end
end
code(j)=1-code(j);
nCard=nCard+2*code(j)-1;
isLast=(nCard==code(n));

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
