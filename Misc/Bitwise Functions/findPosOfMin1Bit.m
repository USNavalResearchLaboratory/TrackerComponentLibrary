function bitPos=findPosOfMin1Bit(n)
%%FINDPOSOFMIN1BIT Given n, which is of an integer data type, determine the
%           position of the least significant bit that is a 1.
%
%INPUTS: n An integer or a matrix of integers. These must be an integer
%          data type (i.e. not a double).
%
%OUTPUTS: bitPos The position of the first 1 in the binary representation
%                of each of the integers in n. Counting starts at 1. A
%                value of 0 means that there is no 1 in the binary
%                representation.
%
%The algorithm just keeps shifting the number (dividing by 2) until it
%the first digit is no longer 0, counting the number of shifts needed.
%
%EXAMPLE:
% n=findPosOfMin1Bit(int32([0,4,1213,2^23]))
%
%The function returns n=[0, 3, 1, 24]. Binary representations of
%the numbers (the least significant bit is to the left) are
%0    0
%4    001
%1213 10111101001
%2^23 000000000000000000000001
%
%January 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(isfloat(n)||~isreal(n))
    error('The values must be real and be an integer data type.')
end

if(ischar(n))
    n=uint8(n); 
end

bitPos=ones(size(n));
numEls=numel(n);

for curEl=1:numEls
    if(n(curEl)~=0)    
        while(bitand(n(curEl),1)==0)
            bitPos(curEl)=bitPos(curEl)+1;
            n(curEl)=bitshift(n(curEl),-1);
        end
    else
        bitPos(curEl)=0;
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
