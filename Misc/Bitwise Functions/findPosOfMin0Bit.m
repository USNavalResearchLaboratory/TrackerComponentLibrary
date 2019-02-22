function bitPos=findPosOfMin0Bit(n)
%%FINDPOSOFMIN1BIT Given a positive integer, determine the position of the
%                  least significant zero bit.
%
%INPUTS: n An integer or a matrix of integers. The integers do not have to
%          be of type uint64 or int32; they can be doubles. However, using
%          a double that is larger than the maximum representable integer
%          without a loss of precision (usually 2^53, determined by the
%          number of bits in the mantissa of the double).
%
%OUTPUTS: bitPos The position(s) of the least significant 0 in the binary
%                representation of each of the integers in n. Counting
%                starts at 1.
%
%This just keeps shifting the numbers (dividng by 2) until twice the last
%number divided by 2 (integer division) equals the last number.
%
%EXAMPLE:
% n1=findPosOfMin0Bit(0)
% n2=findPosOfMin0Bit(4)
% n3=findPosOfMin0Bit(1213)
% n4=findPosOfMin0Bit(2^23)
% n5=findPosOfMin0Bit(2^23-1)
%The function returns n1=1, n2=1, n3=2, n4=1, n5=24. Binary representations
%of the numbers are
%0                             0
%4                           100
%1213                10010111101
%2^23   100000000000000000000000
%2^23-1  11111111111111111111111
%
%November 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(n~=fix(n)&&n~=0||n<0)
    error('The values must be positive integers')
end

bitPos=0;
while(1)
    bitPos=bitPos+1;
    nShift = bitshift(n,-1);
    
    if(n==2*nShift)
        return;
    end
    n=nShift;
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
