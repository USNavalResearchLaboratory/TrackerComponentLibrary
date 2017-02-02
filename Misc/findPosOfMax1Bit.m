function bitPos=findPosOfMax1Bit(n)
%%FINDPOSOFMAX1BIT Given a positive integer, determine the position of the
%                  most significant bit (that is a 1).
%
%INPUTS: n  An integer or a matrix of integers. The integers do not have to
%           be of type uint64 or int32; they can be doubles. However, using
%           a double that is larger than the maximum representable integer
%           without a loss of precision (usually 2^53, determined by the
%           number of bits in the mantissa of the double).
%
%OUTPUTS: bitPos The position(s) of the first 1 in the binary
%            representation of each of the integers in n. Counting starts
%             at 1. A value of 0 means that there is no 1
%            in the binary representation.
%
%The algorithm just keeps shifting the number (dividing by 2) until it
%becomes 0, counting the number of shifts needed.
%
%EXAMPLE:
% n1=findPosOfMax1Bit(0)
% n2=findPosOfMax1Bit(4)
% n3=findPosOfMax1Bit(1213)
% n4=findPosOfMax1Bit(2^23)
%The function returns n1=0, n2=3, n3=11, n4=24. Binary representations of
%the numbers are
%0                           0
%4                         100
%1213              10010111101
%2^23 100000000000000000000000
%
%November 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(n~=fix(n)&&n~=0||n<0)
    error('The values must be positive integers')
end

bitPos=0;
while(n>0)
    bitPos=bitPos+1;
    n = bitshift(n,-1);
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
