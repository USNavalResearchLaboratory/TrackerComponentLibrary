function theString=value2BitString(x)
%%VALUE2BITSTRING This takes a standard numeric type, including floating
%           point types, and creates a string of characters with '0's for
%           the zero bits and '1's for the 1 bits in the data type. Unlike
%           dec2bin, this displays all of the bits in the data type and it
%           does not try to turn floating point types into integers.
%           Rather, it just displays the raw bits. If a complex number is
%           passed, then strings for the real and complex parts will be
%           returned.
%
%INPUTS: x A real or complex scalar numeric value. Supported types are
%          double, single, char, int8, int16, int32, int64, uint8, uint16,
%          uint32, and uint64.
%
%OUTPUTS: theString A 1XnumBits array of '0' and '1' characters
%                   representing x,  if x was real or a 2XnumBits matrix of
%                   characters where the first row is the binary character
%                   representation of the real part of x and the second row
%                   is the binary representation of the imaginary part of
%                   x. The first entry in the string is the first bit in
%                   the number.
%                
%EXAMPLE 1:
% x=Inf;
% value2BitString(x)
%One should get the string
%'0000000000000000000000000000000000000000000000000000111111111110'
%because X is a floating point double by default.
%
%EXAMPLE 2:
% x=uint8(254);
% value2BitString(x)
%One should get the string
%'01111111'
%EXAMPLE 3:
% x=int16(-33+11*1j);
% value2BitString(x)
%One will get the two character strings
% '1111101111111111'
% '1101000000000000'
%for the real and imaginary parts of x.
%
%January 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(~isscalar(x))
    error('x must be scalar.')
end

if(isreal(x))
    numStrings=1;
else
    numStrings=2;
    x=[real(x);imag(x)];
end

switch(class(x))
    case 'double'
        numBits=64;
        x=typecast(x,'uint64');
        oneVal=uint64(1);
    case 'single'
        numBits=32;
        x=typecast(x,'uint32');
        oneVal=uint32(1);
    case 'char'
        numBits=8;
        oneVal=char(1);
    case 'int8'
        numBits=8;
        oneVal=int8(1);
    case 'int16'
        numBits=16;
        oneVal=int16(1);
    case 'int32'
        numBits=32;
        oneVal=int32(1);
    case 'int64'
        numBits=64;
        oneVal=int64(1);
    case 'uint8'
        numBits=8;
        oneVal=uint8(1);
    case 'uint16'
        numBits=16;
        oneVal=uint16(1);
    case 'uint32'
        numBits=32;
        oneVal=uint32(1);
    case 'uint64'
        numBits=64;
        oneVal=uint64(1);
    otherwise
        error('x is an unsupported data type.')
end

theString=repmat('0',[numStrings,numBits]); 

for curString=1:numStrings
    selBit=oneVal;
    for curBit=1:numBits
        if(bitand(x(curString),selBit)~=0)
           theString(curString,curBit)='1'; 
        end

        selBit=bitshift(selBit,1);
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
