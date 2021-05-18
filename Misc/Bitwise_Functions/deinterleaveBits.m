function [x,y]=deinterleaveBits(z)
%%DEINTERLEAVEBITS Given an integer data type, extract a number x from the
%            even bits and another number y from the odd bits. The output
%            has half the number of bits as the input, except it an 8-bit
%            data type is provided. The output is signed/unsigned the same
%            as the input, except the int8 datatype cannot be used on the
%            input, through the utint8 type can. This function is the
%            inverse of the interleaveBits function.
%
%INPUTS: z An integer data type. This can be uint8, int16, uint16, int32,
%          uint32, int64, uint64.
%
%OUTPUTS: x,y The two integers formed by deinterleaving the bits.
%
%This function puts the even and odd bits into x and y. Most of the
%implementation code essentially deals with typecasting since Matlab's
%bitwise operator functions require that all arguments be of the same type.
%
%EXAMPLE:
% [x,y]=deinterleaveBits(interleaveBits(int32(100),-int32(102)))
%One will get x=100 and y=-102, which are the same as the inputs.
%
%January 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

theClass=class(z);
switch(theClass)
    case 'int8'
        error('8-bit datatypes must be unsigned')
    case 'int16'
        numBits=8;
        x=uint16(0);
        y=uint16(0);
        one=uint16(1);
        z=typecast(z,'uint16');
    case 'int32'
        numBits=16;
        x=uint32(0);
        y=uint32(0);
        one=uint32(1);
        z=typecast(z,'uint32');
    case 'int64'
        numBits=32;
        x=uint64(0);
        y=uint64(0);
        one=uint64(1);
        z=typecast(z,'uint64');
    case 'uint8'
        numBits=4;
        x=uint8(0);
        y=uint8(0);
        one=uint8(1);
    case 'uint16'
        numBits=8;
        x=uint16(0);
        y=uint16(0);
        one=uint16(1);
    case 'uint32'
        numBits=16;
        x=uint32(0);
        y=uint32(0);
        one=uint32(1);
    case 'uint64'
        numBits=32;
        x=uint64(0);
        y=uint64(0);
        one=uint64(1);
    otherwise
        error('Unknown data type inputs provided.')
end

for curBit=0:(numBits-1)
    x=bitor(x,bitshift(bitand(z,bitshift(one,2*curBit)),-curBit));
    y=bitor(y,bitshift(bitand(z,bitshift(one,2*curBit+1)),-curBit-1));
end

%Convert to the correct type.
switch(theClass)
    case 'int16'
        x=typecast(uint8(x),'int8');
        y=typecast(uint8(y),'int8');
    case 'int32'
        x=typecast(uint16(x),'int16');
        y=typecast(uint16(y),'int16');
    case 'int64' 
        x=typecast(uint32(x),'int32');
        y=typecast(uint32(y),'int32');
    case 'uint8'
    case 'uint16'
        x=uint8(x);
        y=uint8(y);
    case 'uint32'
        x=uint16(x);
        y=uint16(y);
    case 'uint64'
        x=uint32(x);
        y=uint32(y);
    otherwise
        error('An unknown error occurred')
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
