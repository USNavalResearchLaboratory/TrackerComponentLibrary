function z=interleaveBits(x,y)
%%INTERLEAVEBITS Given integers x and y, interleave their bits. For x and y
%            as integer data types, this function produces an integer
%            having twice as many bits and the same sign (signed unsigned).
%            However, if and or y are 64 bits, they are truncated to 32
%            bits and a 64 bit result is returned. If floating point
%            data types are passed, they must be positive integers less
%            than intmax(uint32) and are treated as unsigned 32 bit
%            integers. If doubles are passed, they are truncated.
%
%INPUTS: x, y Two integers as described above. The data types can be int8,
%             uint8, int16, uint16, int32, uint32, int64, uint64, single,
%             or double. x and y must be of the same type.
%
%OUTPUTS: z The integer obtained by interleaving the bits of x and y with
%           the bits of x in the even positions and those of y in the odd
%           positions. This is an integer data type.
%
%The result z is often referred to as a Morton number. Interlevaing bits
%can be used as a method of turning 2D indexation into 1D indexation of
%items. Morton numbers arise in a number of coding applications.
%
%January 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

theClass=class(x);

if(~strcmp(theClass,class(y)))
    error('x and y must be of the same class.')
end

%Here, we set the class of z to have twice the number of bits of x and y.
%We also recast x and y to be the same class as z so that Matlab's bitwise
%functions do not have errors.
switch(theClass)
    case 'int8'
        numBits=8;
        x=int16(x);
        y=int16(y);
        z=int16(0);
    case 'int16'
        numBits=16;
        x=int32(x);
        y=int32(y);
        z=int32(0);
    case 'int32'
        numBits=32;
        x=int64(x);
        y=int64(y);
        z=int64(0);
    case 'int64'
        numBits=32;%Clip to 32 bits.
        
        if(x>intmax('int32')||x<intmin('int32')||y>intmax('int32')||y<intmin('int32'))
            warning('Output limited to 64 bits. Extra input bits are discarded.')
        end
        
        x=int64(x);
        y=int64(y);
        z=int64(0);
    case 'uint8'
        numBits=8;
        x=uint16(x);
        y=uint16(y);
        z=uint16(0);
    case 'uint16'
        numBits=16;
        x=uint32(x);
        y=uint32(y);
        z=uint32(0);
    case 'uint32'
        numBits=32;
        x=uint64(x);
        y=uint64(y);
        z=uint64(0);
    case 'uint64'
        numBits=32;%Clip to 32 bits.
        
        if(x>intmax('uint32')||y>intmax('uint32'))
            warning('Output limited to 64 bits. Extra input bits are discarded.')
        end
        
        x=int64(x);
        y=int64(y);
        z=int64(0);
    case 'single'
        numBits=log2(flintmax('single'));
        
        if(x~=fix(x)||y~=fix(y))
           error('The inputs are not integers'); 
        end
        
        if(x<0||y<0)
            error('Inputs given as floating point data types must be positive.')
        end

        x=uint64(x);
        y=uint64(y);
        z=uint64(0);
    case 'double'
        numBits=32;%Clip to 32 bits.
        
        if(x~=fix(x)||y~=fix(y))
           error('The inputs are not integers'); 
        end
        
        if(x<0||y<0)
            error('Inputs given as floating point data types must be positive.')
        end
        
        if(x>intmax('uint32')||y>intmax('uint32'))
            warning('Output limited to 64 bits. Extra input bits are discarded.')
        end
        
        x=uint64(x);
        y=uint64(y);
        z=uint64(0);
    otherwise
        error('Unknown data type inputs provided.')
end

for curBit=0:(numBits-1)
    bitVal=bitor(bitshift(bitget(x,curBit+1),2*curBit),bitshift(bitget(y,curBit+1),2*curBit+1));
    z=bitor(z,bitVal);
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
