function theInt=reverseBits(theInt,bitRange)
%%REVERSEBITS Given an integer data type, reverse the ordering of all of
%             the bits or just the bits in a certain range. This can be
%             useful when decoding data with different bit orderings, as
%             can often arise in certain networked data formats.
%
%INPUTS: theInt An integer value having an integer data types, such as
%               uint64.
%      bitRange The range of bits that should be flipped. Counting goes
%               from 0 to the number of bits in the integer-1. If this
%               parameter is omitted or an empty matrix is passed, then all
%               of the bits are reversed.
%
%OUTPUTS: theInt The integer with the ordering of the specified bits
%                reversed.
%
%January 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2)
   bitRange=[]; 
end

if(~isreal(bitRange))
    error('The bit range must be real.')
end

if(isscalar(bitRange))
    return;
end

if(~isempty(bitRange))
    bitStart=bitRange(1);

    if(bitStart<0)
        error('The start of the bit range cannot be less than 0.')
    end
    bitEnd=bitRange(2);

    classType=class(theInt);
    switch(classType)
        case 'int8'
            if(bitEnd>7)
               error('The top of the bit range cannot be larger than 7 for an int8 data type.')
            end
        case 'int16'
            if(bitEnd>15)
               error('The top of the bit range cannot be larger than 15 for an int16 data type.')
            end
        case 'int32'
            if(bitEnd>31)
               error('The top of the bit range cannot be larger than 31 for an int32 data type.')
            end
        case 'int64'
            if(bitEnd>63)
               error('The top of the bit range cannot be larger than 63 for an int64 data type.')
            end
        case 'uint8'
            if(bitEnd>7)
               error('The top of the bit range cannot be larger than 7 for an uint8 data type.')
            end
        case 'uint16'
            if(bitEnd>15)
               error('The top of the bit range cannot be larger than 15 for a uint16 data type.')
            end
        case 'uint32'
            if(bitEnd>31)
               error('The top of the bit range cannot be larger than 31 for a uint32 data type.')
            end
        case 'uint64'
            if(bitEnd>63)
               error('The top of the bit range cannot be larger than 63 for a uint64 data type.')
            end
        otherwise
            error('reverseBits only works on integer data types, such as uint64.')
    end
else
   bitStart=0;
       classType=class(theInt);
    switch(classType)
        case 'int8'
            bitEnd=7;
        case 'int16'
            bitEnd=15;
        case 'int32'
            bitEnd=31;
        case 'int64'
            bitEnd=63;
        case 'uint8'
            bitEnd=7;
        case 'uint16'
            bitEnd=15;
        case 'uint32'
            bitEnd=31;
        case 'uint64'
            bitEnd=63;
        otherwise
            error('reverseBits only works on integer data types, such as uint64.')
    end
end
    
numBits=bitEnd-bitStart+1;

for curBit=0:floor((numBits-1)/2)
   bit1Idx=bitStart+curBit;
   bit2Idx=bitEnd-curBit;
   
   mask1=getBitMask(bit1Idx,classType);
   bit1=(bitand(theInt,mask1)~=0);
   
   mask2=getBitMask(bit2Idx,classType);
   bit2=(bitand(theInt,mask2)~=0);
   
   theInt=bitset(theInt,bit1Idx+1,bit2);
   theInt=bitset(theInt,bit2Idx+1,bit1);
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
