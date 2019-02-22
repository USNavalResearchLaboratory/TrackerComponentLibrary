function theMask=getBitMask(bitRange,dataType)
%%GETBITMASK Get an integer data type with all bit in the specified range
%            equal to ones. This is useful for masking values using bitwise
%            functions (e.g. bit and, bitor, etc.).
%
%INPUTS: bitRange The 2X1 or 1X2 range of bits (starting from 0) that
%                 should be set to 1 in the mask. All other bits are zero.
%                 The top of the range is limited to the length of the
%                 selected data type -1. bitRange(1)<bitRange(2) or the
%                 mask will be all zeros. Alternatively, one can pass a
%                 sscalr in which case only a single bit is set. This is
%                 the same as making the start and end values in bitRange
%                 the same.
%        dataType An optional parameter specifying the data type of the
%                 mask that should be returned. This is a character string.
%                 The default if omitted or an empty matrix is passed is
%                 'uint64'. Possible values are 'int8', 'int16', 'int32',
%                 'int64', 'uint8', 'uint16', 'uint32', and 'uint64'.
%
%OUTPUTS: theMask A mask having 1's in the specified range. This is the
%                 specified data type. Note that Matlab does not always
%                 treat integer data types int he same manner as floating
%                 point types.
%
%EXAMPLE 1:
% theMask=getBitMask(0)
%theMask is just equal to 1.
%
%EXAMPLE 2:
% theMask=getBitMask([52,61])
%theMask is equal to 4607182418800017408, which is the number with bits 52
%to 61 all equal to 1.
%
%January 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(dataType))
    dataType='uint64';
end

if(~isreal(bitRange))
    error('The bit range must be real.')
end

if(isscalar(bitRange))
    bitRange=[bitRange;bitRange];
end

bitStart=bitRange(1);

if(bitStart<0)
    error('The start of the bit range cannot be less than 0.')
end

bitEnd=bitRange(2);

switch(dataType)
    case 'int8'
        if(bitEnd>7)
           error('The top of the bit range cannot be larger than 7 for an int8 data type.')
        end
        
        theMask=int8(0);
    case 'int16'
        if(bitEnd>15)
           error('The top of the bit range cannot be larger than 15 for an int16 data type.')
        end
        
        theMask=int16(0);
    case 'int32'
        if(bitEnd>31)
           error('The top of the bit range cannot be larger than 31 for an int32 data type.')
        end
        
        theMask=int32(0);
    case 'int64'
        if(bitEnd>63)
           error('The top of the bit range cannot be larger than 63 for an int64 data type.')
        end
        
        theMask=int64(0);
    case 'uint8'
        if(bitEnd>7)
           error('The top of the bit range cannot be larger than 7 for an uint8 data type.')
        end
        
        theMask=uint8(0);
    case 'uint16'
        if(bitEnd>15)
           error('The top of the bit range cannot be larger than 15 for a uint16 data type.')
        end
        
        theMask=uint16(0);
    case 'uint32'
        if(bitEnd>31)
           error('The top of the bit range cannot be larger than 31 for a uint32 data type.')
        end
        
        theMask=uint32(0);
    case 'uint64'
        if(bitEnd>63)
           error('The top of the bit range cannot be larger than 63 for a uint64 data type.')
        end
        
        theMask=uint64(0);
    otherwise
        error('Unknown data type specified.')
end

for curBit=bitStart:bitEnd
    bitVal=bitshift(1,curBit); 
    theMask=bitor(theMask,bitVal);
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
