function [mantissa,exponent,signBit]=decomposeFloat(x,decodeExp,addLeading1)
%%DECOMPOSEFLOAT Given a floating point number, extract the sign but, the
%          mantissa and the exponent as three separate integer quantities.
%          Unlike the function [F,E]=log2(x), this provies the raw bits of
%          the mantissa rather than another float as in the F returned by
%          log2. The floating point value is assumed to be in IEEE 754
%          format, which is documented in [1].
%
%INPUTS: x The real, scalar floating point value to decompose. This can be
%          of class 'single' or 'double'.
% decodeExp The raw exponent of the floating point number is biased by 1023
%          for a double and by 127 for a single. The bits are shifted (from
%          the original double) so that the first bit in the return value
%          exponent is the least significant bit in the exponent. If this
%          is true, then the bias is removed. The default if this parameter
%          is omitted or an empty matrix is passed is true. Note that when
%          decoded, the exponent for 1 is 0, but for 0 it is actually -1023
%          for a double and is -127 for a single.
% addLeading1 For mantissas of normalized numbers, if this is true (the
%          default), then the leading 1 that is implied but not expressly
%          stored in the floating point register will be added. No 1 will
%          be added for subnormal numbers, where it is not implied that
%          there is a leading 1. Note that "leading" means that it is added
%          past the last bit in the mantissa.
%
%OUTPUTS: mantissa The mantissa of the float as a uint64 if x is a double
%                  or uint32 if x is a single. If addLeading1=true, then
%                  there are 53 significant bits for double and 24
%                  significant bits for a single. Otherwise, there are
%                  respectively 52 and 23 significant bits for a double or
%                  for a single.
%         exponent The exponent value. If decodeExp is false, then this is
%                  a uint64 if x is a double or a uint32 if x is a single.
%                  If decodeExp is true, then this is respectively an int64
%                  or an int32 for a double or for a single.
%          signBit This is 1 if the sign bit in x is 1 (negative) and 0 if
%                  the sign bit in x is zero. This is a uint64 if x is a
%                  double or a uint32 if x is a single.
%
%The IEEE 754 standard for floating point arithmetic in [1] formats a
%double as a 64-bit floating point data type as
%[mantissa (52 bits)] [exponent (11 bits)] [sign (1 bit)]
%where bit numbering starts at 0 in the mantissa and goes through 63 with
%the sign bit. For a 32-bit floating point single, the bit ordering is
%[mantissa (23 bits)] [exponent (8 bits)] [sign (1 bit)]
%If exponent is all ones and the mantissa is zero, then the value Inf is
%represented with the sign given by the sign bit. If the exponent is all
%ones and the mantissa is any nonzero value, then a NaN is given.
%
%To interpret the mantissa, consider that there is an implied binary point
%after the highest bit and if exponent is nonzero, then an implied 1. If
%exponent is zero, then it is implied that a 0 leads the binary point.
%Thus, the mantissa represents a fraction >=0 and < 0.5. If one reverses
%the order of the bits in the mantissa (because we wrtie numbers in big-
%endian format), then the number represented by the float is of the form
%1.mantissa*2^(exponent-bias) for a nonzero exponent and
%0.mantissa*2^(-bias) for a zero exponent. The bias is 1023 for a double
%and is 127 for a single.
%
%EXAMPLE:
% [mantissa,exponent,signBit]=decomposeFloat(4568.321)
%One gets mantissa=5022922058913284, exponent=12, signBit=0
%
%EXAMPLE:
%When dealing with integers, the mantissa with the leading 1 added
%multiplied by a power of two (shifted) is the integer value. Consider
% [mantissa,exponent,signBit]=decomposeFloat(11)
% %One gets mantissa=6192449487634432, exponent=3, signBit=0
% %Now, shift until the first nonzero bit is at the start.
% mantissaShifted=bitshift(mantissa,-findPosOfMin1Bit(mantissa)+1)
%One will see that the shifted mantissa is the orignal value, 11. The
%actual amount one has to shift the mantissa also depends on the exponent.
%
%REFERENCES:
%[1] IEEE Standard for Floating Point Arithmetic, Institute for Electrical
%    and Electronics Engineers Std. IEEE Std 754-2008, 29 Aug. 2008.
%
%January 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(decodeExp))
    decodeExp=true;
end

if(nargin<3||isempty(addLeading1))
    addLeading1=true;
end

if(~isreal(x)||~isscalar(x))
   error('x must be a real, scalar, floating point value.') 
end

if(isa(x,'double'))
    N=typecast(x,'uint64');
    
    %Get the sign bit.
    theMask=getBitMask(63,'uint64');
    signBit=bitshift(bitand(N,theMask),-63);
    
    %Select the bits that make up the exponent.
    expMask=getBitMask([52,62],'uint64');
    exponent=bitshift(bitand(N,expMask),-52);

    %Select the bits that make up the mantissa.
    theMask=getBitMask([0,51],'uint64');
    mantissa=bitand(N,theMask);
    if(addLeading1&&exponent~=0&&exponent~=expMask)
        %If a leading 1 is implied by the floating point number, then add
        %it. No leading 1 is added if theExp is the value indicating an INf
        %or a NaN.

        mantissa=bitor(mantissa,bitshift(uint64(1),52));
    end

    if(decodeExp)
        %Get rid of the bias. Note that two exceptions are that +/-0 maps
        %to -1023 instead of 0 and +/-Inf and NaN map to 1024.
        exponent=int64(exponent)-1023;
    end
elseif(isa(x,'single'))
    N=typecast(x,'uint32');
    
    %Get the sign bit.
    theMask=getBitMask(31,'uint32');
    signBit=bitshift(bitand(N,theMask),-31);
    
    %Select the bits that make up the exponent.
    expMask=getBitMask([23,30],'uint32');
    exponent=bitshift(bitand(N,expMask),-23);
    
    %Select the bits that make up the mantissa.
    theMask=getBitMask([0,22],'uint32');
    mantissa=bitand(N,theMask);
    if(addLeading1&&exponent~=0&&exponent~=expMask)
        %If a leading 1 is implied by the floating point number, then add
        %it. No leading 1 is added if theExp is the value indicating an INf
        %or a NaN.

        mantissa=bitor(mantissa,bitshift(uint32(1),23));
    end

    if(decodeExp)
        %Get rid of the bias. Note that two exceptions are that +/-0 maps
        %to -127 instead of 0 and +/-Inf and NaN map to 127.
        exponent=int32(exponent)-127;
    end
else
    error('x must be a single or double floating point value.')
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
