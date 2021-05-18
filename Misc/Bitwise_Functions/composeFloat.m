function x=composeFloat(mantissa,exponent,signBit,exponentDecoded)
%%COMPOSEFLOAT Given the mantissa, exponent and sign bit of a floating
%              point number, put them together to get a double or a single.
%              This function is the opposite of decomposeFloat. The
%              floating point standard is in [1].
%
%INPUTS: mantissa A utin64 or int64 (if composing into a double) or a
%                 utin32 or int32 (if composing into a single) value
%                 holding the mantissa of the number. Only the first 53
%                 bits (plus possibly the 54th bit is one added a leading
%                 sign bit) should be nonzero for doubles. For singles,
%                 only the first 23 bits (and possibly the 24th if adding a
%                 leading zero) should be nonzero.
%        exponent The exponent term in the floating point nunmber. This can
%                 be decoded, which means that the bias is subtracted and
%                 exponent must be a signed value (int64 and int32 for
%                 doubles, singles) otherwise it can be a signed or
%                 unsigned integer value that is 64 bits for a double or 32
%                 bits for a single.
%         signBit The binary value (1 or 0) indicating the sign of the
%                 number being composed (1 is negative). This should be a
%                 uint64 or int64 for doubles or a uint32 or int32 for
%                 singles.
% exponentDecoded If true, exponent is a decoded value. Otherwise, it is
%                 not decoed. The default if this parameter is omitted or
%                 an empty matrix is passed is true.
%
%OUTPUTS: x The double or single value composed of the given parts.
%
%See the comments to the function decomposeFloat for more details from [1]
%on the floating point format. 
%
%EXAMPLE:
%Here, we show that if we decompose a number, we can put it back together
%again.
% xOrig=-1.23896e-12;
% [mantissa,exponent,signBit]=decomposeFloat(xOrig);
% x=composeFloat(mantissa,exponent,signBit);
% x==xOrig
%One will see that the statement x==xOrig is 1 (true).
%
%REFERENCES:
%[1] IEEE Standard for Floating Point Arithmetic, Institute for Electrical
%    and Electronics Engineers Std. IEEE Std 754-2008, 29 Aug. 2008.
%
%January 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4||isempty(exponentDecoded))
    exponentDecoded=true;
end

if(isa(mantissa,'uint64')||isa(mantissa,'int64'))
   %If we are constructing a double. 
   
   if(exponentDecoded)
       if(~isa(exponent,'int64'))
           error('When decoded, the exponent of a double should stored in an int64.');
       end
       
       exponent=uint64(exponent+1023);
   elseif(~isa(exponent,'uint64')&&~isa(exponent,'int64'))
       error('When not decoded, the exponent of a double should have its raw bits stored in a uint64 or an int64.');
   end
   
   %Get rid of any possibly leading 1 bit in the mantissa.
    mantissa=bitset(mantissa,53,0);
    
    %Shift everything to the proper positions and put them together.
    x=typecast(bitor(bitor(mantissa,bitshift(exponent,52)),bitshift(signBit,63)),'double');
else
    %If we are constructing a single.
   if(exponentDecoded)
       if(~isa(exponent,'int32'))
           error('When decoded, the exponent of a single should stored in an int32.');
       end
       
       exponent=uint32(exponent+127);
   elseif(~isa(exponent,'uint32')&&~isa(exponent,'int32'))
       error('When not decoded, the exponent of a single should have its raw bits stored in a uint32 or an int32.');
   end
    
    %Get rid of any possibly leading 1 bit in the mantissa.
    mantissa=bitset(mantissa,24,0);
       
    %Shift everything to the proper positions and put them together.
    x=typecast(bitor(bitor(mantissa,bitshift(exponent,23)),bitshift(signBit,31)),'single');
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
