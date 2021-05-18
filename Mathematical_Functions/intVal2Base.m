function [digits,baseVals]=intVal2Base(n,base,firstMostSig,pad2Digits,baseVals)
%%INTVAL2BASE Convert a positive integer value, which is a numeric type
%          such as a double or a uint64, into a different base. The digits
%          of the representation in that base are provided in an array.
%          Numbers are returned, not characters.
%
%INPUTS: n The positive integer value to be converted. This can be any
%          standard integer and floating point type.
%     base The positive base to which the value should be converted.
% firstMostSig If this is true, then the first digit is the most
%          significant. Otherwise, the first digit is the least
%          significant. the default if omitted or an empty matrix is passed
%          is false.
% pad2Digits If the number of digits needed to represent n in the selected
%          base is less than pad2Digits, then the value is padded with
%          zeros to this length. The zeros cover the most significant
%          values and thus their location depends on firstMostSign. If
%          omitted or an empty matrix is passed, pad2Digits=0, which means
%          that nothing is padded. Note that when unpadded, the value 0
%          will take 1 digit.
% baseVals These are powers of base, starting with the zeroth power, which
%          is 1. These are output by the function and can be passed to
%          speed things up if this function is repeatedly called. These
%          should be a uint64 data type. One can just pass back the value
%          returned by this function. If omitted, the values are computed
%          by the function.
%
%OUTPUTS: digits A numDigitsX1 vector containing the representation of n in
%                the given base. This will be cast to the same type as the
%                input.
%       baseVals Powers of base as uint64 values. These are just returned
%                so that they can be passed again on subsequent calls to
%                the function.
%
%The algorithm is similar to that used in the index2NumDim and unrankTuple
%functions. It just uses integer division and subtraction to extract each
%digit.
%
%EXAMPLE:
%Here, we convert the value 7908 to base 10 (thus getting the individual
%digits) and to binary with the first digit the most significant.
% val=7908;
% base10Vals=intVal2Base(val,10,true)'
% base2Vals=intVal2Base(val,2,true)'
%One will get the vector [7,9,0,8] in base 10 and the vector 
%[1,1,1,1,0,1,1,1,0,0,1,0 0] in base 2.
%
%January 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(firstMostSig))
   firstMostSig=false; 
end

if(nargin<4||isempty(pad2Digits))
    pad2Digits=0;
end

nType=class(n);
if(isinteger(n))
    maxInt=intmax(nType);%The value is 2^64-1 if a uint64 is passed.
else%Assume is a floating point value.
    maxInt=flintmax(nType);
end

%If base 2 and we are using 64 bit data types, then this is the most
%possible digits in the converted result.
maxDigits=64;

if(nargin<5||isempty(baseVals))
    baseVals=zeros(maxDigits,1,'uint64');
    baseVals(1)=1;
    for k=2:maxDigits
        baseVals(k)=baseVals(k-1)*base;

        if(baseVals(k)>=maxInt)
            %When it overflows, we have the maximum number of possible
            %digits for this base. In Matlab, integers don't overflow;
            %they just clip to maxInt.
            maxDigits=k-1;
            break;
        end
    end
    baseVals=baseVals(1:maxDigits);
else
    maxDigits=length(baseVals);
end

if(n>maxInt)
    error('This function only supports values of n that can be exactly represented as integers.')
elseif(n==0)%The zero special case
    totalDigits=max(pad2Digits,1);
    digits=zeros(totalDigits,1,nType);
else
    n=uint64(n);
    digits=zeros(maxDigits,1,'uint64');
    foundFirstDigit=false;
    for k=maxDigits:-1:1
        %The number of complete multiples of the base. Note that this is
        %integer division.
        wholeVal=idivide(n,baseVals(k),'fix');
        if(wholeVal~=0)
            if(foundFirstDigit==false)
                foundFirstDigit=true;
                totalDigits=k;
            end
            
            if(firstMostSig)
                digits(totalDigits-k+1)=wholeVal;
            else
                digits(k)=wholeVal;
            end
            n=n-wholeVal*baseVals(k);
            
            if(n==0)
                break;
            end
        end
    end
    
    if(pad2Digits>totalDigits&&firstMostSig)
        %If the first digit is the most significant and the value is padded
        %with zeros, then we need to shift all fo the digits to the end of
        %the array.
        digits=circshift(digits(1:pad2Digits),pad2Digits-totalDigits);
    else
        totalDigits=max(totalDigits,pad2Digits);
        digits=digits(1:totalDigits);
    end
end

%Cast the output to the same type as the input.
digits=cast(digits,nType);

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
