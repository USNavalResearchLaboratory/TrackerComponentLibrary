function charString=integer2CustomCharString(intVal,alphabet)
%%INTEGER2CUSTOMCHARSTRING Given an integer value, convert it into a
%               string. However, rather than using the digits 0-9 for each
%               place, use the characters given in alphabet. The first
%               character refers to the 0's place, the next to the 1's and
%               so on. The number of characters determines the base. For
%               example, for base 2, one could use alphabet='01' or
%               alphabet='ab' if one wanted to use letters instead of
%               numbers.
%
%INPUTS: intVal A positive integer value that one wishes to convert into a
%               custom string. intVal cannot be over the maximum integer
%               representable by a double floating point value without a
%               loss of precision (which is 2^53).
%      alphabet The characters one wishes to use to express the converted
%               integer. The number of characters equals the base of the
%               conversion. For example, to express intVal as a hexadecimal
%               number (base 16), one would use
%               alphabet='0123456789abcdef'. alphabet can be an array of
%               characters or an array of numbers if one wishes to have
%               some sort of a more general conversion.
%
%OUTPUTS: charString  A 1XnumChars character string of intVal expressed in
%                     the given bases/ alphabet. The first element is the
%                     most significant digit. charString is of type char if
%                     alphabet is of type char. otherwise, charString is of
%                     type double.
%
%The conversion is done using the modulo operation with the given base to
%determine which indexed value is in each digit of the solution. Logarithms
%as used to determine how many digits are in the solution.
%
%EXAMPLE 1:
%We want to get the binary representation of 12345 as a string of
%characters
% charString=integer2CustomCharString(12345,'01')
%One gets charString='11000000111001'.
%
%EXAMPLE 2:
%Here, we want to get the binary represetation of 12345 with each element of the output being a (numeric) 0 or 1.
% charString=integer2CustomCharString(12345,[0 1])
%Now one gets charString=[1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1].
%
%EXAMPLE 3:
%One wishes to convert 12345 into hexadecimal:
% charString=integer2CustomCharString(12345,'0123456789abcdef')
%The result is charString='3039'.
%Just to be clear, the same conversion with 15 instead of 12345 would have
%returned 'f'.
%
%March 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

base=length(alphabet);

%Determine the number of digits needed to represent the number.
numDigits=fix(log(intVal)/log(base))+1;

%If intVal==0.
if(numDigits<1)
    charString=alphabet(1);
    return;
end

%Allocate space for the number by creating a vector full of zero digits.
if(ischar(alphabet))
    charString=repmat(alphabet(1),1,numDigits);
else
    charString=alphabet(1)*ones(1,numDigits);
end

for curIdx=1:numDigits
    idx=rem(intVal,base);
    charString(numDigits-curIdx+1)=alphabet(idx+1);
    intVal=fix(intVal/base);
    
    if(intVal==0)
        break;
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
