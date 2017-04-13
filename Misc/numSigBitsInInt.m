function numBits=numSigBitsInInt(intVal)
%%NUMSIGBITSININT Given a positive integer value, return the number of
%                 significant bits in the integer. This is the total number
%                 of bits needed to represent the integer in binary without
%                 leading zeros. For example, numSigBitsInInt(1)=1,
%                 numSigBitsInInt(0)=0 and numSigBitsInInt(8)=4. The
%                 function is presumably only meaningful for positive
%                 integers as negative integers are represented using the
%                 two's complement and thus
%
%INPUTS: intVal   A positve, scalar integer value.
%
%OUTPUTS: numBits The number of bits needed to represent the integer.
%                 Negative inputs just return zero. Non-integer inputs
%                 return an error.
%
%The function just counts the number of right-shifts of the bits in the
%integer that are needed until the result is zero.
%
%For intVal>0, this result is the same as what one would expect from
%evaluating fix(log2(x))+1. However, one would expect a bitshifting
%implementation that works on integers to be faster than an implementation
%that evaluates a logarithm.
%
%December 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    numBits=0;
    while(intVal>0)
        intVal=bitshift(intVal,-1);
        numBits=numBits+1;
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
