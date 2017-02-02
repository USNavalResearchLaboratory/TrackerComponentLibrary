function [hashVal,hexVal]=genHashVal(inputArray,algorithm)
%%GENHASHVAL  Get a hash value of a character string or a 1D array of
%             bytes. That is, integer values from 0 through 255. The hash
%             value can either be a simple integer value, useful for quick
%             comparison of strings, or one can get a cryptographic hash as
%             an an array of integers from 0 to 255 using various
%             algorithms. Cryptographic hashes are also returned has
%             hexadecimal strings.
%
%INPUTS: inputArray A 1-dimensional string of characters or array of
%                   integer values from 0 to 255 for which one desired a
%                   hash code.
%        algorithm  An optional  string specifying the algorithm to use to
%                   generate the hash. Possible values are
%                   'simple' (the default) Use a simple, non-cryptographic
%                            hash function to get a single integer hash for
%                            the input.
%                   'MD5'    Use the MD5 algorithm to get a cryptographic
%                            hash, along with its hexidecimal string.
%                   'SHA-1'  Use the SHA-1 algorithm to get a cryptographic
%                            hash along with its hexidecimal string.
%                 'SHA-256'  Use the SHA-256 algorithm to get a
%                            cryptographic hash along with its hexidecimal
%                            string.
%
%OUTPUTS: hashVal  If algorithmis 'simple', then this is an integer hash
%                  code for the given input. Otherwise, it is an array of
%                  bytes representign the has value.
%         hexVal   A hexidecimal string representing the hash value. 
%
%The simple algorithm is executed using the hashCode method in Java's
%java.lang.String class. The cryptographic hashes are executed using the
%java.security.MessageDigest class. In Matlab, one can directly access Java
%functions, which is what was done to implement this function.
%
%April 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(~ischar(inputArray)&&any(inputArray~=fix(inputArray))||any(inputArray>255)||any(inputArray<0))
    error('Hash values can only be obtained from character arrays or arrays of positive integer values from 0 through 255');
end

if(nargin<2)
    algorithm='simple';
end

%If the simple aglorithm is chosen. otherwise, a cryptographic has is being
%used.
if(strcmp(algorithm,'simple'))
    %Convert the array to a character string.
    if(~ischar(inputArray))
        inputArray=native2unicode(inputArray);
    end

    myString=java.lang.String(inputArray);
    hashVal=myString.hashCode;
    %Java's integers are always 32-bits.
    hexVal=dec2hex(typecast(int32(hashVal),'uint32'),8);
    return;
end

switch(algorithm)
    case 'MD5'
    case 'SHA-1'
    case 'SHA-256'
    otherwise
        error('An unsupported algorithm was given')
end

%The cryptographic has functions work on bytes, not strings.
if(ischar(inputArray))
    inputArray=unicode2native(inputArray);
end

myDigest=java.security.MessageDigest.getInstance(algorithm);

myDigest.update(inputArray);
hashVal=myDigest.digest;
hexVal=javax.xml.bind.DatatypeConverter.printHexBinary(hashVal);
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
