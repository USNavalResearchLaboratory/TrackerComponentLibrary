function [charString,hashValue]=genUUID(UUIDType,data)
%%GENUUID Generate a universally unique identifier (UUID), also known as a
%         globally unique identifier. A UUID is a 128-bit sequence that is
%         generated in such a manner that is will almost certainly be
%         unique across space and time. Thus, entities can create and
%         exchange UUIDs without cantral planning.
%
%INPUTS: UUIDType An optional parameter specifying the UUID type. Possible
%                 values are:
%                 0 (The default if omitted) Create a random UUID. This is
%                   the best if one is creating numerous IDs over time.
%                 1 Create a name-based ID. The second parameter data must
%                   provide the name either as a character string, like
%                   'NAME', or as a vector of integer values representing
%                   the name, like [78;65;77;69], which are the numeric
%                   values corresponding to the characters in 'NAME'. ID's
%                   created with the name name ar identical.
%                 2 Get a UUID given the character string for a UUID given
%                   in data. This will make the output value charString
%                   have the same value as data, but the hashValue will
%                   also be obtained.
%            data This parameter is only needed if UUIDType is not zero. If
%                 UUIDType is one, then this is either the character string
%                 or integer number array of the name. If UUIDType is 2,
%                 then this is the character string of the UUID value.
%
%OUTPUTS:charString A character string representing the UUID. The cahacter
%                  string is of the format <time_low> "-" <time_mid> "-"
%                  <time_high_and_version> "-"  <variant_and_sequence> "-"
%                  <node> For example, f2fdf31c-6ceb-4a67-a888-697dbffe8bd7
%                  is a valid character string. The character string can be
%                  compared to other UUIDs to see if they are the same.
%        hashValue The integer hash code of the UUID.
%
%UUIDs are standardized in the Internet Engineering Task Force's (IETF)
%request for comments (RFC) number 4122, which can be found online at
%https://tools.ietf.org/html/rfc4122
%Code for implementing version 4 and version 5 UUIDs is present in Java's
%util library. In Matlab, one can directly access Java functions, which is
%what was done to implement this function.
%
%April 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<1)
    UUIDType=0;
end

switch(UUIDType)
    case 0%It is a random UUID
        theID=java.util.UUID.randomUUID();
    case 1%It is a name-based UUID.
        %If the name-based UUID is coming from a character string
        if(isa(data,'char'))
            theID=java.util.UUID.nameUUIDFromBytes(unicode2native(data));
        else%Otherwise assume that it is coming from a matrix of numbers.
            theID=java.util.UUID.nameUUIDFromBytes(data);
        end
    case 2%The character string for the UUID is provided.
        theID=java.util.UUID.fromString(data);
    otherwise
        error('Invalid UUID type is provided')
end

charString=char(theID.toString());
hashValue=theID.hashCode();
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
