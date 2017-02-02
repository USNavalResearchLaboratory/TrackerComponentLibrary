function checksum=NMEAChecksum(NMEAString)
%%NMEACHECKSUMVALID  Compute the checksum of a National Maritime
%                    Electronics Association (NMEA) 0183 sentence. This is
%                    the message standard used for Automatic Identification
%                    System (AIS) messages sent by transponders on ships.
%                    This function is useful for validating checksums to
%                    determine whether to discard a message as corrupt.
%
%INPUTS: NMEAString A NMEA 0183 string as a sequence of characters. The
%                   string should start with a ! or $ and everything until
%                   a * goes into the checksum. In such a message, the
%                   checksum is placed after the *. Additional data might
%                   come after that. 
%
%OUTPUTS: checksum A two-character checksum associated with the NMEA
%                  string. Everything after the * is ignored when
%                  this function computes the checksum. If the message
%                  does not have characters surrounded by a ! or $ and
%                  a *, then an empty matrix is returned. The function does
%                  no further checks to make sure that the message is in
%                  the correct format.
%
%The NMEA 0183 standard is expensive and can be purchased from 
%http://www.nmea.org/content/nmea_standards/nmea_0183_v_410.asp
%However, a number of web pages say that the checksum is just the
%hexadecimal representation of xoring the ASCII codes of the characters
%between the ! or $ and the * in the message. 
%
%As an example, calling the function as
%NMEAChecksum('!AIVDM,1,1,,B,177KQJ5000G?tO`K>RA1wUbN0TKH,0*5C')
%where the message is an AIS message would return '5C", indicating that the
%checksum embedded in the message is correct. Alternatively, if one were
%trying to create a checksum to be appended to the message, one would call
%NMEAChecksum('!AIVDM,1,1,,B,177KQJ5000G?tO`K>RA1wUbN0TKH,0*')
%and would get the same result as everything after the * is ignored.
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%If it has an obviously wrong format.
if(isempty(NMEAString)||~ischar(NMEAString)||(~isrow(NMEAString)&&~iscolumn(NMEAString))||(NMEAString(1)~='$'&&NMEAString(1)~='!')||length(NMEAString)<3)
    checksum=[];
    return;
end

%Find the first * in the string.
maxIdx=2;
while(NMEAString(maxIdx)~='*')
    maxIdx=maxIdx+1;
    
    if(maxIdx>length(NMEAString))
        checksum=[];
        return;
    end
end
maxIdx=maxIdx-1;

%The bitxor function will only work with raw UTF8/ASCII codes.
NMEAString=unicode2native(NMEAString);

checksum=0;
for curIdx=2:maxIdx
    checksum=bitxor(checksum,NMEAString(curIdx));
end
checksum=dec2hex(checksum,2);

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
