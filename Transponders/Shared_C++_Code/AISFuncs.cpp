/**AISFUNCS  C++functions that are helpful in parsing NMEA messages for AIS
 *           data. This only deals with message parsing that is not
 *           already handled by functions in AISLib.
 *
 *November 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
 **/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#include "AISFuncs.hpp"

#include <limits>
#include <stdexcept>

void separateMessageAndTimestamp(const std::string &fullMessage,std::string &messagePart,std::string &endPart) {
/**SEPARATEMESSAGEANDTIMESTAMP Given an NMEA string as a string data type,
 *              which might have additional receiver data appended to the
 *              end, separate out the message part from anything that might
 *              come after the message, (anything coming after the message
 *              must be separated from the message by a comma), which often
 *              includes identifiers and timestamps appended by the
 *              receiver. If fullMessage appears to be invalid, then
 *              messagePart and endPart are set to empty strings. If there
 *              is no end part, then endPart is set to an empty string.
 *
 *November 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.*/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/    
        
    //We want to separate everything that comes after the checksum (if
    //anything comes after the checksum). 
    size_t astPos=fullMessage.find_first_of("*",0);
    //The asterisk will only appear in the checksum (and maybe some
    //things after it). Anything that comes after the checksum will
    //come after a comma that is the third character after the
    //asterisk.

    //std::string::npos is the return value if no asterisk is
    //found, in which case there is no checksum and the string
    //should not be valid. Also, if there are not two more
    //characters after the asterisk, then the checksum cannot be
    //valid.
    if(astPos==std::string::npos||astPos+2>fullMessage.length()-1) {
        messagePart=std::string();
        endPart=std::string();
        return;
    }

    //If there are at least three characters after the asterisk and the
    //third one is a comma, then we have to remove that
    //part and process it separately to look for a timestamp or
    //anything else. Otherwise, if there is nothing after the checksum,
    //then set the return parameters accordingly.
    if(astPos+3<=fullMessage.length()-1) {
        //Any other character after the checksum is invalid.
        if(fullMessage[astPos+3]!=',') {
            messagePart=std::string();
            endPart=std::string();
            return;
        }

        //make curString just the part through the end of the checksum.
        messagePart=fullMessage.substr(0,astPos+3);

        //If there is at least one character after the comma, then
        //that is the end part toreturn. Otherwise, there is no end part.
        if(astPos+3<=fullMessage.length()-1) {
            endPart=fullMessage.substr(astPos+4);
        } else {
            endPart=std::string();
        }
    } else {//There is nothing after the checksum.
        messagePart=fullMessage;
        endPart=std::string();
    }
}

double getEndTimestamp(const std::string &theString) {
/**SEPARATEMESSAGEANDTIMESTAMP A string may or may not end with a positive
 *              integer timestamp. The timestamp will either be the entire
 *              message or it will be the final thing at the end of a
 *              message after a comma. This function extracts the
 *              timestamp and returns it as a double, or it returns a NaN
 *              if it cannot extract a timestamp.
 *
 *November 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.*/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/
    
    double timeStamp=std::numeric_limits<double>::quiet_NaN();
    size_t pos=theString.find_last_of(',',std::string::npos);
    
    //If there is no comma in the string, then see if the whole string is
    //the timestamp. Otherwise, extract whatever is after the comma, if
    //anything.
    if(pos==std::string::npos) {
        try {
            timeStamp=static_cast<double>(std::stoul(theString,NULL));
        } catch(const std::invalid_argument& ia) {}
    } else if(pos<=theString.length()-1) {
        try {
            timeStamp=static_cast<double>(std::stoul(theString.substr(pos+1),NULL));
        } catch(const std::invalid_argument& ia) {}
    }
    return timeStamp;
}

/*LICENSE:
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
%OF RECIPIENT IN THE USE OF THE SOFTWARE.*/
