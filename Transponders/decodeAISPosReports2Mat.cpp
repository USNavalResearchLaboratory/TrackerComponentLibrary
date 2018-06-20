/*DECODEAISPOSREPORTS2MAT Given a matrix of Automatic Identification
 *                 System (AIS) messages (or a cell array of AIS messages)
 *                 given as National Marine Electronics Association (NMEA)
 *                 formatted ASCII text messages, extract position reports
 *                 of message types 1, 2, 3, 18, 19, and 27, storing the
 *                 information related to ship location, identity, course,
 *                 and time in different rows of a matrix, where each
 *                 column corresponds to a different report. If a cell
 *                 array is passed, then the matrices for decoding each
 *                 cell array are placed into corresponding cells on
 *                 return. It is also possible to decode timestamps (or
 *                 another integer number) appended as a final field on the
 *                  messages. The manipulation of matrices rather than
 *                 structures is significantly faster in Matlab, which is
 *                 why one might prefer this function to decodeAISString
 *                 when only simple position reports are needed.
 * 
 *INPUTS: NMEAStrings A character string containing one or more Automatic
 *                 Identification System (AIS) messages given as National
 *                 Marine Electronics Association (NMEA) formatted ASCII
 *                 text. Each message is on its own line. Invalid messages
 *                 and messages that are not of types 1, 2, 3, 18, 19, and
 *                 27 are skipped.
 * convert2Metric  If true, all distance units in the result are converted
 *                 to meters, all speed units are converted to meters per
 *                 second, and all angles are converted to radians North of
 *                 East (the mathematical definition) and angular rates are
 *                 converted to radians per second North of East, all
 *                 rather than using typical nautical units. The default if
 *                 omitted is true.
 * 
 *decodedMsgs A 14XN matrix of the N successfully decoded position reports.
 *            The 14 fields that are decoded are stored in the matrix as:
 *                     1) The maritime mobile service identity (MMSI)
 *                        number of the transmitter that can be used to
 *                        identify targets (A 30-bit integer value).
 *                     2) Navigation status. An integer indicating what the
 *                        ship is doing/ whether it is moving. Possible
 *                        valid values are:
 *                        *0 = under way using engine
 *                        *1 = at anchor
 *                        *2 = not under command,
 *                        *3 = restricted maneuverability
 *                        *4 = constrained by her draught
 *                        *5 = moored
 *                        *6 = aground
 *                        *7 = engaged in fishing
 *                        *8 = under way sailing
 *                        *9 = reserved for future amendment of
 *                          navigational status for ships carrying DG, HS,
 *                          or MP, or IMO hazard or pollutant category C,
 *                          high speed craft (HSC).
 *                        *10 = reserved for future amendment of
 *                          navigational status for ships carrying
 *                          dangerous goods (DG), harmful substances (HS)
 *                          or marine pollutants (MP), or IMO hazard or
 *                          pollutant category A, wing in ground (WIG).
 *                        *11 = power- driven vessel towing astern
 *                              (regional use)
 *                        *12 = power-driven vessel pushing ahead or towing
 *                              alongside (regional use)
 *                        *13 = reserved for future use,
 *                        *14 = AIS-SART (active), MOB-AIS, EPIRB-AIS
 *                     3) Turn rate. This is limited to +/- 708 degrees per
 *                        minute (about +/- .2059 radians per second when
 *                        using metric units).
 *                     4) Speed over ground in knots up to 102.2 knots, or
 *                        52.6 meters per second when using metric units.
 *                     5) Position accuracy flag. Possible values are:
 *                        1=high (accuracy beter than 10m)
 *                        0=low (accuracy worse than 10 meters)
 *                     6) Latitude (WGS-84) in degrees North or radians
 *                        North when using the metric units option.
 *                     7) Longitude (WGS 84) in degrees East or radians
 *                        East when using the metric units option.
 *                     8) Course over ground (heading) in degrees East of
 *                        North or radians North of East when using the
 *                        metric units option.
 *                     9) True heading in degrees East of North or radians
 *                        North of East when using the metric units option.
 *                    10) The second number within the universal
 *                        coordinated time (UTC) minute when the ship's
 *                        transmitter sent the position report. NOTE that
 *                        if this is 61 or greater, than this means that the 
 *                        the positioning system is broken or in dead
 *                        reckoning mode, so the return is likely to be
 *                        inaccurate. Also, the standard does not properly
 *                        account for 61 second minutes when a leapsecond
 *                        is added, so it is not clear what will happen
 *                        with leapseconds. This field is not used by
 *                        message 27. Rather, there is no timestamp and a
 *                        position latency indicator is used to indicate
 *                        the extend of the delay of the true time of the
 *                        position from the broadcast time.
 *                    11) Special manoeuvre indicator. This can be
 *                        1 = not engaged in special manoeuvre
 *                        2 = engaged in special manoeuvre
 *                    12) The ID of the message that provided the above
 *                        information. This can be
 *                        1  Class A position report, scheduled.
 *                        2  Class A position report, assigned.
 *                        3  Class A position report, when interrogated.
 *                        18 Standard class B equipment position report.
 *                        19 Extended class B equipment position report.
 *                        27 Long-range automatic identification system
 *                           broadcast message.
 *                    13) A number indicating how many times the message
 *                        has been repeated. This ranged from 0 to 3. 3
 *                        also indicates that the message should not be
 *                        repeated anymore.
 *                    14) Position latency. This is only used with message
 *                        27. This is 0 or 1. 0 means the position latency
 *                        is less than five seconds, one means it is more.
 *          epochTimes If this parameter is requested, then the final
 *                     numbers after the checklsum (separated by a comma)
 *                     are assumed to either be an integer epoch time, or
 *                     anything with a comma separating an integer epoch
 *                     time at the end. Thus, this tries to decode the
 *                     epoch time associated with the message. Note that in
 *                     multi-part messages (message 19 can have multiple
 *                     parts), the epoch time used is that of the LAST
 *                     received part in the sequence. If no time can be
 *                     decoded, then
 *
 *This function relies on libais from 
 *https://github.com/schwehr/libais
 *The library probably does not work on big-endian processors (Intel
 *processors are little-endian).
 *
 *The algorithm can be compiled for use in Matlab using the 
 *CompileCLibraries function.
 *
 *The algorithm is run in Matlab using the command format
 *decodedMsgs=decodeAISPosReports2Mat(NMEAStrings,convert2Metric);
 *or
 *[decodedMsgs,epochTimes]=decodeAISPosReports2Mat(NMEAStrings,convert2Metric);
 *if timestamps are available
 *
 *Example 1:, Decoding a single valid ID=1 string.
 * NMEAStrings='!AIVDM,1,1,,A,13HOI:0P0000VOHLCnHQKwvL05Ip,0*23';
 * decodedMsgs=decodeAISPosReports2Mat(NMEAStrings,false)
 *The decoded message should be:
 * decodedMsgs=[227006760;%MMSI
 *                      0;%Navigation status (underway)
 *                    NaN;%Not used
 *                      0;%Turn rate
 *                      0;%Speed over ground
 *     49.475578308105469;%North latitude
 *      0.131380006670952;%East longitude
 *     36.700000762939453;%Heading East of North
 *                    NaN;%Not used
 *                     14;%Seconds within minute when message sent.
 *                    NaN;%Not used.
 *                      1;%Class A position report, scheduled.
 *                      0;%has not ben rebroadcast.
 *                    NaN];%Not used.
 *
 *Example 2: Decoding multiple messages, one valid string of each type of
 *message. It is most likely that this would arise from reading the
 *messages in from a file, because adding newline characters to strings has
 *to be done using the \n newline character to separate the string in
 *sprintf rather than just typing them. The AIS19 string consists of two
 *parts that are joined.
 * NMEAStrings=sprintf('!AIVDM,1,1,,A,13HOI:0P0000VOHLCnHQKwvL05Ip,0*23\n!AIVDM,1,1,,B,284;UGTdP4301>3L;B@Wk3TnU@A1,0*7C\n!AIVDM,1,1,,B,34hoV<5000Jw95`GWokbFTuf0000,0*6C\n!SAVDM,1,1,4,B,B5NU=J000=l0BD6l590EkwuUoP06,0*61\n!AIVDM,2,1,7,B,C5NMbDQl0NNJC7VNuC<v`7NF4T28V@2g0J6F::00000,0*59\n!AIVDM,2,2,7,B,0J70<RRS0,0*30\n!AIVDM,1,1,,B,K815>P8=5EikdUet,0*6B');
 * decodedMsgs=decodeAISPosReports2Mat(NMEAStrings,false)
 *One will observe that 6 messages are correctly decoded, of types 1, 2, 3,
 *18, 19, and 27.
 *
 *Example 3: Decoding multiple messages, some of them of the wrong type.
 *When messages of the wrong type are used, then they are askipped. Here,
 *two correct messages are read in as well as one incorrect message in
 *between them. The "incorrect" message is of type 26 (not a position report) and is thus ignored.
 *Below are four strings, but only two of them contain AIS position reports
 * NMEAStrings=sprintf('!AIVDM,1,1,,A,13HOI:0P0000VOHLCnHQKwvL05Ip,0*23\n!AIVDM,1,1,,A,H44cj<0DdvlHhuB222222222220,2*46\n!AIVDM,1,1,,B,284;UGTdP4301>3L;B@Wk3TnU@A1,0*7C');
 * decodedMsgs=decodeAISPosReports2Mat(NMEAStrings)
 *The results are in metric units.
 *
 *Example 4: Consider the case where the messages contain some numbers
 *after the checksum, which make up the timestamps of the messages.
 *This is the same as the previous example, but with fake timestamps
 *appended.
 * NMEAStrings=sprintf('!AIVDM,1,1,,A,13HOI:0P0000VOHLCnHQKwvL05Ip,0*23,1241827197\n!AIVDM,1,1,,A,H44cj<0DdvlHhuB222222222220,2*46,1241827198\n!AIVDM,1,1,,B,284;UGTdP4301>3L;B@Wk3TnU@A1,0*7C,1241827199');
 * [decodedMsgs,epochTimes]=decodeAISPosReports2Mat(NMEAStrings)
 *The results are in metric units.
 *
 *Example 5: Some receivers will append aditional information before the
 *timestamp. As long as the timestamp is the last thing in the message,
 *then we can still recover it.
 * NMEAStrings='!AIVDM,1,1,,A,100WhdhP0nJRdiFFHFvm??v00L12,0*13,raishub,1342569600';
 * [decodedMsgs,epochTimes]=decodeAISPosReports2Mat(NMEAStrings)
 *The results are in metric units.
 *
 *November 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.*/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

//This definition is necessary or else the aislib files will cause
//assertions when bad messages are passed.
#ifndef NDEBUG
#define NDEBUG
#endif

//This header is required by Matlab.
#include "mex.h"
//Needed for mxIsClass
#include "matrix.h"

//The AISLib headers.
#include "ais.h"
#include "vdm.h"

//Functions for parsing some parts of the strings that is not handled by
//AISLib.
#include "AISFuncs.hpp"

//Needed to get a NaN value.
#include <limits>

#include "MexValidation.h"
//For fabs
#include <cmath>

//This is for parsing the input into strings that can be provided to the
//VdmStream class to reconstruct the messages.
#include <iostream>
#include <sstream>
#include <string>
//For count
#include <algorithm>

//Multiple a value in degrees by deg2Rad to get its value in radians.
const double deg2Rad=0.0174532925199432957692369076849;
const double halfPi=1.57079632679489661923132169164;
const double min2Sec=60;
//Conversion of knots to meters per hour.
const double knot2Metph=1852;
const double hour2Sec=60*60;

//Prototypes
void decodeCharacterStrings(const char*msgChars,bool hasEndTimeStamps,bool convert2Metric,mxArray *&retMat, mxArray *&timeMat);
bool extractSinglePosReport(std::unique_ptr<libais::AisMsg> &curAisMsg,double *retData);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) { 
    bool hasEndTimestamps=false;
    bool convert2Metric=true;
    
    if(nrhs!=1&&nrhs>2) {
        mexErrMsgTxt("Wrong number of inputs.");
        return;
    }

    if(nlhs>1) {
        hasEndTimestamps=true;
    }
    
    if(nrhs>1) {
        convert2Metric=getBoolFromMatlab(prhs[1]);
    }
            
    if(nlhs>2) {
        mexErrMsgTxt("Wrong number of outputs.");
        return;
    }
    
    if(mxGetN(prhs[0])!=1&&mxGetM(prhs[0])!=1) {
        mexErrMsgTxt("The input has an invalid dimensionality.");
        return;
    }
    
    if(mxIsChar(prhs[0])) {//If a character string is passed.
        mxArray *retMat, *timeMat;
        char *theData;
        
        theData=mxArrayToString(prhs[0]);
        decodeCharacterStrings(theData,hasEndTimestamps,convert2Metric,retMat,timeMat);
        mxFree(theData);
        
        plhs[0]=retMat;
        if(nlhs>1) {
            plhs[1]=timeMat;
        }
    } else {
        mexErrMsgTxt("Invalid data type passed.");
    }
}

void decodeCharacterStrings(const char*msgChars,bool hasEndTimeStamps,bool convert2Metric,mxArray *&retMat, mxArray *&timeMat) {
/*DECODECHARACTERSTRINGS Given a character string, which might contain
 *       multiple NMEA messages, each on a different line, decode them into
 *       a matrix consisting of 14 rows (as described for the return value
 *       at the top of this file).
 *
 *November 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.*/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/
    
    libais::VdmStream decodedAISData;//To decode a set of messages at once.
    
    double *retData, *timeData=nullptr;
    const size_t numRowsInOutput=14;
    size_t maxNumMessages;
    size_t numMsgDecoded;
    std::string allMsg(msgChars);
    std::string curString;
    std::istringstream ss(allMsg);//To parse a bunch of AIS strings
    
    //The maximum possible number of messages to decode.
    maxNumMessages=std::count(allMsg.begin(), allMsg.end(), '\n')+1;
    
    //Allocate space for the maximum possible number of items to return
    retMat=mxCreateDoubleMatrix(numRowsInOutput,maxNumMessages,mxREAL);
    retData=reinterpret_cast<double*>(mxGetData(retMat));
    
    if(hasEndTimeStamps) {
        timeMat=mxCreateDoubleMatrix(maxNumMessages,1,mxREAL);
        timeData=reinterpret_cast<double*>(mxGetData(timeMat));
    }

    //We will now feed all of the lines into the VdmStream object one
    //at a time.
    numMsgDecoded=0;
    while(std::getline(ss, curString,'\n')) {
        //This marks no time stamp as available.
        double timeStamp=std::numeric_limits<double>::quiet_NaN();
        bool pushSucceeded;
        size_t numBeforeAdd=decodedAISData.size();
        std::string messagePart, endPart;
                
        //We want to separate everything that comes after the checksum (if
        //anything comes after the checksum).
        separateMessageAndTimestamp(curString,messagePart,endPart);
        
        //If the message is invalid, then skip it.
        if(messagePart.length()==0) {
            continue;
        }
        
        //Try to extract a timestamp from the end of the end part.
        timeStamp=getEndTimestamp(endPart);

        pushSucceeded=decodedAISData.AddLine(messagePart);

        if(pushSucceeded) {
            //If a complete message was decoded, then try to extact the
            //data and add it.
            if(decodedAISData.size()!=numBeforeAdd) {
                std::unique_ptr<libais::AisMsg> curAisMsg=decodedAISData.PopOldestMessage();
                const size_t outputOffset=numRowsInOutput*numMsgDecoded;
                bool decodingSucceeded;

                decodingSucceeded=extractSinglePosReport(curAisMsg,retData+outputOffset);

                //Get rid of the popped message.
                curAisMsg=nullptr;
                
                if(decodingSucceeded) {
                    if(hasEndTimeStamps) {
                        timeData[numMsgDecoded]=timeStamp;
                    }
                    
                    if(convert2Metric) {
                        //Convert turn rate from degrees per minute East of North to
                        //radians per second North of East
                        retData[outputOffset+2]=-deg2Rad*retData[outputOffset+2]/min2Sec;                    

                        //Convert speed over ground from knots to meters per second.
                        retData[outputOffset+3]=retData[outputOffset+3]*knot2Metph/hour2Sec;
                        //Convert latitude to radians
                        retData[outputOffset+5]=deg2Rad*retData[outputOffset+5];

                        //Convert longitude to radians
                        retData[outputOffset+6]=deg2Rad*retData[outputOffset+6];

                        //Convert course over ground from degrees East of North to
                        //radians North of East
                        retData[outputOffset+7]=halfPi-deg2Rad*retData[outputOffset+7];
                        //Convert true heading from degrees East of North to radians
                        //North of East
                        retData[outputOffset+8]=halfPi-deg2Rad*retData[outputOffset+8];
                    }

                    numMsgDecoded++;
                }
            }
        }
    }
    
    //Resize the matrices to match the actual number of things decoded.
    mxSetN(retMat,numMsgDecoded);
    if(hasEndTimeStamps) {
        mxSetM(timeMat,numMsgDecoded);
    }
}

bool extractSinglePosReport(std::unique_ptr<libais::AisMsg> &curAisMsg,double *retData) {
/**EXTRACTSINGLEMESSAGEDATA Extract the data within a single AIS message
 *  and put it into the array retData. THe formatting is described at the
 *  top of the file. If docoding was a success, then true is returned;
 *  otherwise, false is returned. Failure will be due to the message being
 *  the wrong type.
 *
 *November 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.*/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

    const int msgId=curAisMsg->message_id;

    switch(msgId) {
        case 1://Class A position report, scheduled.
        case 2://Class A position report, assigned.
        case 3://Class A position report, when interrogated.
        {
            libais::Ais1_2_3 *msg=reinterpret_cast<libais::Ais1_2_3*>(curAisMsg.get());

            //MMSI number (An integer).
            retData[0]=static_cast<double>(msg->mmsi);
            //Navigational status
            retData[1]=static_cast<double>(msg->nav_status);

            //Rate of turn in degrees per minute
            if(!msg->rot_over_range) {
                retData[2]=msg->rot;
            } else {
                retData[2]=std::numeric_limits<double>::quiet_NaN();
            }

            //Speed over ground in knots.
            if(msg->sog<=102.2) {
                retData[3]=msg->sog;
            } else {
                retData[3]=std::numeric_limits<double>::quiet_NaN();
            }

            //Position accuracy (0 or 1).
            retData[4]=static_cast<double>(msg->position_accuracy);

            //Latitude (degrees)
            if(fabs(msg->position.lat_deg)<=90) {
                retData[5]=msg->position.lat_deg;
            } else {
                retData[5]=std::numeric_limits<double>::quiet_NaN();
            }

            //Longitude (degrees)
            if(fabs(msg->position.lng_deg)<=180) {
                retData[6]=msg->position.lng_deg;
            } else {
                retData[6]=std::numeric_limits<double>::quiet_NaN();
            }

            //Course over ground (degrees East of North)
            if(msg->cog<360) {
                retData[7]=msg->cog;
            } else {
                retData[7]=std::numeric_limits<double>::quiet_NaN();
            }

            //True heading (degrees East of North)
            if(msg->true_heading<360) {
                retData[8]=static_cast<double>(msg->true_heading);
            } else {
                retData[8]=std::numeric_limits<double>::quiet_NaN();
            }

            //UTC second of minute.
            //The restriction against it being 60 means that an
            //incorrect UTC second will be reported when a
            //leapsecond is added.
            if(msg->timestamp<60) {
                retData[9]=static_cast<double>(msg->timestamp);
            } else {
                retData[9]=std::numeric_limits<double>::quiet_NaN();
            }

            //Special manoeuvre indicator
            if(msg->special_manoeuvre!=0&&msg->special_manoeuvre<3) {
                retData[10]=msg->special_manoeuvre;
            } else {
                retData[10]=std::numeric_limits<double>::quiet_NaN();
            }

            //The message ID.
            retData[11]=static_cast<double>(msg->message_id);

            //The message repeat indicator
            retData[12]=static_cast<double>(msg->repeat_indicator);
            //There is no position latency indicator
            retData[13]=std::numeric_limits<double>::quiet_NaN();
            break;
        }
        case 18://Standard class B equipment position report.
        {
            libais::Ais18* msg=reinterpret_cast<libais::Ais18*>(curAisMsg.get());

            //MMSI number (An integer).
            retData[0]=static_cast<double>(msg->mmsi);
            //There is no navigation status
            retData[1]=std::numeric_limits<double>::quiet_NaN();
            //There is no turn rate
            retData[2]=std::numeric_limits<double>::quiet_NaN();

            //Speed over ground in knots.
            if(msg->sog<=102.2) {
                retData[3]=msg->sog;
            } else {
                retData[3]=std::numeric_limits<double>::quiet_NaN();
            }

            //Position accuracy (0 or 1).
            retData[4]=static_cast<double>(msg->position_accuracy);

            //Latitude (degrees)
            if(fabs(msg->position.lat_deg)<=90) {
                retData[5]=msg->position.lat_deg;
            } else {
                retData[5]=std::numeric_limits<double>::quiet_NaN();
            }

            //Longitude (degrees)
            if(fabs(msg->position.lng_deg)<=180) {
                retData[6]=msg->position.lng_deg;
            } else {
                retData[6]=std::numeric_limits<double>::quiet_NaN();
            }

            //Course over ground (degrees East of North)
            if(msg->cog<360) {
                retData[7]=msg->cog;
            } else {
                retData[7]=std::numeric_limits<double>::quiet_NaN();
            }

            //True heading (degrees East of North)
            if(msg->true_heading<360) {
                retData[8]=static_cast<double>(msg->true_heading);
            } else {
                retData[8]=std::numeric_limits<double>::quiet_NaN();
            }

            //UTC second of minute.
            //The restriction against it being 60 means that an
            //incorrect UTC second will be reported when a
            //leapsecond is added.
            if(msg->timestamp<60) {
                retData[9]=static_cast<double>(msg->timestamp);
            } else {
                retData[9]=std::numeric_limits<double>::quiet_NaN();
            }

            //There is no special manoeuvre.
            retData[10]=std::numeric_limits<double>::quiet_NaN();

            //The message ID.
            retData[11]=static_cast<double>(msg->message_id);

            //The message repeat indicator
            retData[12]=static_cast<double>(msg->repeat_indicator);
            //There is no position latency indicator
            retData[13]=std::numeric_limits<double>::quiet_NaN();
            break;
        }
        case 19://Extended class B equipment position report.
        {
            libais::Ais19*msg=reinterpret_cast<libais::Ais19*>(curAisMsg.get());

            //MMSI number (An integer).
            retData[0]=static_cast<double>(msg->mmsi);
            //There is no navigation status
            retData[1]=std::numeric_limits<double>::quiet_NaN();
            //There is no turn rate
            retData[2]=std::numeric_limits<double>::quiet_NaN();

            //Speed over ground in knots.
            if(msg->sog<=102.2) {
                retData[3]=msg->sog;
            } else {
                retData[3]=std::numeric_limits<double>::quiet_NaN();
            }

            //Position accuracy (0 or 1).
            retData[4]=static_cast<double>(msg->position_accuracy);

            //Latitude (degrees)
            if(fabs(msg->position.lat_deg)<=90) {
                retData[5]=msg->position.lat_deg;
            } else {
                retData[5]=std::numeric_limits<double>::quiet_NaN();
            }

            //Longitude (degrees)
            if(fabs(msg->position.lng_deg)<=180) {
                retData[6]=msg->position.lng_deg;
            } else {
                retData[6]=std::numeric_limits<double>::quiet_NaN();
            }

            //Course over ground (degrees East of North)
            if(msg->cog<360) {
                retData[7]=msg->cog;
            } else {
                retData[7]=std::numeric_limits<double>::quiet_NaN();
            }

            //True heading (degrees East of North)
            if(msg->true_heading<360) {
                retData[8]=static_cast<double>(msg->true_heading);
            } else {
                retData[8]=std::numeric_limits<double>::quiet_NaN();
            }

            //UTC second of minute.
            //The restriction against it being 60 means that an
            //incorrect UTC second will be reported when a
            //leapsecond is added.
            if(msg->timestamp<60) {
                retData[9]=static_cast<double>(msg->timestamp);
            } else {
                retData[9]=std::numeric_limits<double>::quiet_NaN();
            }

            //There is no special manoeuvre.
            retData[10]=std::numeric_limits<double>::quiet_NaN();

            //The message ID.
            retData[11]=static_cast<double>(msg->message_id);

            //The message repeat indicator
            retData[12]=static_cast<double>(msg->repeat_indicator);

            //There is no position latency indicator
            retData[13]=std::numeric_limits<double>::quiet_NaN();
            break;
        }
        case 27://Long-range automatic identification system
                //broadcast message
        {
            libais::Ais27*msg=reinterpret_cast<libais::Ais27*>(curAisMsg.get());

            //MMSI number (An integer).
            retData[0]=static_cast<double>(msg->mmsi);
            //Navigational status
            retData[1]=static_cast<double>(msg->nav_status);
            //There is no turn rate
            retData[2]=std::numeric_limits<double>::quiet_NaN();

            //Speed over ground in knots.
            if(msg->sog<=102.2) {
                retData[3]=msg->sog;
            } else {
                retData[3]=std::numeric_limits<double>::quiet_NaN();
            }

            //Position accuracy (0 or 1).
            retData[4]=static_cast<double>(msg->position_accuracy);

            //Latitude (degrees)
            if(fabs(msg->position.lat_deg)<=90) {
                retData[5]=msg->position.lat_deg;
            } else {
                retData[5]=std::numeric_limits<double>::quiet_NaN();
            }

            //Longitude (degrees)
            if(fabs(msg->position.lng_deg)<=180) {
                retData[6]=msg->position.lng_deg;
            } else {
                retData[6]=std::numeric_limits<double>::quiet_NaN();
            }

            //Course over ground (degrees East of North)
            if(msg->cog<360) {
                retData[7]=msg->cog;
            } else {
                retData[7]=std::numeric_limits<double>::quiet_NaN();
            }

            //There is no true heading information
            retData[8]=std::numeric_limits<double>::quiet_NaN();

            //There is no UTC seconds timestamp in Message 27
            retData[9]=std::numeric_limits<double>::quiet_NaN();
            //There is no special manoeuvre information
            retData[10]=std::numeric_limits<double>::quiet_NaN();

            //The message ID.
            retData[11]=static_cast<double>(msg->message_id);

            //The message repeat indicator
            retData[12]=static_cast<double>(msg->repeat_indicator);

            //The position latency indicator
            retData[13]=static_cast<double>(msg->gnss);
            break;
        }
        default:
        //Unknown message type or it is not a position report
        //supported by this function
            return false;
    }
    return true;
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
