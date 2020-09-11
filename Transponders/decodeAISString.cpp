/**DECODEAISSTRING Decode one or more Automatic Identification System (AIS)
 *                 messages given a full National Marine Electronics
 *                 Association (NMEA) formatted ASCII text messages.
 *                 Multiple part messages can also be decoded. Information
 *                 appended by the receiver after the message checksum is
 *                 ignored.
 *
 *INPUTS: NMEAStrings A character string containing one or more Automatic
 *                 Identification System (AIS) messages given as National
 *                 Marine Electronics Association (NMEA) formatted ASCII
 *                 text. Each message is on its own line.
 *
 *OUTPUTS:decodedMessage A cell array of structures whose components are
 *                    the decoded message elements, one message per cell.
 *                    If an incomplete or invalid message is provided, it
 *                    will be skipped.
 *         reportName A cell array of the names of the general type of
 *                    each message decoded.
 *  reportDescription A cell array of strings describing the type of each
 *                    message decoded.
 *  fieldDescriptions A cell array holding a structure array with the same
 *                    general fields as each cell element of decodedMessage
 *                    but where each field contains a string describing
 *                    what is in the field.
 *
 *AIS messages are used for tracking cooperative ships and for the exchange
 *of information of interest to ships, such as tidal information. The basic
 *standard is [1] with additional messages defined in [2], [3], and at
 *http://www.e-navigation.nl/sites/default/files/asm_files/GN%20Release%202-23MAR15.pdf
 *and in [4]. Note that an older version of the ITU standard as well as [5]
 *contain some message types that have been removed from the newer
 *standards.
 *
 *This function is essentially a Matlab interface for libais from 
 *https://github.com/schwehr/libais
 *The library does not support all possible message types. Messages
 *8_367_22 and 8_366_22 from the library have not been implemented as it
 *appears that the standards are in flux/ or are not fully implemented in
 *the library. Messages 6_0_0 and 8_366_56 will not yet work as
 *decode_body.cpp in the library still lacks a appropriate switch
 *statements.
 *
 *The algorithm can be compiled for use in Matlab using the 
 *CompileCLibraries function.
 *
 *The algorithm is run in Matlab using the command format
 *[decodedMessage,reportName,reportDescription,fieldDescriptions]=decodeAISString(NMEAStrings)
 *
 *Sample NMEA messages from LibAIS for the different message types are
 *Message 1
 *NMEAStrings='!AIVDM,1,1,,A,100WhdhP0nJRdiFFHFvm??v00L12,0*13,raishub,1342569600';
 *Message 2
 *NMEAStrings='!AIVDM,1,1,,B,284;UGTdP4301>3L;B@Wk3TnU@A1,0*7C,raishub,1342572913';
 *Message 3
 *NMEAStrings='!AIVDM,1,1,,B,34hoV<5000Jw95`GWokbFTuf0000,0*6C,raishub,1342569600';
 *Message 4
 *NMEAStrings='!AIVDM,1,1,,B,4h3Owoiuiq000rdhR6G>oQ?020S:,0*10,raishub,1342569600';
 *Message 5 (two fragments put together using sprintf)
 *NMEAStrings=sprintf('!SAVDM,2,1,6,A,55NOvQP1u>QIL@O??SL985`u>0EQ18E=>222221J1p`884i6N344Sll1@m80,0*0C,b003669956,1426118503\n!SAVDM,2,2,6,A,TRA1iH88880,2*6F,b003669956,1426118503');
 *Message 6_0_0
 *NMEAStrings='!AIVDM,1,1,,A,6>l4v:h0006D000000P@0T0,2*56';
 *Message 6_1_0
 *NMEAStrings='!AIVDM,1,1,,B,65Ps:8=:0MjP0420<4U>1@E=B10i>04<fp0,2*23';
 *Message 7
 *NMEAStrings='!AIVDM,1,1,,A,7jvPoD@mMq;U,0*09';
 *Message 8_1_11
 *NMEAStrings='!AIVDM,1,1,,A,8@2<HV@0BkLN:0frqMPaQPtBRRIrwwejwwwwwwwwwwwwwwwwwwwwwwwwwt0,2*34';
 *Message 8_1_22
 *NMEAStrings='!AIVDM,1,1,0,B,803Ovrh0EPM0WB0h2l0MwJUi=6B4G9000aip8<2Bt2Hq2Qhp,0*01,d-084,S1582,t091042.00,T42.19038981,r003669945,1332321042';
 *Message 8_200_10
 *NMEAStrings='!SAVDM,1,1,5,B,85NLn@0j2d<8000000BhI?`50000,0*2A,d,S1241,t002333.00,T33.111314,D08MN-NO-BSABS1,1429316613';
 *Message 8_366_56
 *NMEAStrings='!BSVDM,1,1,,B,853>IhQKf6EQFDdajT?AbaAVhHEWebddhqHC5@?=KwisgP00DWjE,0*6D,b003669701,1429315201';
 *Message 9
 *NMEAStrings='!AIVDM,1,1,,B,9oVAuAI5;rRRv2OqTi?1uoP?=a@1,0*74';
 *Message 10
 *NMEAStrings='!AIVDM,1,1,,A,:5Ovc200B=5H,0*43';
 *Message 11
 *NMEAStrings='!AIVDM,1,1,,B,;028j>iuiq0DoO0ARF@EEmG008Pb,0*25,raishub,1342570856';
 *Message 12 (two fragments put together using sprintf)
 *NMEAStrings=sprintf('!AIVDM,2,1,1,A,<02PeAPpIkF06B?=PB?31P3?>DB?<rP@<51C5P3?>D13DPB?31P3?>DB,0*13\n!AIVDM,2,2,1,A,?<P?>PF86P381>>5<PoqP6?BP=1>41D?BIPB5@?BD@,4*66');
 *Message 14
 *NMEAStrings='!AIVDM,1,1,,A,>>M@rl1<59B1@E=@0000000,2*0D';
 *Message 15
 *NMEAStrings='!SAVDM,1,1,,B,?03OwnB0ACVlD00,2*59';
 *Message 16
 *NMEAStrings='!SAVDO,1,1,,B,@03OwnQ9RgLP3h0000000000,0*32,b003669978,1426173689';
 *Message 17
 *NMEAStrings='!AIVDM,1,1,,A,A6WWW6gP00a3PDlEKLrarOwUr8Mg,0*03';
 *Message 18
 *NMEAStrings='!SAVDM,1,1,4,B,B5NU=J000=l0BD6l590EkwuUoP06,0*61';
 *Message 19 (two fragments put together using sprintf)
 *NMEAStrings=sprintf('!AIVDM,2,1,7,B,C5NMbDQl0NNJC7VNuC<v`7NF4T28V@2g0J6F::00000,0*59\n!AIVDM,2,2,7,B,0J70<RRS0,0*30');
 *Message 20
 *NMEAStrings='!SAVDM,1,1,6,B,Dh3OwjhflnfpLIF9HM1F9HMaF9H,2*3E';
 *Message 21 (two fragments put together using sprintf)
 *NMEAStrings=sprintf('!AIVDM,2,1,9,A,ENk`sO70VQ97aRh1T0W72V@611@=FVj<;V5d@00003v,0*50\n!AIVDM,2,2,9,A,P0<M0,0*3E');
 *Message 22
 *NMEAStrings='!SAVDM,1,1,6,A,F030owj2N2P6Ubib@=4q35b10000,0*74';
 *Message 23
 *NMEAStrings='!AIVDM,1,1,,B,G02:KpP1R`sn@291njF00000900,2*1C';
 *Message 24
 *NMEAStrings='!AIVDM,1,1,,A,H44cj<0DdvlHhuB222222222220,2*46';
 *Message 25
 *NMEAStrings='!AIVDM,1,1,,B,I6S`3Tg@T0a3REBEsjJcT?wSi0fM,0*02';
 *Message 26 (two fragments put together using sprintf)
 *NMEAStrings=sprintf('!AIVDM,2,1,2,B,JfgwlGvNwts9?wUfQswQ<gv9Ow7wCl?nwv0wOi=mwd?,0*03\n!AIVDM,2,2,2,B,oW8uwNg3wNS3tV,5*71');
 *Message 27
 *NMEAStrings='!AIVDM,1,1,,B,K815>P8=5EikdUet,0*6B';
 *
 *REFERENCES:
 *[1] "Technical characteristics for an automatic identification system
 *    using time division multiple access in the VHF maritime mobile
 *    frequency band," International Telecommunication Union,
 *    Radiocommunication Sector Std. ITU-R M.1371-5, Feb. 2014. [Online].
 *    Available: http://www.itu.int/rec/R-REC-M.1371/en
 *[2] "Guidance on the use of AIS application-specific messages,"
 *    International Maritime Organization, 2 Jun. 2013. [Online].
 *    Available: http://www.iho.int/mtg_docs/com_wg/IHOTC/IHOTC_Misc/IMO%20SN_Circ289.pdf
 *[3] "International standard for tracking and tracing on inland waterways
 *    (VTT)," United Nations Economic and Social Council, Economic
 *    Commission for Europe, Inland Transport Committee, New York and
 *    Geneva, 2007. [Online].
 *    Available: http://www.unece.org/fileadmin/DAM/trans/doc/finaldocs/sc3/ECE-TRANS-SC3-176e.pdf
 *[4] "St. Lawrence Seaway AIS Data Messaging Formats and Specifications,"
 *    U.S. Department of Transportation, Revision 4.0A, 9 May 2002.
 *    [Online].
 *    Available: http://www.greatlakes-seaway.com/en/pdf/aisdata.pdf
 *[5] "Guidance on the application of AIS binary messages," International
 *    Maritime Organization, 28 May 2004. [Online].
 *    Available: http://www.imo.org/blast/blastDataHelper.asp?data id=10741&filename=236.pd
 *
 *November 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.*/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

//This definition is necessary or else the aislib files will do
//undesirable things when bad messages are passed.
#ifndef NDEBUG
#define NDEBUG
#endif

/*This header is required by Matlab.*/
#include "mex.h"

//The AISLib headers
#include "ais.h"
#include "vdm.h"

//For fabs
#include<cmath>

//Functions for parsing some parts of the strings that is not handled by
//AISLib.
#include "AISFuncs.hpp"

#include "MexValidation.h"

//This is for parsing the input into strings that can be provided to the
//VdmStream class to reconstruct the messages.
#include <iostream>
#include <sstream>
#include <string>
//For count
#include <algorithm>

static const size_t numReportTypes=28;

//The names from Table 46 in ITU-R M.1371-5 
static const char *AISReportNames[numReportTypes] = {
    "",//There is no zero
    "Position Report",//1
    "Position Report",//2
    "Position Report",//3
    "Base Station Report",//4
    "Static and Voyage Related Data",//5
    "Binary Addressed Message",//6
    "Binary Acknowledgement",//7
    "Binary Broadcast Message",//8
    "Standard search and rescue (SAR) Aircraft Position Report",//9
    "UTC/Date Inquiry",//10
    "UTC/Date Response",//11
    "Addressed Safety Related Message",//12
    "Safety Related Acknowledgement",//13
    "Safety Related Broadcast Message",//14
    "Interrogation",//15
    "Assignment Mode Command",//16
    "DGNSS Broadcast Binary Message",//17
    "Standard Class B Equipment Position Report",//18
    "Extended Class B Equipment Position Report",//19
    "Data Link Management Message",//20
    "Aids-to-Navigation Report",//21
    "Channel management",//22
    "Group Assignment Command",//23
    "Static Data Report",//24
    "Single Slot Binary Message",//25
    "Multiple Slot Binary Message with Communications State",//26
    "Position Report for Long-Range Applications"//27
};

//The descriptions from Table 46 in ITU-R M.1371-5 
static const char *AISReportDescriptions[numReportTypes] = {
    "",//There is no zero
    "Scheduled position report; (Class A shipborne mobile equipment)",//1
    "Assigned scheduled position report; (Class A shipborne mobile equipment)",//2
    "Special position report, response to interrogation; (Class A shipborne mobile equipment)",//3
    "Position, UTC, date and current slot number of base station",//4
    "Scheduled static and voyage related vessel data report; (Class A shipborne mobile equipment)",//5
    "Binary data for addressed communication",//6
    "Acknowledgement of received addressed binary data",//7
    "Binary data for broadcast communication",//8
    "Position report for airborne stations involved in search and rescue (SAR) operations, only",//9
    "Request UTC and date",//10
    "Current UTC and date if available",//11
    "Safety related data for addressed communication",//12
    "Acknowledgement of received addressed safety related message",//13
    "Safety related data for broadcast communication",//14
    "Request for a specific message type (can result in multiple responses from one or several stations)",//15
    "Assignment of a specific report behaviour by competent authority using a Base station",//16
    "DGNSS corrections provided by a base station",//17
    "Standard position report for Class B shipborne mobile equipment to be used instead of Messages 1, 2, 3",//18
    "No longer required; Extended position report for Class B shipborne mobile equipment; contains additional static information",//19
    "Reserve slots for Base station(s)",//20
    "Position and status report for aids-to-navigation",//21
    "Management of channels and transceiver modes by a Base station",//22
    "Assignment of a specific report behaviour by competent authority using a Base station to a specific group of mobiles",//23
    "Additional data assigned to an MMSI Part A: Name Part B: Static Data",//24
    "Short unscheduled binary data transmission (Broadcast or addressed)",//25
    "Scheduled binary data transmission (Broadcast or addressed)",//26
    "Class A and Class B 'SO' shipborne mobile equipment outside base station coverage"//27
};

//Prototypes for functions
bool extractAISMessageData(std::unique_ptr<libais::AisMsg> &aisMsg,mxArray **decodedMessage,mxArray **fieldDescriptions);
void AIS_1_2_3_ToMatlab(libais::Ais1_2_3 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions);
void AIS_4_11_ToMatlab(libais::Ais4_11 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions);
void AIS_5_ToMatlab(libais::Ais5 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions);
void AIS_7_13_ToMatlab(libais::Ais7_13 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions);
void AIS_9_ToMatlab(libais::Ais9 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions);
void AIS_10_ToMatlab(libais::Ais10 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions);
void AIS_12_ToMatlab(libais::Ais12 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions);
void AIS_14_ToMatlab(libais::Ais14 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions);
void AIS_15_ToMatlab(libais::Ais15 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions);
void AIS_16_ToMatlab(libais::Ais16 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions);
void AIS_17_ToMatlab(libais::Ais17 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions);
void AIS_18_ToMatlab(libais::Ais18 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions);
void AIS_19_ToMatlab(libais::Ais19 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions);
void AIS_20_ToMatlab(libais::Ais20 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions);
void AIS_21_ToMatlab(libais::Ais21 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions);
void AIS_22_ToMatlab(libais::Ais22 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions);
void AIS_23_ToMatlab(libais::Ais23 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions);
void AIS_24_ToMatlab(libais::Ais24 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions);
void AIS_25_ToMatlab(libais::Ais25 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions);
void AIS_26_ToMatlab(libais::Ais26 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions);
void AIS_27_ToMatlab(libais::Ais27 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions);
void AIS_6_0_0_ToMatlab(libais::Ais6_0_0 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions);
void AIS_6_1_0_ToMatlab(libais::Ais6_1_0 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions);
void AIS_6_1_1_ToMatlab(libais::Ais6_1_1 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions);
void AIS_6_1_2_ToMatlab(libais::Ais6_1_2 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions);
void AIS_6_1_3_ToMatlab(libais::Ais6_1_3 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions);
void AIS_6_1_4_ToMatlab(libais::Ais6_1_4 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions);
void AIS_6_1_5_ToMatlab(libais::Ais6_1_5 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions);
void AIS_6_1_12_ToMatlab(libais::Ais6_1_12 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions);
void AIS_6_1_14_ToMatlab(libais::Ais6_1_14 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions);
void AIS_6_1_18_ToMatlab(libais::Ais6_1_18 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions);
void AIS_6_1_20_ToMatlab(libais::Ais6_1_20 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions);
void AIS_6_1_25_ToMatlab(libais::Ais6_1_25 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions);
void AIS_6_1_32_ToMatlab(libais::Ais6_1_32 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions);
void AIS_6_1_40_ToMatlab(libais::Ais6_1_40 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions);
void AIS_8_1_0_ToMatlab(libais::Ais8_1_0 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions);
void AIS_8_1_11_ToMatlab(libais::Ais8_1_11 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions);
void AIS_8_1_13_ToMatlab(libais::Ais8_1_13 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions);
void AIS_8_1_15_ToMatlab(libais::Ais8_1_15 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions);
void AIS_8_1_16_ToMatlab(libais::Ais8_1_16 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions);
void AIS_8_1_17_ToMatlab(libais::Ais8_1_17 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions);
void AIS_8_1_19_ToMatlab(libais::Ais8_1_19 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions);
void AIS_8_1_21_ToMatlab(libais::Ais8_1_21 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions);
void AIS_8_1_22_ToMatlab(libais::Ais8_1_22 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions);
void AIS_8_1_24_ToMatlab(libais::Ais8_1_24 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions);
void AIS_8_1_26_ToMatlab(libais::Ais8_1_26 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions);
void AIS_8_1_27_ToMatlab(libais::Ais8_1_27 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions);
void AIS_8_1_29_ToMatlab(libais::Ais8_1_29 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions);
void AIS_8_1_31_ToMatlab(libais::Ais8_1_31 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions);
void AIS_8_200_10_ToMatlab(libais::Ais8_200_10 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions);
void AIS_8_200_23_ToMatlab(libais::Ais8_200_23 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions);
void AIS_8_200_24_ToMatlab(libais::Ais8_200_24 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions);
void AIS_8_200_40_ToMatlab(libais::Ais8_200_40 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions);
void AIS_8_200_55_ToMatlab(libais::Ais8_200_55 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions);
void AIS_8_366_56_ToMatlab(libais::Ais8_366_56 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) { 

    if(nrhs!=1) {
        mexErrMsgTxt("Wrong number of inputs.");
        return;
    }
    
    if(nlhs>4) {
        mexErrMsgTxt("Wrong number of outputs.");
        return;
    }
    
    if(mxGetN(prhs[0])!=1&&mxGetM(prhs[0])!=1) {
        mexErrMsgTxt("The input has an invalid dimensionality.");
        return;
    }
    
    if(mxIsChar(prhs[0])) {
        size_t maxNumMessages,numMessagesDecoded=0;
        std::istringstream ss;
        char *msgChars=mxArrayToString(prhs[0]);
        std::string allMsg(msgChars);
        //This will stitch together a set of messages.
        libais::VdmStream decodedAISData;
        mxArray *decodedMessageCell, *reportNameCell, *reportDescriptionCell, *fieldDescriptionsCell;
        std::string curMsg;
        
        mxFree(msgChars);
        
        //The maximum possible number of messages to decode.
        maxNumMessages=std::count(allMsg.begin(),allMsg.end(),'\n')+1;
     
        //Alocate space for the return variables.
        decodedMessageCell=mxCreateCellMatrix(maxNumMessages,1);
        plhs[0]=decodedMessageCell;
        if(nlhs>1){
            reportNameCell=mxCreateCellMatrix(maxNumMessages,1);
            plhs[1]=reportNameCell;
            if(nlhs>2) {
                reportDescriptionCell=mxCreateCellMatrix(maxNumMessages,1);
                plhs[2]=reportDescriptionCell;
                if(nlhs>3){
                    fieldDescriptionsCell=mxCreateCellMatrix(maxNumMessages,1);
                    plhs[3]=fieldDescriptionsCell;
                }
            }
        }
        
        ss=std::istringstream(allMsg);//To extract one string at a time.
        //Go through all of the strings.
        while(std::getline(ss, curMsg,'\n')) {
            const size_t numBeforeAdd=decodedAISData.size();
            std::string dataPayload, endPart;
            bool pushSucceeded;
            
            //First, remove any additional information that the receiver
            //might have added after the end of the message.
            separateMessageAndTimestamp(curMsg,dataPayload,endPart);
                        
            //If the separation succeeded, then push the
            pushSucceeded=decodedAISData.AddLine(dataPayload);

            //If a complete message was successfully decoded, then decode
            //the message.
            if(pushSucceeded&&numBeforeAdd!=decodedAISData.size()) {
                mxArray *decodedMessage, *fieldDescriptions;
                std::unique_ptr<libais::AisMsg> curAisMsg=decodedAISData.PopOldestMessage();
                bool decodeSuccessful;
                
                if(nlhs<4) {
                    decodeSuccessful=extractAISMessageData(curAisMsg,&decodedMessage,NULL);
                } else {
                    decodeSuccessful=extractAISMessageData(curAisMsg,&decodedMessage,&fieldDescriptions);
                }

                if(decodeSuccessful) {
                    mxSetCell(decodedMessageCell, numMessagesDecoded, decodedMessage);
                    if(nlhs>1) {
                        const char *charString=AISReportNames[curAisMsg->message_id];
                        mxSetCell(reportNameCell, numMessagesDecoded, mxCreateCharMatrixFromStrings(1,&charString));
                        if(nlhs>2) {
                            const char *charString1=AISReportDescriptions[curAisMsg->message_id];
                            mxSetCell(reportDescriptionCell, numMessagesDecoded, mxCreateCharMatrixFromStrings(1,&charString1));
                            if(nlhs>3){
                                mxSetCell(fieldDescriptionsCell, numMessagesDecoded, fieldDescriptions);
                            }
                        }
                    }
                    numMessagesDecoded++;
                }
            }
        }

        //Resize the return values to the actual number of messages
        //decoded.
        mxSetM(decodedMessageCell,numMessagesDecoded);
        if(nlhs>1) {
            mxSetM(reportNameCell,numMessagesDecoded);
            if(nlhs>2) {
                mxSetM(reportDescriptionCell,numMessagesDecoded);
                if(nlhs>3) {
                    mxSetM(fieldDescriptionsCell,numMessagesDecoded);
                }
            }
        }
    } else {
        mexErrMsgTxt("Invalid data type passed."); 
    }
}

bool extractAISMessageData(std::unique_ptr<libais::AisMsg> &aisMsg,mxArray **decodedMessage,mxArray **fieldDescriptions){
/*EXTRACTAISMESSAGEDATA This function takes an NMEA AIS message and returns
 *                  the data in a structure in decodedMessage and if
 *                  fieldDescriptions is not NULL, it also returns a set of
 *                  descriptions of the fields in the structure. If
 *                  decoding is successful, then this function returns
 *                  true. If decoding is not successful, then this function
 *                  returns false and   decodedMessage and
 *                  fieldDescriptions remain unchanged.
 * 
 *This function is essentially just a bunch of switch statements that call
 *the correct function for each message type.
 *
 *November 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.*/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

    switch(aisMsg->message_id){
        case 1://Position report, Class A, scheduled.
        case 2://Position report, Class A, assigned.
        case 3://Position report, Class A, interrogated.
        {
            libais::Ais1_2_3 *msg=reinterpret_cast<libais::Ais1_2_3*>(aisMsg.get());
            
            AIS_1_2_3_ToMatlab(msg,decodedMessage,fieldDescriptions);
            return true;
        }
        case 4://Base station report.
        case 11://UTC/ date response.
        {
            libais::Ais4_11 *msg=reinterpret_cast<libais::Ais4_11*>(aisMsg.get());
            
            AIS_4_11_ToMatlab(msg,decodedMessage,fieldDescriptions);
            return true;
        }
        case 5://Static and voyage related data.
        {
            libais::Ais5 *msg=reinterpret_cast<libais::Ais5*>(aisMsg.get());
            
            AIS_5_ToMatlab(msg,decodedMessage,fieldDescriptions);
            return true;
        }
        case 6://Binary addressed message
        {
            libais::Ais6 *msg=reinterpret_cast<libais::Ais6*>(aisMsg.get());

            switch(msg->dac) {
                case libais::AIS_DAC_0_TEST:
                    if(msg->fi==0) {//Zeni Lite Buoy Co., Ltd buoy status.
                        AIS_6_0_0_ToMatlab(reinterpret_cast<libais::Ais6_0_0*>(msg),decodedMessage,fieldDescriptions);
                        return true;
                    } else {
                        return false;
                    }
                case libais::AIS_DAC_1_INTERNATIONAL:
                    switch(msg->fi) {
                        case 0:// Text message.  ITU 1371-1
                            AIS_6_1_0_ToMatlab(reinterpret_cast<libais::Ais6_1_0*>(msg),decodedMessage,fieldDescriptions);
                            return true;
                        case 1:// Application ack.  ITU 1371-1
                            AIS_6_1_1_ToMatlab(reinterpret_cast<libais::Ais6_1_1*>(msg),decodedMessage,fieldDescriptions);
                            return true;
                        case 2:// Interrogation for a DAC/FI.  ITU 1371-1
                            AIS_6_1_2_ToMatlab(reinterpret_cast<libais::Ais6_1_2*>(msg),decodedMessage,fieldDescriptions);
                            return true;
                        case 3:// Capability interogation.  ITU 1371-1
                            AIS_6_1_3_ToMatlab(reinterpret_cast<libais::Ais6_1_3*>(msg),decodedMessage,fieldDescriptions);
                            return true;
                        case 4:// Capability interogation reply.  ITU 1371-1
                            AIS_6_1_4_ToMatlab(reinterpret_cast<libais::Ais6_1_4*>(msg),decodedMessage,fieldDescriptions);
                            return true;
                        case 5:// International function message 5: Application ack to addr binary message.
                            AIS_6_1_5_ToMatlab(reinterpret_cast<libais::Ais6_1_5*>(msg),decodedMessage,fieldDescriptions);
                            return true;
                        case 12:// IMO Circ 236 Dangerous cargo indication
                            AIS_6_1_12_ToMatlab(reinterpret_cast<libais::Ais6_1_12*>(msg),decodedMessage,fieldDescriptions);
                            return true;
                        case 14:// IMO Circ 236 Tidal window
                            AIS_6_1_14_ToMatlab(reinterpret_cast<libais::Ais6_1_14*>(msg),decodedMessage,fieldDescriptions);
                            return true;
                        case 18:// IMO Circ 289 Clearance time to enter port
                            AIS_6_1_18_ToMatlab(reinterpret_cast<libais::Ais6_1_18*>(msg),decodedMessage,fieldDescriptions);
                            return true;
                        case 20:// IMO Circ 289 Berthing data
                            AIS_6_1_20_ToMatlab(reinterpret_cast<libais::Ais6_1_20*>(msg),decodedMessage,fieldDescriptions);
                            return true;
                        case 25:// IMO Circ 289 Dangerous cargo indication 2
                            AIS_6_1_25_ToMatlab(reinterpret_cast<libais::Ais6_1_25*>(msg),decodedMessage,fieldDescriptions);
                            return true;
                        case 32:// IMO Circ 289 Tidal window
                            AIS_6_1_32_ToMatlab(reinterpret_cast<libais::Ais6_1_32*>(msg),decodedMessage,fieldDescriptions);
                            return true;
                        case 40:// Number of persons on board.  ITU 1371-1
                            AIS_6_1_40_ToMatlab(reinterpret_cast<libais::Ais6_1_40*>(msg),decodedMessage,fieldDescriptions);
                            return true;
                        default:
                            return false;
                    }
                default:
                    return false;
            }
        }
        case 7://Binary acknowledgement
        case 13://Message 13: Safety related acknowledge
        {
            libais::Ais7_13 *msg=reinterpret_cast<libais::Ais7_13*>(aisMsg.get());
            
            AIS_7_13_ToMatlab(msg,decodedMessage,fieldDescriptions);
            return true;
        }
        case 8://Binary broadcast message
        {
            libais::Ais8 *msg=reinterpret_cast<libais::Ais8*>(aisMsg.get());
            
            switch(msg->dac) {
                case libais::AIS_DAC_1_INTERNATIONAL:
                    switch(msg->fi) {
                        case 0:// Text telegram ITU 1371-1
                            AIS_8_1_0_ToMatlab(reinterpret_cast<libais::Ais8_1_0*>(msg),decodedMessage,fieldDescriptions);
                            return true;
                        case 11:// IMO Circ 289 met hydro
                            AIS_8_1_11_ToMatlab(reinterpret_cast<libais::Ais8_1_11*>(msg),decodedMessage,fieldDescriptions);
                            return true;
                        case 13:// IMO Circ 236 Fairway closed
                            AIS_8_1_13_ToMatlab(reinterpret_cast<libais::Ais8_1_13*>(msg),decodedMessage,fieldDescriptions);
                            return true;
                        case 15:// IMO Circ 236 Extended ship static and
                                //voyage data
                            AIS_8_1_15_ToMatlab(reinterpret_cast<libais::Ais8_1_15*>(msg),decodedMessage,fieldDescriptions);
                            return true;
                        case 16:// IMO Circ 236 Number of persons on board
                            AIS_8_1_16_ToMatlab(reinterpret_cast<libais::Ais8_1_16*>(msg),decodedMessage,fieldDescriptions);
                            return true;
                        case 17:// IMO Circ 236 VTS Generated/synthetic
                                //targets
                            AIS_8_1_17_ToMatlab(reinterpret_cast<libais::Ais8_1_17*>(msg),decodedMessage,fieldDescriptions);
                            return true;
                        case 19:// IMO Circ 289 Marine traffic signal
                            AIS_8_1_19_ToMatlab(reinterpret_cast<libais::Ais8_1_19*>(msg),decodedMessage,fieldDescriptions);
                            return true;
                        case 21:// IMO Circ 289 Weather observation report
                                //from ship
                            AIS_8_1_21_ToMatlab(reinterpret_cast<libais::Ais8_1_21*>(msg),decodedMessage,fieldDescriptions);
                            return true;
                        case 22:// Area notice
                            AIS_8_1_22_ToMatlab(reinterpret_cast<libais::Ais8_1_22*>(msg),decodedMessage,fieldDescriptions);
                            return true;
                        case 24:// IMO Circ 289 Extended ship static and
                                //voyage-related
                            AIS_8_1_24_ToMatlab(reinterpret_cast<libais::Ais8_1_24*>(msg),decodedMessage,fieldDescriptions);
                            return true;
                        case 26:// IMO Circ 289 Environmental
                            AIS_8_1_26_ToMatlab(reinterpret_cast<libais::Ais8_1_26*>(msg),decodedMessage,fieldDescriptions);
                            return true;
                        case 27:// IMO Circ 289 Route information
                            AIS_8_1_27_ToMatlab(reinterpret_cast<libais::Ais8_1_27*>(msg),decodedMessage,fieldDescriptions);
                            return true;
                        case 29:// IMO Circ 289 Text description
                            AIS_8_1_29_ToMatlab(reinterpret_cast<libais::Ais8_1_29*>(msg),decodedMessage,fieldDescriptions);
                            return true;
                        case 31:// IMO Circ 289 Meteorological and
                                //Hydrographic data
                            AIS_8_1_31_ToMatlab(reinterpret_cast<libais::Ais8_1_31*>(msg),decodedMessage,fieldDescriptions);
                            return true;
                        default:
                            return false;
                    }
                case libais::AIS_DAC_200_RIS:
                    switch(msg->fi) {
                        case 10:// ECE-TRANS-SC3-2006-10e-RIS.pdf
                                //River Information System, inland ship
                                //static and voyage related data
                            AIS_8_200_10_ToMatlab(reinterpret_cast<libais::Ais8_200_10*>(msg),decodedMessage,fieldDescriptions);
                            return true;
                        case 23://ECE-TRANS-SC3-2006-10e-RIS.pdf 
                                //River Information System
                            AIS_8_200_23_ToMatlab(reinterpret_cast<libais::Ais8_200_23*>(msg),decodedMessage,fieldDescriptions);
                            return true;
                        case 24://ECE-TRANS-SC3-2006-10e-RIS.pdf River
                                //Information System: Water Level
                            AIS_8_200_24_ToMatlab(reinterpret_cast<libais::Ais8_200_24*>(msg),decodedMessage,fieldDescriptions);
                            return true;
                        case 40://ECE-TRANS-SC3-2006-10e-RIS.pdf 
                                //River Information System
                            AIS_8_200_40_ToMatlab(reinterpret_cast<libais::Ais8_200_40*>(msg),decodedMessage,fieldDescriptions);
                            return true;
                        case 55://ECE-TRANS-SC3-2006-10e-RIS.pdf 
                                //River Information System
                            AIS_8_200_55_ToMatlab(reinterpret_cast<libais::Ais8_200_55*>(msg),decodedMessage,fieldDescriptions);
                            return true;
                        default:
                            return false;
                    }
                case libais::AIS_DAC_366_UNITED_STATES_OF_AMERICA:
                    switch(msg->fi){
                        case 56:// USCG Blue Force encrypted message.
                            AIS_8_366_56_ToMatlab(reinterpret_cast<libais::Ais8_366_56*>(msg),decodedMessage,fieldDescriptions);
                            return true;
                        default:
                            return false;
                    }
                default:
                    return false;
            }
        }
        case 9:
        {
            libais::Ais9 *msg=reinterpret_cast<libais::Ais9*>(aisMsg.get());
            
            AIS_9_ToMatlab(msg,decodedMessage,fieldDescriptions);
            return true;
        }
        case 10:
        {
            libais::Ais10 *msg=reinterpret_cast<libais::Ais10*>(aisMsg.get());
            
            AIS_10_ToMatlab(msg,decodedMessage,fieldDescriptions);
            return true;
        }
        case 12:
        {
            libais::Ais12 *msg=reinterpret_cast<libais::Ais12*>(aisMsg.get());
            
            AIS_12_ToMatlab(msg,decodedMessage,fieldDescriptions);
            return true;
        }
        case 14:
        {
            libais::Ais14 *msg=reinterpret_cast<libais::Ais14*>(aisMsg.get());
            
            AIS_14_ToMatlab(msg,decodedMessage,fieldDescriptions);
            return true;
        }
        case 15:
        {
            libais::Ais15 *msg=reinterpret_cast<libais::Ais15*>(aisMsg.get());
            
            AIS_15_ToMatlab(msg,decodedMessage,fieldDescriptions);
            return true;
        }
        case 16:
        {
            libais::Ais16 *msg=reinterpret_cast<libais::Ais16*>(aisMsg.get());
            
            AIS_16_ToMatlab(msg,decodedMessage,fieldDescriptions);
            return true;
        }
        case 17:
        {
            libais::Ais17 *msg=reinterpret_cast<libais::Ais17*>(aisMsg.get());
            
            AIS_17_ToMatlab(msg,decodedMessage,fieldDescriptions);
            return true;
        }
        case 18:
        {
            libais::Ais18 *msg=reinterpret_cast<libais::Ais18*>(aisMsg.get());
            
            AIS_18_ToMatlab(msg,decodedMessage,fieldDescriptions);
            return true;
        }
        case 19:
        {
            libais::Ais19 *msg=reinterpret_cast<libais::Ais19*>(aisMsg.get());
            
            AIS_19_ToMatlab(msg,decodedMessage,fieldDescriptions);
            return true;
        }
        case 20:
        {
            libais::Ais20 *msg=reinterpret_cast<libais::Ais20*>(aisMsg.get());
            
            AIS_20_ToMatlab(msg,decodedMessage,fieldDescriptions);
            return true;
        }
        case 21:
        {
            libais::Ais21 *msg=reinterpret_cast<libais::Ais21*>(aisMsg.get());
            
            AIS_21_ToMatlab(msg,decodedMessage,fieldDescriptions);
            return true;
        }
        case 22:
        {
            libais::Ais22 *msg=reinterpret_cast<libais::Ais22*>(aisMsg.get());
            
            AIS_22_ToMatlab(msg,decodedMessage,fieldDescriptions);
            return true;
        }
        case 23:
        {
            libais::Ais23 *msg=reinterpret_cast<libais::Ais23*>(aisMsg.get());
            
            AIS_23_ToMatlab(msg,decodedMessage,fieldDescriptions);
            return true;
        }
        case 24:
        {
            libais::Ais24 *msg=reinterpret_cast<libais::Ais24*>(aisMsg.get());
            
            AIS_24_ToMatlab(msg,decodedMessage,fieldDescriptions);
            return true;
        }
        case 25:
        {
            libais::Ais25 *msg=reinterpret_cast<libais::Ais25*>(aisMsg.get());
            
            AIS_25_ToMatlab(msg,decodedMessage,fieldDescriptions);
            return true;
        }
        case 26:
        {
            libais::Ais26 *msg=reinterpret_cast<libais::Ais26*>(aisMsg.get());
            
            AIS_26_ToMatlab(msg,decodedMessage,fieldDescriptions);
            return true;
        }
        case 27:
        {
            libais::Ais27 *msg=reinterpret_cast<libais::Ais27*>(aisMsg.get());
            
            AIS_27_ToMatlab(msg,decodedMessage,fieldDescriptions);
            return true;
        }
        default:
            return false;
    }
}

void AIS_1_2_3_ToMatlab(libais::Ais1_2_3 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions) {
    const size_t numberOfFields=26;
    const mwSize dims[2] = {1, 1};
    //These are the names of all of the fields that will be added to the
    //Matlab structre array.
    const char *fieldNames[numberOfFields] = {"message_id", "repeat_indicator",
    "mmsi", "nav_status","rot", "sog",
    "position_accuracy", "lon", "lat", "cog", "true_heading", "timestamp",
    "special_manoeuvre", "spare", "raim", "sync_state", "slot_timeout",
    "received_stations", "slot_number", "utc_hour", "utc_min", "utc_spare",
    "slot_offset", "slot_increment", "slots_to_allocate", "keep_flag"};
    mxArray *theStruct;
    
    //If text descriptions of the fields are desired, then put them into a
    //structure.
    if(fieldDescriptions!=NULL) {
        const char *descriptionStrings[numberOfFields] = {
         "Message identifier.",//0
         "Repeat indicator. Used by the repeater to indicate how many times a message has been repeated. 0 = default; 3 = do not repeat any more.",//1
         "Unique identifier such as a Maritime mobile service identity (MMSI) number",//2
         "0=under way using engine\n1=at anchor,\n2=not under command\n3=restricted maneuverability\n4=constrained by her draught\n5=moored\n6=aground\n7=engaged in fishing\n8=under way sailing\n9=reserved for future amendment of navigational status for ships carrying DG, HS, or MP, or IMO hazard or pollutant category C, high speed craft (HSC)\n10=reserved for future amendment of navigational status for ships carrying dangerous goods (DG), harmful substances (HS) or marine pollutants (MP), or IMO hazard or pollutant category A, wing in ground (WIG)\n11=power-driven vessel towing astern (regional use)\n12=power-driven vessel pushing ahead or towing alongside (regional use)\n13 = reserved for future use\n14=AIS-SART (active), MOB-AIS, EPIRB-AIS\n15=undefined = default (also used by AIS-SART, MOB-AIS and EPIRB- AIS under test)",//3
         "Rate of turn in degrees per minute clockwise",//4
         "Speed over ground in knots. 102.3 means 102.2 knots or higher",//5
         "Position accuracy. 1=high (<=10m), 0=low (>10m), 0=default",//6
         "Longitude East in degrees (WGS-84)",//7
         "Latitude North in degrees (WGS-84)",//8
         "Course over ground in degrees East of North",//9
         "True heading (the direction the ship is pointing/ would be going without wind/currents) in degrees East of North",//10
         "UTC second when report was generated or 61 if positioning system is in manual input mode or 62 if in dead reckining mode or 63 if positioning system is inoperative",//11
         "1=not engaged in special manoeuvre\n2=engaged in special manoeuvre\n(i.e. regional passing arrangement on Inland Waterway)",//12
         "Not used. Should be set to zero. Reserved for future use.",//13
         "Receiver autonomous integrity monitoring (RAIM) flag of electronic position fixing device; 0 = RAIM not in use = default; 1 = RAIM in use.",//14
         "Synchronization state\n0 UTC direct\n1 UTC indirect\n2 Station is synchronized to a base station (base direct)\n3 Station is synchronized to another station based on the highest number of received stations or to another mobile station, which is directly synchronized to a base station",//15
         "(For Message IDs 1 and 2) Specifies frames remaining until a new slot is selected\n0 means that this was the last transmission in this slot\n1-7 means that 1 to 7 frames respectively are left until slot change",//16
         "(For Message IDs 1 and 2) Number of other stations (not own station) which the station currently is receiving (between 0 and 16 383).",//17
         "(For Message IDs 1 and 2) Slot number used for this transmission (between 0 and 2 249).",//18
         "(For Message IDs 1 and 2) Hour (0-23) in UTC time, if available",//19
         "(For Message IDs 1 and 2) Minute (0-23) in UTC time, if available",//20
         "(For Message IDs 1 and 2) Two extra (unused) bits that are associated with the time message",//21
         "(For Message IDs 1 and 2) If the slot time-out value is 0 (zero) then the slot offset should indicate the offset to the slot in which transmission will occur during the next frame. If the slot offset is zero, the slot should be de-allocated after transmission.",//22
         "(For Message ID 3) Offset to next slot to be used, or zero (0) if no more transmissions",//23
         "(For Message ID 3) Number of consecutive slots to allocate.\n0 = 1 slot,\n1 = 2 slots,\n2 = 3 slots,\n3 = 4 slots,\n4 = 5 slots,\n5 = 1 slot; offset = slot increment + 8 192,\n6 = 2 slots; offset = slot increment + 8 192,\n7 = 3 slots; offset = slot increment + 8 192.\nUse of 5 to 7 removes the need for RATDMA broadcast for scheduled transmissions up to 6 min intervals",//24
          "(For Message ID 3) Set to TRUE = 1 if the slot remains allocated for one additional frame"//25
        };
        unsigned int i;
        
        *fieldDescriptions=mxCreateStructArray(2, dims, numberOfFields, fieldNames);
        
        for(i=0;i<numberOfFields;i++) {
            mxSetFieldByNumber(*fieldDescriptions,0,static_cast<int>(i),mxCreateCharMatrixFromStrings(1,&descriptionStrings[i]));
        }
    }
    
    //Create the Matlab structure array.
    theStruct=mxCreateStructArray(2, dims, numberOfFields, fieldNames);
    *decodedMessage=theStruct;
    
    //Fill all of the elements of the structure
    mxSetFieldByNumber(theStruct,0,0,intMat2MatlabDoubles(&msg->message_id,1,1));
    mxSetFieldByNumber(theStruct,0,1,intMat2MatlabDoubles(&msg->repeat_indicator,1,1));
    mxSetFieldByNumber(theStruct,0,2,intMat2MatlabDoubles(&msg->mmsi,1,1));
    mxSetFieldByNumber(theStruct,0,3,intMat2MatlabDoubles(&msg->nav_status,1,1));
    
    if(!msg->rot_over_range) {
       mxSetFieldByNumber(theStruct,0,4,floatMat2MatlabDoubles(&msg->rot,1,1));
    } else {//If the rotation rate is invalid.
       mxSetFieldByNumber(theStruct,0,4,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(msg->sog!=1023) {
        mxSetFieldByNumber(theStruct,0,5,floatMat2MatlabDoubles(&msg->sog,1,1));
    } else {//Speed over ground is not available.
        mxSetFieldByNumber(theStruct,0,5,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    mxSetFieldByNumber(theStruct,0,6,intMat2MatlabDoubles(&msg->position_accuracy,1,1));
    if(fabs(msg->position.lng_deg)<=180) {
        mxSetFieldByNumber(theStruct,0,7,doubleMat2Matlab(&msg->position.lng_deg,1,1));
    } else {//If longitude is not available
        mxSetFieldByNumber(theStruct,0,7,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(fabs(msg->position.lat_deg)<=90) {
        mxSetFieldByNumber(theStruct,0,8,doubleMat2Matlab(&msg->position.lat_deg,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,8,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(msg->cog<360) {
        mxSetFieldByNumber(theStruct,0,9,floatMat2MatlabDoubles(&msg->cog,1,1));
    } else {//Course over ground is not available
        mxSetFieldByNumber(theStruct,0,9,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(msg->true_heading<360) {
        mxSetFieldByNumber(theStruct,0,10,intMat2MatlabDoubles(&msg->true_heading,1,1));
    } else {//True Heading is not available
        mxSetFieldByNumber(theStruct,0,10,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(msg->timestamp!=60) {
        mxSetFieldByNumber(theStruct,0,11,intMat2MatlabDoubles(&msg->timestamp,1,1));  
    } else {//If no time stamp is available.
        mxSetFieldByNumber(theStruct,0,11,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(msg->special_manoeuvre!=0) {
        mxSetFieldByNumber(theStruct,0,12,intMat2MatlabDoubles(&msg->special_manoeuvre,1,1)); 
    } else {//Special maneuver information is not available
        mxSetFieldByNumber(theStruct,0,12,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    mxSetFieldByNumber(theStruct,0,13,intMat2MatlabDoubles(&msg->spare,1,1)); 
    mxSetFieldByNumber(theStruct,0,14,boolMat2Matlab(&msg->raim,1,1));

    //COMM States
    mxSetFieldByNumber(theStruct,0,15,intMat2MatlabDoubles(&msg->sync_state,1,1)); 
    //SODATA
    if(msg->slot_timeout_valid) {
        mxSetFieldByNumber(theStruct,0,16,intMat2MatlabDoubles(&msg->slot_timeout,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,16,mxCreateDoubleMatrix(0,0,mxREAL));
    }

    if(msg->received_stations_valid) {
        mxSetFieldByNumber(theStruct,0,17,intMat2MatlabDoubles(&msg->received_stations,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,17,mxCreateDoubleMatrix(0,0,mxREAL));
    }

    if(msg->slot_number_valid) {
        mxSetFieldByNumber(theStruct,0,18,intMat2MatlabDoubles(&msg->slot_number,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,18,mxCreateDoubleMatrix(0,0,mxREAL));
    }

    if(msg->utc_valid) {
        mxSetFieldByNumber(theStruct,0,19,intMat2MatlabDoubles(&msg->utc_hour,1,1));
        mxSetFieldByNumber(theStruct,0,20,intMat2MatlabDoubles(&msg->utc_min,1,1));
        mxSetFieldByNumber(theStruct,0,21,intMat2MatlabDoubles(&msg->utc_spare,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,19,mxCreateDoubleMatrix(0,0,mxREAL));
        mxSetFieldByNumber(theStruct,0,20,mxCreateDoubleMatrix(0,0,mxREAL));
        mxSetFieldByNumber(theStruct,0,21,mxCreateDoubleMatrix(0,0,mxREAL));
    }

    if(msg->slot_offset_valid) {
        mxSetFieldByNumber(theStruct,0,22,intMat2MatlabDoubles(&msg->slot_offset,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,22,mxCreateDoubleMatrix(0,0,mxREAL));
    }

    //ITDMA
    if(msg->slot_increment_valid) {
        mxSetFieldByNumber(theStruct,0,23,intMat2MatlabDoubles(&msg->slot_increment,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,23,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(msg->slots_to_allocate_valid) {
        mxSetFieldByNumber(theStruct,0,24,intMat2MatlabDoubles(&msg->slots_to_allocate,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,24,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(msg->keep_flag_valid) {
        mxSetFieldByNumber(theStruct,0,25,boolMat2Matlab(&msg->keep_flag,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,25,mxCreateDoubleMatrix(0,0,mxREAL));
    }
}

void AIS_4_11_ToMatlab(libais::Ais4_11 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions) {
    const size_t numberOfFields=24;
    const mwSize dims[2] = {1, 1};
    //These are the names of all of the fields that will be added to the
    //Matlab structre array.
    const char *fieldNames[numberOfFields] = {"message_id", "repeat_indicator",
    "mmsi", "year", "month", "day", "hour", "minute", "second",
    "position_accuracy", "lon", "lat", "fix_type", "transmission_ctl", "spare",
    "raim", "sync_state", "slot_timeout", "received_stations",
    "slot_number", "utc_hour", "utc_min","utc_spare", "slot_offset"};
    mxArray *theStruct;
    
    //If text descrptions of the fields are desired, then put them into a
    //structure.
    if(fieldDescriptions!=NULL) {
        const char *descriptionStrings[numberOfFields] = {
            "Identifier for this message\n4 = UTC and position report from base station\n11 = UTC and position response from mobile station",//0
            "Used by the repeater to indicate how many times a message has been repeated. 0 = default; 3 = do not repeat any more",//1
            "Maritime mobile service identity (MMSI) number",//2
            "UTC year (1-9999)",//3
            "UTC month (1-12)",//4
            "UTC day (1-31)",//5
            "UTC hour (0-23)",//6
            "UTC minute (0-59)",//7
            "UTC second (0-59)",//8
            "Position accuracy. 1=high (<=10m), 0=low (>10m), 0=default",//9
            "Longitude East in degrees (WGS-84)",//10
            "Latitude North in degrees (WGS-84)",//11
            "Type of electronic position fixing device. Use of differential corrections is defined by field position accuracy above:\n1 = global positioning system (GPS)\n2 = GNSS (GLONASS)\n3 = combined GPS/GLONASS\n4 = Loran-C\n5 = Chayka\n6 = integrated navigation system\n7 = surveyed\n8 = Galileo\n9-14 = not used\n15 = internal GNSS",//12
            "Transmission control for long- range broadcast message\n0 = default Class-A AIS station stops transmission of Message 27 within an AIS base station coverage area.\n1 = Request Class-A station to transmit Message 27 within an AIS base station coverage area.",//13
            "Not used. Should be set to zero. Reserved for future use",//14
            "Receiver autonomous integrity monitoring (RAIM) flag of electronic position fixing device; 0 = RAIM not in use = default; 1 = RAIM in use.",//15
            "Synchronization state\n0 UTC direct\n1 UTC indirect\n2 Station is synchronized to a base station (base direct)\n3 Station is synchronized to another station based on the highest number of received stations or to another mobile station, which is directly synchronized to a base station",//16
            "(For Message IDs 1 and 2) Specifies frames remaining until a new slot is selected\n0 means that this was the last transmission in this slot\n1-7 means that 1 to 7 frames respectively are left until slot change",//17
            "Number of other stations (not own station) which the station currently is receiving (between 0 and 16 383).",//18
            "Slot number used for this transmission (between 0 and 2 249).",//19
            "Hour (0-23) in UTC time, if available",//20
            "Minute (0-23) in UTC time, if available",//21
            "Two extra (unused) bits that are associated with the time message",//22
            "If the slot time-out value is 0 (zero) then the slot offset should indicate the offset to the slot in which transmission will occur during the next frame. If the slot offset is zero, the slot should be de-allocated after transmission.",//23
        };
        unsigned int i;

        *fieldDescriptions=mxCreateStructArray(2, dims, numberOfFields, fieldNames);

        for(i=0;i<numberOfFields;i++) {
            mxSetFieldByNumber(*fieldDescriptions,0,static_cast<int>(i),mxCreateCharMatrixFromStrings(1,&descriptionStrings[i]));
        }
    }
    
    //Create the Matlab structure array.
    theStruct=mxCreateStructArray(2, dims, numberOfFields, fieldNames);
    *decodedMessage=theStruct;
    
    //Fill all of the elements of the structure
    mxSetFieldByNumber(theStruct,0,0,intMat2MatlabDoubles(&msg->message_id,1,1));
    mxSetFieldByNumber(theStruct,0,1,intMat2MatlabDoubles(&msg->repeat_indicator,1,1));
    mxSetFieldByNumber(theStruct,0,2,intMat2MatlabDoubles(&msg->mmsi,1,1));
    if(msg->year!=0) {
        mxSetFieldByNumber(theStruct,0,3,intMat2MatlabDoubles(&msg->year,1,1));
    } else {//UTC year not available
        mxSetFieldByNumber(theStruct,0,3,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(msg->month!=0) {
        mxSetFieldByNumber(theStruct,0,4,intMat2MatlabDoubles(&msg->month,1,1));
    } else {//UTC month not available
        mxSetFieldByNumber(theStruct,0,4,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(msg->day!=0) {
        mxSetFieldByNumber(theStruct,0,5,intMat2MatlabDoubles(&msg->day,1,1));
    } else {//UTC day not available
        mxSetFieldByNumber(theStruct,0,5,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(msg->hour<24) {
        mxSetFieldByNumber(theStruct,0,6,intMat2MatlabDoubles(&msg->hour,1,1));
    } else {//UTC hour invalid or not available
        mxSetFieldByNumber(theStruct,0,6,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(msg->minute<60) {
        mxSetFieldByNumber(theStruct,0,7,intMat2MatlabDoubles(&msg->minute,1,1));
    } else {//UTC minute invalid or not available
        mxSetFieldByNumber(theStruct,0,7,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(msg->second<60) {
        mxSetFieldByNumber(theStruct,0,8,intMat2MatlabDoubles(&msg->second,1,1));
    } else {//UTC second invalid or not available
        mxSetFieldByNumber(theStruct,0,8,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    mxSetFieldByNumber(theStruct,0,9,intMat2MatlabDoubles(&msg->position_accuracy,1,1));
    
    if(fabs(msg->position.lng_deg)<=180) {
        mxSetFieldByNumber(theStruct,0,10,doubleMat2Matlab(&msg->position.lng_deg,1,1));
    } else {//If longitude is not available
        mxSetFieldByNumber(theStruct,0,10,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(fabs(msg->position.lat_deg)<=90) {
        mxSetFieldByNumber(theStruct,0,11,doubleMat2Matlab(&msg->position.lat_deg,1,1));
    } else {//Latitude is not available
        mxSetFieldByNumber(theStruct,0,11,mxCreateDoubleMatrix(0,0,mxREAL));
    }

    if(msg->fix_type!=0) {
        mxSetFieldByNumber(theStruct,0,12,intMat2MatlabDoubles(&msg->fix_type,1,1));
    } else {//Fix type is unavailable
        mxSetFieldByNumber(theStruct,0,12,mxCreateDoubleMatrix(0,0,mxREAL));
    }

    mxSetFieldByNumber(theStruct,0,13,intMat2MatlabDoubles(&msg->transmission_ctl,1,1));
    mxSetFieldByNumber(theStruct,0,14,intMat2MatlabDoubles(&msg->spare,1,1));
    mxSetFieldByNumber(theStruct,0,15,boolMat2Matlab(&msg->raim,1,1));
    
    //SOTDMA
    mxSetFieldByNumber(theStruct,0,16,intMat2MatlabDoubles(&msg->sync_state,1,1));
    mxSetFieldByNumber(theStruct,0,17,intMat2MatlabDoubles(&msg->slot_timeout,1,1));

    if(msg->received_stations_valid) {
        mxSetFieldByNumber(theStruct,0,18,intMat2MatlabDoubles(&msg->received_stations,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,18,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(msg->slot_number_valid) {
        mxSetFieldByNumber(theStruct,0,19,intMat2MatlabDoubles(&msg->slot_number,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,19,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(msg->utc_valid) {
        mxSetFieldByNumber(theStruct,0,20,intMat2MatlabDoubles(&msg->utc_hour,1,1));
        mxSetFieldByNumber(theStruct,0,21,intMat2MatlabDoubles(&msg->utc_min,1,1));
        mxSetFieldByNumber(theStruct,0,22,intMat2MatlabDoubles(&msg->utc_spare,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,20,mxCreateDoubleMatrix(0,0,mxREAL)); 
        mxSetFieldByNumber(theStruct,0,21,mxCreateDoubleMatrix(0,0,mxREAL));
        mxSetFieldByNumber(theStruct,0,22,mxCreateDoubleMatrix(0,0,mxREAL));
    }

    if(msg->slot_offset_valid) {
        mxSetFieldByNumber(theStruct,0,22,intMat2MatlabDoubles(&msg->slot_offset,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,23,mxCreateDoubleMatrix(0,0,mxREAL));
    }
}

void AIS_5_ToMatlab(libais::Ais5 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions) {
    const size_t numberOfFields=21;
    const mwSize dims[2] = {1, 1};
    //These are the names of all of the fields that will be added to the
    //Matlab structre array.
    const char *fieldNames[numberOfFields] = {"message_id", "repeat_indicator",
    "mmsi", "ais_version", "imo_num", "callsign", "name", "type_and_cargo",
    "dim_a", "dim_b", "dim_c", "dim_d", "fix_type", "eta_month", "eta_day",
    "eta_hour", "eta_minute", "draught", "destination", "dte", "spare"};
    mxArray *theStruct;
    
    //If text descriptions of the fields are desired, then put them into a
    //structure.
    if(fieldDescriptions!=NULL) {
        const char *descriptionStrings[numberOfFields] = {
             "Identifier for this message\n4 = UTC and position report from base station\n11 = UTC and position response from mobile station",//0
            "Used by the repeater to indicate how many times a message has been repeated. 0 = default; 3 = do not repeat any more",//1
            "Maritime mobile service identity (MMSI) number",//2
            "AIS version indicator\n0 = station compliant with Recommendation ITU-R M.1371-1\n1 = station compliant with Recommendation ITU-R M.1371-3 (or later)\n2 = station compliant with Recommendation ITU-R M.1371-5 (or later)\n3 = station compliant with future editions",//3
            "IMO Number\n0000000001-0000999999 not used\n0001000000-0009999999 = valid IMO number;\n0010000000-1073741823 = official flag state number.",//4
            "Call sign. Craft associated with a parent vessel, should use 'A' followed by the last 6 digits of the MMSI of the parent vessel. Examples of these craft include towed vessels, rescue boats, tenders, lifeboats and liferafts.",//5
            "Name. Maximum 20 characters. The Name should be as shown on the station radio license. For search and rescue (SAR) aircraft, it should be set to 'SAR AIRCRAFT NNNNNNN' where NNNNNNN equals the aircraft registration number.",//6
            "Type of ship and cargo type\n1-99 = as defined in Section 3.3.2 of Annex 8 of ITU-R M.1371-5\n100-199 = reserved, for regional use 200-255 = reserved, for future use\nNot applicable to search and rescue (SAR) aircraft",//7
            "Distance from reference point to bow of ship (meters) 511m means 511m or longer. (If only relative reference point position known=0)",//8
            "Distance from reference point to aft of ship (meters) 511m means 511m or longer. ",//9
            "Distance from reference point to port side of ship (meters) 511m means 511m or longer. (If only relative reference point position known=0)",//10
            "Distance from reference point to starboard side of ship (meters) 511m means 511m or longer.",//11
            "Type of electronic position fixing device:\n1 = global positioning system (GPS)\n2 = GNSS (GLONASS)\n3 = combined GPS/GLONASS\n4 = Loran-C\n5 = Chayka\n6 = integrated navigation system\n7 = surveyed\n8 = Galileo\n9-14 = not used\n15 = internal GNSS",//12
            "Estimated UTC month of arrival (1-12)",//13
            "Estimated UTC day of arrival (1-31)",//14
            "Estimated UTC hour of arrival (0-23)",//15
            "Estimated UTC minute of arrival (0-59)",//16
            "Maximum present static draught in meters. 255m means 255 meters or greater. Not applicable to search and rescue (SAR) aircraft",//17
            "Destination, maximum 20 characters. For search and rescue (SAR) aircraft, the use of this field may be decided by the responsible administration",//18
            "Data terminal equipment (DTE) ready (0 = available, 1 = not available = default)",//19
            "Spare. Not used. Should be set to zero. Reserved for future use"//20
        };
        unsigned int i;

        *fieldDescriptions=mxCreateStructArray(2, dims, numberOfFields, fieldNames);

        for(i=0;i<numberOfFields;i++) {
            mxSetFieldByNumber(*fieldDescriptions,0,static_cast<int>(i),mxCreateCharMatrixFromStrings(1,&descriptionStrings[i]));
        }
    }
    
    //Create the Matlab structure array.
    theStruct=mxCreateStructArray(2, dims, numberOfFields, fieldNames);
    *decodedMessage=theStruct;
    
    //Fill all of the elements of the structure
    mxSetFieldByNumber(theStruct,0,0,intMat2MatlabDoubles(&msg->message_id,1,1));
    mxSetFieldByNumber(theStruct,0,1,intMat2MatlabDoubles(&msg->repeat_indicator,1,1));
    mxSetFieldByNumber(theStruct,0,2,intMat2MatlabDoubles(&msg->mmsi,1,1));
    mxSetFieldByNumber(theStruct,0,3,intMat2MatlabDoubles(&msg->ais_version,1,1));
    
    if(msg->imo_num!=0) {
        mxSetFieldByNumber(theStruct,0,4,intMat2MatlabDoubles(&msg->imo_num,1,1));
    } else {//IMO number not available
        mxSetFieldByNumber(theStruct,0,4,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(msg->callsign.compare("@@@@@@@")!=0) {
        const char *charString=msg->callsign.c_str();
        mxSetFieldByNumber(theStruct,0,5,mxCreateCharMatrixFromStrings(1,&charString));
    } else {//If the call sign is not available
        mxSetFieldByNumber(theStruct,0,5,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(msg->name.compare("@@@@@@@@@@@@@@@@@@@@")!=0) {
        const char *charString=msg->name.c_str();
        mxSetFieldByNumber(theStruct,0,6,mxCreateCharMatrixFromStrings(1,&charString));
    } else {//If the name is not available
        mxSetFieldByNumber(theStruct,0,6,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(msg->type_and_cargo!=0) {
        mxSetFieldByNumber(theStruct,0,7,intMat2MatlabDoubles(&msg->type_and_cargo,1,1));
    } else {//No type and cargo information
        mxSetFieldByNumber(theStruct,0,6,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(msg->dim_a==0&&msg->dim_b==0&&msg->dim_c==0&&msg->dim_d==0) {
        mxSetFieldByNumber(theStruct,0,8,mxCreateDoubleMatrix(0,0,mxREAL));
        mxSetFieldByNumber(theStruct,0,9,mxCreateDoubleMatrix(0,0,mxREAL));
        mxSetFieldByNumber(theStruct,0,10,mxCreateDoubleMatrix(0,0,mxREAL));
        mxSetFieldByNumber(theStruct,0,11,mxCreateDoubleMatrix(0,0,mxREAL));
    } else {//If ship dimensions from reference point are available.
        mxSetFieldByNumber(theStruct,0,8,intMat2MatlabDoubles(&msg->dim_a,1,1));
        mxSetFieldByNumber(theStruct,0,9,intMat2MatlabDoubles(&msg->dim_b,1,1));
        mxSetFieldByNumber(theStruct,0,10,intMat2MatlabDoubles(&msg->dim_c,1,1));
        mxSetFieldByNumber(theStruct,0,11,intMat2MatlabDoubles(&msg->dim_d,1,1));
    }
    
    if(msg->fix_type!=0) {
        mxSetFieldByNumber(theStruct,0,12,intMat2MatlabDoubles(&msg->fix_type,1,1));
    } else {//If the fix type is not available
        mxSetFieldByNumber(theStruct,0,12,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(msg->eta_month!=0) {
        mxSetFieldByNumber(theStruct,0,13,intMat2MatlabDoubles(&msg->eta_month,1,1));
    } else {//No eta month
        mxSetFieldByNumber(theStruct,0,13,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(msg->eta_day!=0) {
        mxSetFieldByNumber(theStruct,0,14,intMat2MatlabDoubles(&msg->eta_day,1,1));
    } else {//No eta day
        mxSetFieldByNumber(theStruct,0,14,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(msg->eta_hour<24) {
        mxSetFieldByNumber(theStruct,0,15,intMat2MatlabDoubles(&msg->eta_hour,1,1));
    } else {//No eta hour
        mxSetFieldByNumber(theStruct,0,15,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(msg->eta_minute<60) {
        mxSetFieldByNumber(theStruct,0,16,intMat2MatlabDoubles(&msg->eta_minute,1,1));
    } else {//No eta minute
        mxSetFieldByNumber(theStruct,0,16,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(msg->draught!=0) {
        mxSetFieldByNumber(theStruct,0,17,floatMat2MatlabDoubles(&msg->draught,1,1));
    } else {//No draught available
        mxSetFieldByNumber(theStruct,0,17,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    if(msg->destination.compare("@@@@@@@@@@@@@@@@@@@@")!=0) {
        const char *charString=msg->destination.c_str();
        mxSetFieldByNumber(theStruct,0,18,mxCreateCharMatrixFromStrings(1,&charString));
    } else {
        mxSetFieldByNumber(theStruct,0,18,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    mxSetFieldByNumber(theStruct,0,19,intMat2MatlabDoubles(&msg->dte,1,1));
    mxSetFieldByNumber(theStruct,0,20,intMat2MatlabDoubles(&msg->spare,1,1));
}


void AIS_7_13_ToMatlab(libais::Ais7_13 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions) {
    const size_t numberOfFields=6;
    const mwSize dims[2] = {1, 1};
    //These are the names of all of the fields that will be added to the
    //Matlab structre array.
    const char *fieldNames[numberOfFields] = {"message_id", "repeat_indicator",
    "mmsi", "dest_mmsi", "seq_num", "spare"};
    mxArray *theStruct;
    
    //If text descriptions of the fields are desired, then put them into a
    //structure.
    if(fieldDescriptions!=NULL) {
        const char *descriptionStrings[numberOfFields] = {
            "Identifier for Messages 7 or 13\n7 = binary acknowledge\n13 = safety related acknowledge",//0
            "Used by the repeater to indicate how many times a message has been repeated. 0 = default; 3 = do not repeat any more",//1
            "Maritime mobile service identity (MMSI) number of source of this acknowledge (ACK)",//2
            "An array of the destination MMSIs of the destinations of this ACK",//3
            "An array of the sequence numbers of the messages to be acknowledged for each destination (0-3)",//4
            "Spare bits. Not used. Should be set to zero. Reserved for future use"//5
        };
        unsigned int i;

        *fieldDescriptions=mxCreateStructArray(2, dims, numberOfFields, fieldNames);

        for(i=0;i<numberOfFields;i++) {
            mxSetFieldByNumber(*fieldDescriptions,0,static_cast<int>(i),mxCreateCharMatrixFromStrings(1,&descriptionStrings[i]));
        }
    }
    
    //Create the Matlab structure array.
    theStruct=mxCreateStructArray(2, dims, numberOfFields, fieldNames);
    *decodedMessage=theStruct;
    
    //Fill all of the elements of the structure
    mxSetFieldByNumber(theStruct,0,0,intMat2MatlabDoubles(&msg->message_id,1,1));
    mxSetFieldByNumber(theStruct,0,1,intMat2MatlabDoubles(&msg->repeat_indicator,1,1));
    mxSetFieldByNumber(theStruct,0,2,intMat2MatlabDoubles(&msg->mmsi,1,1));
    {
        int *theData;
        theData=msg->dest_mmsi.data();
        mxSetFieldByNumber(theStruct,0,3,intMat2MatlabDoubles(theData,msg->dest_mmsi.size(),1));
        theData=msg->seq_num.data();
        mxSetFieldByNumber(theStruct,0,4,intMat2MatlabDoubles(theData,msg->seq_num.size(),1));
    }
    mxSetFieldByNumber(theStruct,0,5,intMat2MatlabDoubles(&msg->spare,1,1)); 
}

void AIS_9_ToMatlab(libais::Ais9 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions) {
    const size_t numberOfFields=28;
    const mwSize dims[2] = {1, 1};
    //These are the names of all of the fields that will be added to the
    //Matlab structre array.
    const char *fieldNames[numberOfFields] = {"message_id", "repeat_indicator",
    "mmsi", "alt", "sog", "position_accuracy", "lon", "lat", "cog",
    "timestamp", "alt_sensor", "spare", "dte", "spare2", "assigned_mode",
    "raim", "commstate_flag", "sync_state", "slot_timeout",
    "received_stations", "slot_number", "utc_hour", "utc_min", "utc_spare",
    "slot_offset", "slot_increment", "slots_to_allocate", "keep_flag"};
    mxArray *theStruct;
    
    //If text descriptions of the fields are desired, then put them into a
    //structure.
    if(fieldDescriptions!=NULL) {
        const char *descriptionStrings[numberOfFields] = {
        	"Message identifier.",//0
            "Repeat indicator. Used by the repeater to indicate how many times a message has been repeated. 0 = default; 3 = do not repeat any more.",//1
            "Maritime mobile service identity (MMSI) number",//2
            "Altitude (derived from GNSS or barometric (see altitude sensor parameter below)) (m) (0-4 094 m) 4094 = 4094 m or higher",//3
            "Speed over ground in 1/10 knots. 1023 means 102.2 knots or higher",//4
            "Position accuracy. 1=high (<=10m), 0=low (>10m), 0=default",//5
            "Longitude East in degrees (WGS-84)",//6
            "Latitude North in degrees (WGS-84)",//7
            "Course over ground in degrees East of North",//8
            "UTC second when report was generated or 61 if positioning system is in manual input mode or 62 if in dead reckining mode or 63 if positioning system is inoperative",//9
            "Altitude sensor. 0 = GNSS, 1 = barometric source",//10
            "Not used. Should be set to zero. Reserved for future use",//11
            "Data terminal ready (0 = available 1 = not available = default)",//12
            "Not used. Should be set to zero. Reserved for future use",//13
            "0 = Station operating in autonomous and continuous mode = default 1 = Station operating in assigned mode",//14
            "Receiver autonomous integrity monitoring (RAIM) flag of electronic position fixing device; 0 = RAIM not in use = default; 1 = RAIM in use.",//15
            "0 = SOTDMA communication state follows 1 = ITDMA communication state follows",//16
            "Synchronization state\n0 UTC direct\n1 UTC indirect\n2 Station is synchronized to a base station (base direct)\n3 Station is synchronized to another station based on the highest number of received stations or to another mobile station, which is directly synchronized to a base station",//17
            "(For SOTDMA) Specifies frames remaining until a new slot is selected\n0 means that this was the last transmission in this slot\n1-7 means that 1 to 7 frames respectively are left until slot change",//18
            "(For SOTDMA) Number of other stations (not own station) which the station currently is receiving (between 0 and 16 383).",//19
            "(For SOTDMA) Slot number used for this transmission (between 0 and 2 249).",//20
            "(For SOTDMA) Hour (0-23) in UTC time, if available",//21
            "(For SOTDMA) Minute (0-23) in UTC time, if available",//22
            "(For SOTDMA) Two extra (unused) bits that are associated with the time message",//23
            "(For SOTDMA) If the slot time-out value is 0 (zero) then the slot offset should indicate the offset to the slot in which transmission will occur during the next frame. If the slot offset is zero, the slot should be de-allocated after transmission.",//24
            "(For ITDMA) Offset to next slot to be used, or zero (0) if no more transmissions",//25
            "(For ITDMA) Number of consecutive slots to allocate.\n0 = 1 slot,\n1 = 2 slots,\n2 = 3 slots,\n3 = 4 slots,\n4 = 5 slots,\n5 = 1 slot; offset = slot increment + 8 192,\n6 = 2 slots; offset = slot increment + 8 192,\n7 = 3 slots; offset = slot increment + 8 192.\nUse of 5 to 7 removes the need for RATDMA broadcast for scheduled transmissions up to 6 min intervals",//26
            "(For ITDMA) Set to TRUE = 1 if the slot remains allocated for one additional frame"//27
        };
        unsigned int i;

        *fieldDescriptions=mxCreateStructArray(2, dims, numberOfFields, fieldNames);

        for(i=0;i<numberOfFields;i++) {
            mxSetFieldByNumber(*fieldDescriptions,0,static_cast<int>(i),mxCreateCharMatrixFromStrings(1,&descriptionStrings[i]));
        }
    }
    
    //Create the Matlab structure array.
    theStruct=mxCreateStructArray(2, dims, numberOfFields, fieldNames);
    *decodedMessage=theStruct;
    
    //Fill all of the elements of the structure
    mxSetFieldByNumber(theStruct,0,0,intMat2MatlabDoubles(&msg->message_id,1,1));
    mxSetFieldByNumber(theStruct,0,1,intMat2MatlabDoubles(&msg->repeat_indicator,1,1));
    mxSetFieldByNumber(theStruct,0,2,intMat2MatlabDoubles(&msg->mmsi,1,1));
    
    if(msg->alt!=4095) {
        mxSetFieldByNumber(theStruct,0,3,intMat2MatlabDoubles(&msg->alt,1,1)); 
    } else {//No altitude information available.
        mxSetFieldByNumber(theStruct,0,3,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(msg->sog!=1023) {
        mxSetFieldByNumber(theStruct,0,4,intMat2MatlabDoubles(&msg->sog,1,1));
    } else {//No speed over ground information available
        mxSetFieldByNumber(theStruct,0,4,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    mxSetFieldByNumber(theStruct,0,5,intMat2MatlabDoubles(&msg->position_accuracy,1,1));
    
    if(fabs(msg->position.lng_deg)<=180) {
        mxSetFieldByNumber(theStruct,0,6,doubleMat2Matlab(&msg->position.lng_deg,1,1));
    } else {//If the longitude is invalid
        mxSetFieldByNumber(theStruct,0,6,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(fabs(msg->position.lat_deg)<=90) {
        mxSetFieldByNumber(theStruct,0,7,doubleMat2Matlab(&msg->position.lat_deg,1,1));
    } else {// If the latitude is invalid
        mxSetFieldByNumber(theStruct,0,7,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(msg->cog<360) {
        mxSetFieldByNumber(theStruct,0,8,floatMat2MatlabDoubles(&msg->cog,1,1));
    } else {//If course over ground is invalid
        mxSetFieldByNumber(theStruct,0,8,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(msg->timestamp!=60) {
        mxSetFieldByNumber(theStruct,0,9,intMat2MatlabDoubles(&msg->timestamp,1,1)); 
    } else {//If the timestamp is invalid
        mxSetFieldByNumber(theStruct,0,9,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    mxSetFieldByNumber(theStruct,0,10,intMat2MatlabDoubles(&msg->alt_sensor,1,1)); 
    mxSetFieldByNumber(theStruct,0,11,intMat2MatlabDoubles(&msg->spare,1,1)); 
    mxSetFieldByNumber(theStruct,0,12,intMat2MatlabDoubles(&msg->dte,1,1)); 
    mxSetFieldByNumber(theStruct,0,13,intMat2MatlabDoubles(&msg->spare2,1,1)); 
    mxSetFieldByNumber(theStruct,0,14,intMat2MatlabDoubles(&msg->assigned_mode,1,1)); 
    mxSetFieldByNumber(theStruct,0,15,boolMat2Matlab(&msg->raim,1,1));
    mxSetFieldByNumber(theStruct,0,16,intMat2MatlabDoubles(&msg->commstate_flag,1,1));    
    mxSetFieldByNumber(theStruct,0,17,intMat2MatlabDoubles(&msg->sync_state,1,1));
    
    if(msg->slot_timeout_valid) {
        mxSetFieldByNumber(theStruct,0,18,intMat2MatlabDoubles(&msg->slot_timeout,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,18,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(msg->received_stations_valid) {
        mxSetFieldByNumber(theStruct,0,19,intMat2MatlabDoubles(&msg->received_stations,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,19,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(msg->slot_number_valid) {
        mxSetFieldByNumber(theStruct,0,20,intMat2MatlabDoubles(&msg->slot_number,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,20,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(msg->utc_valid) {
        mxSetFieldByNumber(theStruct,0,21,intMat2MatlabDoubles(&msg->utc_hour,1,1));
        mxSetFieldByNumber(theStruct,0,22,intMat2MatlabDoubles(&msg->utc_min,1,1));
        mxSetFieldByNumber(theStruct,0,23,intMat2MatlabDoubles(&msg->utc_spare,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,21,mxCreateDoubleMatrix(0,0,mxREAL));
        mxSetFieldByNumber(theStruct,0,22,mxCreateDoubleMatrix(0,0,mxREAL));
        mxSetFieldByNumber(theStruct,0,23,mxCreateDoubleMatrix(0,0,mxREAL));
    }

    if(msg->slot_offset_valid) {
        mxSetFieldByNumber(theStruct,0,24,intMat2MatlabDoubles(&msg->slot_offset,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,24,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(msg->slot_increment_valid) {
        mxSetFieldByNumber(theStruct,0,25,intMat2MatlabDoubles(&msg->slot_increment,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,25,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(msg->slots_to_allocate_valid) {
        mxSetFieldByNumber(theStruct,0,26,intMat2MatlabDoubles(&msg->slots_to_allocate,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,26,mxCreateDoubleMatrix(0,0,mxREAL));
    }

    if(msg->keep_flag_valid) {
        mxSetFieldByNumber(theStruct,0,27,boolMat2Matlab(&msg->keep_flag,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,27,mxCreateDoubleMatrix(0,0,mxREAL));
    }
}

void AIS_10_ToMatlab(libais::Ais10 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions) {
    const size_t numberOfFields=6;
    const mwSize dims[2] = {1, 1};
    //These are the names of all of the fields that will be added to the
    //Matlab structre array.
    const char *fieldNames[numberOfFields] = {"message_id", "repeat_indicator",
    "mmsi", "spare", "dest_mmsi", "spare2"};
    mxArray *theStruct;

    //If text descriptions of the fields are desired, then put them into a
    //structure.
    if(fieldDescriptions!=NULL) {
        const char *descriptionStrings[numberOfFields] = {
            "Identifier for Message 10; always 10.",//0
            "Repeat indicator. Used by the repeater to indicate how many times a message has been repeated. 0 = default; 3 = do not repeat any more.",//1
            "Maritime mobile service identity (MMSI) number of station which inquires UTC",//2
            "Not used. Should be set to zero. Reserved for future use",//3
            "MMSI number of station which is inquired",//4
            "Not used. Should be set to zero. Reserved for future use"//5
        };
        unsigned int i;

        *fieldDescriptions=mxCreateStructArray(2, dims, numberOfFields, fieldNames);

        for(i=0;i<numberOfFields;i++) {
            mxSetFieldByNumber(*fieldDescriptions,0,static_cast<int>(i),mxCreateCharMatrixFromStrings(1,&descriptionStrings[i]));
        }
    }

    //Create the Matlab structure array.
    theStruct=mxCreateStructArray(2, dims, numberOfFields, fieldNames);
    *decodedMessage=theStruct;
    
    //Fill all of the elements of the structure
    mxSetFieldByNumber(theStruct,0,0,intMat2MatlabDoubles(&msg->message_id,1,1));
    mxSetFieldByNumber(theStruct,0,1,intMat2MatlabDoubles(&msg->repeat_indicator,1,1));
    mxSetFieldByNumber(theStruct,0,2,intMat2MatlabDoubles(&msg->mmsi,1,1));
    mxSetFieldByNumber(theStruct,0,3,intMat2MatlabDoubles(&msg->spare,1,1));
    mxSetFieldByNumber(theStruct,0,4,intMat2MatlabDoubles(&msg->dest_mmsi,1,1)); 
    mxSetFieldByNumber(theStruct,0,5,intMat2MatlabDoubles(&msg->spare2,1,1)); 
}

void AIS_12_ToMatlab(libais::Ais12 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions) {
    const size_t numberOfFields=9;
    const mwSize dims[2] = {1, 1};
    //These are the names of all of the fields that will be added to the
    //Matlab structre array.
    const char *fieldNames[numberOfFields] = {"message_id", "repeat_indicator",
    "mmsi", "seq_num", "dest_mmsi", "retransmitted", "spare", "text","spare2"};
    mxArray *theStruct;
    
    //If text descriptions of the fields are desired, then put them into a
    //structure.
    if(fieldDescriptions!=NULL) {
        const char *descriptionStrings[numberOfFields] = {
            "Identifier for Message 12; always 12.",//0
            "Repeat indicator. Used by the repeater to indicate how many times a message has been repeated. 0 = default; 3 = do not repeat any more.",//1
            "Maritime mobile service identity (MMSI) number of station which is the source of the message.",//2
            "Sequence number (0-3) as described in Section 5.3.1, Annex 2 of ITU-R M.1371-5",//3
            "MMSI number of station which is the destination of the message",//4
            "Retransmit flag should be set upon retransmission: 0 = no retransmission = default; 1 = retransmitted",//5
            "Not used. Should be zero. Reserved for future use",//6
            "Safety related text.",//7
            "Extra bits for alignment; not used.",//8
        };
        unsigned int i;

        *fieldDescriptions=mxCreateStructArray(2, dims, numberOfFields, fieldNames);

        for(i=0;i<numberOfFields;i++) {
            mxSetFieldByNumber(*fieldDescriptions,0,static_cast<int>(i),mxCreateCharMatrixFromStrings(1,&descriptionStrings[i]));
        }
    }
    
    //Create the Matlab structure array.
    theStruct=mxCreateStructArray(2, dims, numberOfFields, fieldNames);
    *decodedMessage=theStruct;

    //Fill all of the elements of the structure
    mxSetFieldByNumber(theStruct,0,0,intMat2MatlabDoubles(&msg->message_id,1,1));
    mxSetFieldByNumber(theStruct,0,1,intMat2MatlabDoubles(&msg->repeat_indicator,1,1));
    mxSetFieldByNumber(theStruct,0,2,intMat2MatlabDoubles(&msg->mmsi,1,1));
    mxSetFieldByNumber(theStruct,0,3,intMat2MatlabDoubles(&msg->seq_num,1,1));
    mxSetFieldByNumber(theStruct,0,4,intMat2MatlabDoubles(&msg->dest_mmsi,1,1));
    mxSetFieldByNumber(theStruct,0,5,boolMat2Matlab(&msg->retransmitted,1,1));
    mxSetFieldByNumber(theStruct,0,6,intMat2MatlabDoubles(&msg->spare,1,1));
    {
        const char *charString=msg->text.c_str();
        mxSetFieldByNumber(theStruct,0,7,mxCreateCharMatrixFromStrings(1,&charString));
    }
    mxSetFieldByNumber(theStruct,0,8,intMat2MatlabDoubles(&msg->spare2,1,1));
}

void AIS_14_ToMatlab(libais::Ais14 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions) {
    const size_t numberOfFields=6;
    const mwSize dims[2] = {1, 1};
    //These are the names of all of the fields that will be added to the
    //Matlab structre array.
    const char *fieldNames[numberOfFields] = {"message_id", "repeat_indicator",
    "mmsi", "spare", "text","spare2"};
    mxArray *theStruct;

    //If text descriptions of the fields are desired, then put them into a
    //structure.
    if(fieldDescriptions!=NULL) {
        const char *descriptionStrings[numberOfFields] = {
            "Identifier for Message 14; always 14.",//0
            "Repeat indicator. Used by the repeater to indicate how many times a message has been repeated. 0 = default; 3 = do not repeat any more.",//1
            "Maritime mobile service identity (MMSI) number of source station of message.",//2
            "Not used. Should be set to zero. Reserved for future use",//3
            "Safety related text.",//4
            "Extra bits for alignment; not used.",//5
        };
        unsigned int i;

        *fieldDescriptions=mxCreateStructArray(2, dims, numberOfFields, fieldNames);

        for(i=0;i<numberOfFields;i++) {
            mxSetFieldByNumber(*fieldDescriptions,0,static_cast<int>(i),mxCreateCharMatrixFromStrings(1,&descriptionStrings[i]));
        }
    }
    
    //Create the Matlab structure array.
    theStruct=mxCreateStructArray(2, dims, numberOfFields, fieldNames);
    *decodedMessage=theStruct;
    
    //Fill all of the elements of the structure
    mxSetFieldByNumber(theStruct,0,0,intMat2MatlabDoubles(&msg->message_id,1,1));
    mxSetFieldByNumber(theStruct,0,1,intMat2MatlabDoubles(&msg->repeat_indicator,1,1));
    mxSetFieldByNumber(theStruct,0,2,intMat2MatlabDoubles(&msg->mmsi,1,1));
    mxSetFieldByNumber(theStruct,0,3,intMat2MatlabDoubles(&msg->spare,1,1));
    {
        const char *charString=msg->text.c_str();
        mxSetFieldByNumber(theStruct,0,4,mxCreateCharMatrixFromStrings(1,&charString));
    }
    mxSetFieldByNumber(theStruct,0,5,intMat2MatlabDoubles(&msg->spare2,1,1));
}

void AIS_15_ToMatlab(libais::Ais15 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions) {
    const size_t numberOfFields=15;
    const mwSize dims[2] = {1, 1};
    //These are the names of all of the fields that will be added to the
    //Matlab structre array.
    const char *fieldNames[numberOfFields] = {"message_id", "repeat_indicator",
    "mmsi", "spare", "mmsi_1", "msg_1_1", "slot_offset_1_1", "spare2",
    "dest_msg_1_2","slot_offset_1_2","spare3","mmsi_2","msg_2",
    "slot_offset_2","spare4"};
    mxArray *theStruct;
    
    //If text descriptions of the fields are desired, then put them into a
    //structure.
    if(fieldDescriptions!=NULL) {
        const char *descriptionStrings[numberOfFields] = {
            "Identifier for Message 15; always 15.",//0
            "Repeat indicator. Used by the repeater to indicate how many times a message has been repeated. 0 = default; 3 = do not repeat any more.",//1
            "Maritime mobile service identity (MMSI) number of interrogating station",//2
            "Not used. Should be set to zero. Reserved for future use",//3
            "MMSI number of first interrogated station",//4
            "First requested message type from first interrogated station",//5
            "Response slot offset for first requested message from first interrogated station",//6
            "Not used. Should be set to zero. Reserved for future use",//7
            "Second requested message type from first interrogated station",//8
            "Response slot offset for second requested message from first interrogated station",//9
            "Not used. Should be set to zero. Reserved for future use",//10
            "MMSI number of second interrogated station",//11
            "Requested message type from second interrogated station",//12
            "Response slot offset for requested message from second interrogated station",//13
            "Not used. Should be set to zero. Reserved for future use"//14
        };
        unsigned int i;

        *fieldDescriptions=mxCreateStructArray(2, dims, numberOfFields, fieldNames);

        for(i=0;i<numberOfFields;i++) {
            mxSetFieldByNumber(*fieldDescriptions,0,static_cast<int>(i),mxCreateCharMatrixFromStrings(1,&descriptionStrings[i]));
        }
    }

    //Create the Matlab structure array.
    theStruct=mxCreateStructArray(2, dims, numberOfFields, fieldNames);
    *decodedMessage=theStruct;
    
    //Fill all of the elements of the structure
    mxSetFieldByNumber(theStruct,0,0,intMat2MatlabDoubles(&msg->message_id,1,1));
    mxSetFieldByNumber(theStruct,0,1,intMat2MatlabDoubles(&msg->repeat_indicator,1,1));
    mxSetFieldByNumber(theStruct,0,2,intMat2MatlabDoubles(&msg->mmsi,1,1));
    
    mxSetFieldByNumber(theStruct,0,3,intMat2MatlabDoubles(&msg->spare,1,1));
    mxSetFieldByNumber(theStruct,0,4,intMat2MatlabDoubles(&msg->mmsi_1,1,1));
    mxSetFieldByNumber(theStruct,0,5,intMat2MatlabDoubles(&msg->msg_1_1,1,1));
    mxSetFieldByNumber(theStruct,0,6,intMat2MatlabDoubles(&msg->slot_offset_1_1,1,1));
    mxSetFieldByNumber(theStruct,0,7,intMat2MatlabDoubles(&msg->spare2,1,1));
    mxSetFieldByNumber(theStruct,0,8,intMat2MatlabDoubles(&msg->dest_msg_1_2,1,1));
    mxSetFieldByNumber(theStruct,0,9,intMat2MatlabDoubles(&msg->slot_offset_1_2,1,1));
    mxSetFieldByNumber(theStruct,0,10,intMat2MatlabDoubles(&msg->spare3,1,1));
    mxSetFieldByNumber(theStruct,0,11,intMat2MatlabDoubles(&msg->mmsi_2,1,1));
    mxSetFieldByNumber(theStruct,0,12,intMat2MatlabDoubles(&msg->msg_2,1,1));
    mxSetFieldByNumber(theStruct,0,13,intMat2MatlabDoubles(&msg->slot_offset_2,1,1));
    mxSetFieldByNumber(theStruct,0,14,intMat2MatlabDoubles(&msg->spare4,1,1));
}

void AIS_16_ToMatlab(libais::Ais16 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions) {
    const size_t numberOfFields=11;
    const mwSize dims[2] = {1, 1};
    //These are the names of all of the fields that will be added to the
    //Matlab structre array.
    const char *fieldNames[numberOfFields] = {"message_id", "repeat_indicator",
    "mmsi", "spare", "dest_mmsi_a", "offset_a", "inc_a", "dest_mmsi_b",
    "offset_b", "inc_b", "spare2"};
    mxArray *theStruct;
    
    //If text descriptions of the fields are desired, then put them into a
    //structure.
    if(fieldDescriptions!=NULL) {
        const char *descriptionStrings[numberOfFields] = {
            "Identifier for Message 16; always 16.",//0
            "Repeat indicator. Used by the repeater to indicate how many times a message has been repeated. 0 = default; 3 = do not repeat any more.",//1
            "Maritime mobile service identity (MMSI) number of assigning station",//2
            "Spare. Should be set to zero. Reserved for future use",//3
            "MMSI number. Destination identifier A",//4
            "Offset from current slot to first assigned slot",//5
            "Increment to next assigned slot",//6
            "MMSI number. Destination identifier B. Should be omitted if there is assignment to station A, only",//7
            "Offset from current slot to first assigned slot. Should be omitted if there is assignment to station A, only",//8
            "Increment to next assigned slot. Should be omitted, if there is assignment to station A, only",//9
            "Spare. Not used. Should be set to zero. The number of spare bits, which should be 0 or 4, should be adjusted in order to observe byte boundaries. Reserved for future use"//10
        };
        unsigned int i;

        *fieldDescriptions=mxCreateStructArray(2, dims, numberOfFields, fieldNames);

        for(i=0;i<numberOfFields;i++) {
            mxSetFieldByNumber(*fieldDescriptions,0,static_cast<int>(i),mxCreateCharMatrixFromStrings(1,&descriptionStrings[i]));
        }
    }
    
    //Create the Matlab structure array.
    theStruct=mxCreateStructArray(2, dims, numberOfFields, fieldNames);
    *decodedMessage=theStruct;
    
    //Fill all of the elements of the structure
    mxSetFieldByNumber(theStruct,0,0,intMat2MatlabDoubles(&msg->message_id,1,1));
    mxSetFieldByNumber(theStruct,0,1,intMat2MatlabDoubles(&msg->repeat_indicator,1,1));
    mxSetFieldByNumber(theStruct,0,2,intMat2MatlabDoubles(&msg->mmsi,1,1));
    
    mxSetFieldByNumber(theStruct,0,3,intMat2MatlabDoubles(&msg->spare,1,1));
    mxSetFieldByNumber(theStruct,0,4,intMat2MatlabDoubles(&msg->dest_mmsi_a,1,1));
    mxSetFieldByNumber(theStruct,0,5,intMat2MatlabDoubles(&msg->offset_a,1,1));
    mxSetFieldByNumber(theStruct,0,6,intMat2MatlabDoubles(&msg->inc_a,1,1));
    
    if(msg->dest_mmsi_b!=-1) {
        mxSetFieldByNumber(theStruct,0,7,intMat2MatlabDoubles(&msg->dest_mmsi_b,1,1));
        mxSetFieldByNumber(theStruct,0,8,intMat2MatlabDoubles(&msg->offset_b,1,1));
        mxSetFieldByNumber(theStruct,0,9,intMat2MatlabDoubles(&msg->inc_b,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,7,mxCreateDoubleMatrix(0,0,mxREAL));
        mxSetFieldByNumber(theStruct,0,8,mxCreateDoubleMatrix(0,0,mxREAL));
        mxSetFieldByNumber(theStruct,0,9,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    mxSetFieldByNumber(theStruct,0,10,intMat2MatlabDoubles(&msg->spare2,1,1));
}

void AIS_17_ToMatlab(libais::Ais17 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions) {
    const size_t numberOfFields=12;
    const mwSize dims[2] = {1, 1};
    //These are the names of all of the fields that will be added to the
    //Matlab structre array.
    const char *fieldNames[numberOfFields] = {"message_id", "repeat_indicator",
    "mmsi", "spare", "lon", "lat", "spare2", "gnss_type", "station", "z_cnt", "seq",
    "health"};
    mxArray *theStruct;
        
    //If text descriptions of the fields are desired, then put them into a
    //structure.
    if(fieldDescriptions!=NULL) {
        const char *descriptionStrings[numberOfFields] = {
            "Identifier for Message 17; always 17.",//0
            "Repeat indicator. Used by the repeater to indicate how many times a message has been repeated. 0 = default; 3 = do not repeat any more.",//1
            "Maritime mobile service identity (MMSI) number of the base station",//2
            "Spare. Should be set to zero. Reserved for future use",//3
            "Surveyed longitude of DGNSS reference station in WGS-84 degrees",//4
            "Surveyed latitude of DGNSS reference station in WGS-84 degrees",//5
            "Not used. Should be set to zero. Reserved for future use",//6
            "GNSS Message Type as per ITU-R M.823",//7
            "Recommendation ITU-R M.823 station identifier",//8
            "Modified z-count in units of 0.6 seconds",//9
            "Message sequence number (cyclic 0-7)",//10
            "Reference station health (specified in Recommendation ITU-R M.823)"//11
        };
        unsigned int i;

        *fieldDescriptions=mxCreateStructArray(2, dims, numberOfFields, fieldNames);

        for(i=0;i<numberOfFields;i++) {
            mxSetFieldByNumber(*fieldDescriptions,0,static_cast<int>(i),mxCreateCharMatrixFromStrings(1,&descriptionStrings[i]));
        }
    }
    
    //Create the Matlab structure array.
    theStruct=mxCreateStructArray(2, dims, numberOfFields, fieldNames);
    *decodedMessage=theStruct;
    
    //Fill all of the elements of the structure
    mxSetFieldByNumber(theStruct,0,0,intMat2MatlabDoubles(&msg->message_id,1,1));
    mxSetFieldByNumber(theStruct,0,1,intMat2MatlabDoubles(&msg->repeat_indicator,1,1));
    mxSetFieldByNumber(theStruct,0,2,intMat2MatlabDoubles(&msg->mmsi,1,1));

    mxSetFieldByNumber(theStruct,0,3,intMat2MatlabDoubles(&msg->spare,1,1));
    
    if(fabs(msg->position.lng_deg)<=180) {
        mxSetFieldByNumber(theStruct,0,4,doubleMat2Matlab(&msg->position.lng_deg,1,1));
    } else {//No longitude is available
        mxSetFieldByNumber(theStruct,0,4,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(fabs(msg->position.lat_deg)<=90) {
        mxSetFieldByNumber(theStruct,0,5,doubleMat2Matlab(&msg->position.lat_deg,1,1));
    } else {//No latitude is available
        mxSetFieldByNumber(theStruct,0,5,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    mxSetFieldByNumber(theStruct,0,6,intMat2MatlabDoubles(&msg->spare2,1,1));
    mxSetFieldByNumber(theStruct,0,7,intMat2MatlabDoubles(&msg->gnss_type,1,1));
    mxSetFieldByNumber(theStruct,0,8,intMat2MatlabDoubles(&msg->station,1,1));
    mxSetFieldByNumber(theStruct,0,9,intMat2MatlabDoubles(&msg->z_cnt,1,1));
    mxSetFieldByNumber(theStruct,0,10,intMat2MatlabDoubles(&msg->seq,1,1));
    mxSetFieldByNumber(theStruct,0,11,intMat2MatlabDoubles(&msg->health,1,1));
}

void AIS_18_ToMatlab(libais::Ais18 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions) {
    const size_t numberOfFields=32;
    const mwSize dims[2] = {1, 1};
    //These are the names of all of the fields that will be added to the
    //Matlab structre array.
    const char *fieldNames[numberOfFields] = {"message_id", "repeat_indicator",
    "mmsi", "spare", "sog", "position_accuracy", "lon", "lat", "cog",
    "true_heading", "timestamp", "spare2", "unit_flag", "display_flag",
    "dsc_flag", "band_flag", "m22_flag", "mode_flag", "raim",
    "commstate_flag", "sync_state", "slot_timeout", "received_stations",
    "slot_number", "utc_hour", "utc_min", "utc_spare", "slot_offset",
    "slot_increment", "slots_to_allocate", "keep_flag","commstate_cs_fill"};
    mxArray *theStruct;
    
    //If text descriptions of the fields are desired, then put them into a
    //structure.
    if(fieldDescriptions!=NULL) {
        const char *descriptionStrings[numberOfFields] = {
            "Identifier for Message 18; always 18.",//0
            "Repeat indicator. Used by the repeater to indicate how many times a message has been repeated. 0 = default; 3 = do not repeat any more.",//1
            "Maritime mobile service identity (MMSI) number",//2
            "Not used. Should be set to zero. Reserved for future use",//3
            "Speed over ground in knots. 102.2 means 102.2 knots or higher.",//4
            "Position accuracy. 1=high (<=10m), 0=low (>10m), 0=default",//5
            "Longitude East in degrees (WGS-84)",//6
            "Latitude North in degrees (WGS-84)",//7
            "Course over ground in degrees East of North",//8
            "True heading (the direction the ship is pointing/ would be going without wind/currents) in degrees East of North",//9
            "UTC second when report was generated or 61 if positioning system is in manual input mode or 62 if in dead reckining mode or 63 if positioning system is inoperative",//10
            "Not used. Should be set to zero. Reserved for future use.",//11
            "Class B unit flag\n0 = Class B self-organized time division multiple access (SOTDMA) unit\n1 = Class B carrier sense (CS) unit",//12
            "Class B display flag\n0 = No display available; not capable of displaying Message 12 and 14\n1 = Equipped with integrated display displaying Message 12 and 14",//13
            "Class B digital selective calling (DSC) flag\n0 = Not equipped with DSC function\n1 = Equipped with DSC function (dedicated or time-shared)",//14
            "Class B band flag\n0 = Capable of operating over the upper 525 kHz band of the marine band\n1 = Capable of operating over the whole marine band\n(irrelevant if Class B Message 22 flag is 0)",//15
            "Class B Message 22 flag\n0 = No frequency management via Message 22, operating on AIS 1, AIS 2 only\n1 = Frequency management via Message 22",//16
            "0 = Station operating in autonomous and continuous mode = default 1 = Station operating in assigned mode",//17
            "Receiver autonomous integrity monitoring (RAIM) flag of electronic position fixing device; 0 = RAIM not in use = default; 1 = RAIM in use.",//18
            "0 = SOTDMA communication state follows 1 = ITDMA communication state follows,(always '1' for Class-B CS)",//19
            "Synchronization state\n0 UTC direct\n1 UTC indirect\n2 Station is synchronized to a base station (base direct)\n3 Station is synchronized to another station based on the highest number of received stations or to another mobile station, which is directly synchronized to a base station",//20
            "(For SOTDMA) Specifies frames remaining until a new slot is selected\n0 means that this was the last transmission in this slot\n1-7 means that 1 to 7 frames respectively are left until slot change",//21
            "(For SOTDMA) Number of other stations (not own station) which the station currently is receiving (between 0 and 16 383).",//22
            "(For SOTDMA) Slot number used for this transmission (between 0 and 2 249).",//23
            "(For SOTDMA) Hour (0-23) in UTC time, if available",//24
            "(For SOTDMA) Minute (0-23) in UTC time, if available",//25
            "(For SOTDMA) Two extra (unused) bits that are associated with the time message",//26
            "(For SOTDMA) If the slot time-out value is 0 (zero) then the slot offset should indicate the offset to the slot in which transmission will occur during the next frame. If the slot offset is zero, the slot should be de-allocated after transmission.",//27
            "(For ITDMA) Offset to next slot to be used, or zero (0) if no more transmissions",//28
            "(For ITDMA) Number of consecutive slots to allocate.\n0 = 1 slot,\n1 = 2 slots,\n2 = 3 slots,\n3 = 4 slots,\n4 = 5 slots,\n5 = 1 slot; offset = slot increment + 8 192,\n6 = 2 slots; offset = slot increment + 8 192,\n7 = 3 slots; offset = slot increment + 8 192.\nUse of 5 to 7 removes the need for RATDMA broadcast for scheduled transmissions up to 6 min intervals",//29
            "(For ITDMA) Set to TRUE = 1 if the slot remains allocated for one additional frame",//30
            "The value of the commstate region if commstate is set to 1 for carrier sense. This is supposed to be filled with the binary value 1100000000000000110"//31
        };
        unsigned int i;

        *fieldDescriptions=mxCreateStructArray(2, dims, numberOfFields, fieldNames);

        for(i=0;i<numberOfFields;i++) {
            mxSetFieldByNumber(*fieldDescriptions,0,static_cast<int>(i),mxCreateCharMatrixFromStrings(1,&descriptionStrings[i]));
        }
    }
    
    //Create the Matlab structure array.
    theStruct=mxCreateStructArray(2, dims, numberOfFields, fieldNames);
    *decodedMessage=theStruct;
    
    //Fill all of the elements of the structure
    mxSetFieldByNumber(theStruct,0,0,intMat2MatlabDoubles(&msg->message_id,1,1));
    mxSetFieldByNumber(theStruct,0,1,intMat2MatlabDoubles(&msg->repeat_indicator,1,1));
    mxSetFieldByNumber(theStruct,0,2,intMat2MatlabDoubles(&msg->mmsi,1,1));

    mxSetFieldByNumber(theStruct,0,3,intMat2MatlabDoubles(&msg->spare,1,1));
    
    if(msg->sog<=102.2) {
        mxSetFieldByNumber(theStruct,0,4,floatMat2MatlabDoubles(&msg->sog,1,1));
    } else {//If the speed over ground is unavailable or invalid
        mxSetFieldByNumber(theStruct,0,4,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    mxSetFieldByNumber(theStruct,0,5,intMat2MatlabDoubles(&msg->position_accuracy,1,1));
    
    if(fabs(msg->position.lng_deg)<=180) {
        mxSetFieldByNumber(theStruct,0,6,doubleMat2Matlab(&msg->position.lng_deg,1,1));
    } else {//Longitude is not available
        mxSetFieldByNumber(theStruct,0,6,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(fabs(msg->position.lat_deg)<=90) {
        mxSetFieldByNumber(theStruct,0,7,doubleMat2Matlab(&msg->position.lat_deg,1,1));
    } else {//Latitude is not available
        mxSetFieldByNumber(theStruct,0,7,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(msg->cog<360) {
        mxSetFieldByNumber(theStruct,0,8,floatMat2MatlabDoubles(&msg->cog,1,1));
    } else {//Course over ground is not available
        mxSetFieldByNumber(theStruct,0,8,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(msg->true_heading<360) {
        mxSetFieldByNumber(theStruct,0,9,intMat2MatlabDoubles(&msg->true_heading,1,1));
    } else {//True heading is not available
        mxSetFieldByNumber(theStruct,0,9,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(msg->timestamp!=60) {
        mxSetFieldByNumber(theStruct,0,10,intMat2MatlabDoubles(&msg->timestamp,1,1));
    } else {//Timestamp is not available
        mxSetFieldByNumber(theStruct,0,10,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    mxSetFieldByNumber(theStruct,0,11,intMat2MatlabDoubles(&msg->spare2,1,1));
    mxSetFieldByNumber(theStruct,0,12,intMat2MatlabDoubles(&msg->unit_flag,1,1));
    mxSetFieldByNumber(theStruct,0,13,intMat2MatlabDoubles(&msg->display_flag,1,1));
    mxSetFieldByNumber(theStruct,0,14,intMat2MatlabDoubles(&msg->dsc_flag,1,1));
    mxSetFieldByNumber(theStruct,0,15,intMat2MatlabDoubles(&msg->band_flag,1,1));
    mxSetFieldByNumber(theStruct,0,16,intMat2MatlabDoubles(&msg->m22_flag,1,1));
    mxSetFieldByNumber(theStruct,0,17,intMat2MatlabDoubles(&msg->mode_flag,1,1));
    mxSetFieldByNumber(theStruct,0,18,boolMat2Matlab(&msg->raim,1,1));
    mxSetFieldByNumber(theStruct,0,19,intMat2MatlabDoubles(&msg->commstate_flag,1,1));

    mxSetFieldByNumber(theStruct,0,20,intMat2MatlabDoubles(&msg->sync_state,1,1));
    if(msg->commstate_flag==0) {
        mxSetFieldByNumber(theStruct,0,21,intMat2MatlabDoubles(&msg->slot_timeout,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,21,mxCreateDoubleMatrix(0,0,mxREAL));
    }

    if(msg->received_stations_valid) {
        mxSetFieldByNumber(theStruct,0,22,intMat2MatlabDoubles(&msg->received_stations,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,22,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(msg->slot_number_valid) {
        mxSetFieldByNumber(theStruct,0,23,intMat2MatlabDoubles(&msg->slot_number,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,23,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(msg->utc_valid) {        
        mxSetFieldByNumber(theStruct,0,24,intMat2MatlabDoubles(&msg->utc_hour,1,1));
        mxSetFieldByNumber(theStruct,0,25,intMat2MatlabDoubles(&msg->utc_min,1,1));
        mxSetFieldByNumber(theStruct,0,26,intMat2MatlabDoubles(&msg->utc_spare,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,24,mxCreateDoubleMatrix(0,0,mxREAL));
        mxSetFieldByNumber(theStruct,0,25,mxCreateDoubleMatrix(0,0,mxREAL));
        mxSetFieldByNumber(theStruct,0,26,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(msg->slot_offset_valid) {
        mxSetFieldByNumber(theStruct,0,27,intMat2MatlabDoubles(&msg->slot_offset,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,27,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(msg->slot_increment_valid) {
        mxSetFieldByNumber(theStruct,0,28,intMat2MatlabDoubles(&msg->slot_increment,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,28,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(msg->slots_to_allocate_valid) {
        mxSetFieldByNumber(theStruct,0,29,intMat2MatlabDoubles(&msg->slots_to_allocate,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,29,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(msg->keep_flag_valid) {
        mxSetFieldByNumber(theStruct,0,30,boolMat2Matlab(&msg->keep_flag,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,30,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(msg->commstate_cs_fill_valid) {
        mxSetFieldByNumber(theStruct,0,31,intMat2MatlabDoubles(&msg->commstate_cs_fill,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,31,mxCreateDoubleMatrix(0,0,mxREAL));
    }
}

void AIS_19_ToMatlab(libais::Ais19 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions) {
    const size_t numberOfFields=23;
    const mwSize dims[2] = {1, 1};
    //These are the names of all of the fields that will be added to the
    //Matlab structre array.
    const char *fieldNames[numberOfFields] = {"message_id", "repeat_indicator",
    "mmsi", "spare", "sog", "position_accuracy", "lon", "lat", "cog",
    "true_heading", "timestamp", "spare2", "name", "type_and_cargo",
    "dim_a", "dim_b", "dim_c", "dim_d", "fix_type", "raim", "dte",
    "assigned_mode", "spare3"};
    mxArray *theStruct;
        
    //If text descriptions of the fields are desired, then put them into a
    //structure.
    if(fieldDescriptions!=NULL) {
        const char *descriptionStrings[numberOfFields] = {
            "Identifier for Message 19; always 19.",//0
            "Repeat indicator. Used by the repeater to indicate how many times a message has been repeated. 0 = default; 3 = do not repeat any more.",//1
            "Maritime mobile service identity (MMSI) number",//2
            "Not used. Should be set to zero. Reserved for future use",//3
            "Speed over ground in knots. 102.2 means 102.2 knots or higher.",//4
            "Position accuracy. 1=high (<=10m), 0=low (>10m), 0=default",//5
            "Longitude East in degrees (WGS-84)",//6
            "Latitude North in degrees (WGS-84)",//7
            "Course over ground in degrees East of North",//8
            "True heading (the direction the ship is pointing/ would be going without wind/currents) in degrees East of North",//9
            "UTC second when report was generated or 61 if positioning system is in manual input mode or 62 if in dead reckining mode or 63 if positioning system is inoperative",//10
            "Not used. Should be set to zero. Reserved for future use.",//11
            "Name. Maximum 20 characters.",//12
            "Type of ship and cargo type\n1-99 = as defined in Section 3.3.2 of Annex 8 of ITU-R M.1371-5\n100-199 = reserved, for regional use 200-255 = reserved, for future use\nNot applicable to search and rescue (SAR) aircraft",//13
            "Distance from reference point to bow of ship (meters) 511m means 511m or longer. (If only relative reference point position known=0)",//14
            "Distance from reference point to aft of ship (meters) 511m means 511m or longer. ",//15
            "Distance from reference point to port side of ship (meters) 511m means 511m or longer. (If only relative reference point position known=0)",//16
            "Distance from reference point to starboard side of ship (meters) 511m means 511m or longer.",//17
            "Type of electronic position fixing device:\n1 = global positioning system (GPS)\n2 = GNSS (GLONASS)\n3 = combined GPS/GLONASS\n4 = Loran-C\n5 = Chayka\n6 = integrated navigation system\n7 = surveyed\n8 = Galileo\n9-14 = not used\n15 = internal GNSS",//18
            "Receiver autonomous integrity monitoring (RAIM) flag of electronic position fixing device; 0 = RAIM not in use = default; 1 = RAIM in use.",//19
            "Data terminal equipment (DTE) ready (0 = available, 1 = not available = default)",//20
            "0 = Station operating in autonomous and continuous mode = default 1 = Station operating in assigned mode",//21
            "Not used. Should be set to zero. Reserved for future use"//22
        };
        unsigned int i;

        *fieldDescriptions=mxCreateStructArray(2, dims, numberOfFields, fieldNames);

        for(i=0;i<numberOfFields;i++) {
            mxSetFieldByNumber(*fieldDescriptions,0,static_cast<int>(i),mxCreateCharMatrixFromStrings(1,&descriptionStrings[i]));
        }
    }

    //Create the Matlab structure array.
    theStruct=mxCreateStructArray(2, dims, numberOfFields, fieldNames);
    *decodedMessage=theStruct;
    
    //Fill all of the elements of the structure
    mxSetFieldByNumber(theStruct,0,0,intMat2MatlabDoubles(&msg->message_id,1,1));
    mxSetFieldByNumber(theStruct,0,1,intMat2MatlabDoubles(&msg->repeat_indicator,1,1));
    mxSetFieldByNumber(theStruct,0,2,intMat2MatlabDoubles(&msg->mmsi,1,1));
    mxSetFieldByNumber(theStruct,0,3,intMat2MatlabDoubles(&msg->spare,1,1));
    
     if(msg->sog<=102.2) {
        mxSetFieldByNumber(theStruct,0,4,floatMat2MatlabDoubles(&msg->sog,1,1));
    } else {//If the speed over ground is unavailable or invalid
        mxSetFieldByNumber(theStruct,0,4,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    mxSetFieldByNumber(theStruct,0,5,intMat2MatlabDoubles(&msg->position_accuracy,1,1));
    
    if(fabs(msg->position.lng_deg)<=180) {
        mxSetFieldByNumber(theStruct,0,6,doubleMat2Matlab(&msg->position.lng_deg,1,1));
    } else {//Longitude is not available
        mxSetFieldByNumber(theStruct,0,6,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(fabs(msg->position.lat_deg)<=90) {
        mxSetFieldByNumber(theStruct,0,7,doubleMat2Matlab(&msg->position.lat_deg,1,1));
    } else {//Latitude is not available
        mxSetFieldByNumber(theStruct,0,7,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(msg->cog<360) {
        mxSetFieldByNumber(theStruct,0,8,floatMat2MatlabDoubles(&msg->cog,1,1));
    } else {//Course over ground is not available
        mxSetFieldByNumber(theStruct,0,8,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(msg->true_heading<360) {
        mxSetFieldByNumber(theStruct,0,9,intMat2MatlabDoubles(&msg->true_heading,1,1));
    } else {//True heading is not available
        mxSetFieldByNumber(theStruct,0,9,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(msg->timestamp!=60) {
        mxSetFieldByNumber(theStruct,0,10,intMat2MatlabDoubles(&msg->timestamp,1,1));
    } else {//Timestamp is not available
        mxSetFieldByNumber(theStruct,0,10,mxCreateDoubleMatrix(0,0,mxREAL));
    }

    mxSetFieldByNumber(theStruct,0,11,intMat2MatlabDoubles(&msg->spare2,1,1));

    if(msg->name.compare("@@@@@@@@@@@@@@@@@@@@")!=0) {
        const char *charString=msg->name.c_str();
        mxSetFieldByNumber(theStruct,0,12,mxCreateCharMatrixFromStrings(1,&charString));
    } else {//If the name is not available
        mxSetFieldByNumber(theStruct,0,12,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(msg->type_and_cargo!=0) {
        mxSetFieldByNumber(theStruct,0,13,intMat2MatlabDoubles(&msg->type_and_cargo,1,1));
    } else {//No type and cargo information
        mxSetFieldByNumber(theStruct,0,13,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(msg->dim_a==0&&msg->dim_b==0&&msg->dim_c==0&&msg->dim_d==0) {
        mxSetFieldByNumber(theStruct,0,14,mxCreateDoubleMatrix(0,0,mxREAL));
        mxSetFieldByNumber(theStruct,0,15,mxCreateDoubleMatrix(0,0,mxREAL));
        mxSetFieldByNumber(theStruct,0,16,mxCreateDoubleMatrix(0,0,mxREAL));
        mxSetFieldByNumber(theStruct,0,17,mxCreateDoubleMatrix(0,0,mxREAL));
    } else {//If ship dimensions from reference point are available.
        mxSetFieldByNumber(theStruct,0,14,intMat2MatlabDoubles(&msg->dim_a,1,1));
        mxSetFieldByNumber(theStruct,0,15,intMat2MatlabDoubles(&msg->dim_b,1,1));
        mxSetFieldByNumber(theStruct,0,16,intMat2MatlabDoubles(&msg->dim_c,1,1));
        mxSetFieldByNumber(theStruct,0,17,intMat2MatlabDoubles(&msg->dim_d,1,1));
    }
    
    if(msg->fix_type!=0) {
        mxSetFieldByNumber(theStruct,0,18,intMat2MatlabDoubles(&msg->fix_type,1,1));
    } else {//If the fix type is not available
        mxSetFieldByNumber(theStruct,0,18,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    mxSetFieldByNumber(theStruct,0,19,boolMat2Matlab(&msg->raim,1,1));
    mxSetFieldByNumber(theStruct,0,20,intMat2MatlabDoubles(&msg->dte,1,1));
    mxSetFieldByNumber(theStruct,0,21,intMat2MatlabDoubles(&msg->assigned_mode,1,1));
    mxSetFieldByNumber(theStruct,0,22,intMat2MatlabDoubles(&msg->spare3,1,1));
}

void AIS_20_ToMatlab(libais::Ais20 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions) {
    const size_t numberOfFields=21;
    const mwSize dims[2] = {1, 1};
    //These are the names of all of the fields that will be added to the
    //Matlab structre array.
    const char *fieldNames[numberOfFields] = {"message_id", "repeat_indicator",
    "mmsi", "spare", "offset_1", "num_slots_1", "timeout_1", "incr_1",
    "offset_2", "num_slots_2", "timeout_2", "incr_2", "offset_3",
    "num_slots_3", "timeout_3", "incr_3", "offset_4", "num_slots_4",
    "timeout_4", "incr_4", "spare2"};
    mxArray *theStruct;
      
    //If text descriptions of the fields are desired, then put them into a
    //structure.
    if(fieldDescriptions!=NULL) {
        const char *descriptionStrings[numberOfFields] = {
            "Identifier for Message 20; always 20.",//0
            "Repeat indicator. Used by the repeater to indicate how many times a message has been repeated. 0 = default; 3 = do not repeat any more.",//1
            "Maritime mobile service identity (MMSI) number of base station",//2
            "Not used. Should be set to zero. Reserved for future use",//3
            "Reserved offset number 1",//4
            "Number of reserved consecutive slots 1: 1-15",//5
            "Time-out value 1 in minutes",//6
            "Increment to repeat reservation block 1; 0 = one reservation block per frame",//7
            "Reserved offset number 2 (optional)",//8
            "Number of reserved consecutive slots 2: 1-15; optional",//9
            "Time-out value 2 in minutes (optional)",//10
            "Increment to repeat reservation block 2 (optional)",//11
            "Reserved offset number 3 (optional)",//12
            "Number of reserved consecutive slots 3: 1-15; optional",//13
            "Time-out value 3 in minutes (optional)",//14
            "Increment to repeat reservation block 3 (optional)",//15
            "Reserved offset number 4 (optional)",//16
            "Number of reserved consecutive slots 4: 1-15; optional",//17
            "Time-out value 4 in minutes (optional)",//18
            "Increment to repeat reservation block 4 (optional)",//19
            "Not used. Should be set to zero. The number of spare bits which may be 0, 2, 4 or 6 should be adjusted in order to observe byte boundaries. Reserved for future use"//20
        };
        unsigned int i;

        *fieldDescriptions=mxCreateStructArray(2, dims, numberOfFields, fieldNames);

        for(i=0;i<numberOfFields;i++) {
            mxSetFieldByNumber(*fieldDescriptions,0,static_cast<int>(i),mxCreateCharMatrixFromStrings(1,&descriptionStrings[i]));
        }
    }

    //Create the Matlab structure array.
    theStruct=mxCreateStructArray(2, dims, numberOfFields, fieldNames);
    *decodedMessage=theStruct;
    
    //Fill all of the elements of the structure
    mxSetFieldByNumber(theStruct,0,0,intMat2MatlabDoubles(&msg->message_id,1,1));
    mxSetFieldByNumber(theStruct,0,1,intMat2MatlabDoubles(&msg->repeat_indicator,1,1));
    mxSetFieldByNumber(theStruct,0,2,intMat2MatlabDoubles(&msg->mmsi,1,1));

    mxSetFieldByNumber(theStruct,0,3,intMat2MatlabDoubles(&msg->spare,1,1));
    if(msg->offset_1!=0) {
        mxSetFieldByNumber(theStruct,0,4,intMat2MatlabDoubles(&msg->offset_1,1,1));
    } else {//If the offset number is not available
        mxSetFieldByNumber(theStruct,0,4,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(msg->num_slots_1!=0) {
        mxSetFieldByNumber(theStruct,0,5,intMat2MatlabDoubles(&msg->num_slots_1,1,1));
    } else {//If the number of slots is not available
        mxSetFieldByNumber(theStruct,0,5,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(msg->timeout_1!=0) {
        mxSetFieldByNumber(theStruct,0,6,intMat2MatlabDoubles(&msg->timeout_1,1,1));
    } else {//If the timeout number is not available
        mxSetFieldByNumber(theStruct,0,6,mxCreateDoubleMatrix(0,0,mxREAL));
    }

    mxSetFieldByNumber(theStruct,0,7,intMat2MatlabDoubles(&msg->incr_1,1,1));

    if(msg->group_valid_2) {
        if(msg->offset_2!=0) {
            mxSetFieldByNumber(theStruct,0,8,intMat2MatlabDoubles(&msg->offset_2,1,1));
        } else {
            mxSetFieldByNumber(theStruct,0,8,mxCreateDoubleMatrix(0,0,mxREAL));
        }
        
        if(msg->num_slots_2!=0) {
            mxSetFieldByNumber(theStruct,0,9,intMat2MatlabDoubles(&msg->num_slots_2,1,1));
        } else {
            mxSetFieldByNumber(theStruct,0,9,mxCreateDoubleMatrix(0,0,mxREAL));
        }
        
        if(msg->timeout_2!=0) {
            mxSetFieldByNumber(theStruct,0,10,intMat2MatlabDoubles(&msg->timeout_2,1,1));
        } else {
            mxSetFieldByNumber(theStruct,0,10,mxCreateDoubleMatrix(0,0,mxREAL));
        }
        
        mxSetFieldByNumber(theStruct,0,11,intMat2MatlabDoubles(&msg->incr_2,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,8,mxCreateDoubleMatrix(0,0,mxREAL));
        mxSetFieldByNumber(theStruct,0,9,mxCreateDoubleMatrix(0,0,mxREAL));
        mxSetFieldByNumber(theStruct,0,10,mxCreateDoubleMatrix(0,0,mxREAL));
        mxSetFieldByNumber(theStruct,0,11,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(msg->group_valid_3) {
        if(msg->offset_3!=0) {
            mxSetFieldByNumber(theStruct,0,12,intMat2MatlabDoubles(&msg->offset_3,1,1));
        } else {
            mxSetFieldByNumber(theStruct,0,12,mxCreateDoubleMatrix(0,0,mxREAL));
        }
        
        if(msg->num_slots_3!=0) {
            mxSetFieldByNumber(theStruct,0,13,intMat2MatlabDoubles(&msg->num_slots_3,1,1));
        } else {
            mxSetFieldByNumber(theStruct,0,13,mxCreateDoubleMatrix(0,0,mxREAL));
        }
        
        if(msg->timeout_3!=0) {
            mxSetFieldByNumber(theStruct,0,14,intMat2MatlabDoubles(&msg->timeout_3,1,1));
        } else {
            mxSetFieldByNumber(theStruct,0,14,mxCreateDoubleMatrix(0,0,mxREAL));
        }

        mxSetFieldByNumber(theStruct,0,15,intMat2MatlabDoubles(&msg->incr_3,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,12,mxCreateDoubleMatrix(0,0,mxREAL));
        mxSetFieldByNumber(theStruct,0,13,mxCreateDoubleMatrix(0,0,mxREAL));
        mxSetFieldByNumber(theStruct,0,14,mxCreateDoubleMatrix(0,0,mxREAL));
        mxSetFieldByNumber(theStruct,0,15,mxCreateDoubleMatrix(0,0,mxREAL)); 
    }
    
    if(msg->group_valid_4) {
        if(msg->offset_4!=0) {
            mxSetFieldByNumber(theStruct,0,16,intMat2MatlabDoubles(&msg->offset_4,1,1));
        } else {
            mxSetFieldByNumber(theStruct,0,16,mxCreateDoubleMatrix(0,0,mxREAL));
        }
        
        if(msg->num_slots_4!=0) {
            mxSetFieldByNumber(theStruct,0,17,intMat2MatlabDoubles(&msg->num_slots_4,1,1));
        } else {
            mxSetFieldByNumber(theStruct,0,17,mxCreateDoubleMatrix(0,0,mxREAL));
        }
        
        if(msg->timeout_4!=0) {
            mxSetFieldByNumber(theStruct,0,18,intMat2MatlabDoubles(&msg->timeout_4,1,1));
        } else {
            mxSetFieldByNumber(theStruct,0,18,mxCreateDoubleMatrix(0,0,mxREAL));
        }
        
        mxSetFieldByNumber(theStruct,0,19,intMat2MatlabDoubles(&msg->incr_4,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,16,mxCreateDoubleMatrix(0,0,mxREAL));
        mxSetFieldByNumber(theStruct,0,17,mxCreateDoubleMatrix(0,0,mxREAL));
        mxSetFieldByNumber(theStruct,0,18,mxCreateDoubleMatrix(0,0,mxREAL));
        mxSetFieldByNumber(theStruct,0,19,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    mxSetFieldByNumber(theStruct,0,20,intMat2MatlabDoubles(&msg->spare2,1,1));
}

void AIS_21_ToMatlab(libais::Ais21 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions) {
    const size_t numberOfFields=21;
    const mwSize dims[2] = {1, 1};
    //These are the names of all of the fields that will be added to the
    //Matlab structre array.
    const char *fieldNames[numberOfFields] = {"message_id", "repeat_indicator",
    "mmsi", "aton_type", "name", "position_accuracy", "lon", "lat", "dim_a",
    "dim_b", "dim_c", "dim_d", "fix_type", "timestamp", "off_pos",
    "aton_status", "raim", "virtual_aton", "assigned_mode", "spare",
    "spare2"};
    mxArray *theStruct;
         
    //If text descriptions of the fields are desired, then put them into a
    //structure.
    if(fieldDescriptions!=NULL) {
        const char *descriptionStrings[numberOfFields] = {
            "Identifier for Message 21; always 21.",//0
            "Repeat indicator. Used by the repeater to indicate how many times a message has been repeated. 0 = default; 3 = do not repeat any more.",//1
            "Maritime mobile service identity (MMSI) number (see Article 19 of the RR and Recommendation ITU-R M.585)",//2
            "Type of aids-to-navigation (AtoN); refer to appropriate definition set up by IALA",
            "Name of AtoN, including any extension",
            "Position accuracy. 1=high (<=10m), 0=low (>10m), 0=default",//5
            "Longitude East in degrees (WGS-84)",//6
            "Latitude North in degrees (WGS-84)",//7
            "Distance from reference point to true North of structure (meters) 1 means less than 2 meters",//8
            "Distance from reference point to South of structure (meters) 1 means less than 2 meters",//9
            "Distance from reference point to West side of structure (meters) 1 means less than 2 meters",//10
            "Distance from reference point to East side of structure (meters) 1 means less than 2 meters",//11
            "Type of electronic position fixing device:\n1 = GPS\n2 = GLONASS\n3 = Combined GPS/GLONASS\n4 = Loran-C\n5 = Chayka\n6 = Integrated Navigation System\n7 = surveyed. For fixed AtoN and virtual AtoN, the charted position should be used. The accurate position enhances its function as a radar reference target\n8 = Galileo\n9-14 = not used\n15 = internal GNSS",//12
            "UTC second when report was generated or 61 if positioning system is in manual input mode or 62 if in dead reckining mode or 63 if positioning system is inoperative",//13
            "Off-position indicator. For floating AtoN, only: 0 = on position; 1 = off position.",//14
            "Reserved for the indication of the AtoN status 00000000 = default",//15
            "RAIM (Receiver autonomous integrity monitoring) flag of electronic position fixing device; 0 = RAIM not in use = default; 1 = RAIM in use",//16
            "0 = default = real AtoN at indicated position; 1 = virtual AtoN, does not physically exist.",//17
            "0 = Station operating in autonomous and continuous mode = default 1 = Station operating in assigned mode",//18
            "Spare. Not used. Should be set to zero. Reserved for future use",//19
            "Spare. Should be set to zero. The number of spare bits should be adjusted in order to observe byte boundaries"//20
        };
        unsigned int i;

        *fieldDescriptions=mxCreateStructArray(2, dims, numberOfFields, fieldNames);

        for(i=0;i<numberOfFields;i++) {
            mxSetFieldByNumber(*fieldDescriptions,0,static_cast<int>(i),mxCreateCharMatrixFromStrings(1,&descriptionStrings[i]));
        }
    }
    
    //Create the Matlab structure array.
    theStruct=mxCreateStructArray(2, dims, numberOfFields, fieldNames);
    *decodedMessage=theStruct;
    
    //Fill all of the elements of the structure
    mxSetFieldByNumber(theStruct,0,0,intMat2MatlabDoubles(&msg->message_id,1,1));
    mxSetFieldByNumber(theStruct,0,1,intMat2MatlabDoubles(&msg->repeat_indicator,1,1));
    mxSetFieldByNumber(theStruct,0,2,intMat2MatlabDoubles(&msg->mmsi,1,1));
    
    if(msg->aton_type!=0) {
        mxSetFieldByNumber(theStruct,0,3,intMat2MatlabDoubles(&msg->aton_type,1,1));
    } else {//Aid to navigation type unavailable
        mxSetFieldByNumber(theStruct,0,3,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(msg->name.compare("@@@@@@@@@@@@@@@@@@@@")!=0) {
        const char *charString=msg->name.c_str();
        mxSetFieldByNumber(theStruct,0,4,mxCreateCharMatrixFromStrings(1,&charString));
    } else {//If the name is not available
        mxSetFieldByNumber(theStruct,0,4,mxCreateDoubleMatrix(0,0,mxREAL));
    }

    mxSetFieldByNumber(theStruct,0,5,intMat2MatlabDoubles(&msg->position_accuracy,1,1));
    
    if(fabs(msg->position.lng_deg)<=180) {
        mxSetFieldByNumber(theStruct,0,6,doubleMat2Matlab(&msg->position.lng_deg,1,1));
    } else {//Longitude not available
        mxSetFieldByNumber(theStruct,0,6,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(fabs(msg->position.lat_deg)<=90) {
        mxSetFieldByNumber(theStruct,0,7,doubleMat2Matlab(&msg->position.lat_deg,1,1));
    } else {//Latitude not available
        mxSetFieldByNumber(theStruct,0,7,mxCreateDoubleMatrix(0,0,mxREAL));
    }

    if(msg->dim_a==0&&msg->dim_b==0&&msg->dim_c==0&&msg->dim_d==0) {
        mxSetFieldByNumber(theStruct,0,8,mxCreateDoubleMatrix(0,0,mxREAL));
        mxSetFieldByNumber(theStruct,0,9,mxCreateDoubleMatrix(0,0,mxREAL));
        mxSetFieldByNumber(theStruct,0,10,mxCreateDoubleMatrix(0,0,mxREAL));
        mxSetFieldByNumber(theStruct,0,11,mxCreateDoubleMatrix(0,0,mxREAL));
    } else {//If dimensions are available
        mxSetFieldByNumber(theStruct,0,8,intMat2MatlabDoubles(&msg->dim_a,1,1));
        mxSetFieldByNumber(theStruct,0,9,intMat2MatlabDoubles(&msg->dim_b,1,1));
        mxSetFieldByNumber(theStruct,0,10,intMat2MatlabDoubles(&msg->dim_c,1,1));
        mxSetFieldByNumber(theStruct,0,11,intMat2MatlabDoubles(&msg->dim_d,1,1));
    }
    
    if(msg->fix_type!=0) {
        mxSetFieldByNumber(theStruct,0,12,intMat2MatlabDoubles(&msg->fix_type,1,1));
    } else {//If the fix type is not available
        mxSetFieldByNumber(theStruct,0,12,mxCreateDoubleMatrix(0,0,mxREAL));
    }

    if(msg->timestamp!=60) {
        mxSetFieldByNumber(theStruct,0,13,intMat2MatlabDoubles(&msg->timestamp,1,1));
    } else {//timestamp not available
        mxSetFieldByNumber(theStruct,0,13,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(msg->timestamp<=59) {
        mxSetFieldByNumber(theStruct,0,14,boolMat2Matlab(&msg->off_pos,1,1));
    } else {//Off-position indicator is invalid
        mxSetFieldByNumber(theStruct,0,14,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    mxSetFieldByNumber(theStruct,0,15,intMat2MatlabDoubles(&msg->aton_status,1,1));
    mxSetFieldByNumber(theStruct,0,16,boolMat2Matlab(&msg->raim,1,1));
    mxSetFieldByNumber(theStruct,0,17,boolMat2Matlab(&msg->virtual_aton,1,1));
    mxSetFieldByNumber(theStruct,0,18,boolMat2Matlab(&msg->assigned_mode,1,1));
    mxSetFieldByNumber(theStruct,0,19,intMat2MatlabDoubles(&msg->spare,1,1));
    mxSetFieldByNumber(theStruct,0,20,intMat2MatlabDoubles(&msg->spare2,1,1));
}

void AIS_22_ToMatlab(libais::Ais22 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions) {
    const size_t numberOfFields=18;
    const mwSize dims[2] = {1, 1};
    //These are the names of all of the fields that will be added to the
    //Matlab structre array.
    const char *fieldNames[numberOfFields] = {"message_id", "repeat_indicator",
    "mmsi", "spare", "chan_a", "chan_b", "txrx_mode", "power_low", "lon1",
    "lat1", "lon2", "lat2", "dest_mmsi_1", "dest_mmsi_2", "chan_a_bandwidth",
    "chan_b_bandwidth", "zone_size", "spare2"};
    mxArray *theStruct;
         
    //If text descriptions of the fields are desired, then put them into a
    //structure.
    if(fieldDescriptions!=NULL) {
        const char *descriptionStrings[numberOfFields] = {
            "Identifier for Message 22; always 22.",//0
            "Repeat indicator. Used by the repeater to indicate how many times a message has been repeated. 0 = default; 3 = do not repeat any more.",//1
            "Maritime mobile service identity (MMSI) number of base station",//2
            "Not used. Should be set to zero. Reserved for future use",//3
            "Channel number of 25 kHz simplex or simplex use of 25 kHz duplex in accordance with Recommendation ITU-R M.1084.",//4
            "Channel number of 25 kHz simplex or simplex use of 25 kHz duplex in accordance with Recommendation ITU-R M.1084.",//5
            "Tx/Rx mode\n0 = Tx A/Tx B, Rx A/Rx B (default)\n1 = Tx A, Rx A/Rx B\n2 = Tx B, Rx A/Rx B\n3-15: not used\nWhen the dual channel transmission is suspended by Tx/Rx mode command 1 or 2, the required reporting interval should be maintained using the remaining transmission channel",//6
            "Power: 0 = high (default), 1 = low",//7
            "Longitude (WGS-84) of area to which the assignment applies (degrees); upper right corner (North-East)",//8
            "Latitude (WGS-84) of area to which the assignment applies (degrees); upper right corner (North-East)",//9
            "Longitude (WGS-84) of area to which the assignment applies (degrees); lower left corner (South-West)",//10
            "Latitude (WGS-84) of area to which the assignment applies (degrees); lower left corner (South-West)",//11
            "MMSI of addressed station 1",//12
            "MMSI of addressed station 2",//13
            "Channel A bandwidth:\n0 = default (as specified by channel number);\n1 = spare (formerly 12.5 kHz bandwidth in Recommendation ITU-R M.1371-1 now obsolete)",//14
            "Channel B bandwidth:\n0 = default (as specified by channel number);\n1 = spare (formerly 12.5 kHz bandwidth in Recommendation ITU-R M.1371-1 now obsolete)",//15
            "The transitional zone size in nautical miles should be calculated by adding 1 to this parameter value.",//16
            "Not used. Should be set to zero. Reserved for future use"//17
        };
        unsigned int i;

        *fieldDescriptions=mxCreateStructArray(2, dims, numberOfFields, fieldNames);

        for(i=0;i<numberOfFields;i++) {
            mxSetFieldByNumber(*fieldDescriptions,0,static_cast<int>(i),mxCreateCharMatrixFromStrings(1,&descriptionStrings[i]));
        }
    }
    
    //Create the Matlab structure array.
    theStruct=mxCreateStructArray(2, dims, numberOfFields, fieldNames);
    *decodedMessage=theStruct;
    
    //Fill all of the elements of the structure
    mxSetFieldByNumber(theStruct,0,0,intMat2MatlabDoubles(&msg->message_id,1,1));
    mxSetFieldByNumber(theStruct,0,1,intMat2MatlabDoubles(&msg->repeat_indicator,1,1));
    mxSetFieldByNumber(theStruct,0,2,intMat2MatlabDoubles(&msg->mmsi,1,1));

    mxSetFieldByNumber(theStruct,0,3,intMat2MatlabDoubles(&msg->spare,1,1));
    mxSetFieldByNumber(theStruct,0,4,intMat2MatlabDoubles(&msg->chan_a,1,1));
    mxSetFieldByNumber(theStruct,0,5,intMat2MatlabDoubles(&msg->chan_b,1,1));
    mxSetFieldByNumber(theStruct,0,6,intMat2MatlabDoubles(&msg->txrx_mode,1,1));
    mxSetFieldByNumber(theStruct,0,7,boolMat2Matlab(&msg->power_low,1,1));

    if(msg->pos_valid) {
        if(fabs(msg->position1.lng_deg)<=180) {
            mxSetFieldByNumber(theStruct,0,8,doubleMat2Matlab(&msg->position1.lng_deg,1,1));
        } else {
            mxSetFieldByNumber(theStruct,0,8,mxCreateDoubleMatrix(0,0,mxREAL));
        }
        
        if(fabs(msg->position1.lat_deg)<=90) {
            mxSetFieldByNumber(theStruct,0,9,doubleMat2Matlab(&msg->position1.lat_deg,1,1));
        } else {
            mxSetFieldByNumber(theStruct,0,9,mxCreateDoubleMatrix(0,0,mxREAL));
        }
        
        if(fabs(msg->position2.lng_deg)<=180) {
            mxSetFieldByNumber(theStruct,0,10,doubleMat2Matlab(&msg->position2.lng_deg,1,1));
        } else {
            mxSetFieldByNumber(theStruct,0,10,mxCreateDoubleMatrix(0,0,mxREAL));
        }
        
        if(fabs(msg->position2.lat_deg)<=90) {
            mxSetFieldByNumber(theStruct,0,11,doubleMat2Matlab(&msg->position2.lat_deg,1,1));
        } else {
            mxSetFieldByNumber(theStruct,0,11,mxCreateDoubleMatrix(0,0,mxREAL));
        }
    } else {
        mxSetFieldByNumber(theStruct,0,8,mxCreateDoubleMatrix(0,0,mxREAL));
        mxSetFieldByNumber(theStruct,0,9,mxCreateDoubleMatrix(0,0,mxREAL));
        mxSetFieldByNumber(theStruct,0,10,mxCreateDoubleMatrix(0,0,mxREAL));
        mxSetFieldByNumber(theStruct,0,11,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(msg->dest_valid) {
        mxSetFieldByNumber(theStruct,0,12,intMat2MatlabDoubles(&msg->dest_mmsi_1,1,1));
        mxSetFieldByNumber(theStruct,0,13,intMat2MatlabDoubles(&msg->dest_mmsi_2,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,12,mxCreateDoubleMatrix(0,0,mxREAL));
        mxSetFieldByNumber(theStruct,0,13,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    mxSetFieldByNumber(theStruct,0,14,intMat2MatlabDoubles(&msg->chan_a_bandwidth,1,1));
    mxSetFieldByNumber(theStruct,0,15,intMat2MatlabDoubles(&msg->chan_b_bandwidth,1,1));
    mxSetFieldByNumber(theStruct,0,16,intMat2MatlabDoubles(&msg->zone_size,1,1));
    mxSetFieldByNumber(theStruct,0,17,intMat2MatlabDoubles(&msg->spare2,1,1));
}

void AIS_23_ToMatlab(libais::Ais23 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions) {
    const size_t numberOfFields=15;
    const mwSize dims[2] = {1, 1};
    //These are the names of all of the fields that will be added to the
    //Matlab structre array.
    const char *fieldNames[numberOfFields] = {"message_id", "repeat_indicator",
    "mmsi", "spare", "lon1", "lat1", "lon2", "lat2", "station_type",
    "type_and_cargo", "spare2", "txrx_mode", "interval_raw", "quiet",
    "spare3"};
    mxArray *theStruct;
 
    //If text descriptions of the fields are desired, then put them into a
    //structure.
    if(fieldDescriptions!=NULL) {
        const char *descriptionStrings[numberOfFields] = {
            "Identifier for Message 23; always 23.",//0
            "Repeat indicator. Used by the repeater to indicate how many times a message has been repeated. 0 = default; 3 = do not repeat any more.",//1
            "Maritime mobile service identity (MMSI) number of assigning station",//2
            "Spare. Should be set to zero",//3
            "Longitude (WGS-84) of area to which the assignment applies (degrees); upper right corner (North-East)",//4
            "Latitude (WGS-84) of area to which the assignment applies (degrees); upper right corner (North-East)",//5
            "Longitude (WGS-84) of area to which the assignment applies (degrees); lower left corner (South-West)",//6
            "Latitude (WGS-84) of area to which the assignment applies (degrees); lower left corner (South-West)",//7
            "Station type:\n0 = all types of mobiles (default)\n1 = Class A mobile stations only\n2 = all types of Class B mobile stations\n3 = SAR airborne mobile station\n4 = Class B 'SO' mobile stations only\n5 = Class B 'CS' shipborne mobile station only\n6 = inland waterways\n7 to 9 = regional use\n10 = base station coverage area (see Message 4 and Message 27)\n11 to 15 = for future use",//8
             "Type of ship and cargo type\n1-99 = as defined in Section 3.3.2 of Annex 8 of ITU-R M.1371-5\n100-199 = reserved, for regional use 200-255 = reserved, for future use",//9
            "Not used. Should be set to zero. Reserved for future use",//10
            "This parameter commands the respective stations to one of the following modes:\n0 = TxA/TxB, RxA/RxB (default)\n1 = TxA, RxA/RxB\n2 = TxB, RxA/RxB\n3 = reserved for future use",//11
             "This parameter commands the respective stations to the reporting interval given in Table 77 of ITU-R M.1371-5",//12
             "0 = default = no quiet time commanded; 1-15 = quiet time of 1 to 15 min",//3
             "Not used. Should be set to zero. Reserved for future use"//14
        };
        unsigned int i;

        *fieldDescriptions=mxCreateStructArray(2, dims, numberOfFields, fieldNames);

        for(i=0;i<numberOfFields;i++) {
            mxSetFieldByNumber(*fieldDescriptions,0,static_cast<int>(i),mxCreateCharMatrixFromStrings(1,&descriptionStrings[i]));
        }
    }
    
    //Create the Matlab structure array.
    theStruct=mxCreateStructArray(2, dims, numberOfFields, fieldNames);
    *decodedMessage=theStruct;
    
    //Fill all of the elements of the structure
    mxSetFieldByNumber(theStruct,0,0,intMat2MatlabDoubles(&msg->message_id,1,1));
    mxSetFieldByNumber(theStruct,0,1,intMat2MatlabDoubles(&msg->repeat_indicator,1,1));
    mxSetFieldByNumber(theStruct,0,2,intMat2MatlabDoubles(&msg->mmsi,1,1));
    
    mxSetFieldByNumber(theStruct,0,3,intMat2MatlabDoubles(&msg->spare,1,1));
    mxSetFieldByNumber(theStruct,0,4,doubleMat2Matlab(&msg->position1.lng_deg,1,1));
    mxSetFieldByNumber(theStruct,0,5,doubleMat2Matlab(&msg->position1.lat_deg,1,1));
    mxSetFieldByNumber(theStruct,0,6,doubleMat2Matlab(&msg->position2.lng_deg,1,1));
    mxSetFieldByNumber(theStruct,0,7,doubleMat2Matlab(&msg->position2.lat_deg,1,1));
    mxSetFieldByNumber(theStruct,0,8,intMat2MatlabDoubles(&msg->station_type,1,1));
    
    if(msg->type_and_cargo!=0) {
        mxSetFieldByNumber(theStruct,0,9,intMat2MatlabDoubles(&msg->type_and_cargo,1,1));
    } else {//No type and cargo information
        mxSetFieldByNumber(theStruct,0,9,mxCreateDoubleMatrix(0,0,mxREAL));
    }    
    mxSetFieldByNumber(theStruct,0,10,intMat2MatlabDoubles(&msg->spare2,1,1));

    mxSetFieldByNumber(theStruct,0,11,intMat2MatlabDoubles(&msg->txrx_mode,1,1));
    mxSetFieldByNumber(theStruct,0,12,intMat2MatlabDoubles(&msg->interval_raw,1,1));
    mxSetFieldByNumber(theStruct,0,13,intMat2MatlabDoubles(&msg->quiet,1,1));
    mxSetFieldByNumber(theStruct,0,14,intMat2MatlabDoubles(&msg->spare3,1,1));
}


void AIS_24_ToMatlab(libais::Ais24 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions) {
    const size_t numberOfFields=13;
    const mwSize dims[2] = {1, 1};
    //These are the names of all of the fields that will be added to the
    //Matlab structre array.
    const char *fieldNames[numberOfFields] = {"message_id", "repeat_indicator",
    "mmsi", "part_num", "name", "type_and_cargo", "vendor_id", "callsign",
    "dim_a", "dim_b", "dim_c", "dim_d", "spare"};
    mxArray *theStruct;
                 
    //If text descriptions of the fields are desired, then put them into a
    //structure.
    if(fieldDescriptions!=NULL) {
        const char *descriptionStrings[numberOfFields] = {
            "Identifier for Message 24; always 24.",//0
            "Repeat indicator. Used by the repeater to indicate how many times a message has been repeated. 0 = default; 3 = do not repeat any more.",//1
            "Maritime mobile service identity (MMSI) number",//2
            "Identifier for the message part number; always 0 for Part A, 1 for Part B",//3
            "Name of the MMSI-registered vessel. Maximum 20 characters 6-bit ASCII.",//4
            "Type of ship and cargo type\n1-99 = as defined in Section 3.3.2 of Annex 8 of ITU-R M.1371-5\n100-199 = reserved, for regional use 200-255 = reserved, for future use\nNot applicable to search and rescue (SAR) aircraft",//5
            "Unique identification of the Unit by a number as defined by the manufacturer",//6
            "Call sign. Craft associated with a parent vessel, should use 'A' followed by the last 6 digits of the MMSI of the parent vessel. Examples of these craft include towed vessels, rescue boats, tenders, lifeboats and liferafts.",//7
            "Distance from reference point to bow of ship (meters) 511m means 511m or longer. (If only relative reference point position known=0)",//8
            "Distance from reference point to aft of ship (meters) 511m means 511m or longer. ",//9
            "Distance from reference point to port side of ship (meters) 511m means 511m or longer. (If only relative reference point position known=0)",//10
            "Distance from reference point to starboard side of ship (meters) 511m means 511m or longer.",//11
            "Spare. Should be set to zero",//12     
        };
        unsigned int i;

        *fieldDescriptions=mxCreateStructArray(2, dims, numberOfFields, fieldNames);

        for(i=0;i<numberOfFields;i++) {
            mxSetFieldByNumber(*fieldDescriptions,0,static_cast<int>(i),mxCreateCharMatrixFromStrings(1,&descriptionStrings[i]));
        }
    }

    //Create the Matlab structure array.
    theStruct=mxCreateStructArray(2, dims, numberOfFields, fieldNames);
    *decodedMessage=theStruct;
    
    //Fill all of the elements of the structure
    mxSetFieldByNumber(theStruct,0,0,intMat2MatlabDoubles(&msg->message_id,1,1));
    mxSetFieldByNumber(theStruct,0,1,intMat2MatlabDoubles(&msg->repeat_indicator,1,1));
    mxSetFieldByNumber(theStruct,0,2,intMat2MatlabDoubles(&msg->mmsi,1,1));

    mxSetFieldByNumber(theStruct,0,3,intMat2MatlabDoubles(&msg->part_num,1,1));
    if(msg->part_num==0) {
        //Only a Name
        if(msg->name.compare("@@@@@@@@@@@@@@@@@@@@")!=0) {
            const char *charString=msg->name.c_str();
            mxSetFieldByNumber(theStruct,0,4,mxCreateCharMatrixFromStrings(1,&charString));
        } else {//If the name has actually been set
            mxSetFieldByNumber(theStruct,0,4,mxCreateDoubleMatrix(0,0,mxREAL));
        }
        
        mxSetFieldByNumber(theStruct,0,5,mxCreateDoubleMatrix(0,0,mxREAL));
        mxSetFieldByNumber(theStruct,0,6,mxCreateDoubleMatrix(0,0,mxREAL));
        mxSetFieldByNumber(theStruct,0,7,mxCreateDoubleMatrix(0,0,mxREAL));
        mxSetFieldByNumber(theStruct,0,8,mxCreateDoubleMatrix(0,0,mxREAL));
        mxSetFieldByNumber(theStruct,0,9,mxCreateDoubleMatrix(0,0,mxREAL));
        mxSetFieldByNumber(theStruct,0,10,mxCreateDoubleMatrix(0,0,mxREAL));
        mxSetFieldByNumber(theStruct,0,11,mxCreateDoubleMatrix(0,0,mxREAL));
        mxSetFieldByNumber(theStruct,0,12,mxCreateDoubleMatrix(0,0,mxREAL));
    } else if(msg->part_num==1) {
        //No Name 
        mxSetFieldByNumber(theStruct,0,4,mxCreateDoubleMatrix(0,0,mxREAL));//No name
        
        if(msg->type_and_cargo!=0) {
            mxSetFieldByNumber(theStruct,0,5,intMat2MatlabDoubles(&msg->type_and_cargo,1,1));
        } else {//If the cargo type is not available
            mxSetFieldByNumber(theStruct,0,5,mxCreateDoubleMatrix(0,0,mxREAL));
        }
        
        if(msg->vendor_id.compare("@@@@@@@")!=0) {
            const char *charString=msg->vendor_id.c_str();
            mxSetFieldByNumber(theStruct,0,6,mxCreateCharMatrixFromStrings(1,&charString));
        } else {//If no vendor ID is set
            mxSetFieldByNumber(theStruct,0,6,mxCreateDoubleMatrix(0,0,mxREAL));
        }
        
        if(msg->callsign.compare("@@@@@@@")!=0) {
            const char *charString=msg->callsign.c_str();
            mxSetFieldByNumber(theStruct,0,7,mxCreateCharMatrixFromStrings(1,&charString));
        } else {//If no callsign is set
            mxSetFieldByNumber(theStruct,0,7,mxCreateDoubleMatrix(0,0,mxREAL));
        }
        
        if(msg->dim_a==0&&msg->dim_b==0&&msg->dim_c==0&&msg->dim_d==0) {
            mxSetFieldByNumber(theStruct,0,8,mxCreateDoubleMatrix(0,0,mxREAL));
            mxSetFieldByNumber(theStruct,0,9,mxCreateDoubleMatrix(0,0,mxREAL));
            mxSetFieldByNumber(theStruct,0,10,mxCreateDoubleMatrix(0,0,mxREAL));
            mxSetFieldByNumber(theStruct,0,11,mxCreateDoubleMatrix(0,0,mxREAL));
        } else {//If the dimensions are actually given
            mxSetFieldByNumber(theStruct,0,8,intMat2MatlabDoubles(&msg->dim_a,1,1));
            mxSetFieldByNumber(theStruct,0,9,intMat2MatlabDoubles(&msg->dim_b,1,1));
            mxSetFieldByNumber(theStruct,0,10,intMat2MatlabDoubles(&msg->dim_c,1,1));
            mxSetFieldByNumber(theStruct,0,11,intMat2MatlabDoubles(&msg->dim_d,1,1));
        }
        
        mxSetFieldByNumber(theStruct,0,12,intMat2MatlabDoubles(&msg->spare,1,1));
    } else {
        //Undefined part_num value.
        mxSetFieldByNumber(theStruct,0,4,mxCreateDoubleMatrix(0,0,mxREAL));
        mxSetFieldByNumber(theStruct,0,5,mxCreateDoubleMatrix(0,0,mxREAL));
        mxSetFieldByNumber(theStruct,0,6,mxCreateDoubleMatrix(0,0,mxREAL));
        mxSetFieldByNumber(theStruct,0,7,mxCreateDoubleMatrix(0,0,mxREAL));
        mxSetFieldByNumber(theStruct,0,8,mxCreateDoubleMatrix(0,0,mxREAL));
        mxSetFieldByNumber(theStruct,0,9,mxCreateDoubleMatrix(0,0,mxREAL));
        mxSetFieldByNumber(theStruct,0,10,mxCreateDoubleMatrix(0,0,mxREAL));
        mxSetFieldByNumber(theStruct,0,11,mxCreateDoubleMatrix(0,0,mxREAL));
        mxSetFieldByNumber(theStruct,0,12,mxCreateDoubleMatrix(0,0,mxREAL));
    }
}

void AIS_25_ToMatlab(libais::Ais25 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions) {
    const size_t numberOfFields=6;
    const mwSize dims[2] = {1, 1};
    //These are the names of all of the fields that will be added to the
    //Matlab structre array.
    const char *fieldNames[numberOfFields] = {"message_id", "repeat_indicator",
    "mmsi", "dest_mmsi", "dac", "fi"};
    mxArray *theStruct;
             
    //If text descriptions of the fields are desired, then put them into a
    //structure.
    if(fieldDescriptions!=NULL) {
        const char *descriptionStrings[numberOfFields] = {
            "Identifier for Message 25; always 25.",//0
            "Repeat indicator. Used by the repeater to indicate how many times a message has been repeated. 0 = default; 3 = do not repeat any more.",//1
            "Maritime mobile service identity (MMSI) number of source station",//2
            "Destination MMSI, if used",//3
            "Designated area code",//4
            "Function identifier"//5
        };
        unsigned int i;

        *fieldDescriptions=mxCreateStructArray(2, dims, numberOfFields, fieldNames);

        for(i=0;i<numberOfFields;i++) {
            mxSetFieldByNumber(*fieldDescriptions,0,static_cast<int>(i),mxCreateCharMatrixFromStrings(1,&descriptionStrings[i]));
        }
    }

    //Create the Matlab structure array.
    theStruct=mxCreateStructArray(2, dims, numberOfFields, fieldNames);
    *decodedMessage=theStruct;
    
    //Fill all of the elements of the structure
    mxSetFieldByNumber(theStruct,0,0,intMat2MatlabDoubles(&msg->message_id,1,1));
    mxSetFieldByNumber(theStruct,0,1,intMat2MatlabDoubles(&msg->repeat_indicator,1,1));
    mxSetFieldByNumber(theStruct,0,2,intMat2MatlabDoubles(&msg->mmsi,1,1));

    if(msg->dest_mmsi_valid) {
        mxSetFieldByNumber(theStruct,0,3,intMat2MatlabDoubles(&msg->dest_mmsi,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,3,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(msg->use_app_id) {
        mxSetFieldByNumber(theStruct,0,4,intMat2MatlabDoubles(&msg->dac,1,1));
        mxSetFieldByNumber(theStruct,0,5,intMat2MatlabDoubles(&msg->fi,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,4,mxCreateDoubleMatrix(0,0,mxREAL));
        mxSetFieldByNumber(theStruct,0,5,mxCreateDoubleMatrix(0,0,mxREAL));
    }
}

void AIS_26_ToMatlab(libais::Ais26 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions) {
    const size_t numberOfFields=18;
    const mwSize dims[2] = {1, 1};
    //These are the names of all of the fields that will be added to the
    //Matlab structre array.
    const char *fieldNames[numberOfFields] = {"message_id", "repeat_indicator",
    "mmsi", "dest_mmsi", "dac", "fi", "commstate_flag", "sync_state",
    "slot_timeout", "received_stations", "slot_number", "utc_hour",
    "utc_min", "utc_spare", "slot_offset", "slot_increment",
    "slots_to_allocate", "keep_flag"};
    mxArray *theStruct;
    
    //If text descriptions of the fields are desired, then put them into a
    //structure.
    if(fieldDescriptions!=NULL) {
        const char *descriptionStrings[numberOfFields] = {
            "Identifier for Message 26; always 26.",//0
            "Repeat indicator. Used by the repeater to indicate how many times a message has been repeated. 0 = default; 3 = do not repeat any more.",//1
            "Maritime mobile service identity (MMSI) number of source station",//2
            "Destination MMSI, if used",//3
            "Designated area code",//4
            "Function identifier",//5
            "0 = SOTDMA communication state follows 1 = ITDMA communication state follows",//6
            "Synchronization state\n0 UTC direct\n1 UTC indirect\n2 Station is synchronized to a base station (base direct)\n3 Station is synchronized to another station based on the highest number of received stations or to another mobile station, which is directly synchronized to a base station",//7
            "(For SOTDMA) Specifies frames remaining until a new slot is selected\n0 means that this was the last transmission in this slot\n1-7 means that 1 to 7 frames respectively are left until slot change",//8
            "(For SOTDMA) Number of other stations (not own station) which the station currently is receiving (between 0 and 16 383).",//9
            "(For SOTDMA) Slot number used for this transmission (between 0 and 2 249).",//10
            "(For SOTDMA) Hour (0-23) in UTC time, if available",//11
            "(For SOTDMA) Minute (0-23) in UTC time, if available",//12
            "(For SOTDMA) Two extra (unused) bits that are associated with the time message",//13
            "(For SOTDMA) If the slot time-out value is 0 (zero) then the slot offset should indicate the offset to the slot in which transmission will occur during the next frame. If the slot offset is zero, the slot should be de-allocated after transmission.",//14
            "(For ITDMA) Offset to next slot to be used, or zero (0) if no more transmissions",//15
            "(For ITDMA) Number of consecutive slots to allocate.\n0 = 1 slot,\n1 = 2 slots,\n2 = 3 slots,\n3 = 4 slots,\n4 = 5 slots,\n5 = 1 slot; offset = slot increment + 8 192,\n6 = 2 slots; offset = slot increment + 8 192,\n7 = 3 slots; offset = slot increment + 8 192.\nUse of 5 to 7 removes the need for RATDMA broadcast for scheduled transmissions up to 6 min intervals",//16
            "(For ITDMA) Set to TRUE = 1 if the slot remains allocated for one additional frame"//17
        };
        unsigned int i;

        *fieldDescriptions=mxCreateStructArray(2, dims, numberOfFields, fieldNames);

        for(i=0;i<numberOfFields;i++) {
            mxSetFieldByNumber(*fieldDescriptions,0,static_cast<int>(i),mxCreateCharMatrixFromStrings(1,&descriptionStrings[i]));
        }
    }
    
    //Create the Matlab structure array.
    theStruct=mxCreateStructArray(2, dims, numberOfFields, fieldNames);
    *decodedMessage=theStruct;
    
    //Fill all of the elements of the structure
    mxSetFieldByNumber(theStruct,0,0,intMat2MatlabDoubles(&msg->message_id,1,1));
    mxSetFieldByNumber(theStruct,0,1,intMat2MatlabDoubles(&msg->repeat_indicator,1,1));
    mxSetFieldByNumber(theStruct,0,2,intMat2MatlabDoubles(&msg->mmsi,1,1));
    
    if(msg->dest_mmsi_valid) {
        mxSetFieldByNumber(theStruct,0,3,intMat2MatlabDoubles(&msg->dest_mmsi,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,3,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(msg->use_app_id) {
        mxSetFieldByNumber(theStruct,0,4,intMat2MatlabDoubles(&msg->dac,1,1));
        mxSetFieldByNumber(theStruct,0,5,intMat2MatlabDoubles(&msg->fi,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,4,mxCreateDoubleMatrix(0,0,mxREAL));
        mxSetFieldByNumber(theStruct,0,5,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    mxSetFieldByNumber(theStruct,0,6,intMat2MatlabDoubles(&msg->commstate_flag,1,1));
    mxSetFieldByNumber(theStruct,0,7,intMat2MatlabDoubles(&msg->sync_state,1,1));

    if(msg->slot_timeout_valid) {
        mxSetFieldByNumber(theStruct,0,8,intMat2MatlabDoubles(&msg->slot_timeout,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,8,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(msg->received_stations_valid) {
        mxSetFieldByNumber(theStruct,0,9,intMat2MatlabDoubles(&msg->received_stations,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,9,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(msg->slot_number_valid) {
        mxSetFieldByNumber(theStruct,0,10,intMat2MatlabDoubles(&msg->slot_number,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,10,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(msg->utc_valid) {
        mxSetFieldByNumber(theStruct,0,11,intMat2MatlabDoubles(&msg->utc_hour,1,1));
        mxSetFieldByNumber(theStruct,0,12,intMat2MatlabDoubles(&msg->utc_min,1,1));
        mxSetFieldByNumber(theStruct,0,13,intMat2MatlabDoubles(&msg->utc_spare,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,11,mxCreateDoubleMatrix(0,0,mxREAL));
        mxSetFieldByNumber(theStruct,0,12,mxCreateDoubleMatrix(0,0,mxREAL));
        mxSetFieldByNumber(theStruct,0,13,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(msg->slot_offset_valid) { 
        mxSetFieldByNumber(theStruct,0,14,intMat2MatlabDoubles(&msg->slot_offset,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,14,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(msg->slot_increment_valid) {
        mxSetFieldByNumber(theStruct,0,15,intMat2MatlabDoubles(&msg->slot_increment,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,15,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(msg->slots_to_allocate_valid) {
        mxSetFieldByNumber(theStruct,0,16,intMat2MatlabDoubles(&msg->slots_to_allocate,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,16,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(msg->keep_flag_valid) {
        mxSetFieldByNumber(theStruct,0,17,boolMat2Matlab(&msg->keep_flag,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,17,mxCreateDoubleMatrix(0,0,mxREAL));
    }
}

void AIS_27_ToMatlab(libais::Ais27 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions) {
    const size_t numberOfFields=12;
    const mwSize dims[2] = {1, 1};
    //These are the names of all of the fields that will be added to the
    //Matlab structre array.
    const char *fieldNames[numberOfFields] = {"message_id", "repeat_indicator",
    "mmsi", "position_accuracy", "raim", "nav_status", "lon", "lat", "sog",
    "cog", "gnss", "spare"};
    mxArray *theStruct;
               
    //If text descriptions of the fields are desired, then put them into a
    //structure.
    if(fieldDescriptions!=NULL) {
        const char *descriptionStrings[numberOfFields] = {
            "Identifier for Message 27; always 27.",//0
            "Repeat indicator. Always 3",//1
            "Maritime mobile service identity (MMSI) number of source station",//2
            "Position accuracy. 1=high (<=10m), 0=low (>10m), 0=default",//3
            "Receiver autonomous integrity monitoring (RAIM) flag of electronic position fixing device; 0 = RAIM not in use = default; 1 = RAIM in use.",//4
            "0=under way using engine\n1=at anchor,\n2=not under command\n3=restricted maneuverability\n4=constrained by her draught\n5=moored\n6=aground\n7=engaged in fishing\n8=under way sailing\n9=reserved for future amendment of navigational status for ships carrying DG, HS, or MP, or IMO hazard or pollutant category C, high speed craft (HSC)\n10=reserved for future amendment of navigational status for ships carrying dangerous goods (DG), harmful substances (HS) or marine pollutants (MP), or IMO hazard or pollutant category A, wing in ground (WIG)\n11=power-driven vessel towing astern (regional use)\n12=power-driven vessel pushing ahead or towing alongside (regional use)\n13 = reserved for future use\n14=AIS-SART (active), MOB-AIS, EPIRB-AIS\n15=undefined = default (also used by AIS-SART, MOB-AIS and EPIRB- AIS under test)",//5
            "Longitude East in degrees (WGS-84)",//6
            "Latitude North in degrees (WGS-84)",//7
            "Speed over ground in knots",//8
            "Course over ground in degrees East of North",//9
            "1=Reported position latency is less than 5 seconds, 0=Reported position latency is greater than five seconds",//10
            "Set to zero, to preserve byte boundaries"//11
        };
        unsigned int i;

        *fieldDescriptions=mxCreateStructArray(2, dims, numberOfFields, fieldNames);

        for(i=0;i<numberOfFields;i++) {
            mxSetFieldByNumber(*fieldDescriptions,0,static_cast<int>(i),mxCreateCharMatrixFromStrings(1,&descriptionStrings[i]));
        }
    }
    
    //Create the Matlab structure array.
    theStruct=mxCreateStructArray(2, dims, numberOfFields, fieldNames);
    *decodedMessage=theStruct;
    
    //Fill all of the elements of the structure
    mxSetFieldByNumber(theStruct,0,0,intMat2MatlabDoubles(&msg->message_id,1,1));
    mxSetFieldByNumber(theStruct,0,1,intMat2MatlabDoubles(&msg->repeat_indicator,1,1));
    mxSetFieldByNumber(theStruct,0,2,intMat2MatlabDoubles(&msg->mmsi,1,1));

    mxSetFieldByNumber(theStruct,0,3,intMat2MatlabDoubles(&msg->position_accuracy,1,1));
    mxSetFieldByNumber(theStruct,0,4,boolMat2Matlab(&msg->raim,1,1));
    mxSetFieldByNumber(theStruct,0,5,intMat2MatlabDoubles(&msg->nav_status,1,1));
    
    if(fabs(msg->position.lng_deg)<=180) {
        mxSetFieldByNumber(theStruct,0,6,doubleMat2Matlab(&msg->position.lng_deg,1,1));
    } else {//If longitude not available
        mxSetFieldByNumber(theStruct,0,6,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(msg->position.lat_deg<=90) {
        mxSetFieldByNumber(theStruct,0,7,doubleMat2Matlab(&msg->position.lat_deg,1,1));
    } else {//If latitude not available
        mxSetFieldByNumber(theStruct,0,7,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    mxSetFieldByNumber(theStruct,0,8,intMat2MatlabDoubles(&msg->sog,1,1));
    mxSetFieldByNumber(theStruct,0,9,intMat2MatlabDoubles(&msg->cog,1,1));
    mxSetFieldByNumber(theStruct,0,10,boolMat2Matlab(&msg->gnss,1,1));
    mxSetFieldByNumber(theStruct,0,11,intMat2MatlabDoubles(&msg->spare,1,1));
}

////////////////////////////////////////
///All of the functions for Message 6///
////////////////////////////////////////
void AIS_6_0_0_ToMatlab(libais::Ais6_0_0 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions) {
// This message format is described at
//http://www.e-navigation.nl/content/monitoring-aids-navigation
    const size_t numberOfFields=17;
    const mwSize dims[2] = {1, 1};
    //These are the names of all of the fields that will be added to the
    //Matlab structre array.
    const char *fieldNames[numberOfFields] = {"message_id", "repeat_indicator",
    "mmsi", "seq", "mmsi_dest", "retransmit", "spare", "dac", "fi",
    "sub_id", "voltage", "current", "dc_power_supply", "light_on",
    "battery_low","off_position", "spare2"};
    mxArray *theStruct;
               
    //If text descriptions of the fields are desired, then put them into a
    //structure.
    if(fieldDescriptions!=NULL) {
        const char *descriptionStrings[numberOfFields] = {
            "Identifier for Message 6; always 6.",//0
            "Repeat indicator. Used by the repeater to indicate how many times a message has been repeated. 0 = default; 3 = do not repeat any more.",//1
            "Maritime mobile service identity (MMSI) number of source station",//2
            "Sequence number (0-3); refer to Section 5.3.1 of Annex 2 of ITU-R M.1371-5",//3
            "MMSI number of destination station",//4
            "Retransmit flag should be set upon retransmission: 0 = no retransmission = default; 1 = retransmitted",//5
            "Not used. Should be zero. Reserved for future use",//6
            "Designated area code (DAC)",//7
            "Function identifier (FI)",//8
            "Sub-application ID",//9
            "Voltage Data. Lantern supply voltage data. Max 409.6V",//10
            "Current Data. Lantern drain current data. Max 102.3A",//11
            "Power Supply Type. 0=AC; 1=DC",//12
            "Light Status: 0=Light off; 1=Light on",//13
            "Battery Status: 0=Good; 1=Low Voltage",//14
            "Off Position Status: 0=On position; 1=Off position",//15
            "Set to zero; for future use"//16
        };
        unsigned int i;

        *fieldDescriptions=mxCreateStructArray(2, dims, numberOfFields, fieldNames);

        for(i=0;i<numberOfFields;i++) {
            mxSetFieldByNumber(*fieldDescriptions,0,static_cast<int>(i),mxCreateCharMatrixFromStrings(1,&descriptionStrings[i]));
        }
    }
    
    //Create the Matlab structure array.
    theStruct=mxCreateStructArray(2, dims, numberOfFields, fieldNames);
    *decodedMessage=theStruct;
    
    //Fill all of the elements of the structure
    //This first set is common to all AIS6 messages.
    mxSetFieldByNumber(theStruct,0,0,intMat2MatlabDoubles(&msg->message_id,1,1));
    mxSetFieldByNumber(theStruct,0,1,intMat2MatlabDoubles(&msg->repeat_indicator,1,1));
    mxSetFieldByNumber(theStruct,0,2,intMat2MatlabDoubles(&msg->mmsi,1,1));
    mxSetFieldByNumber(theStruct,0,3,intMat2MatlabDoubles(&msg->seq,1,1));
    mxSetFieldByNumber(theStruct,0,4,intMat2MatlabDoubles(&msg->mmsi_dest,1,1));
    mxSetFieldByNumber(theStruct,0,5,boolMat2Matlab(&msg->retransmit,1,1));
    mxSetFieldByNumber(theStruct,0,6,intMat2MatlabDoubles(&msg->spare,1,1));
    mxSetFieldByNumber(theStruct,0,7,intMat2MatlabDoubles(&msg->dac,1,1));
    mxSetFieldByNumber(theStruct,0,8,intMat2MatlabDoubles(&msg->fi,1,1));
    
    //Fields specific to the Ais6_0_0 type.
    mxSetFieldByNumber(theStruct,0,9,intMat2MatlabDoubles(&msg->sub_id,1,1));
    mxSetFieldByNumber(theStruct,0,10,floatMat2MatlabDoubles(&msg->voltage,1,1));
    mxSetFieldByNumber(theStruct,0,11,floatMat2MatlabDoubles(&msg->current,1,1));
    mxSetFieldByNumber(theStruct,0,12,boolMat2Matlab(&msg->dc_power_supply,1,1));
    mxSetFieldByNumber(theStruct,0,13,boolMat2Matlab(&msg->light_on,1,1));
    mxSetFieldByNumber(theStruct,0,14,boolMat2Matlab(&msg->battery_low,1,1));
    mxSetFieldByNumber(theStruct,0,15,boolMat2Matlab(&msg->off_position,1,1));
    
    mxSetFieldByNumber(theStruct,0,16,intMat2MatlabDoubles(&msg->spare2,1,1));
}

void AIS_6_1_0_ToMatlab(libais::Ais6_1_0 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions){
//This message format is described in Annex 5 of ITU-R M.1371-5.
    const size_t numberOfFields=13;
    const mwSize dims[2] = {1, 1};
    //These are the names of all of the fields that will be added to the
    //Matlab structre array.
    const char *fieldNames[numberOfFields] = {"message_id", "repeat_indicator",
    "mmsi", "seq", "mmsi_dest", "retransmit", "spare", "dac", "fi",
    "ack_required", "msg_seq", "text", "spare2"};
    mxArray *theStruct;

    //If text descriptions of the fields are desired, then put them into a
    //structure.
    if(fieldDescriptions!=NULL) {
        const char *descriptionStrings[numberOfFields] = {
            "Identifier for Message 6; always 6.",//0
            "Repeat indicator. Used by the repeater to indicate how many times a message has been repeated. 0 = default; 3 = do not repeat any more.",//1
            "Maritime mobile service identity (MMSI) number of source station",//2
            "Sequence number (0-3); refer to Section 5.3.1 of Annex 2 of ITU-R M.1371-5",//3
            "MMSI number of destination station",//4
            "Retransmit flag should be set upon retransmission: 0 = no retransmission = default; 1 = retransmitted",//5
            "Not used. Should be zero. Reserved for future use",//6
            "Designated area code (DAC)",//7
            "Function identifier (FI)",//8
            "Acknowledge required flag\n1 = reply is required, optional for addressed binary messages and not used for binary broadcast messages\n0 = reply is not required, optional for an addressed binary message and required for binary broadcast messages",//9
             "Sequence number to be incremented by the application. All zeros indicates that sequence numbers are not being used",//10
             "Text string",//11
             "Not used for data and should be set to zero. The number of bits should be either 0, 2, 4, or 6 to maintain byte boundaries."//12
        };
        unsigned int i;

        *fieldDescriptions=mxCreateStructArray(2, dims, numberOfFields, fieldNames);

        for(i=0;i<numberOfFields;i++) {
            mxSetFieldByNumber(*fieldDescriptions,0,static_cast<int>(i),mxCreateCharMatrixFromStrings(1,&descriptionStrings[i]));
        }
    }
    
    //Create the Matlab structure array.
    theStruct=mxCreateStructArray(2, dims, numberOfFields, fieldNames);
    *decodedMessage=theStruct;
    
    //Fill all of the elements of the structure
    //This first set is common to all AIS6 messages.
    mxSetFieldByNumber(theStruct,0,0,intMat2MatlabDoubles(&msg->message_id,1,1));
    mxSetFieldByNumber(theStruct,0,1,intMat2MatlabDoubles(&msg->repeat_indicator,1,1));
    mxSetFieldByNumber(theStruct,0,2,intMat2MatlabDoubles(&msg->mmsi,1,1));
    mxSetFieldByNumber(theStruct,0,3,intMat2MatlabDoubles(&msg->seq,1,1));
    mxSetFieldByNumber(theStruct,0,4,intMat2MatlabDoubles(&msg->mmsi_dest,1,1));
    mxSetFieldByNumber(theStruct,0,5,boolMat2Matlab(&msg->retransmit,1,1));
    mxSetFieldByNumber(theStruct,0,6,intMat2MatlabDoubles(&msg->spare,1,1));
    mxSetFieldByNumber(theStruct,0,7,intMat2MatlabDoubles(&msg->dac,1,1));
    mxSetFieldByNumber(theStruct,0,8,intMat2MatlabDoubles(&msg->fi,1,1));

    //Fields specific to the Ais6_1_0 type.
    mxSetFieldByNumber(theStruct,0,9,boolMat2Matlab(&msg->ack_required,1,1));
    mxSetFieldByNumber(theStruct,0,10,intMat2MatlabDoubles(&msg->msg_seq,1,1));
    {
        const char *charString=msg->text.c_str();
        mxSetFieldByNumber(theStruct,0,11,mxCreateCharMatrixFromStrings(1,&charString));
    }
    mxSetFieldByNumber(theStruct,0,12,intMat2MatlabDoubles(&msg->spare2,1,1));
}

void AIS_6_1_1_ToMatlab(libais::Ais6_1_1 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions){
//This message format is described in Annex 5 of ITU-R M.1371-1. It is not
//present in ITU-R M.1371-5
    const size_t numberOfFields=12;
    const mwSize dims[2] = {1, 1};
    //These are the names of all of the fields that will be added to the
    //Matlab structre array.
    const char *fieldNames[numberOfFields] = {"message_id", "repeat_indicator",
    "mmsi", "seq", "mmsi_dest", "retransmit", "spare", "dac", "fi",
    "ack_dac", "msg_seq", "spare2"};
    mxArray *theStruct;

    //If text descriptions of the fields are desired, then put them into a
    //structure.
    if(fieldDescriptions!=NULL) {
        const char *descriptionStrings[numberOfFields] = {
            "Identifier for Message 6; always 6.",//0
            "Repeat indicator. Used by the repeater to indicate how many times a message has been repeated. 0 = default; 3 = do not repeat any more.",//1
            "Maritime mobile service identity (MMSI) number of source station",//2
            "Sequence number (0-3); refer to Section 5.3.1 of Annex 2 of ITU-R M.1371-5",//3
            "MMSI number of destination station",//4
            "Retransmit flag should be set upon retransmission: 0 = no retransmission = default; 1 = retransmitted",//5
            "Not used. Should be zero. Reserved for future use",//6
            "Designated area code (DAC)",//7
            "Function identifier (FI)",//8
            "DAC being acknowledged",//9
            "Sequence number of the message being acknowledged",//10
            "Not used. Should be zero."//11
        };
        unsigned int i;

        *fieldDescriptions=mxCreateStructArray(2, dims, numberOfFields, fieldNames);

        for(i=0;i<numberOfFields;i++) {
            mxSetFieldByNumber(*fieldDescriptions,0,static_cast<int>(i),mxCreateCharMatrixFromStrings(1,&descriptionStrings[i]));
        }
    }
    
    //Create the Matlab structure array.
    theStruct=mxCreateStructArray(2, dims, numberOfFields, fieldNames);
    *decodedMessage=theStruct;
    
    //Fill all of the elements of the structure
    //This first set is common to all AIS6 messages.
    mxSetFieldByNumber(theStruct,0,0,intMat2MatlabDoubles(&msg->message_id,1,1));
    mxSetFieldByNumber(theStruct,0,1,intMat2MatlabDoubles(&msg->repeat_indicator,1,1));
    mxSetFieldByNumber(theStruct,0,2,intMat2MatlabDoubles(&msg->mmsi,1,1));
    mxSetFieldByNumber(theStruct,0,3,intMat2MatlabDoubles(&msg->seq,1,1));
    mxSetFieldByNumber(theStruct,0,4,intMat2MatlabDoubles(&msg->mmsi_dest,1,1));
    mxSetFieldByNumber(theStruct,0,5,boolMat2Matlab(&msg->retransmit,1,1));
    mxSetFieldByNumber(theStruct,0,6,intMat2MatlabDoubles(&msg->spare,1,1));
    mxSetFieldByNumber(theStruct,0,7,intMat2MatlabDoubles(&msg->dac,1,1));
    mxSetFieldByNumber(theStruct,0,8,intMat2MatlabDoubles(&msg->fi,1,1));

    //Fields specific to the Ais6_1_1 type.
    mxSetFieldByNumber(theStruct,0,9,intMat2MatlabDoubles(&msg->ack_dac,1,1));
    mxSetFieldByNumber(theStruct,0,10,intMat2MatlabDoubles(&msg->msg_seq,1,1));
    mxSetFieldByNumber(theStruct,0,11,intMat2MatlabDoubles(&msg->spare2,1,1));
}

void AIS_6_1_2_ToMatlab(libais::Ais6_1_2 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions){
//This message format is described in Annex 5 of ITU-R M.1371-5.
    const size_t numberOfFields=11;
    const mwSize dims[2] = {1, 1};
    //These are the names of all of the fields that will be added to the
    //Matlab structre array.
    const char *fieldNames[numberOfFields] = {"message_id", "repeat_indicator",
    "mmsi", "seq", "mmsi_dest", "retransmit", "spare", "dac", "fi",
    "req_dac", "req_fi"};
    mxArray *theStruct;
        
    //If text descriptions of the fields are desired, then put them into a
    //structure.
    if(fieldDescriptions!=NULL) {
        const char *descriptionStrings[numberOfFields] = {
            "Identifier for Message 6; always 6.",//0
            "Repeat indicator. Used by the repeater to indicate how many times a message has been repeated. 0 = default; 3 = do not repeat any more.",//1
            "Maritime mobile service identity (MMSI) number of source station",//2
            "Sequence number (0-3); refer to Section 5.3.1 of Annex 2 of ITU-R M.1371-5",//3
            "MMSI number of destination station",//4
            "Retransmit flag should be set upon retransmission: 0 = no retransmission = default; 1 = retransmitted",//5
            "Not used. Should be zero. Reserved for future use",//6
            "Designated area code (DAC)",//7
            "Function identifier (FI)",//8
            "Requested DAC: International application identifier (IAI), regional application identifier (RAI) or test",//9
            "Requested FI",//10
        };
        unsigned int i;

        *fieldDescriptions=mxCreateStructArray(2, dims, numberOfFields, fieldNames);

        for(i=0;i<numberOfFields;i++) {
            mxSetFieldByNumber(*fieldDescriptions,0,static_cast<int>(i),mxCreateCharMatrixFromStrings(1,&descriptionStrings[i]));
        }
    }
    
    //Create the Matlab structure array.
    theStruct=mxCreateStructArray(2, dims, numberOfFields, fieldNames);
    *decodedMessage=theStruct;
    
    //Fill all of the elements of the structure
    //This first set is common to all AIS6 messages.
    mxSetFieldByNumber(theStruct,0,0,intMat2MatlabDoubles(&msg->message_id,1,1));
    mxSetFieldByNumber(theStruct,0,1,intMat2MatlabDoubles(&msg->repeat_indicator,1,1));
    mxSetFieldByNumber(theStruct,0,2,intMat2MatlabDoubles(&msg->mmsi,1,1));
    mxSetFieldByNumber(theStruct,0,3,intMat2MatlabDoubles(&msg->seq,1,1));
    mxSetFieldByNumber(theStruct,0,4,intMat2MatlabDoubles(&msg->mmsi_dest,1,1));
    mxSetFieldByNumber(theStruct,0,5,boolMat2Matlab(&msg->retransmit,1,1));
    mxSetFieldByNumber(theStruct,0,6,intMat2MatlabDoubles(&msg->spare,1,1));
    mxSetFieldByNumber(theStruct,0,7,intMat2MatlabDoubles(&msg->dac,1,1));
    mxSetFieldByNumber(theStruct,0,8,intMat2MatlabDoubles(&msg->fi,1,1));
    
    //Fields specific to the Ais6_1_2 type.
    mxSetFieldByNumber(theStruct,0,9,intMat2MatlabDoubles(&msg->req_dac,1,1));
    mxSetFieldByNumber(theStruct,0,10,intMat2MatlabDoubles(&msg->req_fi,1,1));
}

void AIS_6_1_3_ToMatlab(libais::Ais6_1_3 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions){
//This message format is described in Annex 5 of ITU-R M.1371-5.
    const size_t numberOfFields=11;
    const mwSize dims[2] = {1, 1};
    //These are the names of all of the fields that will be added to the
    //Matlab structre array.
    const char *fieldNames[numberOfFields] = {"message_id", "repeat_indicator",
    "mmsi", "seq", "mmsi_dest", "retransmit", "spare", "dac", "fi",
    "req_dac", "spare2"};
    mxArray *theStruct;
        
    //If text descriptions of the fields are desired, then put them into a
    //structure.
    if(fieldDescriptions!=NULL) {
        const char *descriptionStrings[numberOfFields] = {
            "Identifier for Message 6; always 6.",//0
            "Repeat indicator. Used by the repeater to indicate how many times a message has been repeated. 0 = default; 3 = do not repeat any more.",//1
            "Maritime mobile service identity (MMSI) number of source station",//2
            "Sequence number (0-3); refer to Section 5.3.1 of Annex 2 of ITU-R M.1371-5",//3
            "MMSI number of destination station",//4
            "Retransmit flag should be set upon retransmission: 0 = no retransmission = default; 1 = retransmitted",//5
            "Not used. Should be zero. Reserved for future use",//6
            "Designated area code (DAC)",//7
            "Function identifier (FI)",//8
            "Requested DAC: International application identifier (IAI), regional application identifier (RAI) or test",//9
            "Not used, should be set to zero, reserved for future use"//10
        };
        unsigned int i;

        *fieldDescriptions=mxCreateStructArray(2, dims, numberOfFields, fieldNames);

        for(i=0;i<numberOfFields;i++) {
            mxSetFieldByNumber(*fieldDescriptions,0,static_cast<int>(i),mxCreateCharMatrixFromStrings(1,&descriptionStrings[i]));
        }
    }
    
    //Create the Matlab structure array.
    theStruct=mxCreateStructArray(2, dims, numberOfFields, fieldNames);
    *decodedMessage=theStruct;
    
    //Fill all of the elements of the structure
    //This first set is common to all AIS6 messages.
    mxSetFieldByNumber(theStruct,0,0,intMat2MatlabDoubles(&msg->message_id,1,1));
    mxSetFieldByNumber(theStruct,0,1,intMat2MatlabDoubles(&msg->repeat_indicator,1,1));
    mxSetFieldByNumber(theStruct,0,2,intMat2MatlabDoubles(&msg->mmsi,1,1));
    mxSetFieldByNumber(theStruct,0,3,intMat2MatlabDoubles(&msg->seq,1,1));
    mxSetFieldByNumber(theStruct,0,4,intMat2MatlabDoubles(&msg->mmsi_dest,1,1));
    mxSetFieldByNumber(theStruct,0,5,boolMat2Matlab(&msg->retransmit,1,1));
    mxSetFieldByNumber(theStruct,0,6,intMat2MatlabDoubles(&msg->spare,1,1));
    mxSetFieldByNumber(theStruct,0,7,intMat2MatlabDoubles(&msg->dac,1,1));
    mxSetFieldByNumber(theStruct,0,8,intMat2MatlabDoubles(&msg->fi,1,1));
    
    //Fields specific to the Ais6_1_3 type.
    mxSetFieldByNumber(theStruct,0,9,intMat2MatlabDoubles(&msg->req_dac,1,1));
    mxSetFieldByNumber(theStruct,0,10,intMat2MatlabDoubles(&msg->spare2,1,1));
}

void AIS_6_1_4_ToMatlab(libais::Ais6_1_4 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions){
//This message format is described in Annex 5 of ITU-R M.1371-5.
    const size_t numberOfFields=13;
    const mwSize dims[2] = {1, 1};
    //These are the names of all of the fields that will be added to the
    //Matlab structre array.
    const char *fieldNames[numberOfFields] = {"message_id", "repeat_indicator",
    "mmsi", "seq", "mmsi_dest", "retransmit", "spare", "dac", "fi",
    "ack_dac", "capabilities", "cap_reserved", "spare2"};
    mxArray *theStruct;

    //If text descriptions of the fields are desired, then put them into a
    //structure.
    if(fieldDescriptions!=NULL) {
        const char *descriptionStrings[numberOfFields] = {
            "Identifier for Message 6; always 6.",//0
            "Repeat indicator. Used by the repeater to indicate how many times a message has been repeated. 0 = default; 3 = do not repeat any more.",//1
            "Maritime mobile service identity (MMSI) number of source station",//2
            "Sequence number (0-3); refer to Section 5.3.1 of Annex 2 of ITU-R M.1371-5",
            "MMSI number of destination station",//4
            "Retransmit flag should be set upon retransmission: 0 = no retransmission = default; 1 = retransmitted",//5
            "Not used. Should be zero. Reserved for future use",//6
            "Designated area code (DAC)",//7
            "Function identifier (FI)",//8
            "DAC code replied",//9
            "A 64X1 array of boolean values, where a 1 indicates whether the FI value for that position is available. In Matlab, capabilities(1) corresponds to FI=0",//10
            "A 64X1 array of values reserved for future use.",//11
            "Not used, should be set to zero, reserved for future use"//12
        };
        unsigned int i;

        *fieldDescriptions=mxCreateStructArray(2, dims, numberOfFields, fieldNames);

        for(i=0;i<numberOfFields;i++) {
            mxSetFieldByNumber(*fieldDescriptions,0,static_cast<int>(i),mxCreateCharMatrixFromStrings(1,&descriptionStrings[i]));
        }
    }
    
    //Create the Matlab structure array.
    theStruct=mxCreateStructArray(2, dims, numberOfFields, fieldNames);
    *decodedMessage=theStruct;
    
    //Fill all of the elements of the structure
    //This first set is common to all AIS6 messages.
    mxSetFieldByNumber(theStruct,0,0,intMat2MatlabDoubles(&msg->message_id,1,1));
    mxSetFieldByNumber(theStruct,0,1,intMat2MatlabDoubles(&msg->repeat_indicator,1,1));
    mxSetFieldByNumber(theStruct,0,2,intMat2MatlabDoubles(&msg->mmsi,1,1));
    mxSetFieldByNumber(theStruct,0,3,intMat2MatlabDoubles(&msg->seq,1,1));
    mxSetFieldByNumber(theStruct,0,4,intMat2MatlabDoubles(&msg->mmsi_dest,1,1));
    mxSetFieldByNumber(theStruct,0,5,boolMat2Matlab(&msg->retransmit,1,1));
    mxSetFieldByNumber(theStruct,0,6,intMat2MatlabDoubles(&msg->spare,1,1));
    mxSetFieldByNumber(theStruct,0,7,intMat2MatlabDoubles(&msg->dac,1,1));
    mxSetFieldByNumber(theStruct,0,8,intMat2MatlabDoubles(&msg->fi,1,1));
    
    //Fields specific to the Ais6_1_4 type.
    mxSetFieldByNumber(theStruct,0,9,intMat2MatlabDoubles(&msg->ack_dac,1,1));
    mxSetFieldByNumber(theStruct,0,10,intMat2MatlabDoubles(msg->capabilities.data(),64,1));
    mxSetFieldByNumber(theStruct,0,11,intMat2MatlabDoubles(msg->cap_reserved.data(),64,1));
    mxSetFieldByNumber(theStruct,0,12,intMat2MatlabDoubles(&msg->spare2,1,1));
}

void AIS_6_1_5_ToMatlab(libais::Ais6_1_5 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions){
//This message format is described in Annex 5 of ITU-R M.1371-5.
    const size_t numberOfFields=15;
    const mwSize dims[2] = {1, 1};
    //These are the names of all of the fields that will be added to the
    //Matlab structre array.
    const char *fieldNames[numberOfFields] = {"message_id", "repeat_indicator",
    "mmsi", "seq", "mmsi_dest", "retransmit", "spare", "dac", "fi",
    "ack_dac", "ack_fi", "seq_num", "ai_available", "ai_response", "spare"
    };
    mxArray *theStruct;

    //If text descriptions of the fields are desired, then put them into a
    //structure.
    if(fieldDescriptions!=NULL) {
        const char *descriptionStrings[numberOfFields] = {
            "Identifier for Message 6; always 6.",//0
            "Repeat indicator. Used by the repeater to indicate how many times a message has been repeated. 0 = default; 3 = do not repeat any more.",//1
            "Maritime mobile service identity (MMSI) number of source station",//2
            "Sequence number (0-3); refer to Section 5.3.1 of Annex 2 of ITU-R M.1371-5",//3
            "MMSI number of destination station",//4
            "Retransmit flag should be set upon retransmission: 0 = no retransmission = default; 1 = retransmitted",//5
            "Not used. Should be zero. Reserved for future use",//6
            "Designated area code (DAC)",//7
            "Function identifier (FI)",//8
            "DAC code of received functional messsage",//9
            "FI code of received functional message",//10
            "Text sequence number: Sequence number in the message being acknowledged as properly received; 0=default",//11
            "AI available: 0=received but AI not available; 1 = AI available",//12
            "AI response: 0=unable to respond\n1 = reception acknowledged\n2 = response to follow\n3 = able to respond but currently inhibited\n4-7 = spare for future use",//13
            "Not used, should be set to zero, reserved for future use"//14
        };
        unsigned int i;

        *fieldDescriptions=mxCreateStructArray(2, dims, numberOfFields, fieldNames);

        for(i=0;i<numberOfFields;i++) {
            mxSetFieldByNumber(*fieldDescriptions,0,static_cast<int>(i),mxCreateCharMatrixFromStrings(1,&descriptionStrings[i]));
        }
    }
    
    //Create the Matlab structure array.
    theStruct=mxCreateStructArray(2, dims, numberOfFields, fieldNames);
    *decodedMessage=theStruct;
    
    //Fill all of the elements of the structure
    //This first set is common to all AIS6 messages.
    mxSetFieldByNumber(theStruct,0,0,intMat2MatlabDoubles(&msg->message_id,1,1));
    mxSetFieldByNumber(theStruct,0,1,intMat2MatlabDoubles(&msg->repeat_indicator,1,1));
    mxSetFieldByNumber(theStruct,0,2,intMat2MatlabDoubles(&msg->mmsi,1,1));
    mxSetFieldByNumber(theStruct,0,3,intMat2MatlabDoubles(&msg->seq,1,1));
    mxSetFieldByNumber(theStruct,0,4,intMat2MatlabDoubles(&msg->mmsi_dest,1,1));
    mxSetFieldByNumber(theStruct,0,5,boolMat2Matlab(&msg->retransmit,1,1));
    mxSetFieldByNumber(theStruct,0,6,intMat2MatlabDoubles(&msg->spare,1,1));
    mxSetFieldByNumber(theStruct,0,7,intMat2MatlabDoubles(&msg->dac,1,1));
    mxSetFieldByNumber(theStruct,0,8,intMat2MatlabDoubles(&msg->fi,1,1));
    
    //Fields specific to the Ais6_1_5 type.
    mxSetFieldByNumber(theStruct,0,9,intMat2MatlabDoubles(&msg->ack_dac,1,1));
    mxSetFieldByNumber(theStruct,0,10,intMat2MatlabDoubles(&msg->ack_fi,1,1));
    mxSetFieldByNumber(theStruct,0,11,intMat2MatlabDoubles(&msg->seq_num,1,1));
    mxSetFieldByNumber(theStruct,0,12,boolMat2Matlab(&msg->ai_available,1,1));
    mxSetFieldByNumber(theStruct,0,13,intMat2MatlabDoubles(&msg->ai_response,1,1));
    mxSetFieldByNumber(theStruct,0,14,intMat2MatlabDoubles(&msg->spare,1,1));
}

void AIS_6_1_12_ToMatlab(libais::Ais6_1_12 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions){
//This message format is described in IMO circular 236, which was revoked
//in IMO circular 289 in 2013.
    const size_t numberOfFields=25;
    const mwSize dims[2] = {1, 1};
    //These are the names of all of the fields that will be added to the
    //Matlab structre array.
    const char *fieldNames[numberOfFields] = {"message_id", "repeat_indicator",
    "mmsi", "seq", "mmsi_dest", "retransmit", "spare", "dac", "fi",
    "last_port", "utc_month_dep", "utc_day_dep", "utc_hour_dep",
    "utc_min_dep", "next_port", "utc_month_next", "utc_day_next",
    "utc_hour_next", "utc_min_next", "main_danger", "imo_cat", "un",
    "value", "value_unit", "spare2"};
    mxArray *theStruct;
    
    //If text descriptions of the fields are desired, then put them into a
    //structure.
    if(fieldDescriptions!=NULL) {
        const char *descriptionStrings[numberOfFields] = {
            "Identifier for Message 6; always 6.",//0
            "Repeat indicator. Used by the repeater to indicate how many times a message has been repeated. 0 = default; 3 = do not repeat any more.",//1
            "Maritime mobile service identity (MMSI) number of source station",//2
            "Sequence number (0-3); refer to Section 5.3.1 of Annex 2 of ITU-R M.1371-5",//3
            "MMSI number of destination station",//4
            "Retransmit flag should be set upon retransmission: 0 = no retransmission = default; 1 = retransmitted",//5
            "Not used. Should be zero. Reserved for future use",//6
            "Designated area code (DAC)",//7
            "Function identifier (FI)",//8
            "The United Nations (UN) locode of the last port of call.",//9
            "UTC month of departure (1-12) from last part of call",//10
            "UTC day of departure (1-31) from last part of call",//11
            "UTC hour of departure (0-23) from last part of call",//12
            "UTC minute of departure (0-59) from last part of call",//13
            "The UN locode of the next port of call.",//14
            "UTC estimated time of arrival (ETA) Month (1-12)",//15
            "UTC ETA day of departure (1-31)",//16
            "UTC ETA hour of departure (0-23)",//17
            "UTC ETA minute of departure (0-59)",//18
            "Main dangeroud good",//19
            "IMD category of main dangerous good",//20
            "UN Number (1-3363) of main dangeroud good",//21
            "Quantity of main dangerous good (1-1023) --Units specified in the next field",//22
            "Units of quantity of main dangerous good\n1 = in kg\n2 = in tons (10e3 kg)\n3 = in 1000 tons (10e6 kg)",//23
            "Not used. Should be set to zero"//24
        };
        unsigned int i;

        *fieldDescriptions=mxCreateStructArray(2, dims, numberOfFields, fieldNames);

        for(i=0;i<numberOfFields;i++) {
            mxSetFieldByNumber(*fieldDescriptions,0,static_cast<int>(i),mxCreateCharMatrixFromStrings(1,&descriptionStrings[i]));
        }
    }
    
    //Create the Matlab structure array.
    theStruct=mxCreateStructArray(2, dims, numberOfFields, fieldNames);
    *decodedMessage=theStruct;
    
    //Fill all of the elements of the structure
    //This first set is common to all AIS6 messages.
    mxSetFieldByNumber(theStruct,0,0,intMat2MatlabDoubles(&msg->message_id,1,1));
    mxSetFieldByNumber(theStruct,0,1,intMat2MatlabDoubles(&msg->repeat_indicator,1,1));
    mxSetFieldByNumber(theStruct,0,2,intMat2MatlabDoubles(&msg->mmsi,1,1));
    mxSetFieldByNumber(theStruct,0,3,intMat2MatlabDoubles(&msg->seq,1,1));
    mxSetFieldByNumber(theStruct,0,4,intMat2MatlabDoubles(&msg->mmsi_dest,1,1));
    mxSetFieldByNumber(theStruct,0,5,boolMat2Matlab(&msg->retransmit,1,1));
    mxSetFieldByNumber(theStruct,0,6,intMat2MatlabDoubles(&msg->spare,1,1));
    mxSetFieldByNumber(theStruct,0,7,intMat2MatlabDoubles(&msg->dac,1,1));
    mxSetFieldByNumber(theStruct,0,8,intMat2MatlabDoubles(&msg->fi,1,1));
    
    //Fields specific to the Ais6_1_12 type.
    if(msg->last_port.compare("@@@@@")!=0) {
        const char *charString=msg->last_port.c_str();
        mxSetFieldByNumber(theStruct,0,9,mxCreateCharMatrixFromStrings(1,&charString));
    } else {//If the last part of call is not available.
        mxSetFieldByNumber(theStruct,0,9,mxCreateDoubleMatrix(0,0,mxREAL));
    }

    if(msg->utc_month_dep!=0) {
        mxSetFieldByNumber(theStruct,0,10,intMat2MatlabDoubles(&msg->utc_month_dep,1,1));
    } else {//Month of departure not available
        mxSetFieldByNumber(theStruct,0,10,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(msg->utc_day_dep!=0) {
        mxSetFieldByNumber(theStruct,0,11,intMat2MatlabDoubles(&msg->utc_day_dep,1,1));
    } else {//Day of departure not available
        mxSetFieldByNumber(theStruct,0,11,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(msg->utc_hour_dep<24) {
        mxSetFieldByNumber(theStruct,0,12,intMat2MatlabDoubles(&msg->utc_hour_dep,1,1));
    } else {//Hour of departure not available
        mxSetFieldByNumber(theStruct,0,12,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(msg->utc_min_dep<60) {
        mxSetFieldByNumber(theStruct,0,13,intMat2MatlabDoubles(&msg->utc_min_dep,1,1));
    } else {//Minute of departure not available
        mxSetFieldByNumber(theStruct,0,13,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(msg->next_port.compare("@@@@@")!=0) {
        const char *charString=msg->next_port.c_str();
        mxSetFieldByNumber(theStruct,0,14,mxCreateCharMatrixFromStrings(1,&charString)); 
    } else {//Next port of call is not available
        mxSetFieldByNumber(theStruct,0,14,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(msg->utc_month_next!=0) {
        mxSetFieldByNumber(theStruct,0,15,intMat2MatlabDoubles(&msg->utc_month_next,1,1));
    } else {//Month of departure not available
        mxSetFieldByNumber(theStruct,0,15,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(msg->utc_day_next!=0) {
        mxSetFieldByNumber(theStruct,0,16,intMat2MatlabDoubles(&msg->utc_day_next,1,1));
    } else {//Day of departure not available
        mxSetFieldByNumber(theStruct,0,16,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(msg->utc_hour_next<24) {
        mxSetFieldByNumber(theStruct,0,17,intMat2MatlabDoubles(&msg->utc_hour_next,1,1));
    } else {//Hour of departure not available
        mxSetFieldByNumber(theStruct,0,17,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(msg->utc_min_next<60) {
        mxSetFieldByNumber(theStruct,0,18,intMat2MatlabDoubles(&msg->utc_min_next,1,1));
    } else {//Minute of departure not available
        mxSetFieldByNumber(theStruct,0,18,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(msg->main_danger.compare("@@@@@@@@@@@@@@@@@@@@")!=0) {
        const char *charString=msg->main_danger.c_str();
        mxSetFieldByNumber(theStruct,0,19,mxCreateCharMatrixFromStrings(1,&charString));
    } else {//Main dangerous good name not available
        mxSetFieldByNumber(theStruct,0,19,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(msg->imo_cat.compare("@@@@@")!=0) {
        const char *charString=msg->imo_cat.c_str();
        mxSetFieldByNumber(theStruct,0,20,mxCreateCharMatrixFromStrings(1,&charString));
    } else {//IMD Category of main dangerous good not available
        mxSetFieldByNumber(theStruct,0,20,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(msg->un!=0) {
        mxSetFieldByNumber(theStruct,0,21,intMat2MatlabDoubles(&msg->un,1,1));
    } else {//UN Number of main dangeroud good not available
        mxSetFieldByNumber(theStruct,0,21,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(msg->value!=0) {
        mxSetFieldByNumber(theStruct,0,22,intMat2MatlabDoubles(&msg->value,1,1));
    } else {//Quantity of main dangerous good not available
        mxSetFieldByNumber(theStruct,0,22,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(msg->value_unit!=0) {
        mxSetFieldByNumber(theStruct,0,23,intMat2MatlabDoubles(&msg->value_unit,1,1));
    } else {//Units of quantity of main dangeroud good not available
        mxSetFieldByNumber(theStruct,0,23,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    mxSetFieldByNumber(theStruct,0,24,intMat2MatlabDoubles(&msg->spare2,1,1));
}

void AIS_6_1_14_ToMatlab(libais::Ais6_1_14 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions){
//This message format is described in IMO circular 236, which was revoked
//in IMO circular 289 in 2013.
    const size_t numberOfFields=12;
    const mwSize dims[2] = {1, 1};
    //These are the names of all of the fields that will be added to the
    //Matlab structre array.
    const char *fieldNames[numberOfFields] = {"message_id", "repeat_indicator",
    "mmsi", "seq", "mmsi_dest", "retransmit", "spare", "dac", "fi",
    "utc_month", "utc_day", "windows"};
    const size_t numberOfWindowFields=8;
    mwSize windowDims[2];
    const char *windowFieldNames[numberOfWindowFields] = {"lon", "lat", "utc_hour_from",
    "utc_min_from", "utc_hour_to", "utc_min_to", "cur_dir", "cur_speed"};
    mxArray *theStruct, *windowArray;
    size_t curWindow;

    //If text descriptions of the fields are desired, then put them into a
    //structure.
    if(fieldDescriptions!=NULL) {
        const char *descriptionStrings[numberOfFields] = {
            "Identifier for Message 6; always 6.",//0
            "Repeat indicator. Used by the repeater to indicate how many times a message has been repeated. 0 = default; 3 = do not repeat any more.",//1
            "Maritime mobile service identity (MMSI) number of source station",//2
            "Sequence number (0-3); refer to Section 5.3.1 of Annex 2 of ITU-R M.1371-5",//3
            "MMSI number of destination station",//4
            "Retransmit flag should be set upon retransmission: 0 = no retransmission = default; 1 = retransmitted",//5
            "Not used. Should be zero. Reserved for future use",//6
            "Designated area code (DAC)",//7
            "Function identifier (FI)",//8
            "UTC Month (1-12)",//9
            "UTC Day (1-31)",//10
            "An array of tidal windows at different times consisting of components:\nLongitude East in degrees (WGS-84)\nLatitude North in degrees (WGS-84)\nUTC hours from (0-23)\nUTC minutes from (0-59)\nUTC hours to (0-23)\nUTC minutes to (0-59)\nDirection of the current in degrees East of North\nSpeed of the current in knots (0-12.6)"//11
        };
        unsigned int i;

        *fieldDescriptions=mxCreateStructArray(2, dims, numberOfFields, fieldNames);

        for(i=0;i<numberOfFields;i++) {
            mxSetFieldByNumber(*fieldDescriptions,0,static_cast<int>(i),mxCreateCharMatrixFromStrings(1,&descriptionStrings[i]));
        }
    }
    
    //Create the Matlab structure array.
    theStruct=mxCreateStructArray(2, dims, numberOfFields, fieldNames);
    *decodedMessage=theStruct;
    
    //Fill all of the elements of the structure
    //This first set is common to all AIS6 messages.
    mxSetFieldByNumber(theStruct,0,0,intMat2MatlabDoubles(&msg->message_id,1,1));
    mxSetFieldByNumber(theStruct,0,1,intMat2MatlabDoubles(&msg->repeat_indicator,1,1));
    mxSetFieldByNumber(theStruct,0,2,intMat2MatlabDoubles(&msg->mmsi,1,1));
    mxSetFieldByNumber(theStruct,0,3,intMat2MatlabDoubles(&msg->seq,1,1));
    mxSetFieldByNumber(theStruct,0,4,intMat2MatlabDoubles(&msg->mmsi_dest,1,1));
    mxSetFieldByNumber(theStruct,0,5,boolMat2Matlab(&msg->retransmit,1,1));
    mxSetFieldByNumber(theStruct,0,6,intMat2MatlabDoubles(&msg->spare,1,1));
    mxSetFieldByNumber(theStruct,0,7,intMat2MatlabDoubles(&msg->dac,1,1));
    mxSetFieldByNumber(theStruct,0,8,intMat2MatlabDoubles(&msg->fi,1,1));
    
    //Fields specific to the Ais6_1_14 type.
    if(msg->utc_month!=0) {
        mxSetFieldByNumber(theStruct,0,9,intMat2MatlabDoubles(&msg->utc_month,1,1));
    } else {//If UTC month is not available
        mxSetFieldByNumber(theStruct,0,9,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(msg->utc_day!=0) {
        mxSetFieldByNumber(theStruct,0,10,intMat2MatlabDoubles(&msg->utc_day,1,1));
    } else {//If UTC day is not available.
        mxSetFieldByNumber(theStruct,0,10,mxCreateDoubleMatrix(0,0,mxREAL));
    }

    windowDims[0]=msg->windows.size();
    windowDims[1]=1;
    windowArray=mxCreateStructArray(2, dims, numberOfWindowFields, windowFieldNames);
    mxSetFieldByNumber(theStruct,0,11,windowArray);

    for (curWindow = 0; curWindow < msg->windows.size(); curWindow++) {
        if(fabs(msg->windows[curWindow].position.lng_deg)<=180) {
            mxSetFieldByNumber(windowArray,curWindow,0,doubleMat2Matlab(&msg->windows[curWindow].position.lng_deg,1,1));
        } else {//Longitude is not available
            mxSetFieldByNumber(windowArray,curWindow,0,mxCreateDoubleMatrix(0,0,mxREAL));
        }
        if(fabs(msg->windows[curWindow].position.lat_deg)<=90) {
            mxSetFieldByNumber(windowArray,curWindow,1,doubleMat2Matlab(&msg->windows[curWindow].position.lat_deg,1,1));
        } else {//Latitude is not available
            mxSetFieldByNumber(windowArray,curWindow,1,mxCreateDoubleMatrix(0,0,mxREAL));
        }
        if(msg->windows[curWindow].utc_hour_from<24) {
            mxSetFieldByNumber(windowArray,curWindow,2,intMat2MatlabDoubles(&msg->windows[curWindow].utc_hour_from,1,1));
        } else {// Hour from is not available
            mxSetFieldByNumber(windowArray,curWindow,2,mxCreateDoubleMatrix(0,0,mxREAL)); 
        }
        if(msg->windows[curWindow].utc_min_from<60) {
            mxSetFieldByNumber(windowArray,curWindow,3,intMat2MatlabDoubles(&msg->windows[curWindow].utc_min_from,1,1));
        } else {//Minute from is not avaiable
            mxSetFieldByNumber(windowArray,curWindow,3,mxCreateDoubleMatrix(0,0,mxREAL)); 
        }
        if(msg->windows[curWindow].utc_hour_to<24) {
            mxSetFieldByNumber(windowArray,curWindow,4,intMat2MatlabDoubles(&msg->windows[curWindow].utc_hour_to,1,1));
        } else {// Hour from is not available
            mxSetFieldByNumber(windowArray,curWindow,4,mxCreateDoubleMatrix(0,0,mxREAL)); 
        }
        if(msg->windows[curWindow].utc_min_to<60) {
            mxSetFieldByNumber(windowArray,curWindow,5,intMat2MatlabDoubles(&msg->windows[curWindow].utc_min_to,1,1));
        } else {//Minute from is not avaiable
            mxSetFieldByNumber(windowArray,curWindow,5,mxCreateDoubleMatrix(0,0,mxREAL)); 
        }
        if(msg->windows[curWindow].cur_dir<360) {
            mxSetFieldByNumber(windowArray,curWindow,6,intMat2MatlabDoubles(&msg->windows[curWindow].cur_dir,1,1));
        } else {//Direction of the current not available
            mxSetFieldByNumber(windowArray,curWindow,6,mxCreateDoubleMatrix(0,0,mxREAL)); 
        }
        if(msg->windows[curWindow].cur_speed!=12.7) {
            mxSetFieldByNumber(windowArray,curWindow,7,floatMat2MatlabDoubles(&msg->windows[curWindow].cur_speed,1,1));
        } else {//Speed of the current not available
            mxSetFieldByNumber(windowArray,curWindow,7,mxCreateDoubleMatrix(0,0,mxREAL)); 
        }
    }
}

void AIS_6_1_18_ToMatlab(libais::Ais6_1_18 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions){
//This message format is described in IMO circular 289
    const size_t numberOfFields=19;
    const mwSize dims[2] = {1, 1};
    //These are the names of all of the fields that will be added to the
    //Matlab structre array.
    const char *fieldNames[numberOfFields] = {"message_id", "repeat_indicator",
    "mmsi", "seq", "mmsi_dest", "retransmit", "spare", "dac", "fi",
    "link_id", "utc_month", "utc_day", "utc_hour", "utc_min", "port_berth",
    "dest", "lon", "lat", "spare2"};
    mxArray *theStruct;
    
    //If text descriptions of the fields are desired, then put them into a
    //structure.
    if(fieldDescriptions!=NULL) {
        const char *descriptionStrings[numberOfFields] = {
            "Identifier for Message 6; always 6.",//0
            "Repeat indicator. Used by the repeater to indicate how many times a message has been repeated. 0 = default; 3 = do not repeat any more.",//1
            "Maritime mobile service identity (MMSI) number of source station",//2
            "Sequence number (0-3); refer to Section 5.3.1 of Annex 2 of ITU-R M.1371-5",//3
            "MMSI number of destination station",//4
            "Retransmit flag should be set upon retransmission: 0 = no retransmission = default; 1 = retransmitted",//5
            "Not used. Should be zero. Reserved for future use",//6
            "Designated area code (DAC)",//7
            "Function identifier (FI)",//8
            "A source specific running number (1-1023), unique across all binary messages equipped with Message Linkage ID. Used to link additional information to the message by a Text Description message. The Message Linkage ID and the first six digits of the source MMSI uniquely identify the sent message.",//9
            "UTC clearance month to enter port (1-12)",//11
            "UTC clearance day to enter port (1-12)",//12
            "UTC clearance hour to enter port (1-12)",//13
            "UTC clearance minute to enter port (1-12)",//14
            "Name of the port and berth",//15
            "UN locode of the destination",//16
            "Longitude East in degrees (WGS-84)",//17
            "Latitude North in degrees (WGS-84)",//18
        };
        unsigned int i;

        *fieldDescriptions=mxCreateStructArray(2, dims, numberOfFields, fieldNames);

        for(i=0;i<numberOfFields;i++) {
            mxSetFieldByNumber(*fieldDescriptions,0,static_cast<int>(i),mxCreateCharMatrixFromStrings(1,&descriptionStrings[i]));
        }
    }
    
    //Create the Matlab structure array.
    theStruct=mxCreateStructArray(2, dims, numberOfFields, fieldNames);
    *decodedMessage=theStruct;
    
    //Fill all of the elements of the structure
    //This first set is common to all AIS6 messages.
    mxSetFieldByNumber(theStruct,0,0,intMat2MatlabDoubles(&msg->message_id,1,1));
    mxSetFieldByNumber(theStruct,0,1,intMat2MatlabDoubles(&msg->repeat_indicator,1,1));
    mxSetFieldByNumber(theStruct,0,2,intMat2MatlabDoubles(&msg->mmsi,1,1));
    mxSetFieldByNumber(theStruct,0,3,intMat2MatlabDoubles(&msg->seq,1,1));
    mxSetFieldByNumber(theStruct,0,4,intMat2MatlabDoubles(&msg->mmsi_dest,1,1));
    mxSetFieldByNumber(theStruct,0,5,boolMat2Matlab(&msg->retransmit,1,1));
    mxSetFieldByNumber(theStruct,0,6,intMat2MatlabDoubles(&msg->spare,1,1));
    mxSetFieldByNumber(theStruct,0,7,intMat2MatlabDoubles(&msg->dac,1,1));
    mxSetFieldByNumber(theStruct,0,8,intMat2MatlabDoubles(&msg->fi,1,1));
    
    //Fields specific to the Ais6_1_18 type.
    if(msg->link_id!=0) {
        mxSetFieldByNumber(theStruct,0,9,intMat2MatlabDoubles(&msg->link_id,1,1));
    } else {//No link ID
        mxSetFieldByNumber(theStruct,0,9,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(msg->utc_month!=0) {
        mxSetFieldByNumber(theStruct,0,10,intMat2MatlabDoubles(&msg->utc_month,1,1));
    } else {//If UTC month is not available
        mxSetFieldByNumber(theStruct,0,10,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(msg->utc_day!=0) {
        mxSetFieldByNumber(theStruct,0,11,intMat2MatlabDoubles(&msg->utc_day,1,1));
    } else {//If UTC day is not available.
        mxSetFieldByNumber(theStruct,0,11,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(msg->utc_hour<24) {
        mxSetFieldByNumber(theStruct,0,12,intMat2MatlabDoubles(&msg->utc_hour,1,1));
    } else {//If UTC hour not available
        mxSetFieldByNumber(theStruct,0,12,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(msg->utc_min<60) {
        mxSetFieldByNumber(theStruct,0,13,intMat2MatlabDoubles(&msg->utc_min,1,1));
    } else {//If UTC minute not available
        mxSetFieldByNumber(theStruct,0,13,mxCreateDoubleMatrix(0,0,mxREAL));
    }

    if(msg->port_berth.compare("@@@@@@@@@@@@@@@@@@@@")!=0) {
        const char *charString=msg->port_berth.c_str();
        mxSetFieldByNumber(theStruct,0,14,mxCreateCharMatrixFromStrings(1,&charString));
    } else {//if name of port and berth is not available
        mxSetFieldByNumber(theStruct,0,14,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(msg->dest.compare("@@@@@")!=0) {
        const char *charString=msg->dest.c_str();
        mxSetFieldByNumber(theStruct,0,15,mxCreateCharMatrixFromStrings(1,&charString));
    } else {//No UN destination locode avaiable
        mxSetFieldByNumber(theStruct,0,15,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    if(fabs(msg->position.lng_deg)<=180) {
        mxSetFieldByNumber(theStruct,0,16,doubleMat2Matlab(&msg->position.lng_deg,1,1));
    } else {//No longitude avaiable
        mxSetFieldByNumber(theStruct,0,16,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    if(fabs(msg->position.lat_deg)<=90) {
        mxSetFieldByNumber(theStruct,0,17,doubleMat2Matlab(&msg->position.lat_deg,1,1));
    } else {//No latitude available
        mxSetFieldByNumber(theStruct,0,17,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    mxSetFieldByNumber(theStruct,0,18,intMat2MatlabDoubles(&msg->spare2[1],2,1));
}

void AIS_6_1_20_ToMatlab(libais::Ais6_1_20 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions){
//This message format is described in IMO circular 289
    const size_t numberOfFields=21;
    const mwSize dims[2] = {1, 1};
    //These are the names of all of the fields that will be added to the
    //Matlab structre array.
    const char *fieldNames[numberOfFields] = {"message_id", "repeat_indicator",
    "mmsi", "seq", "mmsi_dest", "retransmit", "spare", "dac", "fi",
    "link_id", "length", "depth", "mooring_position", "utc_month", "utc_day",
    "utc_hour", "utc_min", "services", "name", "lon", "lat"};
    mxArray *theStruct;

    //If text descriptions of the fields are desired, then put them into a
    //structure.
    if(fieldDescriptions!=NULL) {
        const char *descriptionStrings[numberOfFields] = {
            "Identifier for Message 6; always 6.",//0
            "Repeat indicator. Used by the repeater to indicate how many times a message has been repeated. 0 = default; 3 = do not repeat any more.",//1
            "Maritime mobile service identity (MMSI) number of source station",//2
            "Sequence number (0-3); refer to Section 5.3.1 of Annex 2 of ITU-R M.1371-5",//3
            "MMSI number of destination station",//4
            "Retransmit flag should be set upon retransmission: 0 = no retransmission = default; 1 = retransmitted",//5
            "Not used. Should be zero. Reserved for future use",//6
            "Designated area code (DAC)",//7
            "Function identifier (FI)",//8
            "A source specific running number (1-1023), unique across all binary messages equipped with Message Linkage ID. Used to link additional information to the message by a Text Description message. The Message Linkage ID and the first six digits of the source MMSI uniquely identify the sent message.",//9
            "Berth length in meters. 1 - 510 metres, 511 = 511 metres or greater",//10
            "Water depth at berth in meters. 0.1-25.4 meters. 25.5= 25.5 meters or deeper",//11
            "Mooring Position:\n1 = port-side to\n2 = starboard-side to\n3 = Mediterranean mooring\n4 = mooring buoy\n5 = anchorage\n6 - 7 (reserved for future use)",//12
            "Berth Date/ Time UTC month (1-12)",//13
            "Berth Date/ Time UTC day (1-31)",//14
            "Berth Date/ Time UTC hour (0-23)",//15
            "Berth Date/ Time UTC minute (0-59)",//16
            "A 26X1 vector of values indicating the availability of services specified in IMO circular 289. For each entry:\n0 = service not available or requested = default\n1 = service available\n2 = no data or unknown\n3 = not to be used",//17
            "The name of the berth",//18
            "Central longitude East of the berth in degrees (WGS-84)",//19
            "Central latitude North of the berth in degrees (WGS-84)",//20        
        };
        unsigned int i;

        *fieldDescriptions=mxCreateStructArray(2, dims, numberOfFields, fieldNames);

        for(i=0;i<numberOfFields;i++) {
            mxSetFieldByNumber(*fieldDescriptions,0,static_cast<int>(i),mxCreateCharMatrixFromStrings(1,&descriptionStrings[i]));
        }
    }
    
    //Create the Matlab structure array.
    theStruct=mxCreateStructArray(2, dims, numberOfFields, fieldNames);
    *decodedMessage=theStruct;
    
    //Fill all of the elements of the structure
    //This first set is common to all AIS6 messages.
    mxSetFieldByNumber(theStruct,0,0,intMat2MatlabDoubles(&msg->message_id,1,1));
    mxSetFieldByNumber(theStruct,0,1,intMat2MatlabDoubles(&msg->repeat_indicator,1,1));
    mxSetFieldByNumber(theStruct,0,2,intMat2MatlabDoubles(&msg->mmsi,1,1));
    mxSetFieldByNumber(theStruct,0,3,intMat2MatlabDoubles(&msg->seq,1,1));
    mxSetFieldByNumber(theStruct,0,4,intMat2MatlabDoubles(&msg->mmsi_dest,1,1));
    mxSetFieldByNumber(theStruct,0,5,boolMat2Matlab(&msg->retransmit,1,1));
    mxSetFieldByNumber(theStruct,0,6,intMat2MatlabDoubles(&msg->spare,1,1));
    mxSetFieldByNumber(theStruct,0,7,intMat2MatlabDoubles(&msg->dac,1,1));
    mxSetFieldByNumber(theStruct,0,8,intMat2MatlabDoubles(&msg->fi,1,1));
    
    //Fields specific to the Ais6_1_20 type.
    if(msg->link_id!=0) {
        mxSetFieldByNumber(theStruct,0,9,intMat2MatlabDoubles(&msg->link_id,1,1));
    } else {//No link ID
        mxSetFieldByNumber(theStruct,0,9,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(msg->length!=0) {
        mxSetFieldByNumber(theStruct,0,10,intMat2MatlabDoubles(&msg->length,1,1));
    } else {//No berth length
        mxSetFieldByNumber(theStruct,0,10,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(msg->depth!=0) {
        mxSetFieldByNumber(theStruct,0,11,floatMat2MatlabDoubles(&msg->depth,1,1));
    } else {//No berth depth
        mxSetFieldByNumber(theStruct,0,11,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(msg->mooring_position!=0) {
        mxSetFieldByNumber(theStruct,0,12,intMat2MatlabDoubles(&msg->mooring_position,1,1));
    } else {//No mooring position specified
        mxSetFieldByNumber(theStruct,0,12,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(msg->utc_month!=0) {
        mxSetFieldByNumber(theStruct,0,13,intMat2MatlabDoubles(&msg->utc_month,1,1));
    } else {//If UTC month is not available
        mxSetFieldByNumber(theStruct,0,13,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(msg->utc_day!=0) {
        mxSetFieldByNumber(theStruct,0,14,intMat2MatlabDoubles(&msg->utc_day,1,1));
    } else {//If UTC day is not available.
        mxSetFieldByNumber(theStruct,0,14,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(msg->utc_hour<24) {
        mxSetFieldByNumber(theStruct,0,15,intMat2MatlabDoubles(&msg->utc_hour,1,1));
    } else {//If UTC hour not available
        mxSetFieldByNumber(theStruct,0,15,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(msg->utc_min<60) {
        mxSetFieldByNumber(theStruct,0,16,intMat2MatlabDoubles(&msg->utc_min,1,1));
    } else {//If UTC minute not available
        mxSetFieldByNumber(theStruct,0,16,mxCreateDoubleMatrix(0,0,mxREAL));
    }

    if(msg->services_known) {
        mxSetFieldByNumber(theStruct,0,17,intMat2MatlabDoubles(msg->services.data(),26,1));
    } else {
        mxSetFieldByNumber(theStruct,0,17,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(msg->name.compare("@@@@@@@@@@@@@@@@@@@@")!=0) {
        const char *charString=msg->name.c_str();
        mxSetFieldByNumber(theStruct,0,18,mxCreateCharMatrixFromStrings(1,&charString));
    } else {
        mxSetFieldByNumber(theStruct,0,18,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    if(fabs(msg->position.lng_deg)<=180) {
        mxSetFieldByNumber(theStruct,0,19,doubleMat2Matlab(&msg->position.lng_deg,1,1));
    } else {//No longitude available
        mxSetFieldByNumber(theStruct,0,19,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(fabs(msg->position.lat_deg)<=90){
        mxSetFieldByNumber(theStruct,0,20,doubleMat2Matlab(&msg->position.lat_deg,1,1));
    } else {//No latitude available
        mxSetFieldByNumber(theStruct,0,20,mxCreateDoubleMatrix(0,0,mxREAL));
    }
}

void AIS_6_1_25_ToMatlab(libais::Ais6_1_25 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions){
//This message format is described in IMO circular 289
    const size_t numberOfFields=12;
    const mwSize dims[2] = {1, 1};
    //These are the names of all of the fields that will be added to the
    //Matlab structre array.
    const char *fieldNames[numberOfFields] = {"message_id", "repeat_indicator",
    "mmsi", "seq", "mmsi_dest", "retransmit", "spare", "dac", "fi",
    "amount_unit", "amount", "cargos"};
    const size_t numberOfCargoFields=7;
    mwSize cargoDims[2];
    const char *cargoFieldNames[numberOfCargoFields] = {"code_type",
    "imdg", "spare", "un", "bc", "marpol_oil", "marpol_cat"};
    mxArray *theStruct, *cargoArray;
    size_t curCargo;
         
    //If text descriptions of the fields are desired, then put them into a
    //structure.
    if(fieldDescriptions!=NULL) {
        const char *descriptionStrings[numberOfFields] = {
            "Identifier for Message 6; always 6.",//0
            "Repeat indicator. Used by the repeater to indicate how many times a message has been repeated. 0 = default; 3 = do not repeat any more.",//1
            "Maritime mobile service identity (MMSI) number of source station",//2
            "Sequence number (0-3); refer to Section 5.3.1 of Annex 2 of ITU-R M.1371-5",//3
            "MMSI number of destination station",//4
            "Retransmit flag should be set upon retransmission: 0 = no retransmission = default; 1 = retransmitted",//5
            "Not used. Should be zero. Reserved for future use",//6
            "Designated area code (DAC)",//7
            "Function identifier (FI)",//8
             "Units of quantity of main dangerous good\n1 = in kg\n2 = in tons (10e3 kg)\n3 = in 1000 tons (10e6 kg)",//9
            "Quantity of dangerous cargo (1-1023)",//10
            ""//11 Not used; a structure indicating each of the subelements is saved instead.
        };
        const char *cargoDescriptionStrings[numberOfCargoFields] ={
            "Code under which the cargo is carried. Values are:\n1 = IMDG Code (in packed form)\n2 = IGC Code\n3 = BC Code (from 1.1.2011 IMSBC)\n4 = MARPOL Annex I List of oils (Appendix 1) 5 = MARPOL Annex II IBC Code\n6 = Regional use\n7 - 15 (reserved for future use)",//0
            "IMDG code of the cargo:\n1 - 9 (not used)\n10 - 99 = first digit = main class, second digit = subclass or division (undefined subclasses and divisions should not be used)\n100 - 127 (reserved for future use)",//1
            "Not used. Set to zero.",//2
            "UN Number of cargo, 1 - 3363 = Four digits UN number 3364 - 8191 (reserved for future use)",//3
            "BC code of cargo:\n1=A\n2=B\n3=C\n4=MHB - Material Hazardous in Bulk\n5 - 7 (reserved for future use)",//4
            "Marpol Oil Type:\n1 = asphalt solutions\n2 = oils\n3 = distillates\n4 = gas oil\n5 = gasoline blending stocks\n6 = gasoline\n7 = jet fuels\n8 = naphtha\n9 - 15 (reserved for future use)",//5
            "Marpol category:\n1 = Category X\n2 = Category Y\n3 = Category Z\n4 = other substances\n5 - 7 (reserved for future use)"//6
        };
        unsigned int i;
        mxArray *cargoDescripStruct;

        *fieldDescriptions=mxCreateStructArray(2, dims, numberOfFields, fieldNames);
        for(i=0;i<(numberOfFields-1);i++) {
            mxSetFieldByNumber(*fieldDescriptions,0,static_cast<int>(i),mxCreateCharMatrixFromStrings(1,&descriptionStrings[i]));
        }
        cargoDescripStruct=mxCreateStructArray(2, dims, numberOfCargoFields, cargoFieldNames);
        for(i=0;i<(numberOfCargoFields-1);i++) {
            mxSetFieldByNumber(cargoDescripStruct,0,static_cast<int>(i),mxCreateCharMatrixFromStrings(1,&cargoDescriptionStrings[i]));
        }
        //Set the descriptions for the cargo fields
        mxSetFieldByNumber(*fieldDescriptions,0,(numberOfFields-1),cargoDescripStruct);
    }
    
    //Create the Matlab structure array.
    theStruct=mxCreateStructArray(2, dims, numberOfFields, fieldNames);
    *decodedMessage=theStruct;
    
    //Fill all of the elements of the structure
    //This first set is common to all AIS6 messages.
    mxSetFieldByNumber(theStruct,0,0,intMat2MatlabDoubles(&msg->message_id,1,1));
    mxSetFieldByNumber(theStruct,0,1,intMat2MatlabDoubles(&msg->repeat_indicator,1,1));
    mxSetFieldByNumber(theStruct,0,2,intMat2MatlabDoubles(&msg->mmsi,1,1));
    mxSetFieldByNumber(theStruct,0,3,intMat2MatlabDoubles(&msg->seq,1,1));
    mxSetFieldByNumber(theStruct,0,4,intMat2MatlabDoubles(&msg->mmsi_dest,1,1));
    mxSetFieldByNumber(theStruct,0,5,boolMat2Matlab(&msg->retransmit,1,1));
    mxSetFieldByNumber(theStruct,0,6,intMat2MatlabDoubles(&msg->spare,1,1));
    mxSetFieldByNumber(theStruct,0,7,intMat2MatlabDoubles(&msg->dac,1,1));
    mxSetFieldByNumber(theStruct,0,8,intMat2MatlabDoubles(&msg->fi,1,1));
    
    //Fields specific to the Ais6_1_25 type.
    if(msg->amount_unit!=0) {
        mxSetFieldByNumber(theStruct,0,9,intMat2MatlabDoubles(&msg->amount_unit,1,1));
    } else {//Dangerous cargo units not available
        mxSetFieldByNumber(theStruct,0,9,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    if(msg->amount!=0) {
        mxSetFieldByNumber(theStruct,0,10,intMat2MatlabDoubles(&msg->amount,1,1));
    } else {//Dangerous cargo amount not available
        mxSetFieldByNumber(theStruct,0,10,mxCreateDoubleMatrix(0,0,mxREAL));
    }

    cargoDims[0]=msg->cargos.size();
    cargoDims[1]=1;
    cargoArray=mxCreateStructArray(2, dims, numberOfCargoFields, cargoFieldNames);
    mxSetFieldByNumber(theStruct,0,11,cargoArray);
    
    for (curCargo= 0; curCargo < msg->cargos.size(); curCargo++) {
        if(msg->cargos[curCargo].code_type!=0) {
            mxSetFieldByNumber(cargoArray,curCargo,0,intMat2MatlabDoubles(&msg->cargos[curCargo].code_type,1,1));
        } else {//Code type not available
            mxSetFieldByNumber(cargoArray,curCargo,0,mxCreateDoubleMatrix(0,0,mxREAL));
        }

        if(msg->cargos[curCargo].imdg_valid&&msg->cargos[curCargo].imdg!=0) {
            mxSetFieldByNumber(cargoArray,curCargo,1,intMat2MatlabDoubles(&msg->cargos[curCargo].imdg,1,1));
        } else {
            mxSetFieldByNumber(cargoArray,curCargo,1,mxCreateDoubleMatrix(0,0,mxREAL));
        }
        
        if(msg->cargos[curCargo].spare_valid) {
            mxSetFieldByNumber(cargoArray,curCargo,2,intMat2MatlabDoubles(&msg->cargos[curCargo].spare,1,1));
        } else {
            mxSetFieldByNumber(cargoArray,curCargo,2,mxCreateDoubleMatrix(0,0,mxREAL));
        }
        
        if(msg->cargos[curCargo].un_valid&&msg->cargos[curCargo].un!=0) {
            mxSetFieldByNumber(cargoArray,curCargo,3,intMat2MatlabDoubles(&msg->cargos[curCargo].un,1,1));
        } else {
            mxSetFieldByNumber(cargoArray,curCargo,3,mxCreateDoubleMatrix(0,0,mxREAL));
        }
        
        if(msg->cargos[curCargo].bc_valid&&msg->cargos[curCargo].bc!=0) {
            mxSetFieldByNumber(cargoArray,curCargo,4,intMat2MatlabDoubles(&msg->cargos[curCargo].bc,1,1));
        } else {
            mxSetFieldByNumber(cargoArray,curCargo,4,mxCreateDoubleMatrix(0,0,mxREAL));
        }
        
        if(msg->cargos[curCargo].marpol_oil_valid&&msg->cargos[curCargo].marpol_oil!=0) {
            mxSetFieldByNumber(cargoArray,curCargo,5,intMat2MatlabDoubles(&msg->cargos[curCargo].marpol_oil,1,1));
        } else {
            mxSetFieldByNumber(cargoArray,curCargo,5,mxCreateDoubleMatrix(0,0,mxREAL));
        }
        
        if(msg->cargos[curCargo].marpol_cat_valid&&msg->cargos[curCargo].marpol_cat!=0) {
            mxSetFieldByNumber(cargoArray,curCargo,6,intMat2MatlabDoubles(&msg->cargos[curCargo].marpol_cat,1,1));
        } else {
            mxSetFieldByNumber(cargoArray,curCargo,6,mxCreateDoubleMatrix(0,0,mxREAL));
        }
    }
}

void AIS_6_1_32_ToMatlab(libais::Ais6_1_32 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions){
//This message format is described in IMO circular 289
    const size_t numberOfFields=12;
    const mwSize dims[2] = {1, 1};
    //These are the names of all of the fields that will be added to the
    //Matlab structre array.
    const char *fieldNames[numberOfFields] = {"message_id", "repeat_indicator",
    "mmsi", "seq", "mmsi_dest", "retransmit", "spare", "dac", "fi",
    "utc_month", "utc_day", "windows"};
    const size_t numberOfWindowFields=8;
    mwSize windowDims[2];
    const char *windowFieldNames[numberOfWindowFields] = {"lon", "lat", "from_utc_hour",
    "from_utc_min", "to_utc_hour", "to_utc_min", "cur_dir", "cur_speed"};
    mxArray *theStruct, *windowArray;
    size_t curWindow;
              
    //If text descriptions of the fields are desired, then put them into a
    //structure.
    if(fieldDescriptions!=NULL) {
        const char *descriptionStrings[numberOfFields] = {
            "Identifier for Message 6; always 6.",//0
            "Repeat indicator. Used by the repeater to indicate how many times a message has been repeated. 0 = default; 3 = do not repeat any more.",//1
            "Maritime mobile service identity (MMSI) number of source station",//2
            "Sequence number (0-3); refer to Section 5.3.1 of Annex 2 of ITU-R M.1371-5",
            "MMSI number of destination station",//4
            "Retransmit flag should be set upon retransmission: 0 = no retransmission = default; 1 = retransmitted",//5
            "Not used. Should be zero. Reserved for future use",//6
            "Designated area code (DAC)",//7
            "Function identifier (FI)",//8
            "UTC Month (1-12)",//9
            "UTC Day (1-31)",//10
            "An array of tidal windows at different times consisting of components:\nLongitude East in degrees (WGS-84)\nLatitude North in degrees (WGS-84)\nUTC hours from (0-23)\nUTC minutes from (0-59)\nUTC hours to (0-23)\nUTC minutes to (0-59)\nDirection of the current in degrees East of North\nSpeed of the current in knots, 25.1 means 25.1 knots or greater"//11
        };
        unsigned int i;

        *fieldDescriptions=mxCreateStructArray(2, dims, numberOfFields, fieldNames);

        for(i=0;i<numberOfFields;i++) {
            mxSetFieldByNumber(*fieldDescriptions,0,static_cast<int>(i),mxCreateCharMatrixFromStrings(1,&descriptionStrings[i]));
        }
    }
    
    //Create the Matlab structure array.
    theStruct=mxCreateStructArray(2, dims, numberOfFields, fieldNames);
    *decodedMessage=theStruct;
    
    //Fill all of the elements of the structure
    //This first set is common to all AIS6 messages.
    mxSetFieldByNumber(theStruct,0,0,intMat2MatlabDoubles(&msg->message_id,1,1));
    mxSetFieldByNumber(theStruct,0,1,intMat2MatlabDoubles(&msg->repeat_indicator,1,1));
    mxSetFieldByNumber(theStruct,0,2,intMat2MatlabDoubles(&msg->mmsi,1,1));
    mxSetFieldByNumber(theStruct,0,3,intMat2MatlabDoubles(&msg->seq,1,1));
    mxSetFieldByNumber(theStruct,0,4,intMat2MatlabDoubles(&msg->mmsi_dest,1,1));
    mxSetFieldByNumber(theStruct,0,5,boolMat2Matlab(&msg->retransmit,1,1));
    mxSetFieldByNumber(theStruct,0,6,intMat2MatlabDoubles(&msg->spare,1,1));
    mxSetFieldByNumber(theStruct,0,7,intMat2MatlabDoubles(&msg->dac,1,1));
    mxSetFieldByNumber(theStruct,0,8,intMat2MatlabDoubles(&msg->fi,1,1));
    
    //Fields specific to the Ais6_1_32 type.
    if(msg->utc_month!=0) {
        mxSetFieldByNumber(theStruct,0,9,intMat2MatlabDoubles(&msg->utc_month,1,1));
    } else {//No UTC month given
        mxSetFieldByNumber(theStruct,0,9,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    if(msg->utc_day!=0) {
        mxSetFieldByNumber(theStruct,0,10,intMat2MatlabDoubles(&msg->utc_day,1,1));
    } else {//No UTC day given
        mxSetFieldByNumber(theStruct,0,10,mxCreateDoubleMatrix(0,0,mxREAL));
    }

    windowDims[0]=msg->windows.size();
    windowDims[1]=1;
    windowArray=mxCreateStructArray(2, dims, numberOfWindowFields, windowFieldNames);
    mxSetFieldByNumber(theStruct,0,11,windowArray);

    for (curWindow = 0; curWindow < msg->windows.size(); curWindow++) {
        if(fabs(msg->windows[curWindow].position.lng_deg)<=180) {
            mxSetFieldByNumber(windowArray,curWindow,0,doubleMat2Matlab(&msg->windows[curWindow].position.lng_deg,1,1));
        } else {//No longitude available
            mxSetFieldByNumber(windowArray,curWindow,0,mxCreateDoubleMatrix(0,0,mxREAL));
        }
        
        if(fabs(msg->windows[curWindow].position.lat_deg)<=90) {
            mxSetFieldByNumber(windowArray,curWindow,1,doubleMat2Matlab(&msg->windows[curWindow].position.lat_deg,1,1));
        } else {//No latitude available
            mxSetFieldByNumber(windowArray,curWindow,1,mxCreateDoubleMatrix(0,0,mxREAL));
        }
        
        if(msg->windows[curWindow].from_utc_hour<24) {
            mxSetFieldByNumber(windowArray,curWindow,2,intMat2MatlabDoubles(&msg->windows[curWindow].from_utc_hour,1,1));
        } else {//No UTC hour from available
            mxSetFieldByNumber(windowArray,curWindow,2,mxCreateDoubleMatrix(0,0,mxREAL));
        }
        
        if(msg->windows[curWindow].from_utc_min<60) {
            mxSetFieldByNumber(windowArray,curWindow,3,intMat2MatlabDoubles(&msg->windows[curWindow].from_utc_min,1,1));
        } else {//No UTC minutes from available
            mxSetFieldByNumber(windowArray,curWindow,3,mxCreateDoubleMatrix(0,0,mxREAL));
        }
        
        if(msg->windows[curWindow].to_utc_hour<24) {
            mxSetFieldByNumber(windowArray,curWindow,4,intMat2MatlabDoubles(&msg->windows[curWindow].to_utc_hour,1,1));
        } else {//No UTC hour to available
            mxSetFieldByNumber(windowArray,curWindow,4,mxCreateDoubleMatrix(0,0,mxREAL));
        }
        
        if(msg->windows[curWindow].to_utc_min<60) {
            mxSetFieldByNumber(windowArray,curWindow,5,intMat2MatlabDoubles(&msg->windows[curWindow].to_utc_min,1,1));
        } else {//No UTC minutes to avaiable
            mxSetFieldByNumber(windowArray,curWindow,5,mxCreateDoubleMatrix(0,0,mxREAL));
        }
        
        if(msg->windows[curWindow].cur_dir<360) {
            mxSetFieldByNumber(windowArray,curWindow,6,intMat2MatlabDoubles(&msg->windows[curWindow].cur_dir,1,1));
        } else {//No direction of the current avaiable
            mxSetFieldByNumber(windowArray,curWindow,6,mxCreateDoubleMatrix(0,0,mxREAL));
        }
        
        if(msg->windows[curWindow].cur_speed!=255) {
            mxSetFieldByNumber(windowArray,curWindow,7,floatMat2MatlabDoubles(&msg->windows[curWindow].cur_speed,1,1));
        } else {//No current speed available
            mxSetFieldByNumber(windowArray,curWindow,7,mxCreateDoubleMatrix(0,0,mxREAL));
        }
    }
}

void AIS_6_1_40_ToMatlab(libais::Ais6_1_40 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions){
//This message format is described in Annex 5 of ITU-R M.1371-1. It is not
//present in ITU-R M.1371-5.
    const size_t numberOfFields=11;
    const mwSize dims[2] = {1, 1};
    //These are the names of all of the fields that will be added to the
    //Matlab structre array.
    const char *fieldNames[numberOfFields] = {"message_id", "repeat_indicator",
    "mmsi", "seq", "mmsi_dest", "retransmit", "spare", "dac", "fi",
    "persons", "spare2"};
    mxArray *theStruct;

    //If text descriptions of the fields are desired, then put them into a
    //structure.
    if(fieldDescriptions!=NULL) {
        const char *descriptionStrings[numberOfFields] = {
            "Identifier for Message 6; always 6.",//0
            "Repeat indicator. Used by the repeater to indicate how many times a message has been repeated. 0 = default; 3 = do not repeat any more.",//1
            "Maritime mobile service identity (MMSI) number of source station",//2
            "Sequence number (0-3); refer to Section 5.3.1 of Annex 2 of ITU-R M.1371-5",//3
            "MMSI number of destination station",//4
            "Retransmit flag should be set upon retransmission: 0 = no retransmission = default; 1 = retransmitted",//5
            "Not used. Should be zero. Reserved for future use",//6
            "Designated area code (DAC)",//7
            "Function identifier (FI)",//8
            "Current number of persons on board, including crew members",//9
            "Not used. Should be set to zero"//10
        };
        unsigned int i;

        *fieldDescriptions=mxCreateStructArray(2, dims, numberOfFields, fieldNames);
 
        for(i=0;i<numberOfFields;i++) {
            mxSetFieldByNumber(*fieldDescriptions,0,static_cast<int>(i),mxCreateCharMatrixFromStrings(1,&descriptionStrings[i]));
        }
    }
    
    //Create the Matlab structure array.
    theStruct=mxCreateStructArray(2, dims, numberOfFields, fieldNames);
    *decodedMessage=theStruct;
    
    //Fill all of the elements of the structure
    //This first set is common to all AIS6 messages.
    mxSetFieldByNumber(theStruct,0,0,intMat2MatlabDoubles(&msg->message_id,1,1));
    mxSetFieldByNumber(theStruct,0,1,intMat2MatlabDoubles(&msg->repeat_indicator,1,1));
    mxSetFieldByNumber(theStruct,0,2,intMat2MatlabDoubles(&msg->mmsi,1,1));
    mxSetFieldByNumber(theStruct,0,3,intMat2MatlabDoubles(&msg->seq,1,1));
    mxSetFieldByNumber(theStruct,0,4,intMat2MatlabDoubles(&msg->mmsi_dest,1,1));
    mxSetFieldByNumber(theStruct,0,5,boolMat2Matlab(&msg->retransmit,1,1));
    mxSetFieldByNumber(theStruct,0,6,intMat2MatlabDoubles(&msg->spare,1,1));
    mxSetFieldByNumber(theStruct,0,7,intMat2MatlabDoubles(&msg->dac,1,1));
    mxSetFieldByNumber(theStruct,0,8,intMat2MatlabDoubles(&msg->fi,1,1));
    
    //Fields specific to the Ais6_1_40 type.
    if(msg->persons!=0) {
        mxSetFieldByNumber(theStruct,0,9,intMat2MatlabDoubles(&msg->persons,1,1));
    } else {//The number of persons is not available
        mxSetFieldByNumber(theStruct,0,9,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    mxSetFieldByNumber(theStruct,0,10,intMat2MatlabDoubles(&msg->spare2,1,1));
}

////////////////////////////////////////
///All of the functions for Message 8///
////////////////////////////////////////
void AIS_8_1_0_ToMatlab(libais::Ais8_1_0 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions) {
//This message is described in Annex 5 of ITU-R M.1371-5.
    const size_t numberOfFields=10;
    const mwSize dims[2] = {1, 1};
    //These are the names of all of the fields that will be added to the
    //Matlab structre array.
    const char *fieldNames[numberOfFields] = {"message_id", "repeat_indicator",
    "mmsi", "spare", "dac", "fi", "ack_required", "msg_seq", "text",
    "spare2"};
    mxArray *theStruct;

    //If text descriptions of the fields are desired, then put them into a
    //structure.
    if(fieldDescriptions!=NULL) {
        const char *descriptionStrings[numberOfFields] = {
            "Identifier for Message 8; always 8.",//0
            "Repeat indicator. Used by the repeater to indicate how many times a message has been repeated. 0 = default; 3 = do not repeat any more.",//1
            "Maritime mobile service identity (MMSI) number of source station",//2
            "Not used. Should be zero. Reserved for future use",//3
            "Designated area code (DAC)",//4
            "Function identifier (FI)",//5
            "Acknowledge required flag\n1 = reply is required, optional for addressed binary messages and not used for binary broadcast messages\n0 = reply is not required, optional for an addressed binary message and required for binary broadcast messages",//6
            "Text sequence number to be incremented by the application. Zero indicates that sequence numbers are not being used",//7
            "Text String",//8
            "Not used for data and should be set to zero. The number of bits should be either 0, 2, 4, or 6 to maintain byte boundaries."//9
        };
        unsigned int i;

        *fieldDescriptions=mxCreateStructArray(2, dims, numberOfFields, fieldNames);

        for(i=0;i<numberOfFields;i++) {
            mxSetFieldByNumber(*fieldDescriptions,0,static_cast<int>(i),mxCreateCharMatrixFromStrings(1,&descriptionStrings[i]));
        }
    }
    
    //Create the Matlab structure array.
    theStruct=mxCreateStructArray(2, dims, numberOfFields, fieldNames);
    *decodedMessage=theStruct;
    
    //Fill all of the elements of the structure with the set that is common
    //to all AIS8 messages.
    mxSetFieldByNumber(theStruct,0,0,intMat2MatlabDoubles(&msg->message_id,1,1));
    mxSetFieldByNumber(theStruct,0,1,intMat2MatlabDoubles(&msg->repeat_indicator,1,1));
    mxSetFieldByNumber(theStruct,0,2,intMat2MatlabDoubles(&msg->mmsi,1,1));
    mxSetFieldByNumber(theStruct,0,3,intMat2MatlabDoubles(&msg->spare,1,1));
    mxSetFieldByNumber(theStruct,0,4,intMat2MatlabDoubles(&msg->dac,1,1));
    mxSetFieldByNumber(theStruct,0,5,intMat2MatlabDoubles(&msg->fi,1,1));

    //Fields specific to the Ais8_1_0 type.
    mxSetFieldByNumber(theStruct,0,6,boolMat2Matlab(&msg->ack_required,1,1));
    mxSetFieldByNumber(theStruct,0,7,intMat2MatlabDoubles(&msg->msg_seq,1,1));
    {
        const char *charString=msg->text.c_str();
        mxSetFieldByNumber(theStruct,0,8,mxCreateCharMatrixFromStrings(1,&charString));
    }
    mxSetFieldByNumber(theStruct,0,9,intMat2MatlabDoubles(&msg->spare2,1,1));
}

void AIS_8_1_11_ToMatlab(libais::Ais8_1_11 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions) {
//This message format is described in IMO circular 236, which was revoked
//in IMO circular 289 in 2013. The note says that if information is not
//available, then a field will be set to the max value. Thus, validity
//below is determined by seeing if data is withing the valid range.
//(Message "METEOROLOGICAL AND HYDROLOGICAL DATA")
    const size_t numberOfFields=43;
    const mwSize dims[2] = {1, 1};
    //These are the names of all of the fields that will be added to the
    //Matlab structre array.
    const char *fieldNames[numberOfFields] = {"message_id", "repeat_indicator",
    "mmsi", "spare", "dac", "fi", "lon", "lat", "day", "hour", "minute",
    "wind_ave", "wind_gust", "wind_dir", "wind_gust_dir", "air_temp",
    "rel_humid", "dew_point", "air_pres", "air_pres_trend", "horz_vis",
    "water_level", "water_level_trend", "surf_cur_speed", "surf_cur_dir",
    "cur_speed_2", "cur_dir_2", "cur_depth_2", "cur_speed_3", "cur_dir_3",
    "cur_depth_3", "wave_height", "wave_period", "wave_dir", 
    "swell_height", "swell_period", "swell_dir", "sea_state", "water_temp",
    "precip_type", "salinity", "ice", "spare2"};
    mxArray *theStruct;
    
    //If text descriptions of the fields are desired, then put them into a
    //structure.
    if(fieldDescriptions!=NULL) {
        const char *descriptionStrings[numberOfFields] = {
            "Identifier for Message 8; always 8.",//0
            "Repeat indicator. Used by the repeater to indicate how many times a message has been repeated. 0 = default; 3 = do not repeat any more.",//1
            "Maritime mobile service identity (MMSI) number of source station",//2
            "Not used. Should be zero. Reserved for future use",//3
            "Designated area code (DAC)",//4
            "Function identifier (FI)",//5
            "Longitude East in degrees (WGS-84)",//6
            "Latitude North in degrees (WGS-84)",//7
            "UTC day transmitted (1-31)",//8
            "UTC hour transmitted (0-23)",//9
            "UTC minute transmitted (0-59)",//10
            "Average wind speed for the last 10 minutes in knots (0-120)",//11
            "Maximum wind gust speed during the last 10 minutes in knots (0-120)",//12
            "Wind direction in degrees East of North",//13
            "Wind gust direction in degrees East of North",//14
            "Dry buld air temperature in degrees Celsius (-60 to 60)",//15
            "Relative humidity in percent (0-100)",//16
            "Dew point in degrees Celsius (-20 to 50)",//17
            "Air pressure in hectoPascals (800-1200)",//18
            "Air pressure tendency 0 = steady, 1 = decreasing, 2 = increasing",//19
            "Horizontal visibility in nautical miles (0-25)",//20
            "Water level, including tide, deviation from local datum in meters (-10 to 30)",//21
            "Water level trend 0 = steady, 1 = decreasing, 2 = increasing",//22
            "Surface current speed including tide in knots (0-25)",//23
            "Surface current direction in degrees East of North (0-359)",//24
            "Current 2 speed in knots (0-25)",//25
            "Current 2 speed direction in degrees East of North (0-359)",//26
            "Depth of current 2 in meters (0-30)",//27
            "Current 3 speed in knots (0-25)",//28
            "Current 3 speed direction in degrees East of North (0-359)",//29
            "Depth of current 3 in meters (0-30)",//30
            "Significant wave height in meters (0-25)",//31
            "Wave period in seconds (0-60)",//32
            "Wave direction in degrees East of North (0-359)",//33
            "Swell height in meters (0-25)",//34
            "Swell period in seconds (0-60)",//35
            "Swell direction in degrees East of North (0-359)",//36
            "Sea state according to the Beauford scale (0-12)",//37
            "Water temperature in degrees Celsius (-10 to 50)",//38
            "Precipitation type according to the WMO\n0 = reserved\n1 = rain\n2 = thunderstorm\n3 = freezing rain\n4 = mixed/ice\n5 = snow\n6 = reserved",//39
            "Salinity in parts per thousand (0-50)",//40
            "A boolean value indicating the presence of ice",//41
            "Not used. Should be zero."//42
        };
        unsigned int i;

        *fieldDescriptions=mxCreateStructArray(2, dims, numberOfFields, fieldNames);

        for(i=0;i<numberOfFields;i++) {
            mxSetFieldByNumber(*fieldDescriptions,0,static_cast<int>(i),mxCreateCharMatrixFromStrings(1,&descriptionStrings[i]));
        }
    }
    //Create the Matlab structure array.
    theStruct=mxCreateStructArray(2, dims, numberOfFields, fieldNames);
    *decodedMessage=theStruct;
    
    //Fill all of the elements of the structure with the set that is common
    //to all AIS8 messages.
    mxSetFieldByNumber(theStruct,0,0,intMat2MatlabDoubles(&msg->message_id,1,1));
    mxSetFieldByNumber(theStruct,0,1,intMat2MatlabDoubles(&msg->repeat_indicator,1,1));
    mxSetFieldByNumber(theStruct,0,2,intMat2MatlabDoubles(&msg->mmsi,1,1));
    mxSetFieldByNumber(theStruct,0,3,intMat2MatlabDoubles(&msg->spare,1,1));
    mxSetFieldByNumber(theStruct,0,4,intMat2MatlabDoubles(&msg->dac,1,1));
    mxSetFieldByNumber(theStruct,0,5,intMat2MatlabDoubles(&msg->fi,1,1));

    //Fields specific to the Ais8_1_11 type.
    if(fabs(msg->position.lng_deg)<=180) {
        mxSetFieldByNumber(theStruct,0,6,doubleMat2Matlab(&msg->position.lng_deg,1,1));
    } else {//No longitude
        mxSetFieldByNumber(theStruct,0,6,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    if(fabs(msg->position.lat_deg)<=90) {
        mxSetFieldByNumber(theStruct,0,7,doubleMat2Matlab(&msg->position.lat_deg,1,1));
    } else {//No latitude
        mxSetFieldByNumber(theStruct,0,7,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    if(msg->day!=0&&msg->day<=31) {
        mxSetFieldByNumber(theStruct,0,8,intMat2MatlabDoubles(&msg->day,1,1));
    } else {//No day
        mxSetFieldByNumber(theStruct,0,8,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    if(msg->hour<24) {
        mxSetFieldByNumber(theStruct,0,9,intMat2MatlabDoubles(&msg->hour,1,1));
    } else {//No hour
        mxSetFieldByNumber(theStruct,0,9,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    if(msg->minute<60) {
        mxSetFieldByNumber(theStruct,0,10,intMat2MatlabDoubles(&msg->minute,1,1));
    } else {//No minute
        mxSetFieldByNumber(theStruct,0,10,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    if(msg->wind_ave<=120) {
        mxSetFieldByNumber(theStruct,0,11,intMat2MatlabDoubles(&msg->wind_ave,1,1));
    } else {//No wind average
        mxSetFieldByNumber(theStruct,0,11,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    if(msg->wind_gust<=120) {
        mxSetFieldByNumber(theStruct,0,12,intMat2MatlabDoubles(&msg->wind_gust,1,1));
    } else {//No wind gust speed
        mxSetFieldByNumber(theStruct,0,12,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    if(msg->wind_dir<=359) {
        mxSetFieldByNumber(theStruct,0,13,intMat2MatlabDoubles(&msg->wind_dir,1,1));
    } else {//No wind direction
        mxSetFieldByNumber(theStruct,0,13,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    if(msg->wind_gust_dir<=359) {
        mxSetFieldByNumber(theStruct,0,14,intMat2MatlabDoubles(&msg->wind_gust_dir,1,1));
    } else {//No wind gust direction
        mxSetFieldByNumber(theStruct,0,14,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    if(msg->air_temp<=60&&msg->air_temp>=-60) {
        mxSetFieldByNumber(theStruct,0,15,floatMat2MatlabDoubles(&msg->air_temp,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,15,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    if(msg->rel_humid<=100) {
        mxSetFieldByNumber(theStruct,0,16,intMat2MatlabDoubles(&msg->rel_humid,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,16,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    if(msg->dew_point<=50&&msg->dew_point>=-20) {
        mxSetFieldByNumber(theStruct,0,17,floatMat2MatlabDoubles(&msg->dew_point,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,17,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    if(msg->air_pres<=1200) {
        mxSetFieldByNumber(theStruct,0,18,floatMat2MatlabDoubles(&msg->air_pres,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,18,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    if(msg->air_pres_trend<=2){
        mxSetFieldByNumber(theStruct,0,19,intMat2MatlabDoubles(&msg->air_pres_trend,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,19,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    if(msg->horz_vis<=25) {
        mxSetFieldByNumber(theStruct,0,20,floatMat2MatlabDoubles(&msg->horz_vis,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,20,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    if(msg->water_level<=30&&msg->water_level>=-10) {
        mxSetFieldByNumber(theStruct,0,21,floatMat2MatlabDoubles(&msg->water_level,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,21,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    if(msg->water_level_trend<=2) {
        mxSetFieldByNumber(theStruct,0,22,intMat2MatlabDoubles(&msg->water_level_trend,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,22,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    if(msg->surf_cur_speed<=25) {
        mxSetFieldByNumber(theStruct,0,23,floatMat2MatlabDoubles(&msg->surf_cur_speed,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,23,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    if(msg->surf_cur_dir<=359) {
        mxSetFieldByNumber(theStruct,0,24,intMat2MatlabDoubles(&msg->surf_cur_dir,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,24,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    if(msg->cur_speed_2<=25) {
        mxSetFieldByNumber(theStruct,0,25,floatMat2MatlabDoubles(&msg->cur_speed_2,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,25,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    if(msg->cur_dir_2<=359) {
        mxSetFieldByNumber(theStruct,0,26,intMat2MatlabDoubles(&msg->cur_dir_2,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,26,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    if(msg->cur_depth_2<=30) {
        mxSetFieldByNumber(theStruct,0,27,intMat2MatlabDoubles(&msg->cur_depth_2,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,27,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    if(msg->cur_speed_3<=25) {
        mxSetFieldByNumber(theStruct,0,28,floatMat2MatlabDoubles(&msg->cur_speed_3,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,28,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    if(msg->cur_dir_3<=359) {
        mxSetFieldByNumber(theStruct,0,29,intMat2MatlabDoubles(&msg->cur_dir_3,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,29,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    if(msg->cur_depth_3<=30) {
        mxSetFieldByNumber(theStruct,0,30,intMat2MatlabDoubles(&msg->cur_depth_3,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,30,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    if(msg->wave_height<=25) {
        mxSetFieldByNumber(theStruct,0,31,floatMat2MatlabDoubles(&msg->wave_height,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,31,mxCreateDoubleMatrix(0,0,mxREAL)); 
    }
    if(msg->wave_period<=60) {
        mxSetFieldByNumber(theStruct,0,32,intMat2MatlabDoubles(&msg->wave_period,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,32,mxCreateDoubleMatrix(0,0,mxREAL)); 
    }
    if(msg->wave_dir<=359) {
        mxSetFieldByNumber(theStruct,0,33,intMat2MatlabDoubles(&msg->wave_dir,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,33,mxCreateDoubleMatrix(0,0,mxREAL)); 
    }
    if(msg->swell_height<=25) {
        mxSetFieldByNumber(theStruct,0,34,floatMat2MatlabDoubles(&msg->swell_height,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,34,mxCreateDoubleMatrix(0,0,mxREAL)); 
    }
    if(msg->swell_period<=60) {
        mxSetFieldByNumber(theStruct,0,35,intMat2MatlabDoubles(&msg->swell_period,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,35,mxCreateDoubleMatrix(0,0,mxREAL)); 
    }
    if(msg->swell_dir<=359) {
        mxSetFieldByNumber(theStruct,0,36,intMat2MatlabDoubles(&msg->swell_dir,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,36,mxCreateDoubleMatrix(0,0,mxREAL)); 
    }
    if(msg->sea_state<=12) {
        mxSetFieldByNumber(theStruct,0,37,intMat2MatlabDoubles(&msg->sea_state,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,37,mxCreateDoubleMatrix(0,0,mxREAL)); 
    }
    if(msg->water_temp<=50&&msg->water_temp>=-10) {
        mxSetFieldByNumber(theStruct,0,38,floatMat2MatlabDoubles(&msg->water_temp,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,38,mxCreateDoubleMatrix(0,0,mxREAL)); 
    }
    if(msg->precip_type<7) {
        mxSetFieldByNumber(theStruct,0,39,intMat2MatlabDoubles(&msg->precip_type,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,39,mxCreateDoubleMatrix(0,0,mxREAL)); 
    }
    if(msg->salinity<511) {
        mxSetFieldByNumber(theStruct,0,40,floatMat2MatlabDoubles(&msg->salinity,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,40,mxCreateDoubleMatrix(0,0,mxREAL)); 
    }
    if(msg->ice==0||msg->ice==1) {
        mxSetFieldByNumber(theStruct,0,41,intMat2MatlabDoubles(&msg->ice,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,41,mxCreateDoubleMatrix(0,0,mxREAL)); 
    }

    mxSetFieldByNumber(theStruct,0,42,intMat2MatlabDoubles(&msg->spare2,1,1));
}

void AIS_8_1_13_ToMatlab(libais::Ais8_1_13 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions) {
//This message format is described in IMO circular 236, which was revoked
//in IMO circular 289 in 2013. Message "FAIRWAY CLOSED".
    const size_t numberOfFields=20;
    const mwSize dims[2] = {1, 1};
    //These are the names of all of the fields that will be added to the
    //Matlab structre array.
    const char *fieldNames[numberOfFields] = {"message_id", "repeat_indicator",
    "mmsi", "spare", "dac", "fi", "reason", "location_from", "location_to",
    "radius", "units", "day_from", "month_from", "hour_from",
    "minute_from", "day_to", "month_to", "hour_to", "minute_to", "spare2"};
    mxArray *theStruct;
      
    //If text descriptions of the fields are desired, then put them into a
    //structure.
    if(fieldDescriptions!=NULL) {
        const char *descriptionStrings[numberOfFields] = {
            "Identifier for Message 8; always 8.",//0
            "Repeat indicator. Used by the repeater to indicate how many times a message has been repeated. 0 = default; 3 = do not repeat any more.",//1
            "Maritime mobile service identity (MMSI) number of source station",//2
            "Not used. Should be zero. Reserved for future use",//3
            "Designated area code (DAC)",//4
            "Function identifier (FI)",//5
            "Reason for closing",//6
            "Location of closing from",//7
            "Location of closing to",//8
            "Extension of closed area (radius, 0-1000)",//9
            "Unit of extension value, 0=meters, 1=kilometers, 2=nautical miles, 3=cable length (presumably meaning 0.1 nautical miles)",//10
            "Closing from day (1-31)",//11
            "Closing from month (1-12)",//12
            "From LT hour (0-23) (approx)",//13
            "From LT minute (0-59) (approx)",//14
            "To day (1-31)",//15
            "To month (1-12)",//16
            "To LT hour (0-23) (approx)",//17
            "To LT minute (0-59) (approx)",//18
            "Not used. Should be zero"//19
        };
        unsigned int i;

        *fieldDescriptions=mxCreateStructArray(2, dims, numberOfFields, fieldNames);

        for(i=0;i<numberOfFields;i++) {
            mxSetFieldByNumber(*fieldDescriptions,0,static_cast<int>(i),mxCreateCharMatrixFromStrings(1,&descriptionStrings[i]));
        }
    }

    //Create the Matlab structure array.
    theStruct=mxCreateStructArray(2, dims, numberOfFields, fieldNames);
    *decodedMessage=theStruct;
    
    //Fill all of the elements of the structure with the set that is common
    //to all AIS8 messages.
    mxSetFieldByNumber(theStruct,0,0,intMat2MatlabDoubles(&msg->message_id,1,1));
    mxSetFieldByNumber(theStruct,0,1,intMat2MatlabDoubles(&msg->repeat_indicator,1,1));
    mxSetFieldByNumber(theStruct,0,2,intMat2MatlabDoubles(&msg->mmsi,1,1));
    mxSetFieldByNumber(theStruct,0,3,intMat2MatlabDoubles(&msg->spare,1,1));
    mxSetFieldByNumber(theStruct,0,4,intMat2MatlabDoubles(&msg->dac,1,1));
    mxSetFieldByNumber(theStruct,0,5,intMat2MatlabDoubles(&msg->fi,1,1));

    //Fields specific to the Ais8_1_13 type.
    if(msg->reason.compare("@@@@@@@@@@@@@@@@@@@@")!=0) {
        const char *charString=msg->reason.c_str();
        mxSetFieldByNumber(theStruct,0,6,mxCreateCharMatrixFromStrings(1,&charString));
    } else {//The reason is unavailable
        mxSetFieldByNumber(theStruct,0,6,mxCreateDoubleMatrix(0,0,mxREAL)); 
    }
    if(msg->location_from.compare("@@@@@@@@@@@@@@@@@@@@")!=0) {
        const char *charString=msg->location_from.c_str();
        mxSetFieldByNumber(theStruct,0,7,mxCreateCharMatrixFromStrings(1,&charString));
    } else {//The location of closing from is unavailable
        mxSetFieldByNumber(theStruct,0,7,mxCreateDoubleMatrix(0,0,mxREAL)); 
    }
    if(msg->location_to.compare("@@@@@@@@@@@@@@@@@@@@")!=0) {
        const char *charString=msg->location_to.c_str();
        mxSetFieldByNumber(theStruct,0,8,mxCreateCharMatrixFromStrings(1,&charString));
    } else {//The location of closing to is unavailable
        mxSetFieldByNumber(theStruct,0,8,mxCreateDoubleMatrix(0,0,mxREAL)); 
    }
    if(msg->radius<=1000) {
        mxSetFieldByNumber(theStruct,0,9,intMat2MatlabDoubles(&msg->radius,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,9,mxCreateDoubleMatrix(0,0,mxREAL)); 
    }
    mxSetFieldByNumber(theStruct,0,10,intMat2MatlabDoubles(&msg->units,1,1));
    if(msg->day_from!=0) {
        mxSetFieldByNumber(theStruct,0,11,intMat2MatlabDoubles(&msg->day_from,1,1));
    } else {//Day unavialable
        mxSetFieldByNumber(theStruct,0,11,mxCreateDoubleMatrix(0,0,mxREAL)); 
    }
    if(msg->month_from!=0) {
        mxSetFieldByNumber(theStruct,0,12,intMat2MatlabDoubles(&msg->month_from,1,1));
    } else {//Month unavailable
        mxSetFieldByNumber(theStruct,0,12,mxCreateDoubleMatrix(0,0,mxREAL)); 
    }
    if(msg->hour_from<=23) {
        mxSetFieldByNumber(theStruct,0,13,intMat2MatlabDoubles(&msg->hour_from,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,13,mxCreateDoubleMatrix(0,0,mxREAL)); 
    }
    if(msg->minute_from<=59) {
        mxSetFieldByNumber(theStruct,0,14,intMat2MatlabDoubles(&msg->minute_from,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,14,mxCreateDoubleMatrix(0,0,mxREAL)); 
    }
    if(msg->day_to!=0) {
        mxSetFieldByNumber(theStruct,0,15,intMat2MatlabDoubles(&msg->day_to,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,15,mxCreateDoubleMatrix(0,0,mxREAL)); 
    }
    if(msg->month_to!=0) {
        mxSetFieldByNumber(theStruct,0,16,intMat2MatlabDoubles(&msg->month_to,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,16,mxCreateDoubleMatrix(0,0,mxREAL)); 
    }
    if(msg->hour_to<=23) {
        mxSetFieldByNumber(theStruct,0,17,intMat2MatlabDoubles(&msg->hour_to,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,17,mxCreateDoubleMatrix(0,0,mxREAL)); 
    }
    if(msg->minute_to<=59) {
        mxSetFieldByNumber(theStruct,0,18,intMat2MatlabDoubles(&msg->minute_to,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,18,mxCreateDoubleMatrix(0,0,mxREAL)); 
    }
    mxSetFieldByNumber(theStruct,0,19,intMat2MatlabDoubles(&msg->spare2,1,1));
}

void AIS_8_1_15_ToMatlab(libais::Ais8_1_15 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions) {
//This message format is described in IMO circular 236, which was revoked
//in IMO circular 289 in 2013. Message "EXTENDED SHIP STATIC AND VOYAGE
//RELATED DATA"
    const size_t numberOfFields=8;
    const mwSize dims[2] = {1, 1};
    //These are the names of all of the fields that will be added to the
    //Matlab structre array.
    const char *fieldNames[numberOfFields] = {"message_id", "repeat_indicator",
    "mmsi", "spare", "dac", "fi","air_draught", "spare2"};
    mxArray *theStruct;
    
    //If text descriptions of the fields are desired, then put them into a
    //structure.
    if(fieldDescriptions!=NULL) {
        const char *descriptionStrings[numberOfFields] = {
            "Identifier for Message 8; always 8.",//0
            "Repeat indicator. Used by the repeater to indicate how many times a message has been repeated. 0 = default; 3 = do not repeat any more.",//1
            "Maritime mobile service identity (MMSI) number of source station",//2
            "Not used. Should be zero. Reserved for future use",//3
            "Designated area code (DAC)",//4
            "Function identifier (FI)",//5
            "Air draught, heigh over keel, in meters (0.1-204.7), 204.7 mean 204.7 meters or higher)",//6
            "Not used. Should be set to zero."//7
        };
        unsigned int i;

        *fieldDescriptions=mxCreateStructArray(2, dims, numberOfFields, fieldNames);

        for(i=0;i<numberOfFields;i++) {
            mxSetFieldByNumber(*fieldDescriptions,0,static_cast<int>(i),mxCreateCharMatrixFromStrings(1,&descriptionStrings[i]));
        }
    }
    
    //Create the Matlab structure array.
    theStruct=mxCreateStructArray(2, dims, numberOfFields, fieldNames);
    *decodedMessage=theStruct;
    
    //Fill all of the elements of the structure with the set that is common
    //to all AIS8 messages.
    mxSetFieldByNumber(theStruct,0,0,intMat2MatlabDoubles(&msg->message_id,1,1));
    mxSetFieldByNumber(theStruct,0,1,intMat2MatlabDoubles(&msg->repeat_indicator,1,1));
    mxSetFieldByNumber(theStruct,0,2,intMat2MatlabDoubles(&msg->mmsi,1,1));
    mxSetFieldByNumber(theStruct,0,3,intMat2MatlabDoubles(&msg->spare,1,1));
    mxSetFieldByNumber(theStruct,0,4,intMat2MatlabDoubles(&msg->dac,1,1));
    mxSetFieldByNumber(theStruct,0,5,intMat2MatlabDoubles(&msg->fi,1,1));

    //Fields specific to the Ais8_1_15 type.
    if(msg->air_draught!=0) {
        mxSetFieldByNumber(theStruct,0,6,floatMat2MatlabDoubles(&msg->air_draught,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,6,mxCreateDoubleMatrix(0,0,mxREAL)); 
    }
    mxSetFieldByNumber(theStruct,0,7,intMat2MatlabDoubles(&msg->spare2,1,1));
}

void AIS_8_1_16_ToMatlab(libais::Ais8_1_16 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions) {
//This message format is described in IMO circulars 236 and 289.
//Message "Number of persons on board"
    const size_t numberOfFields=8;
    const mwSize dims[2] = {1, 1};
    //These are the names of all of the fields that will be added to the
    //Matlab structre array.
    const char *fieldNames[numberOfFields] = {"message_id", "repeat_indicator",
    "mmsi", "spare", "dac", "fi","persons", "spare2"};
    mxArray *theStruct;

    //If text descriptions of the fields are desired, then put them into a
    //structure.
    if(fieldDescriptions!=NULL) {
        const char *descriptionStrings[numberOfFields] = {
            "Identifier for Message 8; always 8.",//0
            "Repeat indicator. Used by the repeater to indicate how many times a message has been repeated. 0 = default; 3 = do not repeat any more.",//1
            "Maritime mobile service identity (MMSI) number of source station",//2
            "Not used. Should be zero. Reserved for future use",//3
            "Designated area code (DAC)",//4
            "Function identifier (FI)",//5
            "Number of persons currently on board including crew members (1-8190), 8190 means 8190 or greater",//6
            "Not used. Set to zero."//7
        };
        unsigned int i;

        *fieldDescriptions=mxCreateStructArray(2, dims, numberOfFields, fieldNames);

        for(i=0;i<numberOfFields;i++) {
            mxSetFieldByNumber(*fieldDescriptions,0,static_cast<int>(i),mxCreateCharMatrixFromStrings(1,&descriptionStrings[i]));
        }
    }
    
    //Create the Matlab structure array.
    theStruct=mxCreateStructArray(2, dims, numberOfFields, fieldNames);
    *decodedMessage=theStruct;
    
    //Fill all of the elements of the structure with the set that is common
    //to all AIS8 messages.
    mxSetFieldByNumber(theStruct,0,0,intMat2MatlabDoubles(&msg->message_id,1,1));
    mxSetFieldByNumber(theStruct,0,1,intMat2MatlabDoubles(&msg->repeat_indicator,1,1));
    mxSetFieldByNumber(theStruct,0,2,intMat2MatlabDoubles(&msg->mmsi,1,1));
    mxSetFieldByNumber(theStruct,0,3,intMat2MatlabDoubles(&msg->spare,1,1));
    mxSetFieldByNumber(theStruct,0,4,intMat2MatlabDoubles(&msg->dac,1,1));
    mxSetFieldByNumber(theStruct,0,5,intMat2MatlabDoubles(&msg->fi,1,1));
    
    //Fields specific to the Ais8_1_16 type.
    if(msg->persons!=0) {
        mxSetFieldByNumber(theStruct,0,6,intMat2MatlabDoubles(&msg->persons,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,6,mxCreateDoubleMatrix(0,0,mxREAL)); 
    }
    mxSetFieldByNumber(theStruct,0,7,intMat2MatlabDoubles(&msg->spare2,1,1));
}

void AIS_8_1_17_ToMatlab(libais::Ais8_1_17 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions) {
//This message format is described in IMO circulars 236 and 289.
//Message "VTS-generated/Synthetic targets"
    const size_t numberOfFields=7;
    const mwSize dims[2] = {1, 1};
    //These are the names of all of the fields that will be added to the
    //Matlab structre array.
    const char *fieldNames[numberOfFields] = {"message_id", "repeat_indicator",
    "mmsi", "spare", "dac", "fi","targets"};
    const size_t numberOfTargetFields=8;
    mwSize targetDims[2];
    const char *targetFieldNames[numberOfTargetFields] = {"type", "id",
    "spare", "lon", "lat", "cog", "timestamp" "sog"};
    mxArray *theStruct, *targetArray;
    size_t curTarget;
    
    //If text descriptions of the fields are desired, then put them into a
    //structure.
    if(fieldDescriptions!=NULL) {
        const char *descriptionStrings[numberOfFields] = {
            "Identifier for Message 8; always 8.",//0
            "Repeat indicator. Used by the repeater to indicate how many times a message has been repeated. 0 = default; 3 = do not repeat any more.",//1
            "Maritime mobile service identity (MMSI) number of source station",//2
            "Not used. Should be zero. Reserved for future use",//3
            "Designated area code (DAC)",//4
            "Function identifier (FI)",//5
            ""//6 This is not used. Instead, a strcture describing the target parameters is used.
        };
        const char *targetDescriptions[numberOfTargetFields] = {
            "The type of the target identifier:\n0 = The target identifier is the MMSI number\n1 = The target identifier is the IMO number\n2 = The target identifier is the call sign\n3 = Other (default)",//0
           "Target identifier as a string of characters",//1
           "Spare. Set to zero.",//2
           "Longitude of the target in degrees East (WGS-84)",//3
           "Latitude of the target in degrees North (WGS-84)",//4
           "Course over ground in degrees East of north (0-359)",//5
           "UTC second when the report was generated (0-59)",//6
           "Speed over ground in knots (0-254)"//7
        };
        unsigned int i;
        mxArray *targetFields;

        *fieldDescriptions=mxCreateStructArray(2, dims, numberOfFields, fieldNames);
        for(i=0;i<(numberOfFields-1);i++) {
            mxSetFieldByNumber(*fieldDescriptions,0,static_cast<int>(i),mxCreateCharMatrixFromStrings(1,&descriptionStrings[i]));
        }
        
        targetFields=mxCreateStructArray(2, dims, numberOfTargetFields, targetFieldNames);
        mxSetFieldByNumber(*fieldDescriptions,0,numberOfFields-1,targetFields);

        for(i=0;i<numberOfTargetFields;i++) {
            mxSetFieldByNumber(targetFields,0,static_cast<int>(i),mxCreateCharMatrixFromStrings(1,&targetDescriptions[i]));
        }
    }
    
    //Create the Matlab structure array.
    theStruct=mxCreateStructArray(2, dims, numberOfFields, fieldNames);
    *decodedMessage=theStruct;
    
    //Fill all of the elements of the structure with the set that is common
    //to all AIS8 messages.
    mxSetFieldByNumber(theStruct,0,0,intMat2MatlabDoubles(&msg->message_id,1,1));
    mxSetFieldByNumber(theStruct,0,1,intMat2MatlabDoubles(&msg->repeat_indicator,1,1));
    mxSetFieldByNumber(theStruct,0,2,intMat2MatlabDoubles(&msg->mmsi,1,1));
    mxSetFieldByNumber(theStruct,0,3,intMat2MatlabDoubles(&msg->spare,1,1));
    mxSetFieldByNumber(theStruct,0,4,intMat2MatlabDoubles(&msg->dac,1,1));
    mxSetFieldByNumber(theStruct,0,5,intMat2MatlabDoubles(&msg->fi,1,1));
    
    //Fields specific to the Ais8_1_17 type.
    targetDims[0]=msg->targets.size();
    targetDims[1]=1;
    targetArray=mxCreateStructArray(2, dims, numberOfTargetFields, targetFieldNames);
    mxSetFieldByNumber(theStruct,0,6,targetArray);
    for (curTarget = 0; curTarget < msg->targets.size(); curTarget++) {
        mxSetFieldByNumber(targetArray,curTarget,0,intMat2MatlabDoubles(&msg->targets[curTarget].type,1,1));
        if(msg->targets[curTarget].id.compare("@@@@@@@")!=0) {
            const char *charString=msg->targets[curTarget].id.c_str();
            mxSetFieldByNumber(targetArray,curTarget,1,mxCreateCharMatrixFromStrings(1,&charString));
        } else {//No target ID
            mxSetFieldByNumber(targetArray,curTarget,1,mxCreateDoubleMatrix(0,0,mxREAL)); 
        }
        mxSetFieldByNumber(targetArray,curTarget,2,intMat2MatlabDoubles(&msg->targets[curTarget].spare,1,1));
        if(fabs(msg->targets[curTarget].position.lng_deg)<=180) {
            mxSetFieldByNumber(targetArray,curTarget,3,doubleMat2Matlab(&msg->targets[curTarget].position.lng_deg,1,1));
        } else {//No longitude
            mxSetFieldByNumber(targetArray,curTarget,3,mxCreateDoubleMatrix(0,0,mxREAL)); 
        }
        if(fabs(msg->targets[curTarget].position.lat_deg)<=90) {
            mxSetFieldByNumber(targetArray,curTarget,4,doubleMat2Matlab(&msg->targets[curTarget].position.lat_deg,1,1));
        } else {//No latitude
            mxSetFieldByNumber(targetArray,curTarget,4,mxCreateDoubleMatrix(0,0,mxREAL)); 
        }
        if(msg->targets[curTarget].cog<=359) {
            mxSetFieldByNumber(targetArray,curTarget,5,intMat2MatlabDoubles(&msg->targets[curTarget].cog,1,1));
        } else {//No course over ground
            mxSetFieldByNumber(targetArray,curTarget,5,mxCreateDoubleMatrix(0,0,mxREAL)); 
        }
        if(msg->targets[curTarget].timestamp<=59) {
            mxSetFieldByNumber(targetArray,curTarget,6,intMat2MatlabDoubles(&msg->targets[curTarget].timestamp,1,1));
        } else {//No timestamp
            mxSetFieldByNumber(targetArray,curTarget,6,mxCreateDoubleMatrix(0,0,mxREAL)); 
        }
        if(msg->targets[curTarget].sog<=254) {
            mxSetFieldByNumber(targetArray,curTarget,7,intMat2MatlabDoubles(&msg->targets[curTarget].sog,1,1));
        } else {//No speed over ground
            mxSetFieldByNumber(targetArray,curTarget,7,mxCreateDoubleMatrix(0,0,mxREAL)); 
        }
    }
}

void AIS_8_1_19_ToMatlab(libais::Ais8_1_19 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions) {
//This message format is described in IMO circular 289.
//Message "Marine Traffic Signal"
    const size_t numberOfFields=16;
    const mwSize dims[2] = {1, 1};
    //These are the names of all of the fields that will be added to the
    //Matlab structre array.
    const char *fieldNames[numberOfFields] = {"message_id", "repeat_indicator",
    "mmsi", "spare", "dac", "fi","link_id", "name", "lon", "lat", "status",
    "signal", "utc_hour_next", "utc_min_next", "next_signal", "spare2"};
    mxArray *theStruct;
    
    //If text descriptions of the fields are desired, then put them into a
    //structure.
    if(fieldDescriptions!=NULL) {
        const char *descriptionStrings[numberOfFields] = {
            "Identifier for Message 8; always 8.",//0
            "Repeat indicator. Used by the repeater to indicate how many times a message has been repeated. 0 = default; 3 = do not repeat any more.",//1
            "Maritime mobile service identity (MMSI) number of source station",//2
            "Not used. Should be zero. Reserved for future use",//3
            "Designated area code (DAC)",//4
            "Function identifier (FI)",//5
            "Message Linkage ID: A source specific running number, unique across all binary messages equipped with Message Linkage ID. Used to link additional information to the message by a Text Description message. The Message Linkage ID, source MMSI and message FI uniquely identifies the sent message. (1-1023)",//6
            "Name of Signal Station as defined in ITU-R M. 1371-3, Table 44",//7
            "Longitude of station in degrees East (WGS-84)",//8
            "Latitude of station in degrees North (WGS-84)",//9
            "Status of Signal: 1 = In regular service, 2 = Irregular service, 3 reserved for future use",//10
            "Signal in service message value as defined in Table 8.2 of IMO circular 289",//11
            "UTC hour of next signal shift (0-23)",//12
            "UTC minute of next signal shift (0-59)",//13
            "Expected next signal as defined in Table 8.2 of IMO circular 289"
        };
        unsigned int i;

        *fieldDescriptions=mxCreateStructArray(2, dims, numberOfFields, fieldNames);

        for(i=0;i<numberOfFields;i++) {
            mxSetFieldByNumber(*fieldDescriptions,0,static_cast<int>(i),mxCreateCharMatrixFromStrings(1,&descriptionStrings[i]));
        }
    }
    
    //Create the Matlab structure array.
    theStruct=mxCreateStructArray(2, dims, numberOfFields, fieldNames);
    *decodedMessage=theStruct;
    
    //Fill all of the elements of the structure with the set that is common
    //to all AIS8 messages.
    mxSetFieldByNumber(theStruct,0,0,intMat2MatlabDoubles(&msg->message_id,1,1));
    mxSetFieldByNumber(theStruct,0,1,intMat2MatlabDoubles(&msg->repeat_indicator,1,1));
    mxSetFieldByNumber(theStruct,0,2,intMat2MatlabDoubles(&msg->mmsi,1,1));
    mxSetFieldByNumber(theStruct,0,3,intMat2MatlabDoubles(&msg->spare,1,1));
    mxSetFieldByNumber(theStruct,0,4,intMat2MatlabDoubles(&msg->dac,1,1));
    mxSetFieldByNumber(theStruct,0,5,intMat2MatlabDoubles(&msg->fi,1,1));
    
    //Fields specific to the Ais8_1_19 type.
    if(msg->link_id!=0) {
        mxSetFieldByNumber(theStruct,0,6,intMat2MatlabDoubles(&msg->link_id,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,6,mxCreateDoubleMatrix(0,0,mxREAL)); 
    }
    if(msg->name.compare("@@@@@@@@@@@@@@@@@@@@")!=0) {
        const char *charString=msg->name.c_str();
        mxSetFieldByNumber(theStruct,0,7,mxCreateCharMatrixFromStrings(1,&charString));
    } else {
        mxSetFieldByNumber(theStruct,0,7,mxCreateDoubleMatrix(0,0,mxREAL)); 
    }
    if(fabs(msg->position.lng_deg)<=180) {
        mxSetFieldByNumber(theStruct,0,8,doubleMat2Matlab(&msg->position.lng_deg,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,8,mxCreateDoubleMatrix(0,0,mxREAL)); 
    }
    if(fabs(msg->position.lat_deg)<=90) {
        mxSetFieldByNumber(theStruct,0,9,doubleMat2Matlab(&msg->position.lat_deg,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,9,mxCreateDoubleMatrix(0,0,mxREAL)); 
    }
    if(msg->status!=0) {
        mxSetFieldByNumber(theStruct,0,10,intMat2MatlabDoubles(&msg->status,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,10,mxCreateDoubleMatrix(0,0,mxREAL)); 
    }
    if(msg->signal!=0) {
        mxSetFieldByNumber(theStruct,0,11,intMat2MatlabDoubles(&msg->signal,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,11,mxCreateDoubleMatrix(0,0,mxREAL)); 
    }
    if(msg->utc_hour_next<=23) {
        mxSetFieldByNumber(theStruct,0,12,intMat2MatlabDoubles(&msg->utc_hour_next,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,12,mxCreateDoubleMatrix(0,0,mxREAL)); 
    }
    if(msg->utc_min_next<=59) {
        mxSetFieldByNumber(theStruct,0,13,intMat2MatlabDoubles(&msg->utc_min_next,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,13,mxCreateDoubleMatrix(0,0,mxREAL)); 
    }
    if(msg->next_signal!=0) {
        mxSetFieldByNumber(theStruct,0,14,intMat2MatlabDoubles(&msg->next_signal,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,14,mxCreateDoubleMatrix(0,0,mxREAL)); 
    }
    mxSetFieldByNumber(theStruct,0,15,intMat2MatlabDoubles(msg->spare2.data(),4,1));
}

void AIS_8_1_21_ToMatlab(libais::Ais8_1_21 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions) {
//This message format is described in IMO circular 289.
//Message "Weather observation report from ship"
    const mwSize dims[2] = {1, 1};
    mxArray *theStruct;

    switch(msg->type_wx_report) {
        case 0:
        {
            const size_t numberOfFields=29;
            //These are the names of all of the fields that will be added to the
            //Matlab structre array.
            const char *fieldNames[numberOfFields] = {"message_id",
            "repeat_indicator", "mmsi","spare","dac", "fi",
            "type_wx_report", "location", "lon", "lat", "utc_day", "utc_hour",
            "utc_min", "wx","horz_viz", "humidity", "wind_speed", "wind_dir",
            "pressure", "pressure_tendency", "air_temp", "water_temp",
            "wave_period", "wave_height", "wave_dir", "swell_height",
            "swell_dir", "swell_period", "spare2"};
             if(fieldDescriptions!=NULL) {
                const char *descriptionStrings[numberOfFields] = {
                "Identifier for Message 6; always 6.",//0
                "Repeat indicator. Used by the repeater to indicate how many times a message has been repeated. 0 = default; 3 = do not repeat any more.",//1
                "Maritime mobile service identity (MMSI) number of source station",//2
                "Not used. Should be zero. Reserved for future use",//3
                "Designated area code (DAC)",//4
                "Function identifier (FI)",//5
                "Type of weather report. Here, 1=WMO Weather observation report from ship",//6
                "Geographic location as an ASCII string as defined in ITU-R M. 1371-3, Table 44",//7
                "Longitude of observation in degrees East (WGS-84)",//8
                "Latitude of observation in degrees North (WGS-84)",//9
                "UTC day of observation (1-31)",//10
                "UTC hour of observation (0-24)",//11
                "UTC minute of observation (0-59)",//12
                "Present weather based on WMI Code 45501:\n0 = clear (no clouds at any level)\n1 = cloudy\n2 = rain\n3 = fog\n4 = snow\n5 = typhoon/hurricane\n6 = monsoon\n7 = thunderstorm\n8 = not available (an empty matrix is returned instead)\n9 - 15 (reserved for future use)",//13
                 "Horizontal visibility in nautical miles (0-12.6)",//14 NOTE: The first bit being set makes this an invalid value, and indicates that the equipment was saturated.
                 "Relative humidity in percent (0-100)",//15
                 "Average wind speed in knots (0-126)",//16
                 "Wind direction in degrees East of North",//17
                 "Air pressure in hectoPascals at sea level (799-1201)",//18
                 "Air pressure tendency using WMO FM13 codes for pressure characteristic over the last three hours.",//19
                 "Dry bulb air temperature in degrees Celsius (-60 to 60)",//20
                 "Temperature of the water in degrees Celsius (-10 to 50)",//21
                 "Wave period in seconds",//22
                 "Height of the waves in meters",//23
                 "Wave direction in degrees East of North (0-359)",//24
                 "Swell height in meters",//25
                 "Swell direction in degrees East of North (0-359)",//26
                 "Swell period in seconds",//27
                 "Not used. Set to zero."//28
                };
                unsigned int i;

                *fieldDescriptions=mxCreateStructArray(2, dims, numberOfFields, fieldNames);

                for(i=0;i<numberOfFields;i++) {
                    mxSetFieldByNumber(*fieldDescriptions,0,static_cast<int>(i),mxCreateCharMatrixFromStrings(1,&descriptionStrings[i]));
                }
            };

            //Create the Matlab structure array.
            theStruct=mxCreateStructArray(2, dims, numberOfFields, fieldNames);
            *decodedMessage=theStruct;
            
            //Fill in the components unique to this message.
            if(msg->location.compare("@@@@@@@@@@@@@@@@@@@@")!=0) {
                const char *charString=msg->location.c_str();
                mxSetFieldByNumber(theStruct,0,7,mxCreateCharMatrixFromStrings(1,&charString));
            } else {
                mxSetFieldByNumber(theStruct,0,7,mxCreateDoubleMatrix(0,0,mxREAL)); 
            }
            if(fabs(msg->position.lng_deg)<=180) {
                mxSetFieldByNumber(theStruct,0,8,doubleMat2Matlab(&msg->position.lng_deg,1,1));
            } else {//No longitude
                mxSetFieldByNumber(theStruct,0,8,mxCreateDoubleMatrix(0,0,mxREAL)); 
            }
            if(fabs(msg->position.lat_deg)<=90) {
                mxSetFieldByNumber(theStruct,0,9,doubleMat2Matlab(&msg->position.lat_deg,1,1));
            } else {//No latitude
                mxSetFieldByNumber(theStruct,0,9,mxCreateDoubleMatrix(0,0,mxREAL)); 
            }
            if(msg->utc_day!=0) {
                mxSetFieldByNumber(theStruct,0,10,intMat2MatlabDoubles(&msg->utc_day,1,1));
            } else {//No UTC day
                mxSetFieldByNumber(theStruct,0,10,mxCreateDoubleMatrix(0,0,mxREAL)); 
            }
            if(msg->utc_hour<=23) {
                mxSetFieldByNumber(theStruct,0,11,intMat2MatlabDoubles(&msg->utc_hour,1,1));
            } else{//No UTC hour
                mxSetFieldByNumber(theStruct,0,11,mxCreateDoubleMatrix(0,0,mxREAL)); 
            }
            if(msg->utc_min<=59) {
                mxSetFieldByNumber(theStruct,0,12,intMat2MatlabDoubles(&msg->utc_min,1,1));
            } else {//No UTC minutes
                mxSetFieldByNumber(theStruct,0,12,mxCreateDoubleMatrix(0,0,mxREAL)); 
            }
            mxSetFieldByNumber(theStruct,0,13,intMat2MatlabDoubles(&msg->wx[0],1,1));
            if(msg->horz_viz!=12.7) {
                mxSetFieldByNumber(theStruct,0,14,floatMat2MatlabDoubles(&msg->horz_viz,1,1));
            } else {//No Horizontal visibility
                mxSetFieldByNumber(theStruct,0,14,mxCreateDoubleMatrix(0,0,mxREAL)); 
            }
            if(msg->humidity<=100) {
                mxSetFieldByNumber(theStruct,0,15,intMat2MatlabDoubles(&msg->humidity,1,1));
            } else {//No humidity
                mxSetFieldByNumber(theStruct,0,15,mxCreateDoubleMatrix(0,0,mxREAL));
            }
            if(msg->wind_speed<=126) {
                mxSetFieldByNumber(theStruct,0,16,intMat2MatlabDoubles(&msg->wind_speed,1,1));
            } else {//No wind speed
                mxSetFieldByNumber(theStruct,0,16,mxCreateDoubleMatrix(0,0,mxREAL));
            }
            if(msg->wind_dir<=359) {
                mxSetFieldByNumber(theStruct,0,17,intMat2MatlabDoubles(&msg->wind_dir,1,1));
            } else {
                mxSetFieldByNumber(theStruct,0,17,mxCreateDoubleMatrix(0,0,mxREAL));
            }
            if(msg->pressure<=1201) {
                mxSetFieldByNumber(theStruct,0,18,floatMat2MatlabDoubles(&msg->pressure,1,1));
            } else {
                mxSetFieldByNumber(theStruct,0,18,mxCreateDoubleMatrix(0,0,mxREAL));
            }
            mxSetFieldByNumber(theStruct,0,19,intMat2MatlabDoubles(&msg->pressure_tendency,1,1));
            if(msg->air_temp>=-60&&msg->air_temp<=60) {
                mxSetFieldByNumber(theStruct,0,20,floatMat2MatlabDoubles(&msg->air_temp,1,1));
            } else {
                mxSetFieldByNumber(theStruct,0,20,mxCreateDoubleMatrix(0,0,mxREAL));
            }
            if(msg->water_temp>=-10&&msg->water_temp<=50) {
                mxSetFieldByNumber(theStruct,0,21,floatMat2MatlabDoubles(&msg->water_temp,1,1));
            } else {
                mxSetFieldByNumber(theStruct,0,21,mxCreateDoubleMatrix(0,0,mxREAL));
            }
            if(msg->wave_period<=60) {
                mxSetFieldByNumber(theStruct,0,22,intMat2MatlabDoubles(&msg->wave_period,1,1));
            } else {
                mxSetFieldByNumber(theStruct,0,22,mxCreateDoubleMatrix(0,0,mxREAL));
            }
            if(msg->wave_height<=25.1) {
                mxSetFieldByNumber(theStruct,0,23,floatMat2MatlabDoubles(&msg->wave_height,1,1));
            } else {
                mxSetFieldByNumber(theStruct,0,23,mxCreateDoubleMatrix(0,0,mxREAL));
            }
            if(msg->wave_dir<=359) {
                mxSetFieldByNumber(theStruct,0,24,intMat2MatlabDoubles(&msg->wave_dir,1,1));
            } else {
                mxSetFieldByNumber(theStruct,0,24,mxCreateDoubleMatrix(0,0,mxREAL));
            }
            if(msg->swell_height<=25.1) {
                mxSetFieldByNumber(theStruct,0,25,floatMat2MatlabDoubles(&msg->swell_height,1,1));
            } else {
                mxSetFieldByNumber(theStruct,0,25,mxCreateDoubleMatrix(0,0,mxREAL));
            }
            if(msg->swell_dir<=359) {
                mxSetFieldByNumber(theStruct,0,26,intMat2MatlabDoubles(&msg->swell_dir,1,1));
            } else {
                mxSetFieldByNumber(theStruct,0,26,mxCreateDoubleMatrix(0,0,mxREAL));
            }
            if(msg->swell_period<=60) {
                mxSetFieldByNumber(theStruct,0,27,intMat2MatlabDoubles(&msg->swell_period,1,1));
            } else {
                mxSetFieldByNumber(theStruct,0,27,mxCreateDoubleMatrix(0,0,mxREAL));
            }
            mxSetFieldByNumber(theStruct,0,28,intMat2MatlabDoubles(&msg->spare2,1,1));
        }
            break;
        case 1:
        {
            const size_t numberOfFields=54;
            //These are the names of all of the fields that will be added to the
            //Matlab structre array.
            const char *fieldNames[numberOfFields] = {"message_id",
            "repeat_indicator", "mmsi", "spare", "dac", "fi",
            "type_wx_report", "lon", "lat",
            "utc_month", "utc_day", "utc_hour", "utc_min", "cog", "sog",
            "heading", "pressure", "rel_pressure", "pressure_tendency",
            "wind_dir", "wind_speed_ms", "wind_dir_rel", "wind_speed_rel",
            "wind_gust_speed", "wind_gust_dir", "air_temp_raw", "humidity",
            "water_temp_raw", "horz_viz", "present_weather",
            "past_weather_1", "past_weather_2",
            "cloud_total", "cloud_low", "cloud_low_type",
            "cloud_middle_type", "cloud_high_type",
            "alt_lowest_cloud_base", "wave_period", "wave_height",
            "swell_dir", "swell_period", "swell_height", "swell_dir_2",
            "swell_period_2", "swell_height_2", "ice_thickness",
            "ice_accretion", "ice_accretion_cause",
            "sea_ice_concentration", "amt_type_ice", "ice_situation",
            "ice_devel", "bearing_ice_edge"};
            if(fieldDescriptions!=NULL) {
                const char *descriptionStrings[numberOfFields] = {
                "Identifier for Message 6; always 6.",//0
                "Repeat indicator. Used by the repeater to indicate how many times a message has been repeated. 0 = default; 3 = do not repeat any more.",//1
                "Maritime mobile service identity (MMSI) number of source station",//2
                "Not used. Should be zero. Reserved for future use",//3
                "Designated area code (DAC)",//4
                "Function identifier (FI)",//5
                "Type of weather report. Here, 0=WMO weather observation report from ship",//6
                "Longitude of observation in degrees East (WGS-84)",//7
                "Latitude of observation in degrees North (WGS-84)",//8
                "UTC month of observation (1-12)",//9
                "UTC day of observation (1-31)",//10
                "UTC hour of observation (0-24)",//11
                "UTC minute of observation (0-59)",//12
                "Ship course over ground for the past 10 minutes in degrees East of North (5-360). 0 degrees means stopped",//13
                "Ship speed over ground for the past 10 minutes in meters per second (0-15)",//14
                "Average heading of the ship for the past 10 minutes in degrees East of north (5-360). 0 degrees means stopped",//15
                "Air pressure in hectoPascals reduced to sea level (900-1100)",//16
                "3-hour air pressure change in hectoPascals (-50 to 50)",//17
                "Characteristic of pressure tendency from WMO BUFR table 010063",//18
                "True wind direction averaged over 10 minutes in degrees East of North (5-360). 0 means calm.",//19
                "True wind speed averaged over 10 minutes in meters per second (0-127)",//20
                "Relative wind direction averaged over 10 minutes in degrees East of North (5-360). 0 degrees means calm",//21
                "Relative wind speed averaged over 10 minutes in meters per second (0-127)",//22
                "Maximum wind gust speed in meters per second (0-127)",//23
                "Maximum wind gust direction in degrees East of North. 0 degrees means calm",//24
                "Dry bulb air temperature such that the temperature in degrees Kelvin=air_temp_raw/10+223",//25
                "Relative humidity in percent (0-100)",//26
                "Sea surface temperature such that the temperature in degrees Kelvin=water_temp_raw/10+268",//27
                "Horizontal visibility in meters (0-50000m)",//28
                "Present weather as given by codes 0-510 in WMO BUFR table 020003",//29
                "Past weather 1 as given by codes 0-30 in WMO BUFR table 020004",//30
                "Past weather 2 as given by codes 0-30 in WMO BUFR table 020005",//31
                "Total cloud cover in percent (0-100)",//32
                "Cloud amount low as given by WMO BUFR table 020011",//33
                "Cloud type low as given by WMO BUFR table 020012",//34
                "Cloud type middle as given by WMO BUFR table 020012",//35
                "Cloud type high as given by WMO BUFR table 020012",//36
                "Height of base of lowest cloud in meters (0-2540.16)",//37
                "Period of wind waves in seconds (0-30)",//38
                "Height of wind waves in meters (0-30)",//39
                "Direction from which the first swell is coming in degrees East of North (10-360). 0 means calm",//40
                "Period of the first swell in seconds (0-30)",//41
                "Height of the first swell in meters (0-30)",//42
                "Direction from which the second swell is coming in degrees East of North (10-360). 0 means calm",//43
                "Period of the second swell in seconds (0-30)",//44
                "Height of the second swell in meters (0-30)",//45
                "Ice deposit thickness in centimeters",//46
                "Rate of ice accretion as per WMO BUFR table 020032",//47
                "Cause of ice accretion as per WMO BUFR table 020033",//48
                "Sea ice concentration as per WMO BUFR table 020034",//49
                "Amount and type of ice as per WMO BUFR table 020035",//50
                "Ice situation as per WMO BUFR table 020035",//51
                "Ice development as per WMO BUFR table 020037",//52
                "Bearing of ice edge in degrees East of North (45-360)."//53
                };
                unsigned int i;

                *fieldDescriptions=mxCreateStructArray(2, dims, numberOfFields, fieldNames);

                for(i=0;i<numberOfFields;i++) {
                    mxSetFieldByNumber(*fieldDescriptions,0,static_cast<int>(i),mxCreateCharMatrixFromStrings(1,&descriptionStrings[i]));
                }
             };
            
            //Create the Matlab structure array.
            theStruct=mxCreateStructArray(2, dims, numberOfFields, fieldNames);
            *decodedMessage=theStruct;
            
            if(fabs(msg->position.lng_deg)<=180) {
                mxSetFieldByNumber(theStruct,0,7,doubleMat2Matlab(&msg->position.lng_deg,1,1));
            } else {
                mxSetFieldByNumber(theStruct,0,7,mxCreateDoubleMatrix(0,0,mxREAL));
            }
            if(fabs(msg->position.lat_deg)<=90) {
                mxSetFieldByNumber(theStruct,0,8,doubleMat2Matlab(&msg->position.lat_deg,1,1));
            } else {
                mxSetFieldByNumber(theStruct,0,8,mxCreateDoubleMatrix(0,0,mxREAL));
            }
            if(msg->utc_month!=0) {
                mxSetFieldByNumber(theStruct,0,9,intMat2MatlabDoubles(&msg->utc_month,1,1));
            } else {
                mxSetFieldByNumber(theStruct,0,9,mxCreateDoubleMatrix(0,0,mxREAL));
            }
            if(msg->utc_day!=0) {
                mxSetFieldByNumber(theStruct,0,10,intMat2MatlabDoubles(&msg->utc_day,1,1));
            } else {
                mxSetFieldByNumber(theStruct,0,10,mxCreateDoubleMatrix(0,0,mxREAL));
            }
            if(msg->utc_hour<=23) {
                mxSetFieldByNumber(theStruct,0,11,intMat2MatlabDoubles(&msg->utc_hour,1,1));
            } else {
                mxSetFieldByNumber(theStruct,0,11,mxCreateDoubleMatrix(0,0,mxREAL));
            }
            if(msg->utc_min<=59) {
                mxSetFieldByNumber(theStruct,0,12,intMat2MatlabDoubles(&msg->utc_min,1,1));
            } else {
                mxSetFieldByNumber(theStruct,0,12,mxCreateDoubleMatrix(0,0,mxREAL));
            }
            if(msg->cog<=360) {
                mxSetFieldByNumber(theStruct,0,13,intMat2MatlabDoubles(&msg->cog,1,1));
            } else { 
                mxSetFieldByNumber(theStruct,0,13,mxCreateDoubleMatrix(0,0,mxREAL));
            }
            if(msg->sog<=15) {
                mxSetFieldByNumber(theStruct,0,14,floatMat2MatlabDoubles(&msg->sog,1,1));
            } else {
                mxSetFieldByNumber(theStruct,0,14,mxCreateDoubleMatrix(0,0,mxREAL));
            }
            if(msg->heading<=360){
                mxSetFieldByNumber(theStruct,0,15,intMat2MatlabDoubles(&msg->heading,1,1));
            } else {
                mxSetFieldByNumber(theStruct,0,15,mxCreateDoubleMatrix(0,0,mxREAL));
            }
            if(msg->pressure<=1100) {
                mxSetFieldByNumber(theStruct,0,16,floatMat2MatlabDoubles(&msg->pressure,1,1));
            } else {
                mxSetFieldByNumber(theStruct,0,16,mxCreateDoubleMatrix(0,0,mxREAL));
            }
            if(msg->rel_pressure<=50&&msg->rel_pressure<=50) {
                mxSetFieldByNumber(theStruct,0,17,floatMat2MatlabDoubles(&msg->rel_pressure,1,1));
            } else {
                mxSetFieldByNumber(theStruct,0,17,mxCreateDoubleMatrix(0,0,mxREAL));
            }

            mxSetFieldByNumber(theStruct,0,18,intMat2MatlabDoubles(&msg->pressure_tendency,1,1));
            if(msg->wind_dir<=360) {
                mxSetFieldByNumber(theStruct,0,19,intMat2MatlabDoubles(&msg->wind_dir,1,1));
            } else {
                mxSetFieldByNumber(theStruct,0,19,mxCreateDoubleMatrix(0,0,mxREAL));
            }
            if(msg->wind_speed_ms<=127){
                mxSetFieldByNumber(theStruct,0,20,floatMat2MatlabDoubles(&msg->wind_speed_ms,1,1));
            } else {
                mxSetFieldByNumber(theStruct,0,20,mxCreateDoubleMatrix(0,0,mxREAL));
            }
            if(msg->wind_dir_rel<=360) {
                mxSetFieldByNumber(theStruct,0,21,intMat2MatlabDoubles(&msg->wind_dir_rel,1,1));
            } else {
                mxSetFieldByNumber(theStruct,0,21,mxCreateDoubleMatrix(0,0,mxREAL));
            }
            if(msg->wind_speed_rel<=127) {
                mxSetFieldByNumber(theStruct,0,22,floatMat2MatlabDoubles(&msg->wind_speed_rel,1,1));
            } else {
                mxSetFieldByNumber(theStruct,0,22,mxCreateDoubleMatrix(0,0,mxREAL));
            }
            if(msg->wind_gust_speed<=127) {
                mxSetFieldByNumber(theStruct,0,23,floatMat2MatlabDoubles(&msg->wind_gust_speed,1,1));
            } else {
                mxSetFieldByNumber(theStruct,0,23,mxCreateDoubleMatrix(0,0,mxREAL));
            }
            if(msg->wind_gust_dir<=360) {
                mxSetFieldByNumber(theStruct,0,24,intMat2MatlabDoubles(&msg->wind_gust_dir,1,1));
            } else {
                mxSetFieldByNumber(theStruct,0,24,mxCreateDoubleMatrix(0,0,mxREAL));
            }
            if(msg->air_temp_raw!=1023) {
                mxSetFieldByNumber(theStruct,0,25,intMat2MatlabDoubles(&msg->air_temp_raw,1,1));
            } else {
                mxSetFieldByNumber(theStruct,0,25,mxCreateDoubleMatrix(0,0,mxREAL));
            }
            if(msg->humidity<=100) {
                mxSetFieldByNumber(theStruct,0,26,intMat2MatlabDoubles(&msg->humidity,1,1));
            } else {
                mxSetFieldByNumber(theStruct,0,26,mxCreateDoubleMatrix(0,0,mxREAL));
            }
            if(msg->water_temp_raw!=511) {
                mxSetFieldByNumber(theStruct,0,27,intMat2MatlabDoubles(&msg->water_temp_raw,1,1));
            } else {
                mxSetFieldByNumber(theStruct,0,27,mxCreateDoubleMatrix(0,0,mxREAL));
            } 
            if(msg->horz_viz<=50000) {
                mxSetFieldByNumber(theStruct,0,28,floatMat2MatlabDoubles(&msg->horz_viz,1,1));
            } else {//No Horizontal visibility
                mxSetFieldByNumber(theStruct,0,28,mxCreateDoubleMatrix(0,0,mxREAL)); 
            }
            mxSetFieldByNumber(theStruct,0,29,intMat2MatlabDoubles(&msg->wx[0],3,1));
            mxSetFieldByNumber(theStruct,0,30,intMat2MatlabDoubles(&msg->wx[1],3,1));
            mxSetFieldByNumber(theStruct,0,31,intMat2MatlabDoubles(&msg->wx[2],3,1));
            
            if(msg->cloud_total<=100) {
                mxSetFieldByNumber(theStruct,0,32,intMat2MatlabDoubles(&msg->cloud_total,1,1));
            } else {
                mxSetFieldByNumber(theStruct,0,32,mxCreateDoubleMatrix(0,0,mxREAL));
            }
            if(msg->cloud_low!=15) {
                mxSetFieldByNumber(theStruct,0,33,intMat2MatlabDoubles(&msg->cloud_low,1,1));
            } else {
                mxSetFieldByNumber(theStruct,0,33,mxCreateDoubleMatrix(0,0,mxREAL));
            }
            if(msg->cloud_low_type!=63) {
                mxSetFieldByNumber(theStruct,0,34,intMat2MatlabDoubles(&msg->cloud_low_type,1,1));
            } else {
                mxSetFieldByNumber(theStruct,0,34,mxCreateDoubleMatrix(0,0,mxREAL));
            }
            if(msg->cloud_middle_type!=63) {
                mxSetFieldByNumber(theStruct,0,35,intMat2MatlabDoubles(&msg->cloud_middle_type,1,1));
            } else {
                mxSetFieldByNumber(theStruct,0,35,mxCreateDoubleMatrix(0,0,mxREAL));
            }
            if(msg->cloud_high_type!=63) {
                mxSetFieldByNumber(theStruct,0,36,intMat2MatlabDoubles(&msg->cloud_high_type,1,1));
            } else {
                mxSetFieldByNumber(theStruct,0,36,mxCreateDoubleMatrix(0,0,mxREAL));
            }
            if(msg->alt_lowest_cloud_base<=2540.16) {
                mxSetFieldByNumber(theStruct,0,37,floatMat2MatlabDoubles(&msg->alt_lowest_cloud_base,1,1));
            } else {
                mxSetFieldByNumber(theStruct,0,37,mxCreateDoubleMatrix(0,0,mxREAL));
            }
            if(msg->wave_period<=30) {
                mxSetFieldByNumber(theStruct,0,38,intMat2MatlabDoubles(&msg->wave_period,1,1));
            } else {
                mxSetFieldByNumber(theStruct,0,38,mxCreateDoubleMatrix(0,0,mxREAL));
            }
            if(msg->wave_height<=30) {
                mxSetFieldByNumber(theStruct,0,39,floatMat2MatlabDoubles(&msg->wave_height,1,1));
            } else {
                mxSetFieldByNumber(theStruct,0,39,mxCreateDoubleMatrix(0,0,mxREAL));
            }
            if(msg->swell_dir<=360) {
                mxSetFieldByNumber(theStruct,0,40,intMat2MatlabDoubles(&msg->swell_dir,1,1));
            } else {
                mxSetFieldByNumber(theStruct,0,40,mxCreateDoubleMatrix(0,0,mxREAL));
            }
            if(msg->swell_period<=30) {
                mxSetFieldByNumber(theStruct,0,41,intMat2MatlabDoubles(&msg->swell_period,1,1));
            } else {
                mxSetFieldByNumber(theStruct,0,41,mxCreateDoubleMatrix(0,0,mxREAL));
            }
            if(msg->swell_height<=30) {
                mxSetFieldByNumber(theStruct,0,42,floatMat2MatlabDoubles(&msg->swell_height,1,1));
            } else {
                mxSetFieldByNumber(theStruct,0,42,mxCreateDoubleMatrix(0,0,mxREAL));
            }
            if(msg->swell_dir_2<=360) {
                mxSetFieldByNumber(theStruct,0,43,intMat2MatlabDoubles(&msg->swell_dir_2,1,1));
            } else {
                mxSetFieldByNumber(theStruct,0,43,mxCreateDoubleMatrix(0,0,mxREAL));
            }
            if(msg->swell_period_2<=30) {
                mxSetFieldByNumber(theStruct,0,44,intMat2MatlabDoubles(&msg->swell_period_2,1,1));
            } else {
                mxSetFieldByNumber(theStruct,0,44,mxCreateDoubleMatrix(0,0,mxREAL));
            }
            if(msg->swell_height_2<=30) {
                mxSetFieldByNumber(theStruct,0,45,floatMat2MatlabDoubles(&msg->swell_height_2,1,1));
            } else {
                mxSetFieldByNumber(theStruct,0,45,mxCreateDoubleMatrix(0,0,mxREAL));
            }
            if(msg->ice_thickness<=126) {
                mxSetFieldByNumber(theStruct,0,46,floatMat2MatlabDoubles(&msg->ice_thickness,1,1));
            } else {
                mxSetFieldByNumber(theStruct,0,46,mxCreateDoubleMatrix(0,0,mxREAL));
            }
            if(msg->ice_accretion!=7) {
                mxSetFieldByNumber(theStruct,0,47,intMat2MatlabDoubles(&msg->ice_accretion,1,1));
            } else {
                mxSetFieldByNumber(theStruct,0,47,mxCreateDoubleMatrix(0,0,mxREAL));
            }
            if(msg->ice_accretion_cause!=7) {
                mxSetFieldByNumber(theStruct,0,48,intMat2MatlabDoubles(&msg->ice_accretion_cause,1,1));
            } else {
                mxSetFieldByNumber(theStruct,0,48,mxCreateDoubleMatrix(0,0,mxREAL));
            }
            if(msg->sea_ice_concentration!=31) {
                mxSetFieldByNumber(theStruct,0,49,intMat2MatlabDoubles(&msg->sea_ice_concentration,1,1));
            } else {
                mxSetFieldByNumber(theStruct,0,49,mxCreateDoubleMatrix(0,0,mxREAL));
            }
            if(msg->amt_type_ice!=15) {
                mxSetFieldByNumber(theStruct,0,50,intMat2MatlabDoubles(&msg->amt_type_ice,1,1));
            } else {
                mxSetFieldByNumber(theStruct,0,50,mxCreateDoubleMatrix(0,0,mxREAL));
            }
            if(msg->ice_situation!=31) {
                mxSetFieldByNumber(theStruct,0,51,intMat2MatlabDoubles(&msg->ice_situation,1,1));
            } else {
                mxSetFieldByNumber(theStruct,0,51,mxCreateDoubleMatrix(0,0,mxREAL));
            }
            if(msg->ice_devel!=31) {
                mxSetFieldByNumber(theStruct,0,52,intMat2MatlabDoubles(&msg->ice_devel,1,1));
            } else {
                mxSetFieldByNumber(theStruct,0,52,mxCreateDoubleMatrix(0,0,mxREAL));
            }
            if(msg->bearing_ice_edge!=15) {
                mxSetFieldByNumber(theStruct,0,53,intMat2MatlabDoubles(&msg->bearing_ice_edge,1,1));
            } else {
                mxSetFieldByNumber(theStruct,0,53,mxCreateDoubleMatrix(0,0,mxREAL));
            }
        }
            break;
        default://If the specific type could not be identified, just
        {
            const size_t numberOfFields=7;
            //These are the names of all of the fields that will be added to the
            //Matlab structre array.
            const char *fieldNames[numberOfFields] = {"message_id",
            "repeat_indicator", "mmsi","spare","dac","fi",
            "type_wx_report"};
            //Return the elements common to all types.
            //Create the Matlab structure array.
            if(fieldDescriptions!=NULL) {
                const char *descriptionStrings[numberOfFields] = {
                "Identifier for Message 6; always 6.",//0
                "Repeat indicator. Used by the repeater to indicate how many times a message has been repeated. 0 = default; 3 = do not repeat any more.",//1
                "Maritime mobile service identity (MMSI) number of source station",//2
                "Not used. Should be zero. Reserved for future use",//3
                "Designated area code (DAC)",//4
                "Function identifier (FI)",//5
                "Type of ship weather report"//6
                };
                unsigned int i;

                *fieldDescriptions=mxCreateStructArray(2, dims, numberOfFields, fieldNames);

                for(i=0;i<numberOfFields;i++) {
                    mxSetFieldByNumber(*fieldDescriptions,0,static_cast<int>(i),mxCreateCharMatrixFromStrings(1,&descriptionStrings[i]));
                }
            }
            
            theStruct=mxCreateStructArray(2, dims, numberOfFields, fieldNames);
            *decodedMessage=theStruct;
        }
            break;
    }
        
    //Fill all of the elements of the structure with the set that is common
    //to all AIS8 messages.
    mxSetFieldByNumber(theStruct,0,0,intMat2MatlabDoubles(&msg->message_id,1,1));
    mxSetFieldByNumber(theStruct,0,1,intMat2MatlabDoubles(&msg->repeat_indicator,1,1));
    mxSetFieldByNumber(theStruct,0,2,intMat2MatlabDoubles(&msg->mmsi,1,1));
    mxSetFieldByNumber(theStruct,0,3,intMat2MatlabDoubles(&msg->spare,1,1));
    mxSetFieldByNumber(theStruct,0,4,intMat2MatlabDoubles(&msg->dac,1,1));
    mxSetFieldByNumber(theStruct,0,5,intMat2MatlabDoubles(&msg->fi,1,1));
  
    //Another field specific to this message.
    mxSetFieldByNumber(theStruct,0,6,intMat2MatlabDoubles(&msg->type_wx_report,1,1));
}

void AIS_8_1_22_ToMatlab(libais::Ais8_1_22 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions) {
//This message format is described in IMO circular 289.
//Message "Area notice (broadcast)"
    const size_t numberOfFields=14;
    const mwSize dims[2] = {1, 1};
    mwSize cellArrayDims[2];
    //These are the names of all of the fields that will be added to the
    //Matlab structre array.
    const char *fieldNames[numberOfFields] = {"message_id",
    "repeat_indicator", "mmsi", "spare", "dac", "fi", "link_id",
    "notice_type", "month", "day", "hour", "minute", "duration_minutes",
    "sub_areas"};
    mxArray *theStruct, *subareaCellArray;
    size_t curSubarea;
    
    //If text descriptions of the fields are desired, then put them into a
    //structure.
    if(fieldDescriptions!=NULL) {
        const char *descriptionStrings[numberOfFields] = {
            "Identifier for Message 8; always 8.",//0
            "Repeat indicator. Used by the repeater to indicate how many times a message has been repeated. 0 = default; 3 = do not repeat any more.",//1
            "Maritime mobile service identity (MMSI) number of source station",//2
            "Not used. Should be zero. Reserved for future use",//3
            "Designated area code (DAC)",//4
            "Function identifier (FI)",//5
            "A source specific running number (1-1023), unique across all binary messages equipped with Message Linkage ID. Used to link additional information to the message by a Text Description message. The Message Linkage ID and the first six digits of the source MMSI uniquely identify the sent message.",//6
            "Notice Description, 0-127 as per table 11.11 in IMO circular 289",//7
            "UTC month (1-12)",//8
            "UTC day (1-31)",//9
            "UTC hour (0-23)",//10
            "UTC minute (0-59)",//11
            "Duration in minutes until the end of the area notice measured from the start date and time of the notice. Max is 262142",//12
            "A structure of elements for the subareas as defined in IMO circular 289. The subareas define shapes"//13
        };
        unsigned int i;

        *fieldDescriptions=mxCreateStructArray(2, dims, numberOfFields, fieldNames);

        for(i=0;i<numberOfFields;i++) {
            mxSetFieldByNumber(*fieldDescriptions,0,static_cast<int>(i),mxCreateCharMatrixFromStrings(1,&descriptionStrings[i]));
        }
    }
    
    //Create the Matlab structure array.
    theStruct=mxCreateStructArray(2, dims, numberOfFields, fieldNames);
    *decodedMessage=theStruct;
    
    //Fill all of the elements of the structure with the set that is common
    //to all AIS8 messages.
    mxSetFieldByNumber(theStruct,0,0,intMat2MatlabDoubles(&msg->message_id,1,1));
    mxSetFieldByNumber(theStruct,0,1,intMat2MatlabDoubles(&msg->repeat_indicator,1,1));
    mxSetFieldByNumber(theStruct,0,2,intMat2MatlabDoubles(&msg->mmsi,1,1));
    mxSetFieldByNumber(theStruct,0,3,intMat2MatlabDoubles(&msg->spare,1,1));
    mxSetFieldByNumber(theStruct,0,4,intMat2MatlabDoubles(&msg->dac,1,1));
    mxSetFieldByNumber(theStruct,0,5,intMat2MatlabDoubles(&msg->fi,1,1));
    
    //Fill in the components unique to this message.
    if(msg->link_id!=0) {
        mxSetFieldByNumber(theStruct,0,6,intMat2MatlabDoubles(&msg->link_id,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,6,mxCreateDoubleMatrix(0,0,mxREAL)); 
    }
    mxSetFieldByNumber(theStruct,0,7,intMat2MatlabDoubles(&msg->notice_type,1,1));
    if(msg->month!=0) {
        mxSetFieldByNumber(theStruct,0,8,intMat2MatlabDoubles(&msg->month,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,8,mxCreateDoubleMatrix(0,0,mxREAL)); 
    }
    if(msg->day!=0) {
        mxSetFieldByNumber(theStruct,0,9,intMat2MatlabDoubles(&msg->day,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,9,mxCreateDoubleMatrix(0,0,mxREAL)); 
    }
    if(msg->hour<=23) {
        mxSetFieldByNumber(theStruct,0,10,intMat2MatlabDoubles(&msg->hour,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,10,mxCreateDoubleMatrix(0,0,mxREAL)); 
    }
    if(msg->minute<=59) {
        mxSetFieldByNumber(theStruct,0,11,intMat2MatlabDoubles(&msg->minute,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,11,mxCreateDoubleMatrix(0,0,mxREAL)); 
    }
    if(msg->duration_minutes!=262143) {
        mxSetFieldByNumber(theStruct,0,12,intMat2MatlabDoubles(&msg->duration_minutes,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,12,mxCreateDoubleMatrix(0,0,mxREAL)); 
    }

    //The subareas are put into a cell array.
    cellArrayDims[0]=msg->sub_areas.size();
    cellArrayDims[1]=1;
    //Make space to hold all of the reports in a cell array.
    subareaCellArray=mxCreateCellArray(2, cellArrayDims);
    mxSetFieldByNumber(theStruct,0,13,subareaCellArray);

    //Fill in all of the reports.
    for(curSubarea=0;curSubarea<msg->sub_areas.size();curSubarea++) {
        int area_shape=msg->sub_areas[curSubarea]->getType();
        switch(area_shape) {
            case libais::AIS8_1_22_SHAPE_CIRCLE:
            {
                const size_t numReportFields=6;
                const mwSize reportDims[2] = {1, 1};
                const char *reportFieldNames[numReportFields]={"area_shape",
                "lon", "lat", "precision", "radius_m", "spare"};                      
                mxArray *reportStruct;
                libais::Ais8_1_22_Circle *c =reinterpret_cast<libais::Ais8_1_22_Circle*>(msg->sub_areas[curSubarea]);
                
                //Create the Matlab structure array.
                reportStruct=mxCreateStructArray(2, reportDims, numReportFields, reportFieldNames);
                //Fill in the component common to all subareas.
                mxSetFieldByNumber(reportStruct,0,0,intMat2MatlabDoubles(&area_shape,1,1));

                //Fill in the components for this subarea
                if(fabs(c->position.lng_deg)<=180) {
                    mxSetFieldByNumber(reportStruct,0,1,doubleMat2Matlab(&c->position.lng_deg,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,1,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                if(fabs(c->position.lat_deg)<=90) {
                    mxSetFieldByNumber(reportStruct,0,2,doubleMat2Matlab(&c->position.lat_deg,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,2,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                mxSetFieldByNumber(reportStruct,0,3,intMat2MatlabDoubles(&c->precision,1,1));
                mxSetFieldByNumber(reportStruct,0,4,intMat2MatlabDoubles(&c->radius_m,1,1));
                mxSetFieldByNumber(reportStruct,0,5,mat2MatlabDoubles(&c->spare,1,1));
                
                mxSetCell(subareaCellArray,curSubarea,reportStruct);
            }
                break;
            case libais::AIS8_1_22_SHAPE_RECT:
            {
                const size_t numReportFields=8;
                const mwSize reportDims[2] = {1, 1};
                const char *reportFieldNames[numReportFields]={"area_shape",
                "lon", "lat", "precision", "e_dim_m", "n_dim_m", "orient_deg",
                "spare"};                      
                mxArray *reportStruct;
                libais::Ais8_1_22_Rect *c =reinterpret_cast<libais::Ais8_1_22_Rect*>(msg->sub_areas[curSubarea]);
                
                //Create the Matlab structure array.
                reportStruct=mxCreateStructArray(2, reportDims, numReportFields, reportFieldNames);
                //Fill in the component common to all subareas.
                mxSetFieldByNumber(reportStruct,0,0,intMat2MatlabDoubles(&area_shape,1,1));

                //Fill in the components for this subarea
                if(fabs(c->position.lng_deg)<=180) {
                    mxSetFieldByNumber(reportStruct,0,1,doubleMat2Matlab(&c->position.lng_deg,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,1,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                if(fabs(c->position.lat_deg)<=90) {
                    mxSetFieldByNumber(reportStruct,0,2,doubleMat2Matlab(&c->position.lat_deg,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,2,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                mxSetFieldByNumber(reportStruct,0,3,intMat2MatlabDoubles(&c->precision,1,1));
                mxSetFieldByNumber(reportStruct,0,4,intMat2MatlabDoubles(&c->e_dim_m,1,1));
                mxSetFieldByNumber(reportStruct,0,5,intMat2MatlabDoubles(&c->n_dim_m,1,1));
                mxSetFieldByNumber(reportStruct,0,6,intMat2MatlabDoubles(&c->orient_deg,1,1));
                mxSetFieldByNumber(reportStruct,0,7,mat2MatlabDoubles(&c->spare,1,1));
                
                mxSetCell(subareaCellArray,curSubarea,reportStruct);
            }
                break;
            case libais::AIS8_1_22_SHAPE_SECTOR:
            {
                const size_t numReportFields=7;
                const mwSize reportDims[2] = {1, 1};
                const char *reportFieldNames[numReportFields]={"area_shape",
                "lon", "lat", "precision", "radius_m", "left_bound_deg",
                "right_bound_deg"};                      
                mxArray *reportStruct;
                libais::Ais8_1_22_Sector *c =reinterpret_cast<libais::Ais8_1_22_Sector*>(msg->sub_areas[curSubarea]);
                
                //Create the Matlab structure array.
                reportStruct=mxCreateStructArray(2, reportDims, numReportFields, reportFieldNames);
                //Fill in the component common to all subareas.
                mxSetFieldByNumber(reportStruct,0,0,intMat2MatlabDoubles(&area_shape,1,1));

                //Fill in the components for this subarea
                if(fabs(c->position.lng_deg)<=180) {
                    mxSetFieldByNumber(reportStruct,0,1,doubleMat2Matlab(&c->position.lng_deg,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,1,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                if(fabs(c->position.lat_deg)<=90) {
                    mxSetFieldByNumber(reportStruct,0,2,doubleMat2Matlab(&c->position.lat_deg,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,2,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                mxSetFieldByNumber(reportStruct,0,3,intMat2MatlabDoubles(&c->precision,1,1));
                mxSetFieldByNumber(reportStruct,0,4,intMat2MatlabDoubles(&c->radius_m,1,1));
                mxSetFieldByNumber(reportStruct,0,5,intMat2MatlabDoubles(&c->left_bound_deg,1,1));
                mxSetFieldByNumber(reportStruct,0,6,intMat2MatlabDoubles(&c->right_bound_deg,1,1));

                mxSetCell(subareaCellArray,curSubarea,reportStruct);
            }
                break;
            case libais::AIS8_1_22_SHAPE_POLYLINE:
            {
                const size_t numReportFields=4;
                const mwSize reportDims[2] = {1, 1};
                const char *reportFieldNames[numReportFields]={"area_shape",
                "angles", "dists_m", "spare"};                      
                mxArray *reportStruct, *angleArray, *distArray;
                libais::Ais8_1_22_Polyline *c =reinterpret_cast<libais::Ais8_1_22_Polyline*>(msg->sub_areas[curSubarea]);
                size_t curItem;
                double *data;
                
                //Create the Matlab structure array.
                reportStruct=mxCreateStructArray(2, reportDims, numReportFields, reportFieldNames);
                //Fill in the component common to all subareas.
                mxSetFieldByNumber(reportStruct,0,0,intMat2MatlabDoubles(&area_shape,1,1));

                //Fill in the components for this subarea
                angleArray=mxCreateDoubleMatrix(c->angles.size(),1,mxREAL);
                data=mxGetDoubles(angleArray);
                for(curItem=0;curItem<c->angles.size();curItem++) {
                    data[curItem]=c->angles[curItem];
                }
                mxSetFieldByNumber(reportStruct,0,1,angleArray);
                
                distArray=mxCreateDoubleMatrix(c->dists_m.size(),1,mxREAL);
                data=mxGetDoubles(distArray);
                for(curItem=0;curItem<c->dists_m.size();curItem++) {
                    data[curItem]=c->dists_m[curItem];
                }
                mxSetFieldByNumber(reportStruct,0,2,distArray);
                
                mxSetFieldByNumber(reportStruct,0,3,mat2MatlabDoubles(&c->spare,1,1));

                mxSetCell(subareaCellArray,curSubarea,reportStruct);
            }
                break;
            case libais::AIS8_1_22_SHAPE_POLYGON:
                        {
                const size_t numReportFields=4;
                const mwSize reportDims[2] = {1, 1};
                const char *reportFieldNames[numReportFields]={"area_shape",
                "angles", "dists_m", "spare"};                      
                mxArray *reportStruct,*angleArray, *distArray;
                libais::Ais8_1_22_Polygon *c =reinterpret_cast<libais::Ais8_1_22_Polygon*>(msg->sub_areas[curSubarea]);
                size_t curItem;
                double *data;
                
                //Create the Matlab structure array.
                reportStruct=mxCreateStructArray(2, reportDims, numReportFields, reportFieldNames);
                //Fill in the component common to all subareas.
                mxSetFieldByNumber(reportStruct,0,0,intMat2MatlabDoubles(&area_shape,1,1));

                //Fill in the components for this subarea
                angleArray=mxCreateDoubleMatrix(c->angles.size(),1,mxREAL);
                data=mxGetDoubles(angleArray);
                for(curItem=0;curItem<c->angles.size();curItem++) {
                    data[curItem]=c->angles[curItem];
                }
                mxSetFieldByNumber(reportStruct,0,1,angleArray);
                
                distArray=mxCreateDoubleMatrix(c->dists_m.size(),1,mxREAL);
                data=mxGetDoubles(distArray);
                for(curItem=0;curItem<c->dists_m.size();curItem++) {
                    data[curItem]=c->dists_m[curItem];
                }
                mxSetFieldByNumber(reportStruct,0,2,distArray);
                
                mxSetFieldByNumber(reportStruct,0,3,mat2MatlabDoubles(&c->spare,1,1));

                mxSetCell(subareaCellArray,curSubarea,reportStruct);
            }
                break;
            case libais::AIS8_1_22_SHAPE_TEXT:
            {
                //save the area shape number.
                const size_t numReportFields=2;
                const mwSize reportDims[2] = {1, 1};
                const char *reportFieldNames[numReportFields]={"area_shape",
                "text"};                      
                mxArray *reportStruct;
                libais::Ais8_1_22_Text *c =reinterpret_cast<libais::Ais8_1_22_Text*>(msg->sub_areas[curSubarea]);

                //Create the Matlab structure array.
                reportStruct=mxCreateStructArray(2, reportDims, numReportFields, reportFieldNames);
                //Fill in the component common to all subareas.
                mxSetFieldByNumber(reportStruct,0,0,intMat2MatlabDoubles(&area_shape,1,1));
                
                //Fill in the components for this subarea
                {
                    const char *charString=c->text.c_str();
                    mxSetFieldByNumber(reportStruct,0,1,mxCreateCharMatrixFromStrings(1,&charString));
                }

                mxSetCell(subareaCellArray,curSubarea,reportStruct);
            }
                break;
            default://An error or an unknown shape.
            {
                //save the area shape number.
                const size_t numReportFields=1;
                const mwSize reportDims[2] = {1, 1};
                const char *reportFieldNames[numReportFields]={"area_shape"};                      
                mxArray *reportStruct;
                
                //Create the Matlab structure array.
                reportStruct=mxCreateStructArray(2, reportDims, numReportFields, reportFieldNames);
                //Fill in the component common to all subareas.
                mxSetFieldByNumber(reportStruct,0,0,intMat2MatlabDoubles(&area_shape,1,1));
                
                mxSetCell(subareaCellArray,curSubarea,reportStruct);
            }
                break;
        }
    }
}

void AIS_8_1_24_ToMatlab(libais::Ais8_1_24 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions) {
//This message format is described in IMO circular 289.
//Message "Extended ship static and voyage-related data (broadcast)"
    const size_t numberOfFields=23;
    const mwSize dims[2] = {1, 1};
    //These are the names of all of the fields that will be added to the
    //Matlab structre array.
    const char *fieldNames[numberOfFields] = {"message_id",
    "repeat_indicator", "mmsi", "spare", "dac", "fi", "link_id",
    "air_draught", "last_port", "next_ports", "solas_status", "ice_class",
    "shaft_power", "vhf", "lloyds_ship_type", "gross_tonnage",
    "laden_ballast", "heavy_oil", "light_oil", "diesel", "bunker_oil",
    "persons", "spare2"};
    mxArray *theStruct;

    //If text descriptions of the fields are desired, then put them into a
    //structure.
    if(fieldDescriptions!=NULL) {
        const char *descriptionStrings[numberOfFields] = {
            "Identifier for Message 8; always 8.",//0
            "Repeat indicator. Used by the repeater to indicate how many times a message has been repeated. 0 = default; 3 = do not repeat any more.",//1
            "Maritime mobile service identity (MMSI) number of source station",//2
            "Not used. Should be zero. Reserved for future use",//3
            "Designated area code (DAC)",//4
            "Function identifier (FI)",//5
            "A source specific running number (1-1023), unique across all binary messages equipped with Message Linkage ID. Used to link additional information to the message by a Text Description message. The Message Linkage ID and the first six digits of the source MMSI uniquely identify the sent message.",//6
            "The vertical distance from the waterline to the highest point on the ship in meters. (1-81.9m)",//7
            "UN LOCODE of last port of call",//8
            "UN LOCODEs of next two ports of call",//9
            "SOLAS Equipment status, see IMO circular 289",//10
            "Ice class:\n0 = not classified\n1 = IACS PC 1\n2 = IACS PC 2\n3 = IACS PC 3\n4 = IACS PC 4\n5 = IACS PC 5\n6 = IACS PC 6 / FSICR IA Super / RS Arc5\n7 = IACS PC 7 / FSICR IA / RS Arc4\n8 = FSICR IB / RS Ice3\n9 = FSICR IC / RS Ice2\n10 = RS Ice1\n11 - 14 (reserved for future use)",//11
             "Shaft horse power (0-162142)",//12
             "VHF working channel number as per  ITU-R M.1084",//13
             "Lloyd's ship type, a string",//14
             "Gross tonnage (a measure of size) (0-262142) as per the International Convention on Tonnage Measurement of Ships, 1969",//15
             "1=Laden, 2=Ballast, 3=Not in use",//16
             "Heavy Fuel Oil:\n1 = no\n2 = yes\n3 = not in use",//17
             "Light Fuel Oil:\n1 = no\n2 = yes\n3 = not in use",//18
             "Diesel:\n1 = no\n2 = yes\n3 = not in use",//19
             "Total amount of bunker oil in tonnes",//20
             "Number of persons currently ob board including crew members (1-8191)",//21
             "Not used. Set to zero."//22
        };
        unsigned int i;

        *fieldDescriptions=mxCreateStructArray(2, dims, numberOfFields, fieldNames);

        for(i=0;i<numberOfFields;i++) {
            mxSetFieldByNumber(*fieldDescriptions,0,static_cast<int>(i),mxCreateCharMatrixFromStrings(1,&descriptionStrings[i]));
        }
    }
    
    //Create the Matlab structure array.
    theStruct=mxCreateStructArray(2, dims, numberOfFields, fieldNames);
    *decodedMessage=theStruct;
    
    //Fill all of the elements of the structure with the set that is common
    //to all AIS8 messages.
    mxSetFieldByNumber(theStruct,0,0,intMat2MatlabDoubles(&msg->message_id,1,1));
    mxSetFieldByNumber(theStruct,0,1,intMat2MatlabDoubles(&msg->repeat_indicator,1,1));
    mxSetFieldByNumber(theStruct,0,2,intMat2MatlabDoubles(&msg->mmsi,1,1));
    mxSetFieldByNumber(theStruct,0,3,intMat2MatlabDoubles(&msg->spare,1,1));
    mxSetFieldByNumber(theStruct,0,4,intMat2MatlabDoubles(&msg->dac,1,1));
    mxSetFieldByNumber(theStruct,0,5,intMat2MatlabDoubles(&msg->fi,1,1));
    
    //Fill in the components unique to this message.
    mxSetFieldByNumber(theStruct,0,6,intMat2MatlabDoubles(&msg->link_id,1,1));
    mxSetFieldByNumber(theStruct,0,7,floatMat2MatlabDoubles(&msg->air_draught,1,1));
    if(msg->last_port.compare("@@@@@")!=0) {
        const char *charString=msg->last_port.c_str();
        mxSetFieldByNumber(theStruct,0,8,mxCreateCharMatrixFromStrings(1,&charString));
    } else {
        mxSetFieldByNumber(theStruct,0,8,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    {
        const char *theStrings[2]={msg->next_ports[0].c_str(),msg->next_ports[1].c_str()};

        if(msg->next_ports[0].compare("@@@@@")!=0&&msg->next_ports[1].compare("@@@@@")!=0) {
            mxSetFieldByNumber(theStruct,0,9,mxCreateCharMatrixFromStrings(2,theStrings));
        } else if (msg->next_ports[0].compare("@@@@@")!=0&&msg->next_ports[1].compare("@@@@@")==0) {
            mxSetFieldByNumber(theStruct,0,9,mxCreateCharMatrixFromStrings(1,theStrings));
        } else if (msg->next_ports[0].compare("@@@@@")==0&&msg->next_ports[1].compare("@@@@@")!=0) {
            mxSetFieldByNumber(theStruct,0,9,mxCreateCharMatrixFromStrings(1,&theStrings[1]));
        } else {
            mxSetFieldByNumber(theStruct,0,9,mxCreateDoubleMatrix(0,0,mxREAL));
        }
    }
    
    mxSetFieldByNumber(theStruct,0,10,intMat2MatlabDoubles(msg->solas_status.data(),26,1));

    if(msg->ice_class!=15) {
        mxSetFieldByNumber(theStruct,0,11,intMat2MatlabDoubles(&msg->ice_class,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,11,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    
    if(msg->shaft_power!=262143) {
        mxSetFieldByNumber(theStruct,0,12,intMat2MatlabDoubles(&msg->shaft_power,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,12,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    if(msg->vhf!=0) {
        mxSetFieldByNumber(theStruct,0,13,intMat2MatlabDoubles(&msg->vhf,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,13,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    if(msg->lloyds_ship_type.compare("@@@@@@@")!=0) {
        const char *charString=msg->lloyds_ship_type.c_str();
        mxSetFieldByNumber(theStruct,0,14,mxCreateCharMatrixFromStrings(1,&charString));
    } else {
        mxSetFieldByNumber(theStruct,0,14,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    if(msg->gross_tonnage!=262143){
        mxSetFieldByNumber(theStruct,0,15,intMat2MatlabDoubles(&msg->gross_tonnage,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,15,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    if(msg->laden_ballast!=0) {
        mxSetFieldByNumber(theStruct,0,16,intMat2MatlabDoubles(&msg->laden_ballast,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,16,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    if(msg->heavy_oil!=0) {
        mxSetFieldByNumber(theStruct,0,17,intMat2MatlabDoubles(&msg->heavy_oil,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,17,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    if(msg->light_oil!=0) {
        mxSetFieldByNumber(theStruct,0,18,intMat2MatlabDoubles(&msg->light_oil,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,18,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    if(msg->diesel!=0) {
        mxSetFieldByNumber(theStruct,0,19,intMat2MatlabDoubles(&msg->diesel,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,19,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    if(msg->bunker_oil!=16383) {
        mxSetFieldByNumber(theStruct,0,20,intMat2MatlabDoubles(&msg->bunker_oil,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,20,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    if(msg->persons!=0) {
        mxSetFieldByNumber(theStruct,0,21,intMat2MatlabDoubles(&msg->persons,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,21,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    mxSetFieldByNumber(theStruct,0,22,intMat2MatlabDoubles(&msg->spare2,1,1));
}

void AIS_8_1_26_ToMatlab(libais::Ais8_1_26 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions) {
//This message format is described in IMO circular 289.
//Message "Environmental"
    const size_t numberOfFields=7;
    const mwSize dims[2] = {1, 1};
    //These are the names of all of the fields that will be added to the
    //Matlab structre array.
    const char *fieldNames[numberOfFields] = {"message_id",
    "repeat_indicator", "mmsi", "spare", "dac", "fi", "reports"};
    mwSize cellArrayDims[2];
    mxArray *theStruct, *reportsCellArray;
    size_t curReport;

    //If text descriptions of the fields are desired, then put them into a
    //structure.
    if(fieldDescriptions!=NULL) {
        const char *descriptionStrings[numberOfFields] = {
            "Identifier for Message 8; always 8.",//0
            "Repeat indicator. Used by the repeater to indicate how many times a message has been repeated. 0 = default; 3 = do not repeat any more.",//1
            "Maritime mobile service identity (MMSI) number of source station",//2
            "Not used. Should be zero. Reserved for future use",//3
            "Designated area code (DAC)",//4
            "Function identifier (FI)",//5
            "A cell array containing sensor reports as structures where the structures fields are as given in IMO circular 289"//6
        };
        unsigned int i;

        *fieldDescriptions=mxCreateStructArray(2, dims, numberOfFields, fieldNames);

        for(i=0;i<numberOfFields;i++) {
            mxSetFieldByNumber(*fieldDescriptions,0,static_cast<int>(i),mxCreateCharMatrixFromStrings(1,&descriptionStrings[i]));
        }
    }
    
    //Create the Matlab structure array.
    theStruct=mxCreateStructArray(2, dims, numberOfFields, fieldNames);
    *decodedMessage=theStruct;
    
    //Fill all of the elements of the structure with the set that is common
    //to all AIS8 messages.
    mxSetFieldByNumber(theStruct,0,0,intMat2MatlabDoubles(&msg->message_id,1,1));
    mxSetFieldByNumber(theStruct,0,1,intMat2MatlabDoubles(&msg->repeat_indicator,1,1));
    mxSetFieldByNumber(theStruct,0,2,intMat2MatlabDoubles(&msg->mmsi,1,1));
    mxSetFieldByNumber(theStruct,0,3,intMat2MatlabDoubles(&msg->spare,1,1));
    mxSetFieldByNumber(theStruct,0,4,intMat2MatlabDoubles(&msg->dac,1,1));
    mxSetFieldByNumber(theStruct,0,5,intMat2MatlabDoubles(&msg->fi,1,1));
    
    //Fill in the components unique to this message.
    cellArrayDims[0]=msg->reports.size();
    cellArrayDims[1]=1;
    //Make space to hold all of the reports in a cell array.
    reportsCellArray=mxCreateCellArray(2, cellArrayDims);
    mxSetFieldByNumber(theStruct,0,6,reportsCellArray);

    //Fill in all of the reports.
    for(curReport=0;curReport<msg->reports.size();curReport++) {
        switch (msg->reports[curReport]->report_type) { 
            case libais::AIS8_1_26_SENSOR_LOCATION:
            {
                const size_t numReportFields=11;
                const mwSize reportDims[2] = {1, 1};
                const char *reportFieldNames[numReportFields]={"report_type",
                "utc_day", "utc_hr", "utc_min", "site_id", "lon", "lat", "z",
                "owner", "timeout", "spare"};                      
                libais::Ais8_1_26_Location *rpt = reinterpret_cast<libais::Ais8_1_26_Location *>(msg->reports[curReport]);
                mxArray *reportStruct;
                
                //Create the Matlab structure array.
                reportStruct=mxCreateStructArray(2, reportDims, numReportFields, reportFieldNames);
                //Fill in the components common to all sensor reports.
                mxSetFieldByNumber(reportStruct,0,0,intMat2MatlabDoubles(&rpt->report_type,1,1));
                if(rpt->utc_day!=0) {
                    mxSetFieldByNumber(reportStruct,0,1,intMat2MatlabDoubles(&rpt->utc_day,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,1,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                if(rpt->utc_hr<=23) {
                    mxSetFieldByNumber(reportStruct,0,2,intMat2MatlabDoubles(&rpt->utc_hr,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,2,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                if(rpt->utc_min<=59) {
                    mxSetFieldByNumber(reportStruct,0,3,intMat2MatlabDoubles(&rpt->utc_min,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,3,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                mxSetFieldByNumber(reportStruct,0,4,intMat2MatlabDoubles(&rpt->site_id,1,1));
                
                //Fill in the components for this sensor report type.
                if(fabs(rpt->position.lng_deg)<=180) {
                    mxSetFieldByNumber(reportStruct,0,5,doubleMat2Matlab(&rpt->position.lng_deg,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,5,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                if(fabs(rpt->position.lat_deg)<=90) {
                    mxSetFieldByNumber(reportStruct,0,6,doubleMat2Matlab(&rpt->position.lat_deg,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,6,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                if(rpt->z<=200) {
                    mxSetFieldByNumber(reportStruct,0,7,floatMat2MatlabDoubles(&rpt->z,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,7,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                if(rpt->owner!=0) {
                    mxSetFieldByNumber(reportStruct,0,8,intMat2MatlabDoubles(&rpt->owner,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,8,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                mxSetFieldByNumber(reportStruct,0,9,intMat2MatlabDoubles(&rpt->timeout,1,1));
                mxSetFieldByNumber(reportStruct,0,10,intMat2MatlabDoubles(&rpt->spare,1,1));

                mxSetCell(reportsCellArray,curReport,reportStruct);
            }
                break;
            case libais::AIS8_1_26_SENSOR_STATION:
            {
                const size_t numReportFields=7;
                const mwSize reportDims[2] = {1, 1};
                const char *reportFieldNames[numReportFields]={"report_type",
                "utc_day", "utc_hr", "utc_min", "site_id", "name", "spare"};                      
                libais::Ais8_1_26_Station *rpt = reinterpret_cast<libais::Ais8_1_26_Station *>(msg->reports[curReport]);
                mxArray *reportStruct;
                
                //Create the Matlab structure array.
                reportStruct=mxCreateStructArray(2, reportDims, numReportFields, reportFieldNames);
                //Fill in the components common to all sensor reports.
                mxSetFieldByNumber(reportStruct,0,0,intMat2MatlabDoubles(&rpt->report_type,1,1));
                if(rpt->utc_day!=0) {
                    mxSetFieldByNumber(reportStruct,0,1,intMat2MatlabDoubles(&rpt->utc_day,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,1,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                if(rpt->utc_hr<=23) {
                    mxSetFieldByNumber(reportStruct,0,2,intMat2MatlabDoubles(&rpt->utc_hr,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,2,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                if(rpt->utc_min<=59) {
                    mxSetFieldByNumber(reportStruct,0,3,intMat2MatlabDoubles(&rpt->utc_min,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,3,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                mxSetFieldByNumber(reportStruct,0,4,intMat2MatlabDoubles(&rpt->site_id,1,1));
                
                //Fill in the components for this sensor report type.
                if(rpt->name.compare("@@@@@@@@@@@@@@")!=0) {
                    const char *charString=rpt->name.c_str();
                    mxSetFieldByNumber(reportStruct,0,5,mxCreateCharMatrixFromStrings(1,&charString));
                } else {
                    mxSetFieldByNumber(reportStruct,0,5,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                mxSetFieldByNumber(reportStruct,0,6,intMat2MatlabDoubles(&rpt->spare,1,1));

                mxSetCell(reportsCellArray,curReport,reportStruct);
            }
                break;
            case libais::AIS8_1_26_SENSOR_WIND:
            {
                const size_t numReportFields=18;
                const mwSize reportDims[2] = {1, 1};
                const char *reportFieldNames[numReportFields]={"report_type",
                "utc_day", "utc_hr", "utc_min", "site_id", "wind_speed",
                "wind_gust", "wind_dir", "wind_gust_dir", "sensor_type",
                "wind_forecast", "wind_gust_forecast", "wind_dir_forecast",
                "utc_day_forecast", "utc_hour_forecast",
                "utc_min_forecast", "duration", "spare"};                      
                libais::Ais8_1_26_Wind *rpt = reinterpret_cast<libais::Ais8_1_26_Wind *>(msg->reports[curReport]);
                mxArray *reportStruct;
                
                //Create the Matlab structure array.
                reportStruct=mxCreateStructArray(2, reportDims, numReportFields, reportFieldNames);

                //Fill in the components common to all sensor reports.
                mxSetFieldByNumber(reportStruct,0,0,intMat2MatlabDoubles(&rpt->report_type,1,1));
                if(rpt->utc_day!=0) {
                    mxSetFieldByNumber(reportStruct,0,1,intMat2MatlabDoubles(&rpt->utc_day,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,1,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                if(rpt->utc_hr<=23) {
                    mxSetFieldByNumber(reportStruct,0,2,intMat2MatlabDoubles(&rpt->utc_hr,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,2,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                if(rpt->utc_min<=59) {
                    mxSetFieldByNumber(reportStruct,0,3,intMat2MatlabDoubles(&rpt->utc_min,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,3,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                mxSetFieldByNumber(reportStruct,0,4,intMat2MatlabDoubles(&rpt->site_id,1,1));
                
                //Fill in the components for this sensor report type.
                if(rpt->wind_speed<=121) {
                    mxSetFieldByNumber(reportStruct,0,5,intMat2MatlabDoubles(&rpt->wind_speed,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,5,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                if(rpt->wind_gust<=121) {
                    mxSetFieldByNumber(reportStruct,0,6,intMat2MatlabDoubles(&rpt->wind_gust,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,6,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                if(rpt->wind_dir<=359) {
                    mxSetFieldByNumber(reportStruct,0,7,intMat2MatlabDoubles(&rpt->wind_dir,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,7,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                if(rpt->wind_gust_dir<=359) {
                    mxSetFieldByNumber(reportStruct,0,8,intMat2MatlabDoubles(&rpt->wind_gust_dir,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,8,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                if(rpt->sensor_type!=0) {
                    mxSetFieldByNumber(reportStruct,0,9,intMat2MatlabDoubles(&rpt->sensor_type,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,9,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                if(rpt->wind_forecast<=121) {
                    mxSetFieldByNumber(reportStruct,0,10,intMat2MatlabDoubles(&rpt->wind_forecast,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,10,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                if(rpt->wind_gust_forecast<=121) {
                    mxSetFieldByNumber(reportStruct,0,11,intMat2MatlabDoubles(&rpt->wind_gust_forecast,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,11,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                if(rpt->wind_dir_forecast<=359) {
                    mxSetFieldByNumber(reportStruct,0,12,intMat2MatlabDoubles(&rpt->wind_dir_forecast,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,12,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                if(rpt->utc_day_forecast!=0) {
                    mxSetFieldByNumber(reportStruct,0,13,intMat2MatlabDoubles(&rpt->utc_day_forecast,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,13,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                if(rpt->utc_hour_forecast<=23) {
                    mxSetFieldByNumber(reportStruct,0,14,intMat2MatlabDoubles(&rpt->utc_hour_forecast,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,14,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                if(rpt->utc_min_forecast<=59) {
                    mxSetFieldByNumber(reportStruct,0,15,intMat2MatlabDoubles(&rpt->utc_min_forecast,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,15,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                mxSetFieldByNumber(reportStruct,0,16,intMat2MatlabDoubles(&rpt->duration,1,1));
                mxSetFieldByNumber(reportStruct,0,17,intMat2MatlabDoubles(&rpt->spare,1,1));

                mxSetCell(reportsCellArray,curReport,reportStruct);
            }
                break;
            case libais::AIS8_1_26_SENSOR_WATER_LEVEL:
            {
                const size_t numReportFields=17;
                const mwSize reportDims[2] = {1, 1};
                const char *reportFieldNames[numReportFields]={"report_type",
                "utc_day", "utc_hr", "utc_min", "site_id", "type", "level",
                "trend", "vdatum", "sensor_type", "forecast_type",
                "level_forecast", "utc_day_forecast", "utc_hour_forecast",
                "utc_min_forecast", "duration", "spare"};                      
                libais::Ais8_1_26_WaterLevel *rpt = reinterpret_cast<libais::Ais8_1_26_WaterLevel *>(msg->reports[curReport]);
                mxArray *reportStruct;
                
                //Create the Matlab structure array.
                reportStruct=mxCreateStructArray(2, reportDims, numReportFields, reportFieldNames);
                //Fill in the components common to all sensor reports.
                mxSetFieldByNumber(reportStruct,0,0,intMat2MatlabDoubles(&rpt->report_type,1,1));
                if(rpt->utc_day!=0) {
                    mxSetFieldByNumber(reportStruct,0,1,intMat2MatlabDoubles(&rpt->utc_day,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,1,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                if(rpt->utc_hr<=23) {
                    mxSetFieldByNumber(reportStruct,0,2,intMat2MatlabDoubles(&rpt->utc_hr,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,2,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                if(rpt->utc_min<=59) {
                    mxSetFieldByNumber(reportStruct,0,3,intMat2MatlabDoubles(&rpt->utc_min,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,3,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                mxSetFieldByNumber(reportStruct,0,4,intMat2MatlabDoubles(&rpt->site_id,1,1));
                
                //Fill in the components for this sensor report type.
                mxSetFieldByNumber(reportStruct,0,5,intMat2MatlabDoubles(&rpt->type,1,1));
                if(rpt->level>-327.68) {
                    mxSetFieldByNumber(reportStruct,0,6,floatMat2MatlabDoubles(&rpt->level,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,6,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                if(rpt->trend!=3) {
                    mxSetFieldByNumber(reportStruct,0,7,intMat2MatlabDoubles(&rpt->trend,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,7,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                if(rpt->vdatum!=14) {
                    mxSetFieldByNumber(reportStruct,0,8,intMat2MatlabDoubles(&rpt->vdatum,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,8,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                mxSetFieldByNumber(reportStruct,0,9,intMat2MatlabDoubles(&rpt->sensor_type,1,1));
                mxSetFieldByNumber(reportStruct,0,10,intMat2MatlabDoubles(&rpt->forecast_type,1,1));
                if(rpt->level_forecast>-327.68) {
                    mxSetFieldByNumber(reportStruct,0,11,floatMat2MatlabDoubles(&rpt->level_forecast,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,11,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                if(rpt->utc_day_forecast!=0) {
                    mxSetFieldByNumber(reportStruct,0,12,intMat2MatlabDoubles(&rpt->utc_day_forecast,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,12,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                if(rpt->utc_hour_forecast<=23) {
                    mxSetFieldByNumber(reportStruct,0,13,intMat2MatlabDoubles(&rpt->utc_hour_forecast,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,13,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                if(rpt->utc_min_forecast<=59) {
                    mxSetFieldByNumber(reportStruct,0,14,intMat2MatlabDoubles(&rpt->utc_min_forecast,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,14,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                mxSetFieldByNumber(reportStruct,0,15,intMat2MatlabDoubles(&rpt->duration,1,1));
                mxSetFieldByNumber(reportStruct,0,16,intMat2MatlabDoubles(&rpt->spare,1,1));

                mxSetCell(reportsCellArray,curReport,reportStruct);
            }
                break;
            case libais::AIS8_1_26_SENSOR_CURR_2D:
            {
                const size_t numReportFields=8;
                const mwSize reportDims[2] = {1, 1};
                const char *reportFieldNames[numReportFields]={"report_type",
                "utc_day", "utc_hr", "utc_min", "site_id", "type", "spare",
                "currents"};                      
                libais::Ais8_1_26_Curr2D *rpt = reinterpret_cast<libais::Ais8_1_26_Curr2D *>(msg->reports[curReport]);
                const mwSize currentsDims[2] = {3,1};
                const size_t numberOfCurrentFields=3;
                const char *currentFieldNames[numberOfCurrentFields]={"speed",
                "dir", "depth"};  
                mxArray *reportStruct, *currentsArray;
                size_t curCurrent;
                
                //Create the Matlab structure array.
                reportStruct=mxCreateStructArray(2, reportDims, numReportFields, reportFieldNames);
                //Fill in the components common to all sensor reports.
                mxSetFieldByNumber(reportStruct,0,0,intMat2MatlabDoubles(&rpt->report_type,1,1));
                if(rpt->utc_day!=0) {
                    mxSetFieldByNumber(reportStruct,0,1,intMat2MatlabDoubles(&rpt->utc_day,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,1,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                if(rpt->utc_hr<=23) {
                    mxSetFieldByNumber(reportStruct,0,2,intMat2MatlabDoubles(&rpt->utc_hr,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,2,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                if(rpt->utc_min<=59) {
                    mxSetFieldByNumber(reportStruct,0,3,intMat2MatlabDoubles(&rpt->utc_min,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,3,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                mxSetFieldByNumber(reportStruct,0,4,intMat2MatlabDoubles(&rpt->site_id,1,1));
                
                //Fill in the components for this sensor report type.
                mxSetFieldByNumber(reportStruct,0,5,intMat2MatlabDoubles(&rpt->type,1,1));
                mxSetFieldByNumber(reportStruct,0,6,intMat2MatlabDoubles(&rpt->spare,1,1));

                currentsArray=mxCreateStructArray(2, currentsDims, numberOfCurrentFields, currentFieldNames);
                mxSetFieldByNumber(reportStruct,0,7,currentsArray);
                for(curCurrent=0;curCurrent<3;curCurrent++) {
                    if(rpt->currents[curCurrent].speed<=246) {
                        mxSetFieldByNumber(currentsArray,curCurrent,0,floatMat2MatlabDoubles(&rpt->currents[curCurrent].speed,1,1));
                    } else {
                        mxSetFieldByNumber(currentsArray,curCurrent,0,mxCreateDoubleMatrix(0,0,mxREAL));
                    }
                    if(rpt->currents[curCurrent].dir<=359) {
                        mxSetFieldByNumber(currentsArray,curCurrent,1,intMat2MatlabDoubles(&rpt->currents[curCurrent].dir,1,1));
                    } else {
                        mxSetFieldByNumber(currentsArray,curCurrent,1,mxCreateDoubleMatrix(0,0,mxREAL));
                    }
                    if(rpt->currents[curCurrent].depth<=362) {
                        mxSetFieldByNumber(currentsArray,curCurrent,2,intMat2MatlabDoubles(&rpt->currents[curCurrent].depth,1,1));
                    } else {
                        mxSetFieldByNumber(currentsArray,curCurrent,2,mxCreateDoubleMatrix(0,0,mxREAL));
                    }
                }
                
                mxSetCell(reportsCellArray,curReport,reportStruct);
            }
                break;
            case libais::AIS8_1_26_SENSOR_CURR_3D:
            {
                const size_t numReportFields=8;
                const mwSize reportDims[2] = {1, 1};
                const char *reportFieldNames[numReportFields]={"report_type",
                "utc_day", "utc_hr", "utc_min", "site_id", "type", "spare",
                "currents"};                      
                libais::Ais8_1_26_Curr3D *rpt = reinterpret_cast<libais::Ais8_1_26_Curr3D *>(msg->reports[curReport]);
                const mwSize currentsDims[2] = {2,1};
                const size_t numberOfCurrentFields=4;
                const char *currentFieldNames[numberOfCurrentFields]={"north",
                "east", "up", "depth"};  
                mxArray *reportStruct, *currentsArray;
                size_t curCurrent;
                
                //Create the Matlab structure array.
                reportStruct=mxCreateStructArray(2, reportDims, numReportFields, reportFieldNames);
                //Fill in the components common to all sensor reports.
                mxSetFieldByNumber(reportStruct,0,0,intMat2MatlabDoubles(&rpt->report_type,1,1));
                if(rpt->utc_day!=0) {
                    mxSetFieldByNumber(reportStruct,0,1,intMat2MatlabDoubles(&rpt->utc_day,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,1,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                if(rpt->utc_hr<=23) {
                    mxSetFieldByNumber(reportStruct,0,2,intMat2MatlabDoubles(&rpt->utc_hr,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,2,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                if(rpt->utc_min<=59) {
                    mxSetFieldByNumber(reportStruct,0,3,intMat2MatlabDoubles(&rpt->utc_min,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,3,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                mxSetFieldByNumber(reportStruct,0,4,intMat2MatlabDoubles(&rpt->site_id,1,1));
                
                //Fill in the components for this sensor report type.
                mxSetFieldByNumber(reportStruct,0,5,intMat2MatlabDoubles(&rpt->type,1,1));
                mxSetFieldByNumber(reportStruct,0,6,intMat2MatlabDoubles(&rpt->spare,1,1));
                
                currentsArray=mxCreateStructArray(2, currentsDims, numberOfCurrentFields, currentFieldNames);
                mxSetFieldByNumber(reportStruct,0,7,currentsArray);
                for(curCurrent=0;curCurrent<2;curCurrent++) {
                    if(rpt->currents[curCurrent].north<=24.6) {
                        mxSetFieldByNumber(currentsArray,curCurrent,0,floatMat2MatlabDoubles(&rpt->currents[curCurrent].north,1,1));
                    } else {
                        mxSetFieldByNumber(currentsArray,curCurrent,0,mxCreateDoubleMatrix(0,0,mxREAL));
                    }
                    if(rpt->currents[curCurrent].east<=24.6) {
                        mxSetFieldByNumber(currentsArray,curCurrent,1,floatMat2MatlabDoubles(&rpt->currents[curCurrent].east,1,1));
                    } else {
                        mxSetFieldByNumber(currentsArray,curCurrent,1,mxCreateDoubleMatrix(0,0,mxREAL));
                    }
                    if(rpt->currents[curCurrent].up<=24.6) {
                        mxSetFieldByNumber(currentsArray,curCurrent,2,floatMat2MatlabDoubles(&rpt->currents[curCurrent].up,1,1));
                    } else {
                        mxSetFieldByNumber(currentsArray,curCurrent,2,mxCreateDoubleMatrix(0,0,mxREAL));
                    }
                    mxSetFieldByNumber(currentsArray,curCurrent,3,intMat2MatlabDoubles(&rpt->currents[curCurrent].depth,1,1));
                }
                
                mxSetCell(reportsCellArray,curReport,reportStruct);
            }
                break;
            case libais::AIS8_1_26_SENSOR_HORZ_FLOW:
            {
                const size_t numReportFields=7;
                const mwSize reportDims[2] = {1, 1};
                const char *reportFieldNames[numReportFields]={"report_type",
                "utc_day", "utc_hr", "utc_min", "site_id", "spare",
                "currents"};                      
                libais::Ais8_1_26_HorzFlow *rpt = reinterpret_cast<libais::Ais8_1_26_HorzFlow *>(msg->reports[curReport]);
                const mwSize flowDims[2] = {2,1};
                const size_t numberOfFlowFields=5;
                const char *flowFieldNames[numberOfFlowFields]={"bearing",
                "dist", "speed", "dir", "level"};  
                mxArray *reportStruct, *flowArray;
                size_t curFlow;
                
                //Create the Matlab structure array.
                reportStruct=mxCreateStructArray(2, reportDims, numReportFields, reportFieldNames);
                //Fill in the components common to all sensor reports.
                mxSetFieldByNumber(reportStruct,0,0,intMat2MatlabDoubles(&rpt->report_type,1,1));
                if(rpt->utc_day!=0) {
                    mxSetFieldByNumber(reportStruct,0,1,intMat2MatlabDoubles(&rpt->utc_day,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,1,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                if(rpt->utc_hr<=23) {
                    mxSetFieldByNumber(reportStruct,0,2,intMat2MatlabDoubles(&rpt->utc_hr,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,2,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                if(rpt->utc_min<=59) {
                    mxSetFieldByNumber(reportStruct,0,3,intMat2MatlabDoubles(&rpt->utc_min,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,3,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                mxSetFieldByNumber(reportStruct,0,4,intMat2MatlabDoubles(&rpt->site_id,1,1));
                
                //Fill in the components for this sensor report type.
                mxSetFieldByNumber(reportStruct,0,5,intMat2MatlabDoubles(&rpt->spare,1,1));
                
                flowArray=mxCreateStructArray(2, flowDims, numberOfFlowFields, flowFieldNames);
                mxSetFieldByNumber(reportStruct,0,6,flowArray);
                for(curFlow=0;curFlow<2;curFlow++) {
                    if(rpt->currents[curFlow].bearing!=360) {
                        mxSetFieldByNumber(flowArray,curFlow,0,intMat2MatlabDoubles(&rpt->currents[curFlow].bearing,1,1));
                    } else {
                        mxSetFieldByNumber(flowArray,curFlow,0,mxCreateDoubleMatrix(0,0,mxREAL));
                    }
                    if(rpt->currents[curFlow].dist<=121) {
                        mxSetFieldByNumber(flowArray,curFlow,1,intMat2MatlabDoubles(&rpt->currents[curFlow].dist,1,1));
                    } else {
                        mxSetFieldByNumber(flowArray,curFlow,1,mxCreateDoubleMatrix(0,0,mxREAL));
                    }
                    if(rpt->currents[curFlow].speed<=24.6) {
                        mxSetFieldByNumber(flowArray,curFlow,2,floatMat2MatlabDoubles(&rpt->currents[curFlow].speed,1,1));
                    } else {
                        mxSetFieldByNumber(flowArray,curFlow,2,mxCreateDoubleMatrix(0,0,mxREAL));
                    }
                    if(rpt->currents[curFlow].dir<=359) {
                        mxSetFieldByNumber(flowArray,curFlow,3,intMat2MatlabDoubles(&rpt->currents[curFlow].dir,1,1));
                    } else {
                        mxSetFieldByNumber(flowArray,curFlow,3,mxCreateDoubleMatrix(0,0,mxREAL));
                    }
                    if(rpt->currents[curFlow].level<=361) {
                        mxSetFieldByNumber(flowArray,curFlow,4,intMat2MatlabDoubles(&rpt->currents[curFlow].level,1,1));
                    } else {
                        mxSetFieldByNumber(flowArray,curFlow,4,mxCreateDoubleMatrix(0,0,mxREAL));
                    }
                }
                
                mxSetCell(reportsCellArray,curReport,reportStruct);
            }
                break;
            case libais::AIS8_1_26_SENSOR_SEA_STATE:
            {
                const size_t numReportFields=18;
                const mwSize reportDims[2] = {1, 1};
                const char *reportFieldNames[numReportFields]={"report_type",
                "utc_day", "utc_hr", "utc_min", "site_id", "swell_height",
                "swell_period", "swell_dir", "sea_state", "swell_sensor_type",
                "water_temp", "water_temp_depth", "water_sensor_type",
                "wave_height", "wave_period", "wave_dir", "wave_sensor_type",
                "salinity"};                      
                libais::Ais8_1_26_SeaState *rpt = reinterpret_cast<libais::Ais8_1_26_SeaState *>(msg->reports[curReport]);
                mxArray *reportStruct;
                
                //Create the Matlab structure array.
                reportStruct=mxCreateStructArray(2, reportDims, numReportFields, reportFieldNames);
                //Fill in the components common to all sensor reports.
                mxSetFieldByNumber(reportStruct,0,0,intMat2MatlabDoubles(&rpt->report_type,1,1));
                if(rpt->utc_day!=0) {
                    mxSetFieldByNumber(reportStruct,0,1,intMat2MatlabDoubles(&rpt->utc_day,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,1,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                if(rpt->utc_hr<=23) {
                    mxSetFieldByNumber(reportStruct,0,2,intMat2MatlabDoubles(&rpt->utc_hr,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,2,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                if(rpt->utc_min<=59) {
                    mxSetFieldByNumber(reportStruct,0,3,intMat2MatlabDoubles(&rpt->utc_min,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,3,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                mxSetFieldByNumber(reportStruct,0,4,intMat2MatlabDoubles(&rpt->site_id,1,1));
                
                //Fill in the components for this sensor report type.
                if(rpt->swell_height<=24.6) {
                    mxSetFieldByNumber(reportStruct,0,5,floatMat2MatlabDoubles(&rpt->swell_height,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,5,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                if(rpt->swell_period<=60) {
                    mxSetFieldByNumber(reportStruct,0,6,intMat2MatlabDoubles(&rpt->swell_period,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,6,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                if(rpt->swell_dir<=359) {
                    mxSetFieldByNumber(reportStruct,0,7,intMat2MatlabDoubles(&rpt->swell_dir,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,7,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                if(rpt->sea_state<=12) {
                    mxSetFieldByNumber(reportStruct,0,8,intMat2MatlabDoubles(&rpt->sea_state,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,8,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                if(rpt->swell_sensor_type!=0) {
                    mxSetFieldByNumber(reportStruct,0,9,intMat2MatlabDoubles(&rpt->swell_sensor_type,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,9,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                if(rpt->water_temp>=-10&&rpt->water_temp<=50) {
                    mxSetFieldByNumber(reportStruct,0,10,floatMat2MatlabDoubles(&rpt->water_temp,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,10,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                if(rpt->water_temp_depth<=12.1) {
                    mxSetFieldByNumber(reportStruct,0,11,floatMat2MatlabDoubles(&rpt->water_temp_depth,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,11,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                if(rpt->water_sensor_type!=0) {
                    mxSetFieldByNumber(reportStruct,0,12,intMat2MatlabDoubles(&rpt->water_sensor_type,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,12,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                if(rpt->wave_height<=24.6) {
                    mxSetFieldByNumber(reportStruct,0,13,floatMat2MatlabDoubles(&rpt->wave_height,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,13,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                if(rpt->wave_period<=60) {
                    mxSetFieldByNumber(reportStruct,0,14,intMat2MatlabDoubles(&rpt->wave_period,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,14,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                if(rpt->wave_dir<=359) {
                    mxSetFieldByNumber(reportStruct,0,15,intMat2MatlabDoubles(&rpt->wave_dir,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,15,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                if(rpt->wave_sensor_type!=0) {
                    mxSetFieldByNumber(reportStruct,0,16,intMat2MatlabDoubles(&rpt->wave_sensor_type,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,16,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                if(rpt->salinity<=50.1) {
                    mxSetFieldByNumber(reportStruct,0,17,floatMat2MatlabDoubles(&rpt->salinity,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,17,mxCreateDoubleMatrix(0,0,mxREAL));
                }

                mxSetCell(reportsCellArray,curReport,reportStruct);
            }
                break;
            case libais::AIS8_1_26_SENSOR_SALINITY:
            {
                const size_t numReportFields=12;
                const mwSize reportDims[2] = {1, 1};
                const char *reportFieldNames[numReportFields]={"report_type",
                "utc_day", "utc_hr", "utc_min", "site_id", "water_temp",
                "conductivity", "pressure", "salinity", "salinity_type",
                "sensor_type", "spare"};                      
                libais::Ais8_1_26_Salinity *rpt = reinterpret_cast<libais::Ais8_1_26_Salinity *>(msg->reports[curReport]);
                mxArray *reportStruct;
                
                //Create the Matlab structure array.
                reportStruct=mxCreateStructArray(2, reportDims, numReportFields, reportFieldNames);
                //Fill in the components common to all sensor reports.
                mxSetFieldByNumber(reportStruct,0,0,intMat2MatlabDoubles(&rpt->report_type,1,1));
                if(rpt->utc_day!=0) {
                    mxSetFieldByNumber(reportStruct,0,1,intMat2MatlabDoubles(&rpt->utc_day,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,1,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                if(rpt->utc_hr<=23) {
                    mxSetFieldByNumber(reportStruct,0,2,intMat2MatlabDoubles(&rpt->utc_hr,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,2,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                if(rpt->utc_min<=59) {
                    mxSetFieldByNumber(reportStruct,0,3,intMat2MatlabDoubles(&rpt->utc_min,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,3,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                mxSetFieldByNumber(reportStruct,0,4,intMat2MatlabDoubles(&rpt->site_id,1,1));
                
                //Fill in the components for this sensor report type.
                if(rpt->water_temp<=50&&rpt->water_temp>=-10) {
                    mxSetFieldByNumber(reportStruct,0,5,floatMat2MatlabDoubles(&rpt->water_temp,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,5,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                if(rpt->conductivity<=7) {
                    mxSetFieldByNumber(reportStruct,0,6,floatMat2MatlabDoubles(&rpt->conductivity,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,6,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                if(rpt->pressure<=6000.1) {
                    mxSetFieldByNumber(reportStruct,0,7,floatMat2MatlabDoubles(&rpt->pressure,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,7,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                if(rpt->salinity<=50.1) {
                    mxSetFieldByNumber(reportStruct,0,8,floatMat2MatlabDoubles(&rpt->salinity,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,8,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                mxSetFieldByNumber(reportStruct,0,9,intMat2MatlabDoubles(&rpt->salinity_type,1,1));
                if(rpt->sensor_type!=0) {
                    mxSetFieldByNumber(reportStruct,0,10,intMat2MatlabDoubles(&rpt->sensor_type,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,10,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                mxSetFieldByNumber(reportStruct,0,11,intMat2MatlabDoubles(rpt->spare,2,1));

                mxSetCell(reportsCellArray,curReport,reportStruct);
            }
                break;
            case libais::AIS8_1_26_SENSOR_WX:
            {
                const size_t numReportFields=16;
                const mwSize reportDims[2] = {1, 1};
                const char *reportFieldNames[numReportFields]={"report_type",
                "utc_day", "utc_hr", "utc_min", "site_id", "air_temp",
                "air_temp_sensor_type", "precip", "horz_vis", "dew_point",
                "dew_point_type", "air_pressure", "air_pressure_trend",
                "air_pressor_type", "salinity", "spare"};                      
                libais::Ais8_1_26_Wx *rpt = reinterpret_cast<libais::Ais8_1_26_Wx *>(msg->reports[curReport]);
                mxArray *reportStruct;
                
                //Create the Matlab structure array.
                reportStruct=mxCreateStructArray(2, reportDims, numReportFields, reportFieldNames);
                //Fill in the components common to all sensor reports.
                mxSetFieldByNumber(reportStruct,0,0,intMat2MatlabDoubles(&rpt->report_type,1,1));
                if(rpt->utc_day!=0) {
                    mxSetFieldByNumber(reportStruct,0,1,intMat2MatlabDoubles(&rpt->utc_day,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,1,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                if(rpt->utc_hr<=23) {
                    mxSetFieldByNumber(reportStruct,0,2,intMat2MatlabDoubles(&rpt->utc_hr,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,2,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                if(rpt->utc_min<=59) {
                    mxSetFieldByNumber(reportStruct,0,3,intMat2MatlabDoubles(&rpt->utc_min,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,3,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                mxSetFieldByNumber(reportStruct,0,4,intMat2MatlabDoubles(&rpt->site_id,1,1));
                
                //Fill in the components for this sensor report type.
                if(rpt->air_temp>=-60&&rpt->air_temp<=60) {
                    mxSetFieldByNumber(reportStruct,0,5,floatMat2MatlabDoubles(&rpt->air_temp,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,5,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                if(rpt->air_temp_sensor_type!=0) {
                    mxSetFieldByNumber(reportStruct,0,6,intMat2MatlabDoubles(&rpt->air_temp_sensor_type,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,6,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                mxSetFieldByNumber(reportStruct,0,7,intMat2MatlabDoubles(&rpt->precip,1,1));
                if(rpt->horz_vis<=24.1) {
                    mxSetFieldByNumber(reportStruct,0,8,floatMat2MatlabDoubles(&rpt->horz_vis,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,8,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                if(rpt->dew_point<=50&&rpt->dew_point>=-20) {
                    mxSetFieldByNumber(reportStruct,0,9,floatMat2MatlabDoubles(&rpt->dew_point,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,9,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                if(rpt->dew_point_type!=0) {
                    mxSetFieldByNumber(reportStruct,0,10,intMat2MatlabDoubles(&rpt->dew_point_type,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,10,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                if(rpt->air_pressure<=1201) {
                    mxSetFieldByNumber(reportStruct,0,11,floatMat2MatlabDoubles(&rpt->air_pressure,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,11,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                if(rpt->air_pressure_trend!=3) {
                    mxSetFieldByNumber(reportStruct,0,12,intMat2MatlabDoubles(&rpt->air_pressure_trend,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,12,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                if(rpt->air_pressor_type!=0) {
                    mxSetFieldByNumber(reportStruct,0,13,intMat2MatlabDoubles(&rpt->air_pressor_type,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,13,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                if(rpt->salinity<=50.1) {
                    mxSetFieldByNumber(reportStruct,0,14,floatMat2MatlabDoubles(&rpt->salinity,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,14,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                mxSetFieldByNumber(reportStruct,0,15,intMat2MatlabDoubles(&rpt->spare,1,1));

                mxSetCell(reportsCellArray,curReport,reportStruct);
            }
                break;
            case libais::AIS8_1_26_SENSOR_AIR_DRAUGHT:
            {
                const size_t numReportFields=13;
                const mwSize reportDims[2] = {1, 1};
                const char *reportFieldNames[numReportFields]={"report_type",
                "utc_day", "utc_hr", "utc_min", "site_id", "draught", "gap",
                "forecast_gap", "trend", "utc_day_forecast", "utc_hour_forecast",
                "utc_min_forecast", "spare"};                      
                libais::Ais8_1_26_AirDraught *rpt = reinterpret_cast<libais::Ais8_1_26_AirDraught *>(msg->reports[curReport]);
                mxArray *reportStruct;
                
                //Create the Matlab structure array.
                reportStruct=mxCreateStructArray(2, reportDims, numReportFields, reportFieldNames);

                //Fill in the components common to all sensor reports.
                mxSetFieldByNumber(reportStruct,0,0,intMat2MatlabDoubles(&rpt->report_type,1,1));
                if(rpt->utc_day!=0) {
                    mxSetFieldByNumber(reportStruct,0,1,intMat2MatlabDoubles(&rpt->utc_day,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,1,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                if(rpt->utc_hr<=23) {
                    mxSetFieldByNumber(reportStruct,0,2,intMat2MatlabDoubles(&rpt->utc_hr,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,2,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                if(rpt->utc_min<=59) {
                    mxSetFieldByNumber(reportStruct,0,3,intMat2MatlabDoubles(&rpt->utc_min,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,3,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                mxSetFieldByNumber(reportStruct,0,4,intMat2MatlabDoubles(&rpt->site_id,1,1));
                
                //Fill in the components for this sensor report type.
                if(rpt->draught!=0) {
                    mxSetFieldByNumber(reportStruct,0,5,floatMat2MatlabDoubles(&rpt->draught,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,6,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                if(rpt->gap!=0) {
                    mxSetFieldByNumber(reportStruct,0,6,floatMat2MatlabDoubles(&rpt->gap,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,6,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                if(rpt->forecast_gap!=3) {
                    mxSetFieldByNumber(reportStruct,0,7,floatMat2MatlabDoubles(&rpt->forecast_gap,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,7,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                if(rpt->trend!=3) {
                    mxSetFieldByNumber(reportStruct,0,8,intMat2MatlabDoubles(&rpt->trend,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,8,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                if(rpt->utc_day_forecast!=0) {
                    mxSetFieldByNumber(reportStruct,0,9,intMat2MatlabDoubles(&rpt->utc_day_forecast,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,9,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                if(rpt->utc_hour_forecast<=23) {
                    mxSetFieldByNumber(reportStruct,0,10,intMat2MatlabDoubles(&rpt->utc_hour_forecast,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,10,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                if(rpt->utc_min_forecast<=59) {
                    mxSetFieldByNumber(reportStruct,0,11,intMat2MatlabDoubles(&rpt->utc_min_forecast,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,11,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                mxSetFieldByNumber(reportStruct,0,12,intMat2MatlabDoubles(&rpt->spare,1,1));

                mxSetCell(reportsCellArray,curReport,reportStruct);
            }
                break;
            default:
            {//Error, reserved or unknown report type. Just return the
             //components common to all report types.
                const size_t numReportFields=5;
                const mwSize reportDims[2] = {1, 1};
                const char *reportFieldNames[numReportFields]={"report_type",
                "utc_day", "utc_hr", "utc_min", "site_id"};                      
                libais::Ais8_1_26_SensorReport *rpt = reinterpret_cast<libais::Ais8_1_26_SensorReport*>(msg->reports[curReport]);
                mxArray *reportStruct;
                
                //Create the Matlab structure array.
                reportStruct=mxCreateStructArray(2, reportDims, numReportFields, reportFieldNames);
                //Fill in the components common to all sensor reports.
                mxSetFieldByNumber(reportStruct,0,0,intMat2MatlabDoubles(&rpt->report_type,1,1));
                if(rpt->utc_day!=0) {
                    mxSetFieldByNumber(reportStruct,0,1,intMat2MatlabDoubles(&rpt->utc_day,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,1,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                if(rpt->utc_hr<=23) {
                    mxSetFieldByNumber(reportStruct,0,2,intMat2MatlabDoubles(&rpt->utc_hr,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,2,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                if(rpt->utc_min<=59) {
                    mxSetFieldByNumber(reportStruct,0,3,intMat2MatlabDoubles(&rpt->utc_min,1,1));
                } else {
                    mxSetFieldByNumber(reportStruct,0,3,mxCreateDoubleMatrix(0,0,mxREAL));
                }
                mxSetFieldByNumber(reportStruct,0,4,intMat2MatlabDoubles(&rpt->site_id,1,1));
                
                mxSetCell(reportsCellArray,curReport,reportStruct);
            }
                break;
        }
    }
}

void AIS_8_1_27_ToMatlab(libais::Ais8_1_27 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions) {
//This message format is described in IMO circular 289.
//Message "Route Information (Broadcast)"
    const size_t numberOfFields=15;
    const mwSize dims[2] = {1, 1};
    //These are the names of all of the fields that will be added to the
    //Matlab structre array.
    const char *fieldNames[numberOfFields] = {"message_id",
    "repeat_indicator", "mmsi", "spare", "dac", "fi", "link_id",
    "sender_type", "route_type", "utc_month", "utc_day", "utc_hour",
    "utc_min", "duration", "waypoints"};
    const size_t numberOfWaypointFields=2;
    mwSize waypointDims[2];
    const char *waypointFieldNames[numberOfWaypointFields] = {"lon", "lat"};
    mxArray *theStruct, *waypointArray;
    size_t curWaypoint;
       
    //If text descriptions of the fields are desired, then put them into a
    //structure.
    if(fieldDescriptions!=NULL) {
        const char *descriptionStrings[numberOfFields] = {
            "Identifier for Message 8; always 8.",//0
            "Repeat indicator. Used by the repeater to indicate how many times a message has been repeated. 0 = default; 3 = do not repeat any more.",//1
            "Maritime mobile service identity (MMSI) number of source station",//2
            "Not used. Should be zero. Reserved for future use",//3
            "Designated area code (DAC)",//4
            "Function identifier (FI)",//5
            "A source specific running number (1-1023), unique across all binary messages equipped with Message Linkage ID. Used to link additional information to the message by a Text Description message. The Message Linkage ID and the first six digits of the source MMSI uniquely identify the sent message.",//6
            "Sender Classification\n0 = ship = default\n1 = authority\n2 - 7 (reserved for future use)",//7
            "Route type:\n1 = mandatory route\n2 = recommended route\n3 = alternative route\n4 = recommended route through ice\n5 = ship route plan\n6 - 30 (reserved for future use)\n31 = cancellation (cancel route as identified by Message Linkage ID)",//8
            "UTC Month (1-12)",//9
            "UTC Day (1-31)",//10
            "UTC Hour (0-23)",//11
            "UTC minute (0-59)",//12
            "Minutes until end of validity of the route. Measured from start time of Route Information. 0 = cancel route",//13
            "An array containing WGS-84 longitude East and latitude North coordinates in degrees of the waypoints"//14
        };
        unsigned int i;

        *fieldDescriptions=mxCreateStructArray(2, dims, numberOfFields, fieldNames);

        for(i=0;i<numberOfFields;i++) {
            mxSetFieldByNumber(*fieldDescriptions,0,static_cast<int>(i),mxCreateCharMatrixFromStrings(1,&descriptionStrings[i]));
        }
    }
    
    //Create the Matlab structure array.
    theStruct=mxCreateStructArray(2, dims, numberOfFields, fieldNames);
    *decodedMessage=theStruct;
    
    //Fill all of the elements of the structure with the set that is common
    //to all AIS8 messages.
    mxSetFieldByNumber(theStruct,0,0,intMat2MatlabDoubles(&msg->message_id,1,1));
    mxSetFieldByNumber(theStruct,0,1,intMat2MatlabDoubles(&msg->repeat_indicator,1,1));
    mxSetFieldByNumber(theStruct,0,2,intMat2MatlabDoubles(&msg->mmsi,1,1));
    mxSetFieldByNumber(theStruct,0,3,intMat2MatlabDoubles(&msg->spare,1,1));
    mxSetFieldByNumber(theStruct,0,4,intMat2MatlabDoubles(&msg->dac,1,1));
    mxSetFieldByNumber(theStruct,0,5,intMat2MatlabDoubles(&msg->fi,1,1));
    
    //Fill in the components unique to this message.
    if(msg->link_id!=0) {
        mxSetFieldByNumber(theStruct,0,6,intMat2MatlabDoubles(&msg->link_id,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,6,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    mxSetFieldByNumber(theStruct,0,7,intMat2MatlabDoubles(&msg->sender_type,1,1));
    if(msg->route_type!=0) {
        mxSetFieldByNumber(theStruct,0,8,intMat2MatlabDoubles(&msg->route_type,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,8,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    if(msg->utc_month!=0) {
        mxSetFieldByNumber(theStruct,0,9,intMat2MatlabDoubles(&msg->utc_month,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,9,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    if(msg->utc_day!=0) {
        mxSetFieldByNumber(theStruct,0,10,intMat2MatlabDoubles(&msg->utc_day,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,10,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    if(msg->utc_hour<=23) {
        mxSetFieldByNumber(theStruct,0,11,intMat2MatlabDoubles(&msg->utc_hour,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,11,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    if(msg->utc_min<=59) {
        mxSetFieldByNumber(theStruct,0,12,intMat2MatlabDoubles(&msg->utc_min,1,1));
    } else{
        mxSetFieldByNumber(theStruct,0,12,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    mxSetFieldByNumber(theStruct,0,13,intMat2MatlabDoubles(&msg->duration,1,1));

    waypointDims[0]=msg->waypoints.size();
    waypointDims[1]=1;
    waypointArray=mxCreateStructArray(2, waypointDims, numberOfWaypointFields, waypointFieldNames);
    mxSetFieldByNumber(theStruct,0,14,waypointArray);

    for(curWaypoint= 0; curWaypoint < msg->waypoints.size(); curWaypoint++) {
        if(fabs(msg->waypoints[curWaypoint].lng_deg)<=180) {
            mxSetFieldByNumber(waypointArray,curWaypoint,0,doubleMat2Matlab(&msg->waypoints[curWaypoint].lng_deg,1,1));
        } else {
            mxSetFieldByNumber(waypointArray,curWaypoint,0,mxCreateDoubleMatrix(0,0,mxREAL));
        }
        if(msg->waypoints[curWaypoint].lat_deg<=90) {
            mxSetFieldByNumber(waypointArray,curWaypoint,1,doubleMat2Matlab(&msg->waypoints[curWaypoint].lat_deg,1,1));
        } else {
            mxSetFieldByNumber(waypointArray,curWaypoint,1,mxCreateDoubleMatrix(0,0,mxREAL));
        }
    }
}

void AIS_8_1_29_ToMatlab(libais::Ais8_1_29 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions) {
//This message format is described in IMO circular 289.
//Message "Text description (broadcast)"
    const size_t numberOfFields=9;
    const mwSize dims[2] = {1, 1};
    //These are the names of all of the fields that will be added to the
    //Matlab structre array.
    const char *fieldNames[numberOfFields] = {"message_id",
    "repeat_indicator", "mmsi", "spare", "dac", "fi", "link_id", "text",
    "spare2"};
    mxArray *theStruct;

    //If text descriptions of the fields are desired, then put them into a
    //structure.
    if(fieldDescriptions!=NULL) {
        const char *descriptionStrings[numberOfFields] = {
            "Identifier for Message 8; always 8.",//0
            "Repeat indicator. Used by the repeater to indicate how many times a message has been repeated. 0 = default; 3 = do not repeat any more.",//1
            "Maritime mobile service identity (MMSI) number of source station",//2
            "Not used. Should be zero. Reserved for future use",//3
            "Designated area code (DAC)",//4
            "Function identifier (FI)",//5
            "Message Linkage ID. Used to link the Text Description message with a main message. The Connection ID and source MMSI MID uniquely identifies the main message. (1-1023)",
            "Text String"
        };
        unsigned int i;

        *fieldDescriptions=mxCreateStructArray(2, dims, numberOfFields, fieldNames);

        for(i=0;i<numberOfFields;i++) {
            mxSetFieldByNumber(*fieldDescriptions,0,static_cast<int>(i),mxCreateCharMatrixFromStrings(1,&descriptionStrings[i]));
        }
    }

    //Create the Matlab structure array.
    theStruct=mxCreateStructArray(2, dims, numberOfFields, fieldNames);
    *decodedMessage=theStruct;
    
    //Fill all of the elements of the structure with the set that is common
    //to all AIS8 messages.
    mxSetFieldByNumber(theStruct,0,0,intMat2MatlabDoubles(&msg->message_id,1,1));
    mxSetFieldByNumber(theStruct,0,1,intMat2MatlabDoubles(&msg->repeat_indicator,1,1));
    mxSetFieldByNumber(theStruct,0,2,intMat2MatlabDoubles(&msg->mmsi,1,1));
    mxSetFieldByNumber(theStruct,0,3,intMat2MatlabDoubles(&msg->spare,1,1));
    mxSetFieldByNumber(theStruct,0,4,intMat2MatlabDoubles(&msg->dac,1,1));
    mxSetFieldByNumber(theStruct,0,5,intMat2MatlabDoubles(&msg->fi,1,1));
    
    //Fill in the components unique to this message.
    if(msg->link_id!=0) {
        mxSetFieldByNumber(theStruct,0,6,intMat2MatlabDoubles(&msg->link_id,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,6,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    {
        const char *charString=msg->text.c_str();
        mxSetFieldByNumber(theStruct,0,7,mxCreateCharMatrixFromStrings(1,&charString));
    }
    mxSetFieldByNumber(theStruct,0,8,intMat2MatlabDoubles(&msg->spare2,1,1));
}

void AIS_8_1_31_ToMatlab(libais::Ais8_1_31 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions) {
//This message format is described in IMO circular 289.
//Message "Meteorological and Hydrographic data"
    const size_t numberOfFields=44;
    const mwSize dims[2] = {1, 1};
    //These are the names of all of the fields that will be added to the
    //Matlab structre array.
    const char *fieldNames[numberOfFields] = {"message_id",
    "repeat_indicator", "mmsi", "spare", "dac", "fi", "lon", "lat",
    "position_accuracy", "utc_day", "utc_hour", "utc_min", "wind_ave",
    "wind_gust", "wind_dir", "wind_gust_dir", "air_temp", "rel_humid",
    "dew_point", "air_pres", "air_pres_trend", "horz_vis", "water_level",
    "water_level_trend", "surf_cur_speed", "surf_cur_dir", "cur_speed_2",
    "cur_dir_2", "cur_depth_2", "cur_speed_3", "cur_dir_3", "cur_depth_3",
    "wave_height", "wave_period", "wave_dir", "swell_height", "swell_period",
    "swell_dir", "sea_state", "water_temp", "precip_type", "salinity",
    "ice", "spare2"};
    mxArray *theStruct;

    //If text descriptions of the fields are desired, then put them into a
    //structure.
    if(fieldDescriptions!=NULL) {
        const char *descriptionStrings[numberOfFields] = {
            "Identifier for Message 8; always 8.",//0
            "Repeat indicator. Used by the repeater to indicate how many times a message has been repeated. 0 = default; 3 = do not repeat any more.",//1
            "Maritime mobile service identity (MMSI) number of source station",//2
            "Not used. Should be zero. Reserved for future use",//3
            "Designated area code (DAC)",//4
            "Function identifier (FI)",//5
            "Longitude East in degrees (WGS-84)",//6
            "Latitude North in degrees (WGS-84)",//7
            "Position accuracy. 1=high (<=10m), 0=low (>10m), 0=default",//8
            "UTC day transmitted (1-31)",//9
            "UTC hour transmitted (0-23)",//10
            "UTC minute transmitted (0-59)",//11
            "Average wind speed for the last 10 minutes in knots (0-126)",//12
            "Maximum wind gust speed during the last 10 minutes in knots (0-126)",//13
            "Wind direction in degrees East of North (0-359)",//14
            "Wind gust direction in degrees East of North (0-359)",//15
            "Dry buld air temperature in degrees Celsius (-60 to 60)",//16
            "Relative humidity in percent (0-100)",//17
            "Dew point in degrees Celsius (-20 to 50)",//18
            "Air pressure in hectoPascals (799-1201)",//19
            "Air pressure tendency 0 = steady, 1 = decreasing, 2 = increasing",//20
            "Horizontal visibility in nautical miles (0-12.6)",//21 NOTE: The first bit being set makes this an invalid value, and indicates that the equipment was saturated. We are not using that; an empty matrix is just returned.
            "Water level, including tide, deviation from local datum in meters (-10 to 30)",//22
            "Water level trend 0 = steady, 1 = decreasing, 2 = increasing",//23
            "Surface current speed including tide in knots (0-25.1)",//24
            "Surface current direction in degrees East of North (0-359)",//25
            "Current 2 speed in knots (0-25.1)",//26
            "Current 2 speed direction in degrees East of North (0-359)",//27
            "Depth of current 2 in meters (0-30)",//28
            "Current 3 speed in knots (0-25.1)",//29
            "Current 3 speed direction in degrees East of North (0-359)",//30
            "Depth of current 3 in meters (0-30)",//31
            "Significant wave height in meters (0-25)",//32
            "Wave period in seconds (0-60)",//33
            "Wave direction in degrees East of North (0-359)",//34
            "Swell height in meters (0-25)",//35
            "Swell period in seconds (0-60)",//36
            "Swell direction in degrees East of North (0-359)",//37
            "Sea state according to the Beauford scale (0-12)",//38
            "Water temperature in degrees Celsius (-10 to 50)",//39
            "Precipitation type according to the WMO\n0 = reserved\n1 = rain\n2 = thunderstorm\n3 = freezing rain\n4 = mixed/ice\n5 = snow\n6 = reserved",//40
            "Salinity in parts per thousand (0-50.1)",//41
            "A boolean value indicating the presence of ice",//42
            "Not used. Should be zero."//43
        };
        unsigned int i;

        *fieldDescriptions=mxCreateStructArray(2, dims, numberOfFields, fieldNames);

        for(i=0;i<numberOfFields;i++) {
            mxSetFieldByNumber(*fieldDescriptions,0,static_cast<int>(i),mxCreateCharMatrixFromStrings(1,&descriptionStrings[i]));
        }
    }
    
    //Create the Matlab structure array.
    theStruct=mxCreateStructArray(2, dims, numberOfFields, fieldNames);
    *decodedMessage=theStruct;
    
    //Fill all of the elements of the structure with the set that is common
    //to all AIS8 messages.
    mxSetFieldByNumber(theStruct,0,0,intMat2MatlabDoubles(&msg->message_id,1,1));
    mxSetFieldByNumber(theStruct,0,1,intMat2MatlabDoubles(&msg->repeat_indicator,1,1));
    mxSetFieldByNumber(theStruct,0,2,intMat2MatlabDoubles(&msg->mmsi,1,1));
    mxSetFieldByNumber(theStruct,0,3,intMat2MatlabDoubles(&msg->spare,1,1));
    mxSetFieldByNumber(theStruct,0,4,intMat2MatlabDoubles(&msg->dac,1,1));
    mxSetFieldByNumber(theStruct,0,5,intMat2MatlabDoubles(&msg->fi,1,1));
    
    //Fill in the components unique to this message.
   if(fabs(msg->position.lng_deg)<=180) {
        mxSetFieldByNumber(theStruct,0,6,doubleMat2Matlab(&msg->position.lng_deg,1,1));
    } else {//No longitude
        mxSetFieldByNumber(theStruct,0,6,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    if(fabs(msg->position.lat_deg)<=90) {
        mxSetFieldByNumber(theStruct,0,7,doubleMat2Matlab(&msg->position.lat_deg,1,1));
    } else {//No latitude
        mxSetFieldByNumber(theStruct,0,7,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    mxSetFieldByNumber(theStruct,0,8,intMat2MatlabDoubles(&msg->position_accuracy,1,1));
    
    if(msg->utc_day!=0) {
        mxSetFieldByNumber(theStruct,0,9,intMat2MatlabDoubles(&msg->utc_day,1,1));
    } else {//Day of departure not available
        mxSetFieldByNumber(theStruct,0,9,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    if(msg->utc_hour<24) {
        mxSetFieldByNumber(theStruct,0,10,intMat2MatlabDoubles(&msg->utc_hour,1,1));
    } else {//Hour of departure not available
        mxSetFieldByNumber(theStruct,0,10,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    if(msg->utc_min<60) {
        mxSetFieldByNumber(theStruct,0,11,intMat2MatlabDoubles(&msg->utc_min,1,1));
    } else {//Minute of departure not available
        mxSetFieldByNumber(theStruct,0,11,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    if(msg->wind_ave<=126) {
        mxSetFieldByNumber(theStruct,0,12,intMat2MatlabDoubles(&msg->wind_ave,1,1));
    } else {//No wind average
        mxSetFieldByNumber(theStruct,0,12,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    if(msg->wind_gust<=126) {
        mxSetFieldByNumber(theStruct,0,13,intMat2MatlabDoubles(&msg->wind_gust,1,1));
    } else {//No wind gust speed
        mxSetFieldByNumber(theStruct,0,13,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    if(msg->wind_dir<=359) {
        mxSetFieldByNumber(theStruct,0,14,intMat2MatlabDoubles(&msg->wind_dir,1,1));
    } else {//No wind direction
        mxSetFieldByNumber(theStruct,0,14,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    if(msg->wind_gust_dir<=359) {
        mxSetFieldByNumber(theStruct,0,15,intMat2MatlabDoubles(&msg->wind_gust_dir,1,1));
    } else {//No wind gust direction
        mxSetFieldByNumber(theStruct,0,15,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    if(msg->air_temp<=60&&msg->air_temp>=-60) {
        mxSetFieldByNumber(theStruct,0,16,floatMat2MatlabDoubles(&msg->air_temp,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,16,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    if(msg->rel_humid<=100) {
        mxSetFieldByNumber(theStruct,0,17,intMat2MatlabDoubles(&msg->rel_humid,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,17,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    if(msg->dew_point<=50&&msg->dew_point>=-20) {
        mxSetFieldByNumber(theStruct,0,18,floatMat2MatlabDoubles(&msg->dew_point,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,18,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    if(msg->air_pres<=1201) {
        mxSetFieldByNumber(theStruct,0,19,floatMat2MatlabDoubles(&msg->air_pres,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,19,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    if(msg->air_pres_trend<=2){
        mxSetFieldByNumber(theStruct,0,20,intMat2MatlabDoubles(&msg->air_pres_trend,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,20,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    if(msg->horz_vis<=12.6) {
        mxSetFieldByNumber(theStruct,0,21,floatMat2MatlabDoubles(&msg->horz_vis,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,21,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    if(msg->water_level<=30&&msg->water_level>=-10) {
        mxSetFieldByNumber(theStruct,0,22,floatMat2MatlabDoubles(&msg->water_level,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,22,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    if(msg->water_level_trend<=2) {
        mxSetFieldByNumber(theStruct,0,23,intMat2MatlabDoubles(&msg->water_level_trend,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,23,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    if(msg->surf_cur_speed<=25.1) {
        mxSetFieldByNumber(theStruct,0,24,floatMat2MatlabDoubles(&msg->surf_cur_speed,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,24,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    if(msg->surf_cur_dir<=359) {
        mxSetFieldByNumber(theStruct,0,25,intMat2MatlabDoubles(&msg->surf_cur_dir,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,25,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    if(msg->cur_speed_2<=25.1) {
        mxSetFieldByNumber(theStruct,0,26,floatMat2MatlabDoubles(&msg->cur_speed_2,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,26,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    if(msg->cur_dir_2<=359) {
        mxSetFieldByNumber(theStruct,0,27,intMat2MatlabDoubles(&msg->cur_dir_2,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,27,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    if(msg->cur_depth_2<=30) {
        mxSetFieldByNumber(theStruct,0,28,intMat2MatlabDoubles(&msg->cur_depth_2,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,28,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    if(msg->cur_speed_3<=25) {
        mxSetFieldByNumber(theStruct,0,29,floatMat2MatlabDoubles(&msg->cur_speed_3,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,29,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    if(msg->cur_dir_3<=359) {
        mxSetFieldByNumber(theStruct,0,30,intMat2MatlabDoubles(&msg->cur_dir_3,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,30,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    if(msg->cur_depth_3<=30) {
        mxSetFieldByNumber(theStruct,0,31,intMat2MatlabDoubles(&msg->cur_depth_3,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,31,mxCreateDoubleMatrix(0,0,mxREAL));
    }
    if(msg->wave_height<=25.1) {
        mxSetFieldByNumber(theStruct,0,32,floatMat2MatlabDoubles(&msg->wave_height,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,32,mxCreateDoubleMatrix(0,0,mxREAL)); 
    }
    if(msg->wave_period<=60) {
        mxSetFieldByNumber(theStruct,0,33,intMat2MatlabDoubles(&msg->wave_period,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,33,mxCreateDoubleMatrix(0,0,mxREAL)); 
    }
    if(msg->wave_dir<=359) {
        mxSetFieldByNumber(theStruct,0,34,intMat2MatlabDoubles(&msg->wave_dir,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,34,mxCreateDoubleMatrix(0,0,mxREAL)); 
    }
    if(msg->swell_height<=25.1) {
        mxSetFieldByNumber(theStruct,0,35,floatMat2MatlabDoubles(&msg->swell_height,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,35,mxCreateDoubleMatrix(0,0,mxREAL)); 
    }
    if(msg->swell_period<=60) {
        mxSetFieldByNumber(theStruct,0,36,intMat2MatlabDoubles(&msg->swell_period,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,36,mxCreateDoubleMatrix(0,0,mxREAL)); 
    }
    if(msg->swell_dir<=359) {
        mxSetFieldByNumber(theStruct,0,37,intMat2MatlabDoubles(&msg->swell_dir,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,37,mxCreateDoubleMatrix(0,0,mxREAL)); 
    }
    if(msg->sea_state<=12) {
        mxSetFieldByNumber(theStruct,0,38,intMat2MatlabDoubles(&msg->sea_state,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,38,mxCreateDoubleMatrix(0,0,mxREAL)); 
    }
    if(msg->water_temp<=50&&msg->water_temp>=-10) {
        mxSetFieldByNumber(theStruct,0,39,floatMat2MatlabDoubles(&msg->water_temp,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,39,mxCreateDoubleMatrix(0,0,mxREAL)); 
    }
    if(msg->precip_type<7) {
        mxSetFieldByNumber(theStruct,0,40,intMat2MatlabDoubles(&msg->precip_type,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,40,mxCreateDoubleMatrix(0,0,mxREAL)); 
    }
    if(msg->salinity<=50) {
        mxSetFieldByNumber(theStruct,0,41,floatMat2MatlabDoubles(&msg->salinity,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,41,mxCreateDoubleMatrix(0,0,mxREAL)); 
    }
    if(msg->ice==0||msg->ice==1) {
        mxSetFieldByNumber(theStruct,0,42,intMat2MatlabDoubles(&msg->ice,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,42,mxCreateDoubleMatrix(0,0,mxREAL)); 
    }

    mxSetFieldByNumber(theStruct,0,43,intMat2MatlabDoubles(&msg->spare2,1,1));
}

void AIS_8_200_10_ToMatlab(libais::Ais8_200_10 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions) {
//This message is defined in ECE-TRANS-SC3-176e.pdf
//Message "Inland ship static and voyage related data"
    const size_t numberOfFields=17;
    const mwSize dims[2] = {1, 1};
    //These are the names of all of the fields that will be added to the
    //Matlab structre array.
    const char *fieldNames[numberOfFields] = {"message_id",
    "repeat_indicator", "mmsi", "spare", "dac", "fi", "eu_id", "length",
    "beam", "ship_type", "haz_cargo", "draught", "loaded", "speed_qual",
    "course_qual", "heading_qual", "spare2"};
    mxArray *theStruct;

    //If text descriptions of the fields are desired, then put them into a
    //structure.
    if(fieldDescriptions!=NULL) {
        const char *descriptionStrings[numberOfFields] = {
            "Identifier for Message 8; always 8.",//0
            "Repeat indicator. Used by the repeater to indicate how many times a message has been repeated. 0 = default; 3 = do not repeat any more.",//1
            "Maritime mobile service identity (MMSI) number of source station",//2
            "Not used. Should be zero. Reserved for future use",//3
            "Designated area code (DAC)",//4
            "Function identifier (FI)",//5
            "European Identifier string",//6
            "Length of ship in meters (0=default)",//7
            "Beam of ship in meters (0=default)",//8
            "Numeric ERI Classification for ship or combination type",//9
            "Hazardous Cargo, number of blue cones/ lights0-3, 4=B-flag",//10
            "Draught in meters",//11
            "1=loaded, 2=unloaded, 3 should not be used",//12
            "Quality of speed information1=high, 0=low/GNSS",//13
            "Quality of course information1=high, 0=low/GNSS",//14
            "Quality of heading information1=high, 0=low",//15
            "Not used, should be set to zero. Reserved for future use."//16
        };
        unsigned int i;

        *fieldDescriptions=mxCreateStructArray(2, dims, numberOfFields, fieldNames);

        for(i=0;i<numberOfFields;i++) {
            mxSetFieldByNumber(*fieldDescriptions,0,static_cast<int>(i),mxCreateCharMatrixFromStrings(1,&descriptionStrings[i]));
        }
    }
    
    //Create the Matlab structure array.
    theStruct=mxCreateStructArray(2, dims, numberOfFields, fieldNames);
    *decodedMessage=theStruct;
    
    //Fill all of the elements of the structure with the set that is common
    //to all AIS8 messages.
    mxSetFieldByNumber(theStruct,0,0,intMat2MatlabDoubles(&msg->message_id,1,1));
    mxSetFieldByNumber(theStruct,0,1,intMat2MatlabDoubles(&msg->repeat_indicator,1,1));
    mxSetFieldByNumber(theStruct,0,2,intMat2MatlabDoubles(&msg->mmsi,1,1));
    mxSetFieldByNumber(theStruct,0,3,intMat2MatlabDoubles(&msg->spare,1,1));
    mxSetFieldByNumber(theStruct,0,4,intMat2MatlabDoubles(&msg->dac,1,1));
    mxSetFieldByNumber(theStruct,0,5,intMat2MatlabDoubles(&msg->fi,1,1));
    
    //Fill in the components unique to this message.
    {
        const char *charString=msg->eu_id.c_str();
        mxSetFieldByNumber(theStruct,0,6,mxCreateCharMatrixFromStrings(1,&charString));
    }
    mxSetFieldByNumber(theStruct,0,7,floatMat2MatlabDoubles(&msg->length,1,1));
    mxSetFieldByNumber(theStruct,0,8,floatMat2MatlabDoubles(&msg->beam,1,1));
    mxSetFieldByNumber(theStruct,0,9,intMat2MatlabDoubles(&msg->ship_type,1,1));
    if(msg->haz_cargo<=4) {
        mxSetFieldByNumber(theStruct,0,10,intMat2MatlabDoubles(&msg->haz_cargo,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,10,mxCreateDoubleMatrix(0,0,mxREAL)); 
    }
    if(msg->draught!=0) {
        mxSetFieldByNumber(theStruct,0,11,floatMat2MatlabDoubles(&msg->draught,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,11,mxCreateDoubleMatrix(0,0,mxREAL)); 
    }
    if(msg->loaded!=0) {
        mxSetFieldByNumber(theStruct,0,12,intMat2MatlabDoubles(&msg->loaded,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,12,mxCreateDoubleMatrix(0,0,mxREAL)); 
    }
    mxSetFieldByNumber(theStruct,0,13,intMat2MatlabDoubles(&msg->speed_qual,1,1));
    mxSetFieldByNumber(theStruct,0,14,intMat2MatlabDoubles(&msg->course_qual,1,1));
    mxSetFieldByNumber(theStruct,0,15,intMat2MatlabDoubles(&msg->heading_qual,1,1));
    mxSetFieldByNumber(theStruct,0,16,intMat2MatlabDoubles(&msg->spare2,1,1));
}

void AIS_8_200_23_ToMatlab(libais::Ais8_200_23 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions) {
//This message is defined in ECE-TRANS-SC3-176e.pdf
//Message "EMMA warning"
    const size_t numberOfFields=26;
    const mwSize dims[2] = {1, 1};
    //These are the names of all of the fields that will be added to the
    //Matlab structre array.
    const char *fieldNames[numberOfFields] = {"message_id",
    "repeat_indicator", "mmsi", "spare", "dac", "fi", "utc_year_start",
    "utc_month_start", "utc_day_start", "utc_hour_start", "utc_min_start"
    "utc_year_end", "utc_month_end", "utc_day_end", "utc_hour_end",
    "utc_min_end", "lon1", "lat1", "lon2", "lat2", "type", "min", "max",
    "classification", "wind_dir", "spare2"};
    mxArray *theStruct;
    
    //If text descriptions of the fields are desired, then put them into a
    //structure.
    if(fieldDescriptions!=NULL) {
        const char *descriptionStrings[numberOfFields] = {
            "Identifier for Message 8; always 8.",//0
            "Repeat indicator. Used by the repeater to indicate how many times a message has been repeated. 0 = default; 3 = do not repeat any more.",//1
            "Maritime mobile service identity (MMSI) number of source station",//2
            "Not used. Should be zero. Reserved for future use",//3
            "Designated area code (DAC)",//4
            "Function identifier (FI)",//5
            "UTC year of the starting validity period given as an offset from the year 2000 (1-255)",//6
            "UTC month of the starting validity period (1-12)",//7
            "UTC day of the starting validity period (1-31)",//8
            "UTC hour of the starting validity period (0-23)",//9
            "UTC minute of the starting validity period (0-59)",//10
            "UTC year of the ending validity period given as an offset from the year 2000 (1-255)",//11
            "UTC month of the ending validity period (1-12)",//12
            "UTC day of the ending validity period (1-31)",//13
            "UTC hour of the ending validity period (0-23)",//14
            "UTC minute of the ending validity period (0-59)",//15
            "Longitude East in degrees for the beginning of the fairway section (WGS-84)",//16
            "Latitude North in degrees for the beginning of the fairway section (WGS-84)",//17
            "Longitude East in degrees for the end of the fairway section (WGS-84)",//18
            "Latitude North in degrees for the end of the fairway section (WGS-84)",//19
            "Type of weather warning. Possible values are:\n1 Wind\n2 Rain\n3 Snow and Ice\n4 Thunderstorm\n5 Fog\n6 Low Temperature\n7 High Temperature\n8 Flood\n9 Fire in the Forests",//20
            "Min Value",//21
            "Max Value",//22
            "Classification of warning\n1  slight\n2  medium\n3 strong/heavy",//23
            "Wind direction such that\n1 North\n2 North East\n3 East\n4 South East\n5 South\n6 South West\n7 West\n8 North West",//24
            "Not used, should be set to zero. Reserved for future use."//25
        };
        unsigned int i;

        *fieldDescriptions=mxCreateStructArray(2, dims, numberOfFields, fieldNames);

        for(i=0;i<numberOfFields;i++) {
            mxSetFieldByNumber(*fieldDescriptions,0,static_cast<int>(i),mxCreateCharMatrixFromStrings(1,&descriptionStrings[i]));
        }
    }
    
    //Create the Matlab structure array.
    theStruct=mxCreateStructArray(2, dims, numberOfFields, fieldNames);
    *decodedMessage=theStruct;
    
    //Fill all of the elements of the structure with the set that is common
    //to all AIS8 messages.
    mxSetFieldByNumber(theStruct,0,0,intMat2MatlabDoubles(&msg->message_id,1,1));
    mxSetFieldByNumber(theStruct,0,1,intMat2MatlabDoubles(&msg->repeat_indicator,1,1));
    mxSetFieldByNumber(theStruct,0,2,intMat2MatlabDoubles(&msg->mmsi,1,1));
    mxSetFieldByNumber(theStruct,0,3,intMat2MatlabDoubles(&msg->spare,1,1));
    mxSetFieldByNumber(theStruct,0,4,intMat2MatlabDoubles(&msg->dac,1,1));
    mxSetFieldByNumber(theStruct,0,5,intMat2MatlabDoubles(&msg->fi,1,1));
    
    //Fill in the components unique to this message.
    if(msg->utc_year_start!=0) {
        mxSetFieldByNumber(theStruct,0,6,intMat2MatlabDoubles(&msg->utc_year_start,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,6,mxCreateDoubleMatrix(0,0,mxREAL)); 
    }
    if(msg->utc_month_start!=0) {
        mxSetFieldByNumber(theStruct,0,7,intMat2MatlabDoubles(&msg->utc_month_start,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,7,mxCreateDoubleMatrix(0,0,mxREAL)); 
    }
    if(msg->utc_day_start!=0) {
        mxSetFieldByNumber(theStruct,0,8,intMat2MatlabDoubles(&msg->utc_day_start,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,8,mxCreateDoubleMatrix(0,0,mxREAL)); 
    }
    if(msg->utc_hour_start<=23) {
        mxSetFieldByNumber(theStruct,0,9,intMat2MatlabDoubles(&msg->utc_hour_start,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,9,mxCreateDoubleMatrix(0,0,mxREAL)); 
    }
    if(msg->utc_min_start<=59) {
        mxSetFieldByNumber(theStruct,0,10,intMat2MatlabDoubles(&msg->utc_min_start,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,10,mxCreateDoubleMatrix(0,0,mxREAL)); 
    }
    if(msg->utc_year_end!=0) {
        mxSetFieldByNumber(theStruct,0,11,intMat2MatlabDoubles(&msg->utc_year_end,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,11,mxCreateDoubleMatrix(0,0,mxREAL)); 
    }
    if(msg->utc_month_end!=0) {
        mxSetFieldByNumber(theStruct,0,12,intMat2MatlabDoubles(&msg->utc_month_end,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,12,mxCreateDoubleMatrix(0,0,mxREAL)); 
    }
    if(msg->utc_day_end!=0) {
        mxSetFieldByNumber(theStruct,0,13,intMat2MatlabDoubles(&msg->utc_day_end,1,1)); 
    } else {
        mxSetFieldByNumber(theStruct,0,13,mxCreateDoubleMatrix(0,0,mxREAL)); 
    }
    if(msg->utc_hour_end<=23) {
        mxSetFieldByNumber(theStruct,0,14,intMat2MatlabDoubles(&msg->utc_hour_end,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,14,mxCreateDoubleMatrix(0,0,mxREAL)); 
    }
    if(msg->utc_min_end<=59) {
        mxSetFieldByNumber(theStruct,0,15,intMat2MatlabDoubles(&msg->utc_min_end,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,15,mxCreateDoubleMatrix(0,0,mxREAL)); 
    }
    mxSetFieldByNumber(theStruct,0,16,doubleMat2Matlab(&msg->position1.lng_deg,1,1));
    mxSetFieldByNumber(theStruct,0,17,doubleMat2Matlab(&msg->position1.lat_deg,1,1));
    mxSetFieldByNumber(theStruct,0,18,doubleMat2Matlab(&msg->position2.lng_deg,1,1));
    mxSetFieldByNumber(theStruct,0,19,doubleMat2Matlab(&msg->position2.lat_deg,1,1));
    if(msg->type!=0){
        mxSetFieldByNumber(theStruct,0,20,intMat2MatlabDoubles(&msg->type,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,20,mxCreateDoubleMatrix(0,0,mxREAL)); 
    }
    mxSetFieldByNumber(theStruct,0,21,intMat2MatlabDoubles(&msg->min,1,1));
    mxSetFieldByNumber(theStruct,0,22,intMat2MatlabDoubles(&msg->max,1,1));
    if(msg->classification!=0) {
        mxSetFieldByNumber(theStruct,0,23,intMat2MatlabDoubles(&msg->classification,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,23,mxCreateDoubleMatrix(0,0,mxREAL)); 
    }
    if(msg->wind_dir!=0) {
        mxSetFieldByNumber(theStruct,0,24,intMat2MatlabDoubles(&msg->wind_dir,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,24,mxCreateDoubleMatrix(0,0,mxREAL)); 
    }
    mxSetFieldByNumber(theStruct,0,25,intMat2MatlabDoubles(&msg->spare2,1,1));
}

void AIS_8_200_24_ToMatlab(libais::Ais8_200_24 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions) {
//This message is defined in ECE-TRANS-SC3-176e.pdf
//Message "Water levels"
    const size_t numberOfFields=9;
    const mwSize dims[2] = {1, 1};
    //These are the names of all of the fields that will be added to the
    //Matlab structre array.
    const char *fieldNames[numberOfFields] = {"message_id",
    "repeat_indicator", "mmsi", "spare", "dac", "fi", "country",
    "gauge_ids", "levels"};
    mxArray *theStruct;
    
    //If text descriptions of the fields are desired, then put them into a
    //structure.
    if(fieldDescriptions!=NULL) {
        const char *descriptionStrings[numberOfFields] = {
            "Identifier for Message 8; always 8.",//0
            "Repeat indicator. Used by the repeater to indicate how many times a message has been repeated. 0 = default; 3 = do not repeat any more.",//1
            "Maritime mobile service identity (MMSI) number of source station",//2
            "Not used. Should be zero. Reserved for future use",//3
            "Designated area code (DAC)",//4
            "Function identifier (FI)",//5
            "United Nations country code string according to ERI specification",//6
            "An array of four gauge IDs, (0-2047) 0 means unknown.",//7
            "An array of water levels in meters (0-204.7), 0 means unknown",//8
        };
        unsigned int i;

        *fieldDescriptions=mxCreateStructArray(2, dims, numberOfFields, fieldNames);

        for(i=0;i<numberOfFields;i++) {
            mxSetFieldByNumber(*fieldDescriptions,0,static_cast<int>(i),mxCreateCharMatrixFromStrings(1,&descriptionStrings[i]));
        }
    }
    
    //Create the Matlab structure array.
    theStruct=mxCreateStructArray(2, dims, numberOfFields, fieldNames);
    *decodedMessage=theStruct;
    
    //Fill all of the elements of the structure with the set that is common
    //to all AIS8 messages.
    mxSetFieldByNumber(theStruct,0,0,intMat2MatlabDoubles(&msg->message_id,1,1));
    mxSetFieldByNumber(theStruct,0,1,intMat2MatlabDoubles(&msg->repeat_indicator,1,1));
    mxSetFieldByNumber(theStruct,0,2,intMat2MatlabDoubles(&msg->mmsi,1,1));
    mxSetFieldByNumber(theStruct,0,3,intMat2MatlabDoubles(&msg->spare,1,1));
    mxSetFieldByNumber(theStruct,0,4,intMat2MatlabDoubles(&msg->dac,1,1));
    mxSetFieldByNumber(theStruct,0,5,intMat2MatlabDoubles(&msg->fi,1,1));
    
    //Fill in the components unique to this message.
    {
        const char *charString=msg->country.c_str();
        mxSetFieldByNumber(theStruct,0,6,mxCreateCharMatrixFromStrings(1,&charString));
    }
    mxSetFieldByNumber(theStruct,0,7,intMat2MatlabDoubles(msg->gauge_ids.data(),4,1));
    mxSetFieldByNumber(theStruct,0,8,floatMat2MatlabDoubles(msg->levels.data(),4,1));
}

void AIS_8_200_40_ToMatlab(libais::Ais8_200_40 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions) {
//This message is defined in ECE-TRANS-SC3-176e.pdf
//Message "Signal Status"
    const size_t numberOfFields=13;
    const mwSize dims[2] = {1, 1};
    //These are the names of all of the fields that will be added to the
    //Matlab structre array.
    const char *fieldNames[numberOfFields] = {"message_id",
    "repeat_indicator", "mmsi", "spare", "dac", "fi", "lon", "lat", "form",
    "dir", "stream_dir", "status_raw", "spare2"};
    mxArray *theStruct;
    
    //If text descriptions of the fields are desired, then put them into a
    //structure.
    if(fieldDescriptions!=NULL) {
        const char *descriptionStrings[numberOfFields] = {
            "Identifier for Message 8; always 8.",//0
            "Repeat indicator. Used by the repeater to indicate how many times a message has been repeated. 0 = default; 3 = do not repeat any more.",//1
            "Maritime mobile service identity (MMSI) number of source station",//2
            "Not used. Should be zero. Reserved for future use",//3
            "Designated area code (DAC)",//4
            "Function identifier (FI)",//5
            "Signal position East longitue in degrees (WGS-84)",//6
            "Signal position North latitude in degrees (WGS-84)",//7
            "Signal form according to Annex C of ECE-TRANS-SC3-2006-10e.pdf",//8
            "Orientation of signal in degrees East of North (0-359)",//9
            "Direction of impact:\n1 Upstream\n2 Downstream\n3 To the left bank\n4 To the right bank",//10
            "The light status as per ECE-TRANS-SC3-2006-10e.pdf",//11
            "Not used, should be set to zero. Reserved for future use."//12
        };
        unsigned int i;

        *fieldDescriptions=mxCreateStructArray(2, dims, numberOfFields, fieldNames);

        for(i=0;i<numberOfFields;i++) {
            mxSetFieldByNumber(*fieldDescriptions,0,static_cast<int>(i),mxCreateCharMatrixFromStrings(1,&descriptionStrings[i]));
        }
    }
    
    //Create the Matlab structure array.
    theStruct=mxCreateStructArray(2, dims, numberOfFields, fieldNames);
    *decodedMessage=theStruct;
    
    //Fill all of the elements of the structure with the set that is common
    //to all AIS8 messages.
    mxSetFieldByNumber(theStruct,0,0,intMat2MatlabDoubles(&msg->message_id,1,1));
    mxSetFieldByNumber(theStruct,0,1,intMat2MatlabDoubles(&msg->repeat_indicator,1,1));
    mxSetFieldByNumber(theStruct,0,2,intMat2MatlabDoubles(&msg->mmsi,1,1));
    mxSetFieldByNumber(theStruct,0,3,intMat2MatlabDoubles(&msg->spare,1,1));
    mxSetFieldByNumber(theStruct,0,4,intMat2MatlabDoubles(&msg->dac,1,1));
    mxSetFieldByNumber(theStruct,0,5,intMat2MatlabDoubles(&msg->fi,1,1));
    
    //Fill in the components unique to this message.
    if(fabs(msg->position.lng_deg)<=180) {
        mxSetFieldByNumber(theStruct,0,6,doubleMat2Matlab(&msg->position.lng_deg,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,6,mxCreateDoubleMatrix(0,0,mxREAL)); 
    }
    if(fabs(msg->position.lat_deg)<=90) {
        mxSetFieldByNumber(theStruct,0,7,doubleMat2Matlab(&msg->position.lat_deg,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,7,mxCreateDoubleMatrix(0,0,mxREAL)); 
    }
    if(msg->form!=0&&msg->form!=15) {
        mxSetFieldByNumber(theStruct,0,8,intMat2MatlabDoubles(&msg->form,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,8,mxCreateDoubleMatrix(0,0,mxREAL)); 
    }
    if(msg->dir<=359) {
        mxSetFieldByNumber(theStruct,0,9,intMat2MatlabDoubles(&msg->dir,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,9,mxCreateDoubleMatrix(0,0,mxREAL)); 
    }
    if(msg->stream_dir!=0) {
        mxSetFieldByNumber(theStruct,0,10,intMat2MatlabDoubles(&msg->stream_dir,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,10,mxCreateDoubleMatrix(0,0,mxREAL)); 
    }
    mxSetFieldByNumber(theStruct,0,11,intMat2MatlabDoubles(&msg->status_raw,1,1));
    mxSetFieldByNumber(theStruct,0,12,intMat2MatlabDoubles(&msg->spare2,1,1));
}

void AIS_8_200_55_ToMatlab(libais::Ais8_200_55 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions) {
//This message is defined in ECE-TRANS-SC3-176e.pdf
//Message "number of persons on board"
    const size_t numberOfFields=10;
    const mwSize dims[2] = {1, 1};
    //These are the names of all of the fields that will be added to the
    //Matlab structre array.
    const char *fieldNames[numberOfFields] = {"message_id",
    "repeat_indicator", "mmsi", "spare", "dac", "fi", "crew", "passengers",
    "yet_more_personnel", "spare2"};
    mxArray *theStruct;
    
    //If text descriptions of the fields are desired, then put them into a
    //structure.
    if(fieldDescriptions!=NULL) {
        const char *descriptionStrings[numberOfFields] = {
            "Identifier for Message 8; always 8.",//0
            "Repeat indicator. Used by the repeater to indicate how many times a message has been repeated. 0 = default; 3 = do not repeat any more.",//1
            "Maritime mobile service identity (MMSI) number of source station",//2
            "Not used. Should be zero. Reserved for future use",//3
            "Designated area code (DAC)",//4
            "Function identifier (FI)",//5
            "Numer of crew members on board (0-254)",//6
            "Number of passengers on board (0-8190)",//7
            "Number of shipboard personnel on board (0-254)",//8
            "Not used, should be set to zero. Reserved for future use."//9
        };
        unsigned int i;

        *fieldDescriptions=mxCreateStructArray(2, dims, numberOfFields, fieldNames);

        for(i=0;i<numberOfFields;i++) {
            mxSetFieldByNumber(*fieldDescriptions,0,static_cast<int>(i),mxCreateCharMatrixFromStrings(1,&descriptionStrings[i]));
        }
    }
    
    //Create the Matlab structure array.
    theStruct=mxCreateStructArray(2, dims, numberOfFields, fieldNames);
    *decodedMessage=theStruct;
    
    //Fill all of the elements of the structure with the set that is common
    //to all AIS8 messages.
    mxSetFieldByNumber(theStruct,0,0,intMat2MatlabDoubles(&msg->message_id,1,1));
    mxSetFieldByNumber(theStruct,0,1,intMat2MatlabDoubles(&msg->repeat_indicator,1,1));
    mxSetFieldByNumber(theStruct,0,2,intMat2MatlabDoubles(&msg->mmsi,1,1));
    mxSetFieldByNumber(theStruct,0,3,intMat2MatlabDoubles(&msg->spare,1,1));
    mxSetFieldByNumber(theStruct,0,4,intMat2MatlabDoubles(&msg->dac,1,1));
    mxSetFieldByNumber(theStruct,0,5,intMat2MatlabDoubles(&msg->fi,1,1));
    
    //Fill in the components unique to this message.
    if(msg->crew!=255) {
        mxSetFieldByNumber(theStruct,0,6,intMat2MatlabDoubles(&msg->crew,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,6,mxCreateDoubleMatrix(0,0,mxREAL)); 
    }
    if(msg->passengers!=8191) {
        mxSetFieldByNumber(theStruct,0,7,intMat2MatlabDoubles(&msg->passengers,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,7,mxCreateDoubleMatrix(0,0,mxREAL)); 
    }
    if(msg->yet_more_personnel!=255) {
        mxSetFieldByNumber(theStruct,0,8,intMat2MatlabDoubles(&msg->yet_more_personnel,1,1));
    } else {
        mxSetFieldByNumber(theStruct,0,8,mxCreateDoubleMatrix(0,0,mxREAL)); 
    }
    mxSetFieldByNumber(theStruct,0,9,intMat2MatlabDoubles(msg->spare2.data(),3,1));
}

void AIS_8_366_56_ToMatlab(libais::Ais8_366_56 *msg, mxArray **decodedMessage,mxArray **fieldDescriptions) {
    const size_t numberOfFields=7;
    const mwSize dims[2] = {1, 1};
    //These are the names of all of the fields that will be added to the
    //Matlab structre array.
    const char *fieldNames[numberOfFields] = {"message_id",
    "repeat_indicator", "mmsi", "spare", "dac", "fi", "encrypted"};
    mxArray *theStruct;
    
    //If text descriptions of the fields are desired, then put them into a
    //structure.
    if(fieldDescriptions!=NULL) {
        const char *descriptionStrings[numberOfFields] = {
            "Identifier for Message 8; always 8.",//0
            "Repeat indicator. Used by the repeater to indicate how many times a message has been repeated. 0 = default; 3 = do not repeat any more.",//1
            "Maritime mobile service identity (MMSI) number of source station",//2
            "Not used. Should be zero. Reserved for future use",//3
            "Designated area code (DAC)",//4
            "Function identifier (FI)",//5
            "Encrypted data string as unsigned characters"//6
        };
        unsigned int i;

        *fieldDescriptions=mxCreateStructArray(2, dims, numberOfFields, fieldNames);

        for(i=0;i<numberOfFields;i++) {
            mxSetFieldByNumber(*fieldDescriptions,0,static_cast<int>(i),mxCreateCharMatrixFromStrings(1,&descriptionStrings[i]));
        }
    }
    
    //Create the Matlab structure array.
    theStruct=mxCreateStructArray(2, dims, numberOfFields, fieldNames);
    *decodedMessage=theStruct;
    
    //Fill all of the elements of the structure with the set that is common
    //to all AIS8 messages.
    mxSetFieldByNumber(theStruct,0,0,intMat2MatlabDoubles(&msg->message_id,1,1));
    mxSetFieldByNumber(theStruct,0,1,intMat2MatlabDoubles(&msg->repeat_indicator,1,1));
    mxSetFieldByNumber(theStruct,0,2,intMat2MatlabDoubles(&msg->mmsi,1,1));
    mxSetFieldByNumber(theStruct,0,3,intMat2MatlabDoubles(&msg->spare,1,1));
    mxSetFieldByNumber(theStruct,0,4,intMat2MatlabDoubles(&msg->dac,1,1));
    mxSetFieldByNumber(theStruct,0,5,intMat2MatlabDoubles(&msg->fi,1,1));
    
    //Now for the data associated with this message.
    mxSetFieldByNumber(theStruct,0,5,unsignedCharMat2Matlab(msg->encrypted.data(),msg->encrypted.size(),1));
}