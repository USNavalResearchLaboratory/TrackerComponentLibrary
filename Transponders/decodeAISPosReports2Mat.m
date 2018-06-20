function [decodedMsgs,epochTimes]=decodeAISPosReports2Mat(NMEAStrings,convert2Metric)
%%DECODEAISPOSREPORTS2MAT Given a matrix of Automatic Identification
%                 System (AIS) messages (or a cell array of AIS messages)
%                 given as National Marine Electronics Association (NMEA)
%                 formatted ASCII text messages, extract position reports
%                 of message types 1, 2, 3, 18, 19, and 27, storing the
%                 information related to ship location, identity, course,
%                 and time in different rows of a matrix, where each
%                 column corresponds to a different report. If a cell
%                 array is passed, then the matrices for decoding each
%                 cell array are placed into corresponding cells on
%                 return. It is also possible to decode timestamps (or
%                 another integer number) appended as a final field on the
%                  messages. The manipulation of matrices rather than
%                 structures is significantly faster in Matlab, which is
%                 why one might prefer this function to decodeAISString
%                 when only simple position reports are needed.
% 
%INPUTS: NMEAStrings A character string containing one or more Automatic
%                 Identification System (AIS) messages given as National
%                 Marine Electronics Association (NMEA) formatted ASCII
%                 text. Each message is on its own line. Invalid messages
%                 and messages that are not of types 1, 2, 3, 18, 19, and
%                 27 are skipped.
% convert2Metric  If true, all distance units in the result are converted
%                 to meters, all speed units are converted to meters per
%                 second, and all angles are converted to radians North of
%                 East (the mathematical definition) and angular rates are
%                 converted to radians per second North of East, all
%                 rather than using typical nautical units. The default if
%                 omitted is true.
% 
%decodedMsgs A 14XN matrix of the N successfully decoded position reports.
%            The 14 fields that are decoded are stored in the matrix as:
%                     1) The maritime mobile service identity (MMSI)
%                        number of the transmitter that can be used to
%                        identify targets (A 30-bit integer value).
%                     2) Navigation status. An integer indicating what the
%                        ship is doing/ whether it is moving. Possible
%                        valid values are:
%                       *0 = under way using engine
%                       *1 = at anchor
%                       *2 = not under command,
%                       *3 = restricted maneuverability
%                       *4 = constrained by her draught
%                       *5 = moored
%                       *6 = aground
%                       *7 = engaged in fishing
%                       *8 = under way sailing
%                       *9 = reserved for future amendment of
%                          navigational status for ships carrying DG, HS,
%                          or MP, or IMO hazard or pollutant category C,
%                          high speed craft (HSC).
%                       *10 = reserved for future amendment of
%                          navigational status for ships carrying
%                          dangerous goods (DG), harmful substances (HS)
%                          or marine pollutants (MP), or IMO hazard or
%                          pollutant category A, wing in ground (WIG).
%                       *11 = power- driven vessel towing astern
%                              (regional use)
%                       *12 = power-driven vessel pushing ahead or towing
%                              alongside (regional use)
%                       *13 = reserved for future use,
%                       *14 = AIS-SART (active), MOB-AIS, EPIRB-AIS
%                     3) Turn rate. This is limited to +/- 708 degrees per
%                        minute (about +/- .2059 radians per second when
%                        using metric units).
%                     4) Speed over ground in knots up to 102.2 knots, or
%                        52.6 meters per second when using metric units.
%                     5) Position accuracy flag. Possible values are:
%                        1=high (accuracy beter than 10m)
%                        0=low (accuracy worse than 10 meters)
%                     6) Latitude (WGS-84) in degrees North or radians
%                        North when using the metric units option.
%                     7) Longitude (WGS 84) in degrees East or radians
%                        East when using the metric units option.
%                     8) Course over ground (heading) in degrees East of
%                        North or radians North of East when using the
%                        metric units option.
%                     9) True heading in degrees East of North or radians
%                        North of East when using the metric units option.
%                    10) The second number within the universal
%                        coordinated time (UTC) minute when the ship's
%                        transmitter sent the position report. NOTE that
%                        if this is 61 or greater, than this means that the 
%                        the positioning system is broken or in dead
%                        reckoning mode, so the return is likely to be
%                        inaccurate. Also, the standard does not properly
%                        account for 61 second minutes when a leapsecond
%                        is added, so it is not clear what will happen
%                        with leapseconds. This field is not used by
%                        message 27. Rather, there is no timestamp and a
%                        position latency indicator is used to indicate
%                        the extend of the delay of the true time of the
%                        position from the broadcast time.
%                    11) Special manoeuvre indicator. This can be
%                        1 = not engaged in special manoeuvre
%                        2 = engaged in special manoeuvre
%                    12) The ID of the message that provided the above
%                        information. This can be
%                        1  Class A position report, scheduled.
%                        2  Class A position report, assigned.
%                        3  Class A position report, when interrogated.
%                        18 Standard class B equipment position report.
%                        19 Extended class B equipment position report.
%                        27 Long-range automatic identification system
%                           broadcast message.
%                    13) A number indicating how many times the message
%                        has been repeated. This ranged from 0 to 3. 3
%                        also indicates that the message should not be
%                        repeated anymore.
%                    14) Position latency. This is only used with message
%                        27. This is 0 or 1. 0 means the position latency
%                        is less than five seconds, one means it is more.
%          epochTimes If this parameter is requested, then the final
%                     numbers after the checklsum (separated by a comma)
%                     are assumed to either be an integer epoch time, or
%                     anything with a comma separating an integer epoch
%                     time at the end. Thus, this tries to decode the
%                     epoch time associated with the message. Note that in
%                     multi-part messages (message 19 can have multiple
%                     parts), the epoch time used is that of the LAST
%                     received part in the sequence. If no time can be
%                     decoded, then
%
%This function relies on libais from 
%https://github.com/schwehr/libais
%The library probably does not work on big-endian processors (Intel
%processors are little-endian).
%
%The algorithm can be compiled for use in Matlab using the 
%CompileCLibraries function.
%
%The algorithm is run in Matlab using the command format
%decodedMsgs=decodeAISPosReports2Mat(NMEAStrings,convert2Metric);
%or
%[decodedMsgs,epochTimes]=decodeAISPosReports2Mat(NMEAStrings,convert2Metric);
%if timestamps are available
%
%The algorithm is run in Matlab using the command format
%decodedMsgs=decodeAISPosReports2Mat(NMEAStrings,convert2Metric);
%or
%[decodedMsgs,epochTimes]=decodeAISPosReports2Mat(NMEAStrings,convert2Metric);
%if timestamps are available
%
%Example 1:, Decoding a single valid ID=1 string.
% NMEAStrings='!AIVDM,1,1,,A,13HOI:0P0000VOHLCnHQKwvL05Ip,0*23';
% decodedMsgs=decodeAISPosReports2Mat(NMEAStrings,false)
%The decoded message should be:
% decodedMsgs=[227006760;%MMSI
%                      0;%Navigation status (underway)
%                    NaN;%Not used
%                      0;%Turn rate
%                      0;%Speed over ground
%     49.475578308105469;%North latitude
%      0.131380006670952;%East longitude
%     36.700000762939453;%Heading East of North
%                    NaN;%Not used
%                     14;%Seconds within minute when message sent.
%                    NaN;%Not used.
%                      1;%Class A position report, scheduled.
%                      0;%has not ben rebroadcast.
%                    NaN];%Not used.
%
%Example 2: Decoding multiple messages, one valid string of each type of
%message. It is most likely that this would arise from reading the
%messages in from a file, because adding newline characters to strings has
%to be done using the \n newline character to separate the string in
%sprintf rather than just typing them. The AIS19 string consists of two
%parts that are joined.
% NMEAStrings=sprintf('!AIVDM,1,1,,A,13HOI:0P0000VOHLCnHQKwvL05Ip,0*23\n!AIVDM,1,1,,B,284;UGTdP4301>3L;B@Wk3TnU@A1,0*7C\n!AIVDM,1,1,,B,34hoV<5000Jw95`GWokbFTuf0000,0*6C\n!SAVDM,1,1,4,B,B5NU=J000=l0BD6l590EkwuUoP06,0*61\n!AIVDM,2,1,7,B,C5NMbDQl0NNJC7VNuC<v`7NF4T28V@2g0J6F::00000,0*59\n!AIVDM,2,2,7,B,0J70<RRS0,0*30\n!AIVDM,1,1,,B,K815>P8=5EikdUet,0*6B');
% decodedMsgs=decodeAISPosReports2Mat(NMEAStrings,false)
%One will observe that 6 messages are correctly decoded, of types 1, 2, 3,
%18, 19, and 27.
%
%Example 3: Decoding multiple messages, some of them of the wrong type.
%When messages of the wrong type are used, then they are askipped. Here,
%two correct messages are read in as well as one incorrect message in
%between them. The "incorrect" message is of type 26 (not a position report) and is thus ignored.
%Below are four strings, but only two of them contain AIS position reports
% NMEAStrings=sprintf('!AIVDM,1,1,,A,13HOI:0P0000VOHLCnHQKwvL05Ip,0*23\n!AIVDM,1,1,,A,H44cj<0DdvlHhuB222222222220,2*46\n!AIVDM,1,1,,B,284;UGTdP4301>3L;B@Wk3TnU@A1,0*7C');
% decodedMsgs=decodeAISPosReports2Mat(NMEAStrings)
%The results are in metric units.
%
%Example 4: Consider the case where the messages contain some numbers
%after the checksum, which make up the timestamps of the messages.
%This is the same as the previous example, but with fake timestamps
%appended.
% NMEAStrings=sprintf('!AIVDM,1,1,,A,13HOI:0P0000VOHLCnHQKwvL05Ip,0*23,1241827197\n!AIVDM,1,1,,A,H44cj<0DdvlHhuB222222222220,2*46,1241827198\n!AIVDM,1,1,,B,284;UGTdP4301>3L;B@Wk3TnU@A1,0*7C,1241827199');
% [decodedMsgs,epochTimes]=decodeAISPosReports2Mat(NMEAStrings)
%The results are in metric units.
%
%Example 5: Some receivers will append aditional information before the
%timestamp. As long as the timestamp is the last thing in the message,
%then we can still recover it.
% NMEAStrings='!AIVDM,1,1,,A,100WhdhP0nJRdiFFHFvm??v00L12,0*13,raishub,1342569600';
% [decodedMsgs,epochTimes]=decodeAISPosReports2Mat(NMEAStrings)
%The results are in metric units.
%
%November 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

error('This function is only implemented as a mexed C or C++ function. Please run CompileCLibraries.m to compile the function for use.')

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
