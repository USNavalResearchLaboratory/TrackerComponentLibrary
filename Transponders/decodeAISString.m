function [decodedMessage,reportName,reportDescription,fieldDescriptions]=decodeAISString(NMEAStrings)
%%DECODEAISSTRING Decode one or more Automatic Identification System (AIS)
%                 messages given a full National Marine Electronics
%                 Association (NMEA) formatted ASCII text messages.
%                 Multiple part messages can also be decoded. Information
%                 appended by the receiver after the message checksum is
%                 ignored.
%
%INPUTS: NMEAStrings A character string containing one or more Automatic
%                 Identification System (AIS) messages given as National
%                 Marine Electronics Association (NMEA) formatted ASCII
%                 text. Each message is on its own line.
%
%OUTPUTS: decodedMessage A cell array of structures whose components are
%                    the decoded message elements, one message per cell.
%                    If an incomplete or invalid message is provided, it
%                    will be skipped.
%         reportName A cell array of the names of the general type of
%                    each message decoded.
%  reportDescription A cell array of strings describing the type of each
%                    message decoded.
%  fieldDescriptions A cell array holding a structure array with the same
%                    general fields as each cell element of decodedMessage
%                    but where each field contains a string describing
%                    what is in the field.
%
%AIS messages are used for tracking cooperative ships and for the exchange
%of information of interest to ships, such as tidal information. The basic
%standard is [1] with additional messages defined in [2], [3], and at
%http://www.e-navigation.nl/sites/default/files/asm_files/GN%20Release%202-23MAR15.pdf
%and in [4]. Note that an older version of the ITU standard as well as [5]
%contain some message types that have been removed from the newer
%standards.
%
%This function is essentially a Matlab interface for libais from 
%https://github.com/schwehr/libais
%The library does not support all possible message types. Messages
%8_367_22 and 8_366_22 from the library have not been implemented as it
%appears that the standards are in flux/ or are not fully implemented in
%the library. Messages 6_0_0 and 8_366_56 will not yet work as
%decode_body.cpp in the library still lacks a appropriate switch
%statements.
%
%The algorithm can be compiled for use in Matlab using the 
%CompileCLibraries function.
%
%The algorithm is run in Matlab using the command format
%[decodedMessage,reportName,reportDescription,fieldDescriptions]=decodeAISString(NMEAStrings)
%
%Sample NMEA messages from LibAIS for the different message types are
%Message 1
%NMEAStrings='!AIVDM,1,1,,A,100WhdhP0nJRdiFFHFvm??v00L12,0*13,raishub,1342569600';
%Message 2
%NMEAStrings='!AIVDM,1,1,,B,284;UGTdP4301>3L;B@Wk3TnU@A1,0*7C,raishub,1342572913';
%Message 3
%NMEAStrings='!AIVDM,1,1,,B,34hoV<5000Jw95`GWokbFTuf0000,0*6C,raishub,1342569600';
%Message 4
%NMEAStrings='!AIVDM,1,1,,B,4h3Owoiuiq000rdhR6G>oQ?020S:,0*10,raishub,1342569600';
%Message 5 (two fragments put together using sprintf)
%NMEAStrings=sprintf('!SAVDM,2,1,6,A,55NOvQP1u>QIL@O??SL985`u>0EQ18E=>222221J1p`884i6N344Sll1@m80,0*0C,b003669956,1426118503\n!SAVDM,2,2,6,A,TRA1iH88880,2*6F,b003669956,1426118503');
%Message 6_0_0
%NMEAStrings='!AIVDM,1,1,,A,6>l4v:h0006D000000P@0T0,2*56';
%Message 6_1_0
%NMEAStrings='!AIVDM,1,1,,B,65Ps:8=:0MjP0420<4U>1@E=B10i>04<fp0,2*23';
%Message 7
%NMEAStrings='!AIVDM,1,1,,A,7jvPoD@mMq;U,0*09';
%Message 8_1_11
%NMEAStrings='!AIVDM,1,1,,A,8@2<HV@0BkLN:0frqMPaQPtBRRIrwwejwwwwwwwwwwwwwwwwwwwwwwwwwt0,2*34';
%Message 8_1_22
%NMEAStrings='!AIVDM,1,1,0,B,803Ovrh0EPM0WB0h2l0MwJUi=6B4G9000aip8<2Bt2Hq2Qhp,0*01,d-084,S1582,t091042.00,T42.19038981,r003669945,1332321042';
%Message 8_200_10
%NMEAStrings='!SAVDM,1,1,5,B,85NLn@0j2d<8000000BhI?`50000,0*2A,d,S1241,t002333.00,T33.111314,D08MN-NO-BSABS1,1429316613';
%Message 8_366_56
%NMEAStrings='!BSVDM,1,1,,B,853>IhQKf6EQFDdajT?AbaAVhHEWebddhqHC5@?=KwisgP00DWjE,0*6D,b003669701,1429315201';
%Message 9
%NMEAStrings='!AIVDM,1,1,,B,9oVAuAI5;rRRv2OqTi?1uoP?=a@1,0*74';
%Message 10
%NMEAStrings='!AIVDM,1,1,,A,:5Ovc200B=5H,0*43';
%Message 11
%NMEAStrings='!AIVDM,1,1,,B,;028j>iuiq0DoO0ARF@EEmG008Pb,0*25,raishub,1342570856';
%Message 12 (two fragments put together using sprintf)
%NMEAStrings=sprintf('!AIVDM,2,1,1,A,<02PeAPpIkF06B?=PB?31P3?>DB?<rP@<51C5P3?>D13DPB?31P3?>DB,0*13\n!AIVDM,2,2,1,A,?<P?>PF86P381>>5<PoqP6?BP=1>41D?BIPB5@?BD@,4*66');
%Message 14
%NMEAStrings='!AIVDM,1,1,,A,>>M@rl1<59B1@E=@0000000,2*0D';
%Message 15
%NMEAStrings='!SAVDM,1,1,,B,?03OwnB0ACVlD00,2*59';
%Message 16
%NMEAStrings='!SAVDO,1,1,,B,@03OwnQ9RgLP3h0000000000,0*32,b003669978,1426173689';
%Message 17
%NMEAStrings='!AIVDM,1,1,,A,A6WWW6gP00a3PDlEKLrarOwUr8Mg,0*03';
%Message 18
%NMEAStrings='!SAVDM,1,1,4,B,B5NU=J000=l0BD6l590EkwuUoP06,0*61';
%Message 19 (two fragments put together using sprintf)
%NMEAStrings=sprintf('!AIVDM,2,1,7,B,C5NMbDQl0NNJC7VNuC<v`7NF4T28V@2g0J6F::00000,0*59\n!AIVDM,2,2,7,B,0J70<RRS0,0*30');
%Message 20
%NMEAStrings='!SAVDM,1,1,6,B,Dh3OwjhflnfpLIF9HM1F9HMaF9H,2*3E';
%Message 21 (two fragments put together using sprintf)
%NMEAStrings=sprintf('!AIVDM,2,1,9,A,ENk`sO70VQ97aRh1T0W72V@611@=FVj<;V5d@00003v,0*50\n!AIVDM,2,2,9,A,P0<M0,0*3E');
%Message 22
%NMEAStrings='!SAVDM,1,1,6,A,F030owj2N2P6Ubib@=4q35b10000,0*74';
%Message 23
%NMEAStrings='!AIVDM,1,1,,B,G02:KpP1R`sn@291njF00000900,2*1C';
%Message 24
%NMEAStrings='!AIVDM,1,1,,A,H44cj<0DdvlHhuB222222222220,2*46';
%Message 25
%NMEAStrings='!AIVDM,1,1,,B,I6S`3Tg@T0a3REBEsjJcT?wSi0fM,0*02';
%Message 26 (two fragments put together using sprintf)
%NMEAStrings=sprintf('!AIVDM,2,1,2,B,JfgwlGvNwts9?wUfQswQ<gv9Ow7wCl?nwv0wOi=mwd?,0*03\n!AIVDM,2,2,2,B,oW8uwNg3wNS3tV,5*71');
%Message 27
%NMEAStrings='!AIVDM,1,1,,B,K815>P8=5EikdUet,0*6B';
%
%REFERENCES:
%[1] "Technical characteristics for an automatic identification system
%    using time division multiple access in the VHF maritime mobile
%    frequency band," International Telecommunication Union,
%    Radiocommunication Sector Std. ITU-R M.1371-5, Feb. 2014. [Online].
%    Available: http://www.itu.int/rec/R-REC-M.1371/en
%[2] "Guidance on the use of AIS application-specific messages,"
%    International Maritime Organization, 2 Jun. 2013. [Online].
%    Available: http://www.iho.int/mtg_docs/com_wg/IHOTC/IHOTC_Misc/IMO%20SN_Circ289.pdf
%[3] "International standard for tracking and tracing on inland waterways
%    (VTT)," United Nations Economic and Social Council, Economic
%    Commission for Europe, Inland Transport Committee, New York and
%    Geneva, 2007. [Online].
%    Available: http://www.unece.org/fileadmin/DAM/trans/doc/finaldocs/sc3/ECE-TRANS-SC3-176e.pdf
%[4] "St. Lawrence Seaway AIS Data Messaging Formats and Specifications,"
%    U.S. Department of Transportation, Revision 4.0A, 9 May 2002.
%    [Online].
%    Available: http://www.greatlakes-seaway.com/en/pdf/aisdata.pdf
%[5] "Guidance on the application of AIS binary messages," International
%    Maritime Organization, 28 May 2004. [Online].
%    Available: http://www.imo.org/blast/blastDataHelper.asp?data id=10741&filename=236.pd
%
%November 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
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
