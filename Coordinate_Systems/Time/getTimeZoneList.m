function [timeZoneNames,timeZoneProperties]=getTimeZoneList()
%%GETTIMEZONELIST Get a list of available time zones along with general
%                 information about the time zones. This information comes
%                 from the values built into the standard Java library,
%                 which can be accessed from within Matlab.
%
%INPUTS: None
%
%OUTPUTS: timeZoneNames A 2XnumTimeZones cell array containing the names of
%              the time zones. timeZoneNames{1,i} is the short name of the
%              ith time zone, such as one can pass to the getTimeZoneInfo
%              function, and timeZoneNames{2,i} is a longer version of the
%              name.
%         timeZoneProperties A 2XnumTimeZones list of properties of the
%              time zones. timeZoneProperties(1,i) is the offset of the ith
%              time zone from UTC in seconds (without daylight savings
%              time). timeZoneProperties(2,i) is true if the time zone is
%              listed as having ever supported daylight savings time or if
%              it is listed as ever supporting it in the future.
%
%This function just calls the appropriate methods in java.util.TimeZone.
%
%May 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    timeZones=java.util.TimeZone.getAvailableIDs();
    numZones=length(timeZones);
    
    timeZoneNames=cell(2,numZones);
    timeZoneProperties=zeros(2,numZones);
    for curZone=1:numZones
        timeZone=java.util.TimeZone.getTimeZone(timeZones(curZone));
        timeZoneNames{2,curZone}=char(timeZone.getDisplayName());
        timeZoneProperties(1,curZone)=timeZone.getRawOffset()/1e3;%1e3 put it in seconds
        timeZoneProperties(2,curZone)=timeZone.observesDaylightTime() || timeZone.useDaylightTime();
    end
    timeZoneNames(1,:)=cell(timeZones);
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
