function boolVal=fileIsInvisible(path)
%%FILEISINVISIBLE  Given a full path to a file or folder, indicate whether
%                  the path appears to go through invisible fodlers or
%                  whether the folder/ file is invisible based on the names
%                  of the folders in the path and the filename.
%
%INPUTS:  path    The path to the file. For example, '/test.txt'. The path
%                 should be UNIX-style with '/' to indicate folders, not
%                 '\'.
%
%OUTPUTS: boolVal True if the file is invisible or is a folder based on the
%                 path.
%
%This function does not check the metadata of the file/ of the folders in
%the path. Thus, a file marked as invisible using HFS+ flags will not
%appear to be invisible. Rather, this function looks for common indicators
%in the path and filename that indicate that the file is invisible. If a
%foldername or the filename begin with __ or ., then the file is assumed to
%be invisible. The beginning of a folder name/ filename that is not at the
%beginning of the path is identified by a / character that has not been
%escaped by a \ character.
%
%This function is useful when decompressing zip files that were compressed
%under Mac OS X, where a bunch of .DS_Store files indicating window
%placement and a __MACOSX folder holds metadata associated with the files.
%Often, one does not care about any of that and would prefer to just not
%decompress thoe files at all.
%
%November 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    numChars=length(path);
    boolVal=false;
    
    if(numChars==0)
        return;
    end
    
    %If the first folder/ file begins with a period, then mark it as being
    %invisible.
    if(path(1)=='.')
       boolVal=true;
       return;
    end
    
    if(numChars==1)
        return;
    end
    
    if(path(1)=='_'&&path(2)=='_')
       boolVal=true;
       return;
    end
    
    for curChar=3:numChars
        %If a period is encountered at the beginning of a file or folder
        %name (which is identified by a / in the path that has not been
        %escaped by using a \ before it.
        if(path(curChar)=='.'&&path(curChar-1)=='/'&&path(curChar-1)~='\')
            boolVal=true;
            return;
        end
        
        %See if a __ is present at the beginning of a file name, whereby
        %the file after the previous folder is identified by a non-escaped
        %/ character.
        if(path(curChar)=='_'&&path(curChar-1)=='_'&&path(curChar-2)=='/'&&((curChar==3)||path(curChar-3)~='\'))
            boolVal=true;
            return;
        end
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
