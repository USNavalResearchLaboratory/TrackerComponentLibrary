function fileNames=fileNamesInZipArchive(path2File, includeFolders, omitInvisibles)
%%NUMFILESINZIPARCHIVE Get the names of the files in a zip archive. These
%                      are the full path strings of the files such each
%                      string can be passed to the readZipArchive function
%                      to decode the file. Thus is is the full path in the
%                      zip archive, for example folder1/folder2/test.txt
%                      not just test.txt. One can optionally include
%                      folder names as separate entries, such as
%                      folder1/folder2/ and one can optionally include the
%                      names of files that would notmally be invisible.
%
%INPUTS:  path2File A Matlab character string representing the path to the
%                   zip archive. The path should be UNIX-style with '/' to
%                   indicate folders, not '\'.
%    includeFolders If true, folders are counted. If false, folders are not
%                   counted. The default if omitted is false.
%    omitInvisibles Omit folders and files that are invisible. The default
%                   if omitted is true.  See the function fileIsInvisible
%                   for more on how files are determined invisible.
%
%OUTPUTS fileNames A numNamesX1 cell array containing the names of the
%                  files in the archive, as per the includeFolders and
%                  omitInvisibles flags. The names include the path within
%                  the zip archive. For example, a file called file.txt in
%                  a folder called folder in the zip archive would be
%                  returned as 'folder/file.txt'. The Matlab function
%                  fileparts can be used on fileNames if one only wants the
%                  name without the folder.
%
%One can call functions in Java directly from Matlab. Java's util library
%contains a zip class for manipulating zipped files. This function just
%uses the appropriate routines from Java's library. 
%
%November 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3)
    omitInvisibles=true;
end

if(nargin<2)
    includeFolders=false;
end

%Open the zip file
zipObject=java.util.zip.ZipFile(path2File);

fileNames=cell(zipObject.size(),1);%Allocate space
curFile=0;

entryEnum=zipObject.entries();
while(entryEnum.hasMoreElements)
    nextEl=entryEnum.nextElement;
    
    fileName=char(nextEl.getName);
    numChars=length(fileName);
    
    %Folders are determined by a / not escaped with a \ at the end of the
    %name.
    isFile=(fileName(end)~='/'||(fileName(end)=='/'&&(numChars==1||fileName(end-1)=='\')));
    
    if((isFile||includeFolders)&&(~omitInvisibles||~fileIsInvisible(fileName)))
        curFile=curFile+1;
        fileNames{curFile}=fileName;
    end
end

%Resize it to fit the actual number of file names found, which is
%different from zipObject.size() if invisible files are omitted.
fileNames=fileNames(1:curFile);

%Close the zip file.
zipObject.close();

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
