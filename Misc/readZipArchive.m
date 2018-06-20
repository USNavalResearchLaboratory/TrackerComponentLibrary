function unzippedData=readZipArchive(path2File,param2)
%%READZIPARCHIVE  Read and decompress a single gzip compressed file or
%                 one or more files from a zip archive into memory without
%                 saving them to disk. Gzip archives only store a single
%                 file, zip archives store multiple files. When reading
%                 multiple files from a zip, one has the option of omitting 
%                 invisible files. Note that decompressing a large file 
%                 into memory might cause the Java virtual machine (JVM) to
%                 run out of memory (and thus, the decompression fails), as
%                 its memory usage is capped in Matlab. The memory for the
%                 JVM must be a few times larger than the size of the file
%                 being decompressed. The memory that Matlab allows for the
%                 JVM is independent of the amount of RAM in the computer.
%                 One might have to increase the memory allotment for the
%                 JVM to extract large files. This is done via
%                 Preferences->General->Java Heap Memory.
%                 
%INPUTS:  path2File  A Matlab character string representing the path to the
%                    zip archive. For example './text.zip'. The path should
%                    be OS X/UNIX-style with '/' to indicate folders, not '\'.
%                    If the file extension is .gz, then it is assumed to be
%                    a gzip file. Otherwise, it is assumed to be a zip
%                    archive.
%         param2     An optional parameter that is only relevant when
%                    decoding a zip file, not a gzip file. If this is
%                    omitted, then all of the files in the archive are
%                    decoded and invisibles files are omitted. If param2 is
%                    a character string, then it is fileName:
%                    An optional string representing the name of a single
%                    file in the archive that is to be decompressed. This
%                    must be the full path within the zip file. For
%                    example, if someone compressed test.txt with other
%                    files in the folder test_folder, then the path to that
%                    file would be 'test_folder/test.txt'. If fileName is
%                    provided but is not in the archive, then an error will
%                    occur.
%                    Otherwise, one can set param2 to a boolean value
%                    decodeInvis, whereby all files are decoded such that:
%                    If this is true, then invisible files are also
%                    decoded. if this is false, then invisible files are
%                    not decoded. The default if omitted is false. See the
%                    function fileIsInvisible for more on how files are
%                    determined invisible.
%
%OUTPUTS: unzippedData A numFiles cell array containing the unzipped
%                      data and the names of the unzipped files.
%                      unzippedData{1,i} contains the name of the ith 
%                      files and unzippedData{2,i} contains the contents
%                      of the file as an array of uint8 values. If a zipped
%                      text file is being read, then char(unzippedData{2,i})
%                      will return a Matlab string of the contents of the
%                      file. When decoding a gzip archive,
%                      unzippedData{1,1} will just be the name of the
%                      archive without the .gz extension as gzip is a pure
%                      compression routine and does not provide information
%                      on what was compressed.
%
%One can call functions in Java directly from Matlab. Java's util library
%contains a zip class for manipulating zipped files. This function just
%uses the appropriate routines from Java's library. It is not as fast as it
%could be, though, since the data must be read one byte at a time due to
%limitations in Matlab.
%
%Matlab has a built-in function, unzip. However, that function unzips the
%data onto the disk, not into memory. This function makes use of the
%undocumented com.mathworks.mlwidgets.io.InterruptibleStreamCopier Java
%class so that one does not have to copy data very slowly byte by byte to
%decompress it.
%
%This function can be called as
%unzippedData=readZipArchive(path2File,fileName);
%To decode a single file or as
%unzippedData=readZipArchive(path2File,decodeInvis);
%to decode the whole archive, with the boolean decodeInvis to indicate
%whether invisible files should be decoded or as
%unzippedData=readZipArchive(path2File);
%where the whole archive without invisible files is decoded into memory.
%
%January 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    %Undocumented Matlab classes to copy byte arrays instead of having to
    %do it one byte at a time on an input stream, because Matlab normally
    %does everything by value, so one cannot normally just pass it a byte
    %array and tell it to fill it.
    streamCopier = com.mathworks.mlwidgets.io.InterruptibleStreamCopier.getInterruptibleStreamCopier;
    byteOutputSteam = java.io.ByteArrayOutputStream;

    %If the file is a gzip file, then  there can only be one thing in it
    %and this will just decompress that one thing into unzippedData. The
    %gzip file does not contain any information on the name of the file.
    [~,filename,ext] = fileparts(path2File);
    if(strcmp(ext,'.gz')==1)
        if(nargin>1)
           error('Too many inputs provided with a gzipped file'); 
        end
        
        unzippedData=cell(2,1);
        unzippedData{1,1}=filename;
        
        myFileInputStream=java.io.FileInputStream(path2File);
        gzipInputStream=java.util.zip.GZIPInputStream(myFileInputStream);
        
        %Perform garbage collection before trying to decompress the file so
        %as to lessen the likelihood that the Java virtual machine in
        %Matlab runs out of memory.
        java.lang.System.gc();
        %Decompress the data using Matlab's undocumented method to copy
        %the data stream (as opposed to having to use a loop to copy
        %everything one byte at a time, very slowly).
        streamCopier.copyStream(gzipInputStream,byteOutputSteam);
        
        unzippedData{2,1}=typecast(byteOutputSteam.toByteArray,'uint8')';

        %CLose the gzip file.
        myFileInputStream.close();
        return;
    end
    
    %Open the zip file.
    zipObject=java.util.zip.ZipFile(path2File);
    
    %If one just wants to decode a specific file
    if(nargin>1&&isa(param2,'char'))
        fileName=param2;
        unzippedData=cell(2,1);
        
        unzippedData{1,1}=fileName;
        theEntry=zipObject.getEntry(fileName);
        inputStream=zipObject.getInputStream(theEntry);
        
        %Perform garbage collection before trying to decompress the file so
        %as to lessen the likelihood that the Java virtual machine in
        %Matlab runs out of memory.
        java.lang.System.gc();
        
        %Decompress the data using Matlab's undocumented method to copy
        %the data stream (as opposed to having to use a loop to copy
        %everything one byte at a time, very slowly).
        streamCopier.copyStream(inputStream,byteOutputSteam);
        unzippedData{2,1}=typecast(byteOutputSteam.toByteArray,'uint8')';

        %Close the zip file.
        zipObject.close();
        return;
    elseif(nargin>1)
        decodeInvis=param2;
    else
        decodeInvis=false;
    end

    %If one wants to decode all of the files.
    entryEnum=zipObject.entries();
    curFile=0;
    unzippedData=cell(2,zipObject.size());%Allocate space
    while(entryEnum.hasMoreElements)
        nextEl=entryEnum.nextElement;
        fileName=char(nextEl.getName);
        numChars=length(fileName);
        
        %If no part of the path is invisible (an invisible folder or file)
        %and if the file is not actually a folder (which is indicated by a
        %/ not escaped with a \), then decode it.
        if((decodeInvis||~fileIsInvisible(fileName))&&(fileName(end)~='/'||(fileName(end)=='/'&&(numChars==1||fileName(end-1)=='\'))))
            curFile=curFile+1;
            %Save the name of the file.
            unzippedData{1,curFile}=fileName;
            inputStream=zipObject.getInputStream(nextEl);
            
            %In this instace, garbage collection is not performed before
            %decompressing each file lest it be slow.
            
            %Decompress the data using Matlab's undocumented method to copy
            %the data stream (as opposed to having to use a loop to copy
            %everything one byte at a time, very slowly).
            streamCopier.copyStream(inputStream,byteOutputSteam);
            unzippedData{2,curFile}=typecast(byteOutputSteam.toByteArray,'uint8')';
        end
    end
    
    %Resize it to fit the actual number of files decoded, which is
    %different from zipObject.size() if invisible files are omitted.
    unzippedData=unzippedData(:,1:curFile);
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

