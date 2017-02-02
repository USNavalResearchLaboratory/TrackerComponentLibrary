function numBits=systemNumberOfBits()
%%SYSTEMNUMBEROFBITS Return the number of bits of the system. For example,
%                    a 32-bit or 64,bit system.
%
%INPUTS:  None
%
%OUTPUTS: numBits   The number of bits of the system. This will typically
%                   be 32 or 64.
%
%This function just uses the fact that the last two characters of the file
%extension of a mex file are the number of bits of the system. If systems
%with more bits arise, for example 128-bit computers, then this function
%will not properly report the number of bits as it only looks at the last
%two characters.
%
%January 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numBits=mexext();
numBits=str2double(numBits(end-1:end));
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
