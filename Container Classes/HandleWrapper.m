classdef HandleWrapper < handle
%%HANDLEWRAPPER In Matlab, everything that is not derived from the handle
%               class is passed by reference. This can make it difficult or
%               for recursive functions or functions that traverse trees to
%               efficiently fill data into an array. Thus, this class is
%               here just so that one can put any type of data into the
%               data property and let different levels of recursion modify
%               the original data without making copies of it.
%
%Note that when this function is deleted, it does not call delete on any
%data that might be in it.
%
%September 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    properties
        data
    end
methods
    function theWrapper=HandleWrapper(dataInit)
    %%HANDLEWRAPPER The constructor method of the handleWrapper class class.
    %
    %INPUTS: dataInit The initial data to put in the class. If this is
    %                 omitted, then data is set to an empty matrix.
    %
    %OUTPUTS: theTree  A new HandleWrapper object.

        if(nargin<1)
            theWrapper.data=[];
        else
            theWrapper.data=dataInit;
        end
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
