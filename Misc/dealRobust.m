function varargout=dealRobust(varargin)
%%DEALROBUST The deal function in Matlab takes a variable number of inputs
%            in and then outputs them as the same number of outputs. The
%            function is useful when creating anonymous functions to pass
%            to other functions. However, the deal function will throw an
%            error if the function calling the anonymous function asks for
%            fewer outputs than are passed. For example, this mgiht arise
%            if a routine wants a passed function handle to return a value
%            and its derivative, but the calling function does not need the
%            derivative at some point and only requests the value. This
%            function is the same as matlab's deal function, except it does
%            not throw an error if fewer outputs are requested than inputs
%            passed.
%
%INPUTS: A variable number of items that must be at least as large as the
%       number of outputs requested.
%
%OUTPUTS: The outputs are the same as the inputs.
%
%October 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    varargout=varargin(1:nargout);
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
