function theOutput=getOutputN(f,N,varargin)
%%GETOUTPUTN Given the handle to a function f and parameters to pass to
%       that function, this function calls f and only returns the nth
%       parameter. For functions that check how many outputs were
%       requested, exactly N are requested when calling f. This function
%       can be useful when creating anonymout functions where one wishes to
%       use something other than this first output of the function.
%
%INPUTS: f A handle to the function to call.
%        N The output of f that is desired.
%        A variable number of other items may be passed to this function.
%        They will all be passed to f when called.
%
%OUTPUTS: The outputs if the Nth output of f(varargin{:}).
%
%EXAMPLE:
%Here, we call the function directly and can see that it correctly returns
%the third output.
% N=3;
% f=@(X)eig(X);
% A=magic(4);
% W=getOutputN(f,N,A)
% [~,~,W]=eig(A)
%
%July 2023 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

switch(N)
    case 1
        theOutput=f(varargin{:});
    case 2
        [~,theOutput]=f(varargin{:});
    case 3
        [~,~,theOutput]=f(varargin{:});
    case 4
        [~,~,~,theOutput]=f(varargin{:});
    case 5
        [~,~,~,~,theOutput]=f(varargin{:});
    case 6
        [~,~,~,~,~,theOutput]=f(varargin{:});
    case 7
        [~,~,~,~,~,~,theOutput]=f(varargin{:});
    case 8
        [~,~,~,~,~,~,~,theOutput]=f(varargin{:});
    case 9
        [~,~,~,~,~,~,~,~,theOutput]=f(varargin{:});
    case 10
        [~,~,~,~,~,~,~,~,~,theOutput]=f(varargin{:});
    otherwise
        error('N must be a real integer from 1 to 10.')
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
