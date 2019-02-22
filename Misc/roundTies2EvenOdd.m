function y=roundTies2EvenOdd(x,mode)
%%ROUNDTIES2EVENODD Round the values in x to the nearest integer. Ties will
%                   either round to even or odd integers.
%
%INPUTS:   x A scalar, vector, matrix of values to round to the nearest
%            integer.
%       mode A parameter specifying how ties (when the value is an integer
%           +0.5) are handled. Possible values are
%           0 (The default if omitted or an empty matrix is passed) Round
%             ties to even integers.
%           1 Round ties to odd integers.
%
%OUTPUTS: y x rounded with the selected handling of ties.
%
%EXAMPLE:
%As an example of rounding to even versus odd integers, consider:
% x=[-1.5,-2.5,-3.5,-1.1,-1.9,-2.1,-2.9,0,1.5,2.5,3.5,1.1,1.9,2.1,2.9];
% y1=roundTies2EvenOdd(x,0)
% y2=roundTies2EvenOdd(x,1)
%One will find that
%y1=[-2,-2,-4,-1,-2,-2,-3,0,2,2,4,1,2,2,3]
%y2=[-1,-3,-3,-1,-2,-2,-3,0,1,3,3,1,2,2,3]
%
%November 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    if(nargin<2||isempty(mode))
        mode=0;
    end
  
    switch(mode)
        case 0%Round with ties going to even
            y=round(x);
            sel1=abs(x-fix(x))==0.5 & mod(fix(x),2)==0;
            y(sel1)=fix(x(sel1));
        case 1%Round with ties going to odd.
            y=round(x);
            sel1=abs(x-fix(x))==0.5 & mod(fix(x),2)==1;
            y(sel1)=fix(x(sel1));
        otherwise
            error('Wrong Mode specified')
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
