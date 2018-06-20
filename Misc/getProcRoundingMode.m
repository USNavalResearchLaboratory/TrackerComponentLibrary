function roundMode=getProcRoundingMode()
%%GETPROCROUNDIGNMODE Get the rounding mode used in the processor. In
%          Matlab, the fegetround function in C99 does not work. Thus,
%          rather than using a mex file to read the processor rounding mode
%          directly, we must infer the rounding mode based on how things
%          actually are rounded. Control of the rounding mode is useful
%          when implementing routines utilizing interval algebra.
%
%INPUT: None
%
%OUTPUTS: roundMode An integer specifying the rounding mode to use.
%                   Possible values are
%                   0 Rounding is done towards negative infinity.
%                   1 Rounding is done towards zero.
%                   2 Rounding is done to the nearest value.
%                   3 Rounding is done towards positive infinity. 
%
%The function will have an error if the rounding mode is inconsistent with
%all of the above modes. To set the processor rounding mode, use the
%function setProcRoundingMode.
%
%January 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

test1=1+eps(0)>1;
test2=1-eps(0)<1;
test3=-1+eps(0)>-1;
test4=-1-eps(0)<-1;

if(test1==0&&test2==1&&test3==0&&test4==1)
    roundMode=0;
elseif(test1==0&&test2==1&&test3==1&&test4==0)
    roundMode=1;
elseif(test1==0&&test2==0&&test3==0&&test4==0)
    roundMode=2;
elseif(test1==1&&test2==0&&test3==1&&test4==0)
    roundMode=3;
else
    error('Unknown Rounding Mode Encountered.')
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
