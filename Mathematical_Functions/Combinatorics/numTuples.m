function numTup=numTuples(maxVals)
%%NUMTUPLES Determine the number of tuples of a certain number of digits
%           (indices).  Each index of the tuple counts from 0 to the
%           corresponding value in maxVals. tuples can be obtained using
%           the genAllTuples function.
%
%INPUTS: maxVals An NX1 or 1XN vector of the bases of each of the N digits
%                in the tuples minus 1. All elements must be >=0.
%
%OUTPUTS: The number of possible tuples.
%
%The number of tuples is just prod(maxVals+1).
%
%EXAMPLE:
% maxVals=[4;4;8;2];
% numTup=numTuples(maxVals)
%One will get numTup=675.
%
%September 2020 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    numTup=prod(maxVals+1);

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
