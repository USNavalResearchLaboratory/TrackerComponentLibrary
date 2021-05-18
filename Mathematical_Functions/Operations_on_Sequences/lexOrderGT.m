function compRes=lexOrderGT(a,b,firstElMostSig)
%%LEXORDERGT Given two sequences of numbers, determine whether sequence a
%            >b in lexicographic ordering. The ordering can be determined
%            either taking the first element as the most significant, or
%            taking the last element as the most significant.
%
%INPUTS      a, b Two length n vectors whose elements can be compared using
%                 < and > operations.
%  firstElMostSig An optional boolean value indicating whether or not the
%                 first element of a,b is the most signficiant (or whetehr
%                 the last element is most significant). The default if
%                 this parameter is omitted or an empty matrix is passed is
%                 false.
%
%OUTPUTS compRes A value indicating the relation between sets a and b.
%                Possible values are
%                1 a>b in lexicographic ordering
%                0 a=b
%               -1 a<b in lexicographic ordering
%
%Lexicographic ordering is akin to alphabetical ordering, but is often
%applied to combinations and other items where each element takes a
%restricted set of values. For example, combinations of 6 choose 5 given in
%lexicographic order with the first digit being the most significant is
% 210
% 310
% 320
% 321
% 410
% 420
% 421
% 430
% 431
% 432
% 510
% 520
% 521
% 530
% 531
% 532
% 540
% 541
% 542
% 543
%Consequently,
% compRes=lexOrderGT([4;1;0],[3;2;1],true)
%Returns 1.
%
%The ordering is true is a(i)>b(i) and if a(i)=b(i_, then one compare
%additional elements. For example
% compRes=lexOrderGT([1;1;1;4],[1;1;2;3],true)
%returns -1.
%
%January 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(firstElMostSig))
    firstElMostSig=false;
end

numIdx=length(a);
compRes=0;

%If the first element is the most significant.
if(firstElMostSig)
    for curIdx=1:numIdx
        if(a(curIdx)>b(curIdx))
            compRes=1;
            break;
        elseif(a(curIdx)<b(curIdx))
            compRes=-1;
            break;
        end
        %The loop continues if the values are equal.
    end
else
    for curIdx=numIdx:-1:1
        if(a(curIdx)>b(curIdx))
            compRes=1;
            break;
        elseif(a(curIdx)<b(curIdx))
            compRes=-1;
            break;
        end
        %The loop continues if the values are equal.
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
