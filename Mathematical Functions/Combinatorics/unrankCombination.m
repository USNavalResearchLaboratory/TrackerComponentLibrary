function combo=unrankCombination(rank,n,m,firstElMostSig)
%%UNRANKCOMBIANTION Return the combination of the given rank in the
%                   lexicographic ordering of combinations consisting of
%                   m elements chosen from a total set of n elements,
%                   where the first element of combo is the least
%                   significant, unless otherwise specified.
%
%INPUTS: rank The order of the desired combination in lexicographic order.
%             Note that 0<=rank<binomial(n,k).
%           n The number of items from which m items are chosen for the
%             ranked combinations.
%           m The number of items chosen.
% firstElMostSig An optional parameter specifying whether the first element
%             is the most or least significant. The default if omitted or
%             an empty matrix is passed is false (the first element is the
%             least significant).
%
%OUTPUTS: combo An mX1 vector containing the combination with values in
%               INCREASING order. The lowest item is indexed zero. If a
%               rank equal to or greater than the total number of unique
%               combinations is given, then an empty matrix is returned.
%
%The rank represents the  combinatorial number system of degree m,
%where m is the length of the vector comb. The system is described in
%Chapter 7.2.1.3 of [1]. Under that system, it can be observed that the
%first time the highest-value element of a combination vector acquaires a
%new value, all of the other terms in the binomial sum needed to find the
%rank are zero. Thus, one can pick off values in a ranked combinations
%similarly to how one might use logarithms repeatedly to extract the digits
%of a number.
%
%REFERENCES:
%[1] D. E. Knuth, The Art of Computer Programming. Vol. 4, Fascicle 3:
%    Generating all Combinations and Partitions, Upper Saddle River, NJ:
%    Addison-Wesley, 2009.
%
%September 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4||isempty(firstElMostSig))
    firstElMostSig=false;
end

combo=zeros(m,1);

if(rank>=binomial(n,m))
    combo=[];
    return
end

%Indices in the combination start from zero, so the maximum number that can
%be in the combination is one less than the total number of things from
%which one can choose.
cap=n-1;
curFloor=rank;
for i=0:(m-1)
    %If all of the remaining values count down to zero.
    if(curFloor==0) 
        for k=i:(m-1)
            combo(k+1)=m-k-1;
        end
        combo=fliplr(combo')';
        return;
    end
    
    j=0;
    curBinom=binomial(cap-j,m-i);
    while(curBinom>curFloor)
        j=j+1;
        curBinom=binomial(cap-j,m-i);
    end
    
    combo(i+1)=cap-j;
    cap=cap-j-1;
    curFloor=curFloor-curBinom;
end

if(firstElMostSig==false)
    combo=flipud(combo(:));
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
