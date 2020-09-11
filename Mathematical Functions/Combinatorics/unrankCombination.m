function combo=unrankCombination(theRank,n,m,firstElMostSig,startVal,algorithm)
%%UNRANKCOMBINATION Return the combination of the given rank in a specified
%                ordering of combinations consisting of m elements chosen
%                from a total set of n elements, where the first element of
%                combo is the least significant, unless otherwise
%                specified.
%
%INPUTS: theRank The order of the desired combination in colexicographic
%             order. Note that 0<=rank<binomial(n,m).
%           n The number of items from which m items are chosen for the
%             ranked combinations.
%           m The number of items chosen.
% firstElMostSig An optional parameter specifying whether the first element
%             is the most or least significant. The default if omitted or
%             an empty matrix is passed is false (the first element is the
%             least significant).
%    startVal This is zero or 1, indicating which value the value at which
%             the elements in combo can start. The default if omitted or an
%             empty matrix is passed is 0.
%   algorithm An optional parameter, which specifies the algorithm to use
%             and thus also defines the ordering used. Possible values are:
%             0 (The default if omitted or an empty matrix is passed) Use
%               colexicographic ordering via the combinatorial number
%               system that is described in Chapter 7.2.1.3 of [3].
%             1 Use lexicographic ordering via the algorithm of [1] with
%               the corrections of [2].
%
%OUTPUTS: combo An mX1 vector containing the combination with values in
%               INCREASING order if firstElMostSig=false. The lowest item
%               is indexed startVal. If a rank equal to or greater than the
%               total number of unique combinations is given, then an empty
%               matrix is returned.
%
%EXAMPLE:
%We demonstrate that the output of this function agrees with that of the
%equivalent ordering in genAllCombinations when irstElMostSig=false.
% n=11;
% m=7;
% startVal=0;
% algorithm=1;
% firstElMostSig=false;
% theCombos=genAllCombinations(n,m,startVal,algorithm);
% 
% numCombos=binomial(n,m);
% theUnrankedCombos=zeros(m,numCombos);
% for curRank=0:(numCombos-1)
%     theUnrankedCombos(:,curRank+1)=unrankCombination(curRank,n,m,firstElMostSig,startVal,algorithm);
% end
% assert(all(theCombos(:)==theUnrankedCombos(:)))
%There will be no error, because the assertion is true and both outputs are
%the same (even if we switch algorithm from 0 to 1.
%
%REFERENCES:
%[1] B. P. Buckles and M. Lybanon, "Algorithm 515: Generation of a vector
%    from the lexicographic index [G6]," ACM Transactions on Mathematical
%    Software, vol. 3, no. 2, pp. 180-182, Jun. 1977.
%[2] D. F. Crouse, "Remark on Algorithm 515: Generation of a Vector from
%    the lexicographic index combinations," ACM Transactions on
%    Mathematical Software, vol. 33, no. 2, Article 12, Jun. 2006
%[3] D. E. Knuth, The Art of Computer Programming. Vol. 4, Fascicle 3:
%    Generating all Combinations and Partitions, Upper Saddle River, NJ:
%    Addison-Wesley, 2009.
%
%July 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<6||isempty(algorithm))
    algorithm=0;
end

if(nargin<5||isempty(startVal))
    startVal=0;
end

if(nargin<4||isempty(firstElMostSig))
    firstElMostSig=false;
end

switch(algorithm)
    case 0%Colexicographic order.
        combo=unrankColexCombo(theRank,n,m,firstElMostSig,startVal);
    case 1%Lexicographic order.
        theRank=theRank+1;%The function unrankLexCombo starts numbering its
        %lexicographic combinations starting at 1, so we adjust the rank.
        combo=unrankLexCombo(theRank,n,m,firstElMostSig,startVal);
    otherwise
        error('Unknown algorithm selected.')
end
end

function c=unrankLexCombo(theRank,n,m,firstElMostSig,startVal)
%%UNRANKLEXCOMBO Return the combination of the given rank in the
%                lexicographic ordering of combinations consisting of m
%                elements chosen from a total set of n elements, where the
%                first element of combo is the least significant, unless
%                otherwise specified.
%
%INPUTS: rank The order of the desired combination in colexicographic
%             order. Note that 1<=rank<=binomial(n,m).
%           n The number of items from which m items are chosen for the
%             ranked combinations.
%           m The number of items chosen.
% firstElMostSig An optional parameter specifying whether the first element
%             is the most or least significant. The default if omitted or
%             an empty matrix is passed is false (the first element is the
%             least significant).
%    startVal This is zero or 1, indicating which value the value at which
%             the elements in combo can start. The default if omitted or an
%             empty matrix is passed is 0.
%
%OUTPUTS: combo An mX1 vector containing the combination with values in
%               INCREASING order if firstElMostSig=false. The lowest item
%               is indexed startVal. If a rank equal to or greater than the
%               total number of unique combinations is given, then an empty
%               matrix is returned.
%
%The algorithm of [1] is used with the corrections of [2].
%
%REFERENCES:
%[1] B. P. Buckles and M. Lybanon, "Algorithm 515: Generation of a vector
%    from the lexicographic index [G6]," ACM Transactions on Mathematical
%    Software, vol. 3, no. 2, pp. 180-182, Jun. 1977.
%[2] D. F. Crouse, "Remark on Algorithm 515: Generation of a Vector from
%    the lexicographic index combinations," ACM Transactions on
%    Mathematical Software, vol. 33, no. 2, Article 12, Jun. 2007.
%
%July 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.

if(theRank>binomial(n,m))
    c=[];
    return
end

if(m==1)
    c=theRank;
    return;
end

c=zeros(m,1);

k=0;%The lower bound index.

p1=m-1;
c(1)=0;%Initialize the first element.
for i=1:p1
    if(i~=1)
        c(i)=c(i-1);
    end
   
    while(1)
        c(i)=c(i)+1;
        r=binomial(n-c(i),m-i);
        k=k+r;
        if(k>=theRank)
            break;
        end
    end
    k=k-r;
end

c(m)=c(p1)+theRank-k;

c=c+(startVal-1);
if(firstElMostSig)
    c=flipud(c(:));
end
end

function combo=unrankColexCombo(theRank,n,m,firstElMostSig,startVal)
%%UNRANKCOLEXCOMBIANTION Return the combination of the given rank in the
%                   colexicographic ordering of combinations consisting of
%                   m elements chosen from a total set of n elements,
%                   where the first element of combo is the least
%                   significant, unless otherwise specified.
%
%INPUTS: theRank The order of the desired combination in colexicographic
%             order. Note that 0<=rank<binomial(n,m).
%           n The number of items from which m items are chosen for the
%             ranked combinations.
%           m The number of items chosen.
% firstElMostSig An optional parameter specifying whether the first element
%             is the most or least significant. The default if omitted or
%             an empty matrix is passed is false (the first element is the
%             least significant).
%    startVal This is zero or 1, indicating which value the value at which
%             the elements in combo can start. The default if omitted or an
%             empty matrix is passed is 0.
%
%OUTPUTS: combo An mX1 vector containing the combination with values in
%               INCREASING order if firstElMostSig=false. The lowest item
%               is indexed startVal. If a rank equal to or greater than the
%               total number of unique combinations is given, then an empty
%               matrix is returned.
%
%The rank represents the combinatorial number system of degree m,
%where m is the length of the vector comb. The system is described in
%Chapter 7.2.1.3 of [1]. Under that system, it can be observed that the
%first time the highest-value element of a combination vector acquires a
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

if(nargin<5||isempty(startVal))
    startVal=0;
end

if(nargin<4||isempty(firstElMostSig))
    firstElMostSig=false;
end

combo=zeros(m,1);

if(theRank>=binomial(n,m))
    combo=[];
    return
end

%Indices in the combination start from zero, so the maximum number that can
%be in the combination is one less than the total number of things from
%which one can choose.
cap=n-1;
curFloor=theRank;
for i=0:(m-1)
    %If all of the remaining values count down to zero.
    if(curFloor==0) 
        for k=i:(m-1)
            combo(k+1)=m-k-1;
        end
        break;
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

combo=combo+startVal;
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
