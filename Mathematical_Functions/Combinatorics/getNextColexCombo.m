function [theCombo,recurData]=getNextColexCombo(param1,t,startVal)
%%GETNEXTCOLEXCOMBO Return the next combination in colexicographic order
%              given the current combination. If the final combination in
%              the sequence has been reached, then an empty matrix is
%              returned. This function is similar to getNextCombo, except
%              the sequence is in colexicographic order, not lexicographic
%              order, and the parameters for the recursion are different.
%
%INPUTS: param1 If the first combination in the sequence is desired, then
%               more than one input must be given and param1 is n, the
%               total number of items from which combinations are taken. If
%               a subsequent value in the sequence is desired, then param1
%               is recurData and no other inputs can be provided.
%             t The number of items to choose from the total of n. This
%               must be given if the first combination is desired and this
%               must not be given when going on to subsequent combinations.
%      startVal This is zero or 1, indicating which value the value at
%               which the elements in I can start. The default if omitted
%               or an empty matrix is passed is 0. This parameter can only
%               be given when requesting the first combination. It must not
%               be provided when requesting subsequent combinations (it is
%               saved in recurData).
%
%OUTPUTS: theCombo The next combination in the colexicographic sequence, or
%                  an emptycmatrix if the final combination in the
%                  colexicographic order has been passed.
%
%Algorithm T of Section 7.2.1.3 of [1] is used.
%
%EXAMPLE:
%Here, we show that the sequence produced by this function is consistent
%with that given by the correct algorithmic choice in genAllCombinations and
%in unrankCombination.
% n=10;
% m=7;
% startVal=1;
% firstElMostSig=false;
% numCombos=binomial(n,m);
% 
% allCombosBatch=genAllCombinations(n,m,startVal,0);
% allRecurCombos=zeros(m,numCombos);
% allUnrankedCombos=zeros(m,numCombos);
% numCombos=binomial(n,m);
% 
% [allRecurCombos(:,1),recurData]=getNextColexCombo(n,m,startVal);
% allUnrankedCombos(:,1)=unrankCombination(0,n,m,firstElMostSig,startVal);
% for k=2:numCombos
%     allUnrankedCombos(:,k)=unrankCombination(k-1,n,m,firstElMostSig,startVal);
%     [allRecurCombos(:,k),recurData]=getNextColexCombo(recurData);
% end
% assert(isempty(getNextColexCombo(recurData)))%End of sequence
% assert(all(allRecurCombos(:)==allUnrankedCombos(:)))
% assert(all(allRecurCombos(:)==allCombosBatch(:)))
%All of the assertions will be true, indicating no error, which signifies
%the equivalence of the sequence produced by this function and the other
%functions.
%
%REFERENCES:
%[1] D. E. Knuth, The Art of Computer Programming. Vol. 4A: Combinatorial
%    Algorithms, Part I, Boston: Addison-Wesley, 2011.
%
%July 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%The first combo is 0:n.
if(nargin>1)
    %If the first combination in the sequence is desired.
    if(nargin<3||isempty(startVal))
       startVal=0; 
    end

    n=param1;
    if(n==t)
        theCombo=(0:(n-1)).'+startVal;
        
        recurData=[];
        return
    end

    c=zeros(t+2,1);

    %Step T1
    for j=1:t
        c(j)=j-1;
    end
    c(t+1)=n;
    c(t+2)=0;
    j=t;
    
    theCombo=c(1:t)+startVal;
    recurData.t=t;
    recurData.c=c;
    recurData.j=j;
    recurData.x=[];
    recurData.startVal=startVal;
    return
else
    recurData=param1;
    if(isempty(recurData))
        theCombo=[];
        return
    end
    t=recurData.t;
    c=recurData.c;
    j=recurData.j;
    x=recurData.x;
    startVal=recurData.startVal;
end

%The while loop is used so that the break commands all put everything to
%the step where we record the combination.
while(1)    
    if(j>0)
        x=j;
        %Step T6
        c(j)=x;
        j=j-1;
        break; 
    end
    
    %Step T3
    if(c(1)+1<c(2))
        c(1)=c(1)+1;
        break;
    end
    j=2;
    
    %Step T4
    while(1)
        c(j-1)=j-2;
        x=c(j)+1;
        if(x~=c(j+1))
            break;
        end
        j=j+1;
    end
    
    %Step T5
    if(j>t)
        theCombo=[];
        return;
    end
    
    %Step T6
    c(j)=x;
    j=j-1;
    break;
end

theCombo=c(1:t)+startVal;
recurData.t=t;
recurData.c=c;
recurData.j=j;
recurData.x=x;
recurData.startVal=startVal;
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
