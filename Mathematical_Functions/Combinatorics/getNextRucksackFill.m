function [cCombo,dataRecur]=getNextRucksackFill(param1,N)
%%GETNEXTRUCKSACKFILL Get the next combination of items with given positive
%               weights that will fit into a rucksack (knapsack, backpack)
%               with a capacity of N. This goes through all feasible
%               combinations (starting with putting nothing in the bag).
%               Going through all feasible combinations is not an efficient
%               approach to determining the set of items that can be placed
%               into the bag to maximize its weight.
%
%INPUTS: param1 To generate the first combination (which will be the empty
%               set - put nothing in the bag) this is w, an nX1 vector
%               holding the weights of each of the items that can be placed
%               in the bag. All weights>0, but need not be integer. 
%               Otherwise, this is the dadaRecur output from the previous
%               call to this function.
%             N The capacity of the rucksack. N>0.
%
%OUTPUTS: cCombo The tX1 vector of indices of the items in w that form the
%                current combination of elements that fit into the
%                rucksack. An empty matrix is the first set (out nothing
%                in) and is also returned when the final value has been
%                passed.
%      dataRecur A data structure that can be passed to this function on
%                subsequent calls to get the next value in the sequence. If
%                an empty matrix is reurned, then the final value in the
%                sequence has been passed.
%
%This function implements Algorithm F of Section 7.2.1.3 of [1].
%
%EXAMPLE:
%This displays all of the combination indices for a simple problem.
% w=[1;1.5;5;3;1.5];
% N=6;
% [cCombo,dataRecur]=getNextRucksackFill(w,N);
% while(~isempty(dataRecur))
%     cCombo
%     [cCombo,dataRecur]=getNextRucksackFill(dataRecur);
% end
%
%REFERENCES:
%[1] D. E. Knuth, The Art of Computer Programming. Vol. 4A: Combinatorial
%    Algorithms, Part I, Boston: Addison-Wesley, 2011.
%
%October 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%If we are getting the first filling possibility (the empty sack).
if(nargin==2)
    w=param1;%The set of all-positive weights.
    [w,idx]=sort(w,'ascend');

    n=length(w);
    deltaJ=w(2:n)-w(1:(n-1));

    %The extra one at the beginning is for unassigned weights.
    c=zeros(n+1,1);

    %Step F1, Initialize
    t=0;
    c(0+1)=n;
    r=N;
    %The initial combination is the empty set.
    cCombo=[];
    dataRecur.w=w;
    dataRecur.idx=idx;
    dataRecur.deltaJ=deltaJ;
    dataRecur.c=c;
    dataRecur.r=r;
    dataRecur.t=t;
    return;
else
    dataRecur=param1;
    w=dataRecur.w;
    idx=dataRecur.idx;
    deltaJ=dataRecur.deltaJ;
    c=dataRecur.c;
    r=dataRecur.r;
    t=dataRecur.t;
end

%Step F3, try to add w(0).
if(c(t+1)>0&&r>=w(0+1))
    t=t+1;
    c(t+1)=0;
    r=r-w(0+1);

    %Step F2, visit the current combination.
    cCombo=idx(c(2:(t+1))+1);
    dataRecur.t=t;
    dataRecur.c=c;
    dataRecur.r=r;
    return;
end

while(1)
    if(t==0)%We have passed the last combination.
        cCombo=[];
        dataRecur=[];
        return;
    end

    if(c(t-1+1)>c(t+1)+1&&r>=deltaJ(c(t+1)+1))
        %Step F4
        c(t+1)=c(t+1)+1;
        r=r-deltaJ(c(t+1));
        %Step F2, visit the current combination.
        cCombo=idx(c(2:(t+1))+1);
        dataRecur.t=t;
        dataRecur.c=c;
        dataRecur.r=r;
        return
    else
        %Step F5
        r=r+w(c(t+1)+1);
        t=t-1;
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
