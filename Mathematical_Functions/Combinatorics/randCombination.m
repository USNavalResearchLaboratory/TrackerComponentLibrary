function combo=randCombination(n,m,algorithm)
%%RANDCOMBINATION Generate a random combination of m elements drawn from a
%                 set of n elements. Element values start at 0.
%
%INPUTS: n The number of items from which m items are chosen for the
%          combination. n>=1.
%        m The number of items chosen. 1<=m<=n.
% algorithm An optional parameter specifying the algorithm to use. Possible
%          values are:
%          0) Use Floyd's algorithm as described in [1]. This algorithm
%             will work even if binomial(n,m) overflows.
%          1) Generate a random value from 0 to binomial(n,m) and then
%             unrank the corresponding combination. with the
%             unrankCombination function. This will not perform properly if
%             binomial(n,m)>flintmax().
%
%OUTPUTS: combo An mX1 vector containing the combination with values in
%               INCREASING order. The lowest item is indexed zero.
%
%EXAMPLE:
%Here, one can see that a histogram of the indices of the ranked
%combinations from algorithm 0 is approximately uniform.
% numRuns=10000;
% xIdx=zeros(numRuns,1);
% for curRun=1:numRuns
%     combo=randCombination(10,8);
%     xIdx(curRun)=rankCombination(combo);
% end
% figure(1)
% clf
% histogram(xIdx)
%
%REFERENCES:
%[1] J. Bentley and B. Floyd, "Programming pearls: A sample of brilliance,"
%   Communications of the ACM, vol. 30, no. 9, pp. 754-757, Sep. 1987.
%
%October 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(algorithm))
    algorithm=0;
end

switch(algorithm)
     case 0%Floyd's algorithm
        combo=zeros(m,1);
        curIdx=1;
        for j=(n-m+1):n
            t=randi(j);
            %When curIdx==1, this comparison is false.
            if(any(t==combo(1:(curIdx-1))))
                combo(curIdx)=j;
            else
                combo(curIdx)=t;
            end
            curIdx=curIdx+1;
        end

        combo=sort(combo-1,'ascend');
    case 1
        totalCombo=binomial(n,m);
        %The min is for the (presumably zero probability) case that the random
        %variable is 1.
        rank=min(fix(rand(1)*totalCombo),totalCombo-1);

        combo=unrankCombination(rank,n,m);

    otherwise
        error('Unknown algorithm specified.')
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
