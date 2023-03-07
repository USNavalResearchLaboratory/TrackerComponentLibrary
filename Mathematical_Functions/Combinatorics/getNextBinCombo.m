function [curCombo,recurData]=getNextBinCombo(param1,k)
%%GETNEXTBINCOMBO Get the next all combinations of k items chosen from a
%              set of n items, presenting the results as binary strings.
%              A Gray code algorithm is used and subsequent combinations
%              differ by 2 bits.
%
%INPUTS: param1 To get the first combination, param1, should be n, The
%               number of items from which k items are chosen. Otherwise,
%               to get the next combination in the sequence, k should not
%               be provided and param1 should be the recurData output from
%               the previous call to this function.
%             k The number of items chosen, k<=n. This should only be
%               provided when the first combination in the sequence is
%               desired.
%
%OUTPUTS: curCombo The nX1 array of n bits choosing k of them to be 1.
%        recurData A data structure that can be passed back to this
%                  function to get subsequent binary combinations.
%
%Algorithm 2 in [1] is used. There is a total of binomial(n,k)
%combinations.
%
%EXAMPLE:
%We demonstrate that this produces the same results as
%genAllBinCombinations when the appropriate algorithm is selected.
% n=7;
% k=6;
% 
% numCombos=binomial(n,k);
% allRecurCombos=zeros(n,numCombos);
% 
% [curCombo,recurData]=getNextBinCombo(n,k);
% comboIdx=1;
% while(~isempty(curCombo))
%     allRecurCombos(:,comboIdx)=curCombo;
%     comboIdx=comboIdx+1;
%     [curCombo,recurData]=getNextBinCombo(recurData);
% end
% allCombos=genAllBinCombinations(n,k,1);
% assert(all(allCombos(:)==allRecurCombos(:)))
%The assertion should be without error, because the values should be equal.
%
%REFERENCES:
%[1] J. R. Bitner, G. Ehrlich, and E. M. Reingold, "Efficient generation of
%    the binary reflected gray code and its applications," Communications
%    of the ACM, vol. 19, no. 9, pp. 517-521, Sep. 1976.
%
%July 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin==2)
    %Return the first combination.
    n=param1;
    g=[ones(k,1);zeros(n-k+2,1)];
    tau=2:(n+2);
    t=k;
    tau(1)=k+1;
    i=0;
    
    %Special case.
    if(k==0)
        i=n+1;
        curCombo=g(1:n);
        recurData.g=g;
        recurData.tau=tau;
        recurData.t=t;
        recurData.i=i;
        recurData.n=n;
        return;
    end
else
    g=param1.g;
    tau=param1.tau;
    t=param1.t;
    i=param1.i;
    n=param1.n;
end

%If we have passed the final combination.
if(i>n)
    curCombo=[];
    recurData=param1;
    return
end

curCombo=g(1:n);

i=tau(1);
tau(1)=tau(i);
tau(i)=i+1;

if(g(i)==1)
    if(t~=0)
        g(t)=1-g(t);
    else
        g(i-1)=1-g(i-1);
    end
    t=t+1;
else
    if(t~=1)
        g(t-1)=1-g(t-1);
    else
        g(i-1)=1-g(i-1);
    end
    t=t-1;
end

g(i)=1-g(i);

if(t==i-1||t==0)
    t=t+1;
else
    t=t-g(i-1);
    tau(i-1)=tau(1);

    if(t==0)
        tau(1)=i-1;
    else
        tau(1)=t+1;
    end
end

recurData.g=g;
recurData.tau=tau;
recurData.t=t;
recurData.i=i;
recurData.n=n;
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
