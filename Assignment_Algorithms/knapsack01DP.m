function [valuesSelected,costVal]=knapsack01DP(v,w,W)
%%KNAPSACK01DP Solve the 0-1 knapsack problem using dynamic programming.
%        The 0-1 knapsack problem is to find a subset sel of the items in v
%        that solves the optimization problem
%        max_{sel} sum(v(sel))
%        such that sum(w(sel))<=W
%        This can be thought of as fitting the most things in a knapsack if
%        the v's are all 1's, the w's are the sizes of things and W is the
%        capacity of the knapsack.
%
%INPUTS: v An nX1 or 1Xn vector of positive real values to maximize.
%        w An nX1 or 1Xn vector of positive integer capacities required for
%          each value.
%        W The positive integer maximum allowable capacity that can be
%          filled.
%
%OUTPUTS: valuesSelected The indices of the items in v (and w) that are
%                  selected. If the problem is infeasible, then an empty
%                  matrix is returned.
%          costVal The sum of the selected items in v. 0 is returned if the
%                  problem is infeasible.
%
%The knapsack problem is NP-hard. The complexity of this algorithm
%scaled with W. It is O(n*W) The algorithm is sovled via dynamic
%programming. A development of the recursion in such a dynamic programming
%approach is given in Chapter 5.4 of [1]. However, this is a very common
%method for solving this problem and thus solutions can be found in a
%number of basic textbooks on algorithms.
%
%EXAMPLE:
% W=7;
% v=[2.1;3.125;66;4.3;6.2];
% w=[1;2;3;4;5];
% [valuesSelected,costVal]=knapsack01DP(v,w,W)
%One will get vales 3, 2, and 1 selected. 
%
%REFERENCES:
%[1] L. A. Wolsey, Integer Programming. New York: Wiley-Interscience, 1998.
%
%January 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.
    
n=length(v);

V=zeros(n+1,W+1);
V(0+1,:)=0;
keep=zeros(n,W+1);

for i=1:n
    for wIdx=0:W
        if(wIdx-w(i)>=0)
            VTest=V(i-1+1,wIdx-w(i)+1);
        else
            VTest=-Inf;
        end
        
        if((w(i)<=wIdx)&&(v(i)+VTest)>V(i-1+1,wIdx+1))
            V(i+1,wIdx+1)=v(i)+VTest;
            keep(i,wIdx+1)=1;
        else
            V(i+1,wIdx+1)=V(i-1+1,wIdx+1);
            keep(i,wIdx+1)=0;
        end
    end
end

K=W;

valuesSelected=zeros(n,1);
numSelected=0;

for i=n:-1:1
    if(keep(i,K+1)==1)
        numSelected=numSelected+1;
        valuesSelected(numSelected)=i;
        K=K-w(i); 
    end
end

valuesSelected=valuesSelected(1:numSelected);

costVal=V(n+1,W+1);
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
