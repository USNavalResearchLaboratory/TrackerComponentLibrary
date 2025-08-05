function aList=genAllBinCombinations(n,t,algorithm)
%%GENALLBINARYCOMBINATIONS Generate all combinations of t items chosen from
%              a set of n items, presenting the results as binary strings.
%
%INPUTS: n The number of items from which t items are chosen .
%        t The number of items chosen, t<=n.
% algorithm An optional parameter affecting which algorithm is used (and
%          thus the ordering of the codes). Possible values are:
%          0 (The default if omitted or an empty matrix is passed). Use
%            Algorithm C in Chapter 7.2.1.3 of [1]. The combinations are
%            given in the order of Chase's sequence.
%          1 Use Algorithm 2 in [2]. The combinations are given as a gray
%            code. Consecutive combinations differ in only 2 bits.
%
%OUTPUTS: aList An nXnumCombos list of all of the binary strings of n bits
%               choosing t of them to be 1. 
%
%There is a total of numCombos=binomial(n,t) combinations.
%
%EXAMPLE:
% aList0=genAllBinCombinations(4,2,0)
% aList1=genAllBinCombinations(4,2,1)
%One will find that
% aList0=[0, 1, 0, 0, 1, 1;
%         0, 0, 1, 1, 0, 1;
%         1, 0, 0, 1, 1, 0;
%         1, 1, 1, 0, 0, 0];
% aList1=[1, 0, 1, 0, 0, 1;
%         1, 1, 0, 0, 1, 0;
%         0, 1, 1, 1, 0, 0;
%         0, 0, 0, 1, 1, 1];
%
%REFERENCES:
%[1] D. E. Knuth, The Art of Computer Programming. Vol. 4A: Combinatorial
%    Algorithms, Part I, Boston: Addison-Wesley, 2011.
%[2] J. R. Bitner, G. Ehrlich, and E. M. Reingold, "Efficient generation of
%    the binary reflected gray code and its applications," Communications
%    of the ACM, vol. 19, no. 9, pp. 517-521, Sep. 1976.
%
%October 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(algorithm))
    algorithm=0;
end

switch(algorithm)
    case 0
        aList=genAllBinCombinationsKnuth(n,t);
    case 1
        aList=genAllBinCombinationsGray(n,t);
    otherwise
        error('Unknown algorithm specified.')
end

end

function aList=genAllBinCombinationsKnuth(n,t)
%%GENALLBINARYCOMBINATIONS Generate all combinations of t items chosen from
%              a set of n items, presenting the results as binary strings.
%
%INPUTS: n The number of items from which t items are chosen.
%        t The number of items chosen, t<=n.
%
%OUTPUTS: aList An nXnumCombos list of all of the binary strings of n bits
%               choosing t of them to be 1. 
%
%There is a total of numCombos=binomial(n,t) combinations. Algorithm C in
%Chapter 7.2.1.3 of [1] is used to generate the strings in the order of
%Chase's sequence.
%
%EXAMPLE:
% aList=genAllBinCombinations(4,2)
%One will find that
% aList=[0,1,0,0,1,1;
%        0,0,1,1,0,1;
%        1,0,0,1,1,0;
%        1,1,1,0,0,0];
%
%REFERENCES:
%[1] D. E. Knuth, The Art of Computer Programming. Vol. 4A: Combinatorial
%    Algorithms, Part I, Boston: Addison-Wesley, 2011.
%
%October 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.

if(t==n)
    aList=ones(n,1);
    return
elseif(t==0)
    aList=zeros(n,1);
    return
end

s=n-t;

numCombos=binomial(n,t);
aList=zeros(n,numCombos);

%Step C1. a and w are indexed from 0 in [1].
a=[zeros(s,1);ones(t,1)];%
w=ones(n+1,1);
r=s;

for curCombo=1:numCombos
    %Step C1
    aList(:,curCombo)=a;
   
    %Step C3
    j=r;
    while(w(j+1)==0)
        w(j+1)=1;
        j=j+1;
        if(j==n)
            return;
        end
    end
    
    w(j+1)=0;
    if(mod(j,2)==1)%j is odd
        if(a(j+1)~=0)
            %Step C4 
            a(j-1+1)=1;
            a(j+1)=0;
            if(r==j&&j>1)
                r=j-1;
            elseif(r==j-1)
                r=j;
            end
        else
            %Step C7
            if(a(j-1+1)~=0)
                %Step C6
                a(j+1)=1;
                a(j-1+1)=0;
                if(r==j&&j>1)
                    r=j-1;
                elseif(r==j-1)
                    r=j;
                end
            else
               a(j+1)=1;
               a(j-2+1)=0;
               if(r==j-2)
                   r=j;
               elseif(r==j-1)
                   r=j-2;
               end
            end
        end
    else%j is even
        if(a(j+1)~=0)
            %Step C5
            if(a(j-2+1)~=0)
                %Step C4
                a(j-1+1)=1;
                a(j+1)=0;
                if(r==j&&j>1)
                    r=j-1;
                elseif(r==j-1)
                    r=j;
                end
            else
                a(j-2+1)=1;
                a(j+1)=0;
                if(r==j)
                    r=max(j-2,1); 
                elseif(r==j-2)
                    r=j-1;
                end
            end
        else
            %Step C6
            a(j+1)=1;
            a(j-1+1)=0;
            if(r==j&&j>1)
                r=j-1;
            elseif(r==j-1)
                r=j;
            end
        end
    end

end

end

function theCodes=genAllBinCombinationsGray(n,k)
%%GENALLBINCOMBINATIONSGRAY Generate all combinations of k items chosen
%              from a set of n items, presenting the results as binary
%              strings. Subsequent combinations differ by 2 bits.
%
%INPUTS: n The number of items from which k items are chosen.
%        t The number of items chosen, t<=n.
%
%OUTPUTS: aList An nXnumCombos list of all of the binary strings of n bits
%               choosing k of them to be 1. 
%
%There is a total of numCombos=binomial(n,t) combinations. Algorithm 2 in
%[1] is used to generate all combinations are a gray code. Subsequent
%combinations differ in only 2 bits.
%
%REFERENCES:
%[1] J. R. Bitner, G. Ehrlich, and E. M. Reingold, "Efficient generation of
%    the binary reflected gray code and its applications," Communications
%    of the ACM, vol. 19, no. 9, pp. 517-521, Sep. 1976.
%
%July 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.

if(k==0)
    theCodes=zeros(n,1);
    return;
end

numCodes=binomial(n,k);
theCodes=zeros(n,numCodes);

g=[ones(k,1);zeros(n-k+2,1)];
tau=2:(n+2);
t=k;
tau(1)=k+1;
i=0;

curGrayCode=1;
while(i<n+1)
    theCodes(:,curGrayCode)=g(1:n);
    curGrayCode=curGrayCode+1;
    
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
