function theCombos=genAllCombinations(n,t,startVal,algorithm)
%%GENALLCOMBINAIONS Generate all combinations of length t chosen from n
%              elements. There are binomial(n,t) possible combinations.
%
%INPUTS: n The number of items from which t items are chosen.
%        t The number of items chosen, t<=n.
% startVal This is zero or 1, indicating which value the value at which the
%          elements in the combinations can start. The default if omitted
%          or an empty matrix is passed is 0;
% algorithm An optional parameter specifying the algorithm to use. Possible
%          values are:
%          0 (The default if omitted or an empty matrix is passed) Use
%            Algorithm T of Section 7.2.1.3 of [1]. The combinations are in
%            colexicographic order.
%          1 Use the algorithm of [1]. The combinations are in
%            lexicographic order.
%
%OUTPUTS: theCombos A tXnumCombos matrix containing all possible
%                   combinations of values.
%
%REFERENCES:
%[1] D. E. Knuth, The Art of Computer Programming. Vol. 4A: Combinatorial
%    Algorithms, Part I, Boston: Addison-Wesley, 2011.
%[1] C. J. Mifsud, "Algorithm 154: Combination in lexicographical order," 
%    Communications of the ACM, vol. 6, no. 3 pp. 103, Mar. 1963.
%    modified to start from zero instead of one.
%
%October 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(startVal))
    startVal=0;
end

if(nargin<4||isempty(algorithm))
    algorithm=0;
end

switch(algorithm)
    case 0
        theCombos=genAllCombinations0(n,t,startVal);
    case 1
        theCombos=genAllCombinations1(n,t,startVal);
    otherwise
        error('Unknown algorithm specified.')
end
end

function theCombos=genAllCombinations0(n,t,startVal)
%%GENALLCOMBINATIONS0 Generate all combinations of t items chosen from n in
%                    colexicographic order.
%
%This function implements Algorithm T in Section 7.2.1.3 of [1].
%
%REFERENCES:
%[1] D. E. Knuth, The Art of Computer Programming. Vol. 4A: Combinatorial
%    Algorithms, Part I, Boston: Addison-Wesley, 2011.
%
%October 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.

if(nargin<3||isempty(startVal))
   startVal=0; 
end

if(n==t)
    if(startVal==0)
        theCombos=(0:(n-1)).';
    else
        theCombos=(1:n).';
    end
    return
end

numCombos=binomial(n,t);
theCombos=zeros(t,numCombos);

c=zeros(t+2,1);

%Step T1
for j=1:t
    c(j)=j-1;
end
c(t+1)=n;
c(t+2)=0;
j=t;

for curCombo=1:numCombos
    %Step T2
    theCombos(:,curCombo)=c(1:t);
    
    if(j>0)
        x=j;
        %Step T6
        c(j)=x;
        j=j-1;
        continue; 
    end
    
    %Step T3
    if(c(1)+1<c(2))
        c(1)=c(1)+1;
        continue;
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
        break;
    end
    
    %Step T6
    c(j)=x;
    j=j-1;
end

if(startVal~=0)
    theCombos=theCombos+startVal; 
end
end

function theCombos=genAllCombinations1(n,r,startVal)
%%GENALLCOMBINATIONS Generate all combinations of t items chosen from n in
%                    colexicographic order.
%
%This function implements the algorithm of [1].
%
%REFERENCES:
%[1] C. J. Mifsud, "Algorithm 154: Combination in lexicographical order," 
%    Communications of the ACM, vol. 6, no. 3 pp. 103, Mar. 1963.
%    modified to start from zero instead of one.
%
%October 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.

if(nargin<3||isempty(startVal))
    startVal=0;
end

numCombos=binomial(n,r);
theCombos=zeros(r,numCombos);

if(startVal==0)
    I=(0:(r-1)).';
    theCombos(:,1)=I;
    for curCombo=2:numCombos
        if(I(r)<n-1)
            I(r)=I(r)+1;
            theCombos(:,curCombo)=I;
            continue;
        else
            for j=r:-1:2
               if(I(j-1)<n-r+j-2)
                   I(j-1)=I(j-1)+1;
                   for s=j:1:r
                      I(s)=I(j-1)+s-(j-1); 
                   end
                   theCombos(:,curCombo)=I;
                   break;
               end
            end
        end
    end
else
    I=(1:r).';
    theCombos(:,1)=I;
    for curCombo=2:numCombos
        if(I(r)<n)
            I(r)=I(r)+1;
            theCombos(:,curCombo)=I;
            continue;
        else
            for j=r:-1:2
               if(I(j-1)<n-r+j-1)
                   I(j-1)=I(j-1)+1;
                   for s=j:1:r
                      I(s)=I(j-1)+s-(j-1); 
                   end
                   theCombos(:,curCombo)=I;
                   break;
               end
            end
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
