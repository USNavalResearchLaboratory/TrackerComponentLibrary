function theCombos=genAllRDCombinations(n,t,startVal)
%%GENALLRDCOMBINATIONS Generate all length t combinations of n elements
%                   (n>=t>1) in a resolving door gray-code order. This is
%                   such that the difference between each successiuve
%                   combination involves changing a single number.
%
%INPUTS: n The number of items from which t items are chosen; n>0.
%        t The number of items chosen, t<=n.
% startVal This is zero or 1, indicating which value the value at which the
%          elements in the combinations can start. The default if omitted
%          or an empty matrix is passed is 0;
%
%OUTPUTS: theCombos A tXnumCombos matrix containing all possible
%                   combinations of values in the "revolving door" order.
%
%This function implements algorithm R of Chapter 7.2.1.3 of [1].
%
%REFERENCES:
%[1] D. E. Knuth, The Art of Computer Programming. Vol. 4A: Combinatorial
%    Algorithms, Part I, Boston: Addison-Wesley, 2011.
%
%October 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(startVal))
    startVal=0; 
end

if(n==t&&t==1)
    theCombos=startVal;
    return
end

numCombos=binomial(n,t);
theCombos=zeros(t,numCombos);

%Step R1, Initialization
c=[(0:(t-1)).';n];

for curCombo=1:numCombos
    %Step R2, visit the combination.
    theCombos(:,curCombo)=c(1:t);

    %Step R3
    if(mod(t,2)==1)
        if(c(1)+1<c(2))
            c(1)=c(1)+1;
            continue;
        else
            j=2;
            skip4=false;
        end
    else
        if(c(1)>0)
            c(1)=c(1)-1;
            continue;
        else
           j=2;
           skip4=true;
        end
    end
    
    while(1)
        if(skip4==false)
            %Step R4
            if(c(j)>=j)
                c(j)=c(j-1);
                c(j-1)=j-2;
                break;
            else
                j=j+1;
            end
        end
        skip4=false;
        
        %Step R5.
        if(c(j)+1<c(j+1))
            c(j-1)=c(j);
            c(j)=c(j)+1;
            break;
        end
        
        j=j+1;
        
        if(j<=t)
            continue;
        else
            if(startVal~=0)
                theCombos=theCombos+1;
            end
 
            return;
        end
    end
end

if(startVal~=0)
    theCombos=theCombos+1;
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
