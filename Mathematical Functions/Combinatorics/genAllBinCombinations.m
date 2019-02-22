function aList=genAllBinCombinations(n,t)
%%GENALLBINARYCOMBINATIONS Generate all combinations of t items chosen from
%              a set pf n items, presenting the results as binary strings.
%
%INPUTS: n The number of items from which t items are chosen .
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
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

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
