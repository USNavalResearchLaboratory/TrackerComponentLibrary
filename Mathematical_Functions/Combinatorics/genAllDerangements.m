function derangements=genAllDerangements(n)
%%GENALLDERANGEMENTS Generate all derangements of length n. A derangement
%             is a permutation (of 1:n) such that no value remains in its
%             original spot. There is a total of subfactorial(n)
%             derangements.
%
%INPUTS: n The positive integer dimensionality of the derangements desired.
%
%OUTPUTS: derangements An nXnumDerange matrix of all derangements of n
%                  items. For n<=1, an empty matrix is returned as there
%                  are no derangements.
%
%The algorithm labeled DERANGE2 in [1] is used. However, the recursion has
%been removed.
%
%REFERENCES:
%[1] S. G. Akl, "A new algorithm for generating derangements," BIT
%    Numerical Mathematics, vol. 20, no. 1, pp. 2-7, 1980.
%
%November 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(n<=1)
    derangements=[];
    return;
end

numDerange=subfactorial(n);

derangements=zeros(n,numDerange);

curDerange=1;

p=1:n;
jVec=zeros(n,1);
lVec=zeros(n,1);

curLevel=n;
goingUp=false;
while(curLevel<=n)
    m=curLevel;
    
    if(goingUp)
        %Swap p(jVec(m)) and p(m)
        temp=p(jVec(m));
        p(jVec(m))=p(m);
        p(m)=temp;
        jVec(m)=jVec(m)-1;
        if(jVec(m)<1);
            %Keep going up
            curLevel=curLevel+1;
        else
            %Swap p(jVec(m)) and p(m)
            temp=p(m);
            p(m)=p(jVec(m));
            p(jVec(m))=temp;
            goingUp=false;
            %Go down
            curLevel=curLevel-1;
        end
    else
        if(m~=1)
            if(p(m)==m)
                lVec(m)=m-1;
            else
                lVec(m)=m;
            end
            jVec(m)=lVec(m);
            %Swap p(m) and p(jVec(m))
            temp=p(m);
            p(m)=p(jVec(m));
            p(jVec(m))=temp;
            %Go down
            curLevel=curLevel-1;
        elseif(p(1)~=1)
            derangements(:,curDerange)=p;
            curDerange=curDerange+1;
            goingUp=true;
            curLevel=curLevel+1;
        else
            goingUp=true;
            curLevel=curLevel+1;
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
