function p=randTComposition(n,t)
%%RANDTCOMPOSITION Generate a random composition of n unlableled items into
%                 t labeled slots.
%
%INPUTS: n The number of items that are composed into slots; n>=1, n>=t.
%        t The number of slots that can hold items; t>=1.
%        
%OUTPUTS: p A tX1 vector holding the random length t composition, whose
%           elements sum to n. Each element is the number of "balls" in
%           that slot.
%
%A random combination is generated using the randCombination function, and
%then the function relation described in Chapter 7.2.1.3 of [1] for mapping
%combinations to compositions is used to obtain a composition.
%
%October 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

n=n-1;
nCombo=n;
tCombo=t-1;

if(tCombo>0)
    c=randCombination(nCombo,tCombo);
    if(isempty(c))
        p=[];
        return;
    end

    %Transform the combination into a valid composition.
    p=zeros(t,1);
    p(1)=c(1)+1;
    for curIdx=2:tCombo
        p(curIdx)=c(curIdx)-c(curIdx-1);
    end
    p(t)=n-c(tCombo);
else%The t=1 case.
    p=n+1;
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
