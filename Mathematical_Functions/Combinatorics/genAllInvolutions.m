function piVals=genAllInvolutions(n)
%%GENALLINVOLUTIONS Generate all involutions of the numbers 1,2,...n. There
%  are numInvolutions(n) involutions. Involutions are permutations where
%  the permutation is its own inverse permutation --so rearranging the
%  elements of a vector by the same permutation twice results in the
%  original permutation being produced.
%
%INPUTS: n The length of the involution n>=0. If 0, an empty matrix is
%          returned.
%
%OUTPUTS: piVals An nXnumInvolutions(n) matrix holding all of the
%                involutions --except for n=0, in which case an empty
%                matrix is returned.
%
%This function implements Algorithm 3 of [1]. However, the implementation
%differs from [1] in that the recursion is avoided. As noted in [2],
%permutation matrices formed from involutions are symmetric.
%
%EXAMPLE:
%This generates all the involutions of length n and then verifies that they
%are actually involutions using the isInvolution function.
% n=11;
% p=genAllInvolutions(n);
% numInv=numInvolutions(n);
% for k=1:numInv
%     assert(isInvolution(p(:,k))==true)
% end
%
%REFERENCES:
%[1] V. Vajnovszki, "Generating involutions, derangements, and relatives
%    by ECO," Discrete Mathematics and Theoretical Computer Science,
%    vol. 12, no. 1, pp. 109-122, Jan. 2010. [Online]. Available:
%    https://dmtcs.episciences.org/479/pdf
%[2] Weisstein, Eric W. "Permutation Involution." From MathWorld--A Wolfram
%    Web Resource. https://mathworld.wolfram.com/PermutationInvolution.html
%
%September 2024 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(n==0)
    piVals=zeros(n,0);
    return
end

numInv=numInvolutions(n);
piVals=zeros(n,numInv);

F=zeros(n,n);
F(1,1)=1;
numInF=zeros(n,1);
numInF(1)=1;
%The index of the loop on each level.
i=zeros(n,1);
isPastLoop=false(n,1);
piVal=(1:n).';
t=1;%The current recursion level.
goingDown=true;
curSol=0;
isPastLoop(1)=false;
while(t>0)
    if(goingDown)
        if(t==n)
            curSol=curSol+1;
            piVals(:,curSol)=piVal;
            goingDown=false;
            t=t-1;
            continue;
        else
            %We just entered this level.
            i(t)=1;
            if(numInF(t)>0)
                %Swap i(t) and t+1;
                idx=F(i(t),t);
                temp=piVal(idx);
                piVal(idx)=piVal(t+1);
                piVal(t+1)=temp;
    
                %Remove the i(t)th element in F for the next level.
                F(1:(numInF(t)-1),t+1)=F([1:(i(t)-1),(i(t)+1):numInF(t)],t);
                numInF(t+1)=numInF(t)-1;
                isPastLoop(t+1)=false;
                t=t+1;
                continue;
            else
                %F is empty. Augment F and go down again.
                isPastLoop(t)=true;
                F(1:numInF(t),t+1)=F(1:numInF(t),t);
                F(numInF(t)+1,t+1)=t+1;
                numInF(t+1)=numInF(t)+1;

                isPastLoop(t+1)=false;
                t=t+1;
                continue;
            end
        end
    else
        %We are going back up a level, so we have to deal with adding back 
        if(isPastLoop(t)==false)
            %We are returning in the loop and thus have to undo the swap. F
            %for this level is unchanged.
            idx=F(i(t),t);
            temp=piVal(idx);
            piVal(idx)=piVal(t+1);
            piVal(t+1)=temp;

            i(t)=i(t)+1;
            if(i(t)<=numInF(t))
                %Redo the loop and go down another level.
                idx=F(i(t),t);
                temp=piVal(idx);
                piVal(idx)=piVal(t+1);
                piVal(t+1)=temp;
    
                %Remove the i(t)th element in F for the next level.
                F(1:(numInF(t)-1),t+1)=F([1:(i(t)-1),(i(t)+1):numInF(t)],t);
                numInF(t+1)=numInF(t)-1;
                isPastLoop(t+1)=false;
                t=t+1;
                goingDown=true;
                continue;
            else
                %We just got past the loop. Augment F and go down again.
                isPastLoop(t)=true;
                F(1:numInF(t),t+1)=F(1:numInF(t),t);
                F(numInF(t)+1,t+1)=t+1;
                numInF(t+1)=numInF(t)+1;

                isPastLoop(t+1)=false;
                t=t+1;
                goingDown=true;
                continue;
            end
        end
        %We are past the loop and already went down a level, so go back up
        %a level.
        t=t-1;
        continue;
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
