function piVals=genAllPureInvolutions(n)
%GENALLPUREINVOLUTIONS Generate all pure involutions of the numbers
%  1,2,...n. A pure involution is one having no fixed points. That is none
%  of the values in the involution is in its original place. There are
%  numPureInvolutions(n)  pure involutions. Involutions are permutations
%  where the permutation is its own inverse permutation --so rearranging
%  the elements of a vector by the same permutation twice results in the
%  original permutation being produced. Pure involutions (fixed-point free
%   involutions) are both involutions and derangements and they arise in
%   cryptanalysis, among other areas.
%
%INPUTS: n The length of the involution n>=0. If 0, an empty matrix is
%          returned.
%
%OUTPUTS: piVals An nXnumPureInvolutions(n) matrix holding all of the
%                involutions. Note that if n is odd an empty matrix is
%                returned, because there are no pure involutions for an odd
%                n.
%
%This function implements Algorithm 5 of [1].
%
%EXAMPLE:
%Here, all pure involutions for a particular n are generated and they are
%verified as being pure involutions by calling the isPureInvolution on all
%of them. 
% n=12;
% allPureInv=genAllPureInvolutions(n);
% for k=1:numPureInvolutions(n)
%     assert(isPureInvolution(allPureInv(:,k)))
% end
%
%REFERENCES:
%[1] V. Vajnovszki, "Generating involutions, derangements, and relatives
%    by ECO," Discrete Mathematics and Theoretical Computer Science,
%    vol. 12, no. 1, pp. 109-122, Jan. 2010. [Online]. Available:
%    https://dmtcs.episciences.org/479/pdf
%
%September 2024 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(mod(n,2)==1)
    %There are no pure involutions when n is odd.
    piVals=zeros(n,0);
    return;
end

numPureInv=numPureInvolutions(n);
piVals=zeros(n,numPureInv);
curInv=0;
t=0;
piVal=(1:n).';
i=zeros(n,1);
sigmaInv=zeros(n,n);

isDescending=true;
while(t>=0)
    if(t==n)
        curInv=curInv+1;
        piVals(:,curInv)=piVal;
        isDescending=false;
        t=t-2;
        continue;
    else
        if(isDescending)
            %We have just entered this level going down.
            i(t+1)=0;
            t1=t+1;
            t2=t+2;
            %Swap t+1 and t+2.
            temp=piVal(t2);
            piVal(t2)=piVal(t1);
            piVal(t1)=temp;
            t=t+2;
            continue;
        else
            %If we are going up in the levels.
            if(i(t+1)==0)
                %Undo the previous swap.
                t1=t+1;
                t2=t+2;
                %Swap t+1 and t+2.
                temp=piVal(t2);
                piVal(t2)=piVal(t1);
                piVal(t1)=temp;
            else
                %Undo the permutation with sigmaInv.
                piVal=piVal(sigmaInv(:,t));
            end
            i(t+1)=i(t+1)+1;

            if(i(t+1)>t)
                %Keep ascending.
                t=t-2;
                continue;
            else
                sigmaCur=1:n;
                %Swap i and t+1
                t1=i(t+1);
                t2=t+1;
                temp=sigmaCur(t2);
                sigmaCur(t2)=sigmaCur(t1);
                sigmaCur(t1)=temp;
                %Swap piVal(i) and t+2
                t1=piVal(i(t+1));
                t2=t+2;
                temp=sigmaCur(t2);
                sigmaCur(t2)=sigmaCur(t1);
                sigmaCur(t1)=temp;
                %Swap t+1 and t+2
                t1=t+1;
                t2=t+2;
                temp=sigmaCur(t2);
                sigmaCur(t2)=sigmaCur(t1);
                sigmaCur(t1)=temp;
                piVal=piVal(sigmaCur);
                sigmaInv(:,t)=inversePermutation(sigmaCur);
                t=t+2;
                isDescending=true;
                continue;
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
