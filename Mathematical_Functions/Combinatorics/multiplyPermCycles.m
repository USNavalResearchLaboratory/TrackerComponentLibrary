function T=multiplyPermCycles(cycles,numInCycle,n)
%%MULTIPLYPERMUTATIONCYCLES Given a list of cycles making up a permutation,
%               combine them together to get the permutation. The cycle
%               notation does not need to be minimal. All cycles will be
%               combined. The cycles have the same format as used by the
%               getPermutationCycles function.
%
%INPUTS: cycles A maxCycleLengthXnumCycles matrix where each cycle is
%               given. The number of elements in each cycle is given by the
%               numInCycle input. Elements that are not part of a cycle
%               should be set to zero. For a length-n permutation, a cycle
%               consists of numbers from 1 to n. Singleton cycles (elements
%               remain in the same place), can be omitted if n is provided.
%    numInCycle A numCyclesX1 or 1XnumCycles vector listing the number of
%               elements in the cycles.
%             n The length of the permutation. If this parameter is omitted
%               or an empty matrix is passed, then the length is inferred
%               by the highest value in a cycles. If the last element(s) in
%               the permutation is(are) a singleton cycle (i.e. it remains
%               in the same spot), and the singleton cycle is not given in
%               cycles, then it will be omitted in the output.
%
%OUTPUTS: T The nX1 permutation corresponding to the cycles.
%
%This function implements Algorithm B in Section 1.3.3 of [1]. It has been
%modifed to handle cycles provided as an array rather than as a string with
%parentheses. See the comments to getPermutationCycles for more information
%on permutation cycles.
%
%EXAMPLE:
%This is n example in Section 1.3.3 of [1]. Rather than using letters,
%numbers are used here.
% cycles=[1,3,6,7,0;
%         2,3,4,0,0;
%         1,5,4,0,0;
%         6,1,4,5,0;
%         2,7,6,1,5]';
% numInCycle=[4;3;3;4;5];
% n=7;
% T=multiplyPermCycles(cycles,numInCycle,n)
%One will find that T=[4;3;5;7;2;6;1], which corresponds to the letters in
%the first column of Table 2 in Section 1.3.3 of [1].
%
%REFERENCES:
%[1] D. E. Knuth, The Art of Computer Programming. Vol. 1: Fundamental
%    Algorithms, 3rd Edition, Boston: Addison-Wesley, 2011.
%
%November 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

j=0;
numCycles=length(numInCycle);
curCycle=numCycles;
cycleIdx=numInCycle(curCycle);

nIn=max(cycles(:));

if(nargin<3)
    n=nIn;
end

if(n<nIn)
    error('n is inconsistent with the values in cycles.')
end

%Step B1
T=(1:n).';

%Step B2
while(1)
    if(numInCycle(curCycle)==cycleIdx)
        %Closed parenthesis (starting a new cycle).
        Z=0;
    end

%Step B3
    i=cycles(cycleIdx,curCycle);
    temp=Z;
    Z=T(i);
    T(i)=temp;
    if(T(i)==0)
        j=i;
        
        cycleIdx=cycleIdx-1;
        if(cycleIdx==0)
            T(j)=Z;%Finished a cycle (open parenthesis).
            curCycle=curCycle-1;
            if(curCycle==0)
                break;
            end
            cycleIdx=numInCycle(curCycle);
        end
        continue;
    end
    
%Step B4
    T(j)=Z;
    cycleIdx=cycleIdx-1;
    if(cycleIdx==0)
        %Finished a cycle (open parenthesis).
        curCycle=curCycle-1;
        if(curCycle==0)
            break;
        end
        cycleIdx=numInCycle(curCycle);
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
