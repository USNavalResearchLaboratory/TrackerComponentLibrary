function [optAlign,score]=alignSequences(s,t,scores,algorithm)
%%ALIGNSEQUENCES Given two sequences, find the best alignment between them.
%                This means lining up the sequences taking into account the
%                possibility of insertions, deletions and changing values.
%                The result is the set of the two sequences aligned, with
%                insertions and deletions flagged by NaNs. Different cost
%                functions can be specified for penalizing insertions/
%                deletions.
%
%INPUTS: s An mX1 or 1Xm sequence to align with t. This cannot contain
%          NaNs.
%        t An nX1 or 1Xn sequence to align with s. This cannot contain
%          NaNs.
%   scores An optional function handle or 4X1 vector that affects how
%          mismatches are penalized. If a vector is passe, then this
%          contains the following values when comparing two characters in
%          the sequences or a character in the sequence with an empty
%          character (always a NaN):
%          scores(1) For x==y, this is usually positive.
%          scores(2) for x==NaN, this is usually negative.
%          scores(3) for y==NaN, this is usually negative.
%          scores(4) for x~=y, this is usually negative.
%          Alternatively, if scores is a function handle, it should take
%          two arguments and return the score. The ability to pass a
%          function handle allows one to make a function that scores
%          mismatches differently depending on what the mismatched
%          characters are. If this parameter is omitted or an empty matrix
%          is passed, the default is to use scores=[2;-2;-2;-1];
%  algorithm An optional parameter specifying the algorithm to use to align
%          the sequences. Possible values are:
%          0 (The default if omitted or an empty matrix is passed) Use the
%            Needleman-Wunsch algorithm of [1]. The complexity and memory
%            requirements scale as O(m*n).
%          1 Use the Hirschberg algorithm. The ideas from the algorithm
%            come from [2], though there is not a clear description. A
%            simple description is given in [3], though the steps for the
%            algorithm contain a few mistakes. The implementation used here
%            is recursive.
%
%OUTPUTS: optAlign The optimal alignment of sequences. optAlign(1,:) is s
%            with gaps inserted for where t has a character that s does
%            not. optAlign(2,:) is t with similar gaps inserted. When
%            considering editing s, the gaps in optAlign(1,:) can be viewed
%            as insertions (add the character form t), and the gaps in
%            optAlign(2,:) can be viewed as deletions (remove the character
%            from s. Instances where the characters in s and t are
%            different represent changing characters to edit one string
%            into the other. optAlign consists of doubles, so gaps are
%            NaNs. If s and t were character strings, one can do
%            char(optAlign) to make the result character strings, and the
%            NaNs become ASCII zero characters, which Matlab displays as
%            spaces.
%      score The score of the alignment. This is the sum of the costs of
%            all of the matches, mismatches, insertions and deletions.
%
%Both of the algorithms used are forms of dynamic programming. Compare the
%type of score produced by this function to the edit distance in
%findEditDistance.
%
%EXAMPLE 1:
% s='acgtgtcaacgt';
% t='acgtcgtagcta';
% [optAlign,score]=alignSequences(s,t,[],0)
% [optAlign1,score1]=alignSequences(s,t,[],1)
% %To display the modified strings as text again, use
% optAlign=char(optAlign)
% optAlign1=char(optAlign1)
%Using the two different algorithms, one gets different optimally aligned
%sequences, but the scores are the same (9). optAlign is
%['acgt gtcaacgt';'acgtcgt agcta'], where instead of spaces a null
%character is used, and optAlign1 is ['acgt gtcaacgt ';'acgtcgt agc ta'];
%
%EXAMPLE 2:
% s='agtacgca'
% t='tatgc'
% [optAlign,score]=alignSequences(s,t,[],0)
% [optAlign1,score1]=alignSequences(s,t,[],1)
% %To display the modified strings as text again, use
% optAlign=char(optAlign)
% optAlign1=char(optAlign1)
%In this instace, both of the solutions are the same with a score of 1. The
%solutions are ['agtacgca';'  tatgc '];
%
%REFERENCES:
%[1] S. B. Needleman and C. D. Wunsch, "A general method applicable to the
%    search for similarities in the amino acid sequenced of two proteins,"
%    Journal of Molecular Biology, vol. 48, no. 3, pp. 443-453, 28 Mar.
%    1970.
%[2] D. S. Hirschberg, "A linear space algorithm for computing maximal
%    common subsequences," Communications of the ACM, vol. 18, no. 6, pp.
%    341-343, Jun. 1975.
%[3] E. Basic, "Divide-and-conquer sequence alignment within a specific
%    band," 2005, thesis for a diploma in computer science, Lund
%    University, Sweden. [Online].
%    Available: http://lup.lub.lu.se/student-papers/record/1326454
%
%December 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4||isempty(algorithm))
   algorithm=0; 
end

if(nargin>2&&~isempty(scores))
    if(isa(scores,'function_handle'))
        scoreFun=scores;
    else%Assume is a 3X1 or 1X3 matrix of scores.
        scoreFun=@(x,y)twoElScore(x,y,scores);
    end
else%Default score function.
    costList(1)=2;
    costList(2)=-2;%Must be negative or zero
    costList(3)=-2;%Must be negative or zero
    costList(4)=-1;%Must be negative or zero
    %costList(1) for x==y
	%costList(2) for x==NaN
    %costList(3) for y==NaN
    %costList(4) for x~=y

    scoreFun=@(x,y)twoElScore(x,y,costList);
end

switch(algorithm)
    case 0%The Needleman-Wunsch algorithm
        [optAlign,score]=NeedlemanWunschAlign(s,t,scoreFun);
    case 1%The Hirschberg algorithm
        m=length(s);
        n=length(t);
        maxLen=m+n+1;

        %To hold the optimal alignment.
        optAlign=zeros(2,maxLen);
        %prefScores and suffScores are temporary buffers
        prefScores=zeros(n+1,1);
        suffScores=zeros(n+1,1);
        [optAlign,score,pos]=HirschbergAlign(s,t,1,m,1,n,prefScores,suffScores,optAlign,0,0,scoreFun);
        %Shrink to fit.
        optAlign=optAlign(:,1:pos);
    otherwise
        error('Unknown algorithm specified.')
end
end

function [optAlign,score]=NeedlemanWunschAlign(s,t,scoreFun)
%Non-recursive version
    A=makeSimilarityMat(s,t,scoreFun);
    score=A(end,end);
    
    m=length(s);
    n=length(t);
    maxLen=m+n+1;
    %To hold the optimal alignment.
    optAlign=zeros(2,maxLen);
    
    %These level values are used to unwind the recursion in the
    %Needleman-Wunsch algorithm. This is necessary, because very large
    %sequences can make the recursion go too deep.
    iLevels=zeros(maxLen,1);
    jLevels=zeros(maxLen,1);
    %This will hold a number indicating which area the function was in when
    %it went down in the recursion so that the correct option can be chosen
    %when coming up from the recursion.
    caseLevels=zeros(maxLen,1);
    
    pos=0;
    descending=true;
    curLevel=1;
    iLevels(1)=m;
    jLevels(1)=n;
    
    while(curLevel>0)
        i=iLevels(curLevel);
        j=jLevels(curLevel);
    
        if(descending)
            %If going down in the recursion
            if(i==0&&j==0)
                %Innermost recursion; the algorithm is finished. It is time
                %to ascend the levels.
                descending=false;
                curLevel=curLevel-1;
                continue;
            elseif(i>0&&A(i+1,j+1)==A(i+1-1,j+1)+scoreFun(s(i),NaN))
                caseLevels(curLevel)=0;
                iLevels(curLevel+1)=i-1;
                jLevels(curLevel+1)=j;
            elseif(i>0&&j>0&&A(i+1,j+1)==A(i+1-1,j+1-1)+scoreFun(s(i),t(j)))
                caseLevels(curLevel)=1;
                iLevels(curLevel+1)=i-1;
                jLevels(curLevel+1)=j-1;
            else%Must be j>0 and A(i+1,j+1)==A(i+1,j+1-1)+scoreFun(NaN,t(j))
                caseLevels(curLevel)=2;
                iLevels(curLevel+1)=i;
                jLevels(curLevel+1)=j-1;
            end
            curLevel=curLevel+1;
        else
            %If coming up from the recursion.
            if(caseLevels(curLevel)==0)%(i>0&&A(i+1,j+1)==A(i+1-1,j+1)+scoreFun(s(i),NaN))
                pos=pos+1;
                optAlign(1,pos)=s(i);
                optAlign(2,pos)=NaN;
            elseif(caseLevels(curLevel)==1)%(i>0&&j>0&&A(i+1,j+1)==A(i+1-1,j+1-1)+scoreFun(s(i),t(j)))
            	pos=pos+1;
                optAlign(1,pos)=s(i);
                optAlign(2,pos)=t(j);
            else%j>0 and A(i+1,j+1)==A(i+1,j+1-1)+scoreFun(NaN,t(j))
            	pos=pos+1;
                optAlign(1,pos)=NaN;
                optAlign(2,pos)=t(j);
            end
            curLevel=curLevel-1;
        end
    end
    
    %Shrink to fit.
    optAlign=optAlign(:,1:pos);
end

function A=makeSimilarityMat(s,t,scoreFun)
%%Create the score matrix for the Needleman-Wunsch algorithm.
    m=length(s);
    n=length(t);
    A=zeros(m+1,n+1);
    
    A(1,1)=0;
    for j=1:n
        A(1,j+1)=A(1,j+1-1)+scoreFun(NaN,t(j));
    end

    for i=1:m
        A(i+1,1)=A(i+1-1,1)+scoreFun(s(i),NaN);
        
        for j=1:n
            val1=A(i+1-1,j+1)+scoreFun(s(i),NaN);%Vertical
            val2=A(i+1-1,j+1-1)+scoreFun(s(i),t(j));%Diagonal
            val3=A(i+1,j+1-1)+scoreFun(NaN,t(j));%Horizontal
            
            A(i+1,j+1)=max([val1;val2;val3]);
        end
    end
end

function A=HirschbergScore(s,t,sIdxSpan,tIdxSpan,isReversed,scoreFun)
%%Compute the final score row for the Hirschberg algorithm. A(end) is the
% score. sIdxSpan and tIdxSpan are so that one can work with subsets of the
% strings s and t.
%
%The isReversed flag indicates whether the values in sidxSpan and tIdxSpan
%of s and t are backwards, meaning we should index them differently. This
%does not change the final score, but it does change the values in A.

    m=sIdxSpan(2)-sIdxSpan(1)+1;
    n=tIdxSpan(2)-tIdxSpan(1)+1;
    A=zeros(n+1,1);
    
    A(1)=0;
    
    %jIdx is to deal with using only a subset of t.
    if(isReversed==false)
        jIdx=tIdxSpan(1)-1;
    else
        jIdx=tIdxSpan(2)+1;
    end
    for j=1:n
        if(isReversed==false)
            jIdx=jIdx+1;
        else
            jIdx=jIdx-1;
        end
        A(j+1)=A(j+1-1)+scoreFun(NaN,t(jIdx));
    end

    if(isReversed==false)
        iIdx=sIdxSpan(1)-1;
    else
        iIdx=sIdxSpan(2)+1;
    end
    for i=1:m
        if(isReversed==false)
            iIdx=iIdx+1;
        else
            iIdx=iIdx-1;
        end
        
        oldVal=A(1);
        %iIdx is to deal with using only a subset of s.
        
        A(1)=A(1)+scoreFun(s(iIdx),NaN);
        
        if(isReversed==false)
            jIdx=tIdxSpan(1)-1;
        else
            jIdx=tIdxSpan(2)+1;
        end
        for j=1:n
            if(isReversed==false)
                jIdx=jIdx+1;
            else
                jIdx=jIdx-1;
            end
            
            temp=A(j+1);
            
            val1=A(j+1)+scoreFun(s(iIdx),NaN);%Vertical
            val2=oldVal+scoreFun(s(iIdx),t(jIdx));%Diagonal
            val3=A(j+1-1)+scoreFun(NaN,t(jIdx));%Horizontal
            
            A(j+1)=max([val1;val2;val3]);
            
            oldVal=temp;
        end
    end
end

function [optAlign,score,pos]=HirschbergAlign(s,t,a,b,c,d,prefScores,suffScores,optAlign,score,pos,scoreFun)

    if(b-a+1==0)
        jIdx=c-1;
        for j=1:(d-c+1)
            jIdx=jIdx+1;
            pos=pos+1;
            optAlign(1,pos)=NaN;
            optAlign(2,pos)=t(jIdx);
            score=score+scoreFun(NaN,t(jIdx));
        end
    elseif(d-c+1==0)
        iIdx=a-1;
        for i=1:(b-a+1)
            iIdx=iIdx+1;
            pos=pos+1;
            optAlign(1,pos)=s(iIdx);
            optAlign(2,pos)=NaN;
            score=score+scoreFun(s(iIdx),NaN);
        end
    elseif(b-a+1==1||d-c+1==1)
        [optAlignCur,scoreCur]=NeedlemanWunschAlign(s(a:b),t(c:d),scoreFun);
        num2Add=size(optAlignCur,2);
        optAlign(1,(pos+1):(pos+num2Add))=optAlignCur(1,:);
        optAlign(2,(pos+1):(pos+num2Add))=optAlignCur(2,:);
        pos=pos+num2Add;
        score=score+scoreCur;
    else
        i=a+fix((b-a)/2);%Middle index in s
        %Compute prefix and suffix scores
        prefScores((c+1-1):(d+1))=HirschbergScore(s,t,[a;i],[c;d],false,scoreFun);
        suffScores((c+1-1):(d+1))=HirschbergScore(s,t,[i+1;b],[c;d],true,scoreFun);

        %Determine the optimal index posMax in t to partition. In this
        %instance, we are going through suffScores in reversed order. That
        %makes the indexation a bit tricky.
        posMax=c-1;
        vMax=prefScores(c+1-1)+suffScores(d+1);
        jIdxRev=d+1;
        for j=c:d
            jIdxRev=jIdxRev-1;
            vCur=prefScores(j+1)+suffScores(jIdxRev);
            if(vCur>vMax)
               vMax=vCur;
               posMax=j;
            end
        end

        %Recursively solve the two subproblems.
        [optAlign,score,pos]=HirschbergAlign(s,t,a,i,c,posMax,prefScores,suffScores,optAlign,score,pos,scoreFun);
        [optAlign,score,pos]=HirschbergAlign(s,t,i+1,b,posMax+1,d,prefScores,suffScores,optAlign,score,pos,scoreFun);
    end
end

function val=twoElScore(x,y,costList)
%TWOELSCORE A very basic generic score function.
    if(x==y)
        val=costList(1);
    elseif(isnan(x))
        val=costList(2);
    elseif(isnan(y))
        val=costList(3);
    else
        val=costList(4);
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
