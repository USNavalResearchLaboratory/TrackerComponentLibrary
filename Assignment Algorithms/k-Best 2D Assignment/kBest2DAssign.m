function [col4rowBest,row4colBest,gainBest]=kBest2DAssign(C,k,maximize)
%%KBEST2DASSIGN Find the k lowest (or highest) cost 2D assignments for the
%               two-dimensional assignment problem with a rectangular cost
%               matrix C.
%
%INPUTS: C A numRowXnumCol cost matrix that does not contain any NaNs and
%          where the largest finite element minus the smallest element is a
%          finite quantity (does not overflow) when performing minimization
%          and where the smallest finite element minus the largest element
%          is finite when performing maximization. Forbidden assignments
%          can be given costs of +Inf for minimization and -Inf for
%          maximization.
%        k The number >=1 of hypotheses to generate. If k is less than the
%          total number of unique hypotheses, then all possible hypotheses
%          will be returned.
% maximize If true, the minimization problem is transformed into a
%          maximization problem. The default if this parameter is omitted
%          is false.
%
%OUTPUTS: col4rowBest A numRowXk vector where the entry in each element
%                     is an assignment of the element in that row to a
%                     column. 0 entries signify unassigned rows.
%         row4colbest A numColXk vector where the entry in each element
%                     is an assignment of the element in that column to a
%                     row. 0 entries signify unassigned columns.
%            gainBest A kX1 vector containing the sum of the values of the
%                     assigned elements in C for all of the hypotheses.
%
%DEPENDENCIES: BinaryHeap.m
%              KeyVal.m
%              MurtyData.m
%
%This is an implementation of Murty's method, which is described in [1].
%The algorithm relies on the existence of a 2D assignment algorithm. The 2D
%assignment algorithm of [2] is used. Additionally, the dual variable
%inheritance methods described in [3] is used to reduce the computational
%complexity of the technique.
%
%Murty's algorithm runs 2D assignment algorithms a number of times with an
%increasing number of constraints. Much of the assignment code is in the
%handle subclass MurtyData. Instances of MurtyData are stored in an ordered
%list implemented using the BinaryHeap class.
%
%EXAMPLE:
%This is the example used in [1]. We show how to turn the col4row and
%row4col outputs into the tuples as shown in the paper.
% C=[7,   51,  52,  87,  38,  60,  74,  66,   0   20;
%    50,  12,   0,  64,   8,  53,   0,  46,  76,  42;
%    27,  77,   0,  18,  22,  48,  44,  13,   0,  57;
%    62,   0,   3,   8,   5,   6,  14,   0   26,  39;
%     0,  97,   0,   5,  13,   0,  41,  31,  62,  48;
%    79,  68,   0,   0,  15,  12,  17,  47,  35,  43;
%    76,  99,  48,  27,  34,   0,   0,   0,  28,   0;
%    0,  20,   9,  27,  46,  15,  84,  19,   3,  24;
%    56,  10,  45,  39,   0,  93,  67,  79,  19,  38;
%    27,   0,  39,  53,  46,  24,  69,  46,  23,   1];
% numSol=2;
% [col4row, row4col, gain]=kBest2DAssign(C,numSol);
% N=size(C,1);%It is a square matrix.
% tuples=zeros(2,N,numSol);
% for curSol=1:numSol
%     for curRow=1:N
%         tuples(1,curRow,curSol)=curRow;
%         tuples(2,curRow,curSol)=col4row(curRow,curSol);
%     end
% end
% %The costs of the two solutions are
% gain
% %The tuples of each solution are:
% tuples(:,:,1)
% tuples(:,:,2)
% %One will see that the gain and tuples match the paper.
%
%REFERENCES:
%[1] K. G. Murty, "An algorithm for ranking all the assignments in order of
%    increasing cost," Operations Research, vol. 16, no. 3, pp. 682-687,
%    May-Jun. 1968.
%[2] D. F. Crouse, "On Implementing 2D Rectangular Assignment Algorithms,"
%    IEEE Transactions on Aerospace and Electronic Systems, vol. 52, no. 4,
%    pp. 1679-1696, Aug. 2016.
%[3] M. L. Miller, H. S. Stone, and J. Cox, Ingemar, "Optimizing Murty's
%    ranked assignment method," IEEE Transactions on Aerospace and
%    Electronic Systems, vol. 33, no. 3, pp. 851-862, Jul. 1997.
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    if(nargin<3)
        maximize=false;
    end

    numRow=size(C,1);
    numCol=size(C,2);
    
%The cost matrix must have all non-negative elements for the assignment
%algorithm to work. This forces all of the elements to be positive. The
%delta is added back in when computing the gain in the end.
    if(maximize==true)
        CDelta=max(max(C));
        C=-C+CDelta;
    else
        CDelta=min(min(C));
        C=C-CDelta;
    end

    didFlip=false;
    if(numRow>numCol)
        C=C';
        temp=numRow;
        numRow=numCol;
        numCol=temp;
        didFlip=true;
    end
    
    col4rowBest=zeros(numRow,k);%The assignment of rows to columns in each of the k-hypotheses.
    row4colBest=zeros(numCol,k);%The assignment of rows to columns in each of the k-hypotheses.
    gainBest=zeros(k,1);%The costs of each of the k-hypotheses.
    
%Step 1: Find the best assignment. It will become the root from which
%everything else splits.
    %MUST PAD IT TO MAKE IT SQUARE, because we are using suboptimal
    %dual variables.
    numPad=numCol-numRow;
    C=[C;zeros(numPad,numCol)];

    LCHyp=MurtyData(C,numRow);

%Check for feasibility.
    if(LCHyp.gainFull==-1)
        col4rowBest=[];
        row4colBest=[];
        gainBest=-1;
        return;
    end
    
    col4rowBest(:,1)=LCHyp.col4rowLCFull(1:numRow);
    row4colBest(:,1)=LCHyp.row4colLCFull;
    gainBest(1)=LCHyp.gainFull;
    
    HypList=BinaryHeap(10,false);
    HypList.insert(LCHyp,[]);
    for curSweep=2:k
        %We have to successively split the LC hypothesis either k times or
        %until splitting is no longer possible and get rid of the
        %hypothesis from the list.
        smallestSol=HypList.deleteTop();
        
        smallestSol.key.split(HypList);
        
        %Save the ordered best solutions.
        smallestSol=HypList.getTop();
        if(HypList.heapSize~=0)
            col4rowBest(:,curSweep)=smallestSol.key.col4rowLCFull(1:numRow);
            row4colBest(:,curSweep)=smallestSol.key.row4colLCFull;
            gainBest(curSweep)=smallestSol.key.gainFull;
        else
           %Reduce the matrix to the number that are actually found.
            col4rowBest=col4rowBest(:,1:(curSweep-1));
            row4colBest=row4colBest(:,1:(curSweep-1));
            gainBest=gainBest(1:(curSweep-1));
            
            break;%All possible hypotheses have been enumerated.
        end
    end
    
    %Free the list.
    HypList.delete();
    
    %Set the indices corresponding to padded values to zero.
    if(numPad>0)
        sel=row4colBest>numRow;
        row4colBest(sel)=0;
    end
    
    %Adjust the gain for the initial offset of the cost matrix.
    if(maximize==true)
        gainBest=-gainBest+CDelta*numRow;
    else
        gainBest=gainBest+CDelta*numRow;
    end
    
    if(didFlip==true)
        temp=row4colBest;
        row4colBest=col4rowBest;
        col4rowBest=temp;
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
