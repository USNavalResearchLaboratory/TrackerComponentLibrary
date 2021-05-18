function [tuples,gain,u,v]=assign2DFull(C,maximize,algorithm)
%%ASSIGN2DFULL This functions performs the type of 2D assignment that often
%       arises as a subproblem when solving an S-dimensional assignment
%       problem for target tracking. The cost function being minimized (or
%       maximized) is
%       \sum_{i=1}^{numRow}\sum_{j=1}^{numCol}C_{i,j}*x_{i,j}
%       subject to
%       \sum_{j=1}^{numCol}x_{i,j}=1 for i=2:numRow
%       \sum_{i=1}^{numRow}x_{i,j}=1 for j=2:numCol
%       x_{i,j}=0 or 1.
%       Note that the first row and the first column can be multiply
%       assigned and every row and every column aside from the first ones
%       must be assigned to exactly one value.
%
%INPUTS: C A (numRow)X(numCol) cost matrix where the first row is the
%          cost of not assigning each column (from number 2) and the first
%          column is the cost of not assigning each row (from number 2) and
%          C(1,1) is an assignment cost that is unconstrained and only
%          assigned if it improves the cost function (given above). The
%          matrix cannot contain any NaNs and if maximizing, cannot
%          contain any Inf values and if minimizing cannot contain any -Inf
%          values. If algorithm=1 is selected, then it is not allowed that
%          both C(2:end,1) and C(1,2:end) contain non-finite terms.
%          However, either one individually can contain non-finite terms.
%          When algorithm=1 and both C(2:end,1) and C(1,2:end) contain non-
%          finite terms, the algorithm will report that the problem is
%          infeasible, because it cannot solve such problems. 
% maximize If true, the minimization problem is transformed into a
%          maximization problem. The default if this parameter is omitted
%          or an empty matrix is passed is false.
% algorithm A value selecting which algorithm will be used. The basis for
%          the algorithms are discussed in more detail below. Possible
%          values are:
%          0 (The default if omitted or an empty matrix is passed) Use an
%            algorithm that effectively solves an optimization problem on a
%            (numRow-1)+(numCol-1) augmented matrix, though the matrix is
%            not explicitly formed.
%          1 Use an algorithm that removes the first column and modifies
%            the cost matrix so that the optimization problem can be
%            solved using the assign2DMissedDetect function, which
%            effectively solves a problem on a
%            (numRow+numCol-1)X(numCol-1) augmented matrix. Unlike
%            algorithm 0, this algorithm requires that there not be
%            non-finite terms in both the first row and the first column of
%            C. Otherwise, the problem will be labeled as infeasible.
%
%OUTPUTS: tuples A 2XnumTuples set of assignment values. This is ordered
%               [Row Index;Column Index]. the number of tuples ranges from
%               max(numRow-1,numCol-1), for the case where as many things
%               as possible are not assigned to the uncontrained values and
%               C(1,1) is not assigned, to numRow+numCol-1 for the case
%               where all rows and columns are assigned to the
%               unconstrained hypothesis and C(1,1) is provided. If the
%               problem is infeasible, this is an empty matrix.
%          gain This is the value of the cost. This is the sum of the
%               values in C corresponding to the tuples. If
%               the problem is infeasible, this is -1.
%           u,v These are the dual variables for the constrained columns
%               and for constrained rows of C. These will satisfy the
%               complementary slackness condition of the problem as
%               transformed. The transformation differes dependong on
%               whether algorithm 0 or algorithm 1 is chosen. See the
%               examples below.
%
%Algorithm 0 utilizes the equivalence of the minimization for the problem
%with two unconstrained indices to one with no unconstrained indices. The
%problem described above, removing the C(1,1) element, is equivalent to the
%optimization problem 
%           minimize (or maximize)
%          \sum_{i=1}^{numAug}\sum_{j=1}^{numAug}D_{i,j}*x_{i,j}
%          subject to
%          \sum_{j=1}^{numAug}x_{i,j}=1 for all i
%          \sum_{i=1}^{numAug}x_{i,j}=1 for all j
%          x_{i,j}=0 or 1.
%where the matrix D is a numAugXnumAug matrix with numAug=(numRow+numCol-2)
%that can be formed as follows:
% D=Inf(numAug,numAug);
% D(1:(numRow-1),1:(numCol-1))=C(2:end,2:end);
% D(numRow:totalAug,numCol:totalAug)=0;
% for curRow=1:(numRow-1)
%     missedCol=(numCol-1)+curRow;
%     D(curRow,missedCol)=C(curRow+1,1);
% end
% for curCol=1:(numCol-1)
%     missedRow=(numRow-1)+curCol;
%     D(missedRow,curCol)=C(1,curCol+1);
% end
%Basically, C(2:end,2:end) is an upper-left submatrix. The upper-right and
%lower-left submatrix consist of Inf terms with diagonals that are the
%missed detection costs for the rows and columns. The lower-right submatrix
%is all zeros. This function solves the problem by implementing the 2D
%Jonker-Volgenant algorithm of [1] for the augmented problem, but the
%matrix C itself is never augmented.
%
%Algorithm 1 is based on the fact that since every unconstrained column
%must be assigned, one can subtract the cost of the unconstrained row from
%all columns and then perform regular 2D assignment using
%assign2DMissedDetect. One can subtract it, because an assignment must be
%made for each constrained column. Subtracting the unconstrained row costs
%changes the gain, but it does not change the optimal assignment of the
%contrained columns. We find the optimal column assignment with a standard
%2D assignment algorithm. The unassigned rows then just need to be added
%back in as assigned to the unconstrained hypothesis and the gain must be
%adjusted. In general, the cost matrix will not be all-positive, so a
%constant value is added to it to make all of the elements positive. The
%hypothesis C(1,1) (from the unmodified matrix) is only added in if it
%improves the cost function. If the unconstrained column contains non-
%finite terms, then the matrix is transposed and the transposed problem is
%solved, after which the results are flipped back.
%
%EXAMPLE 1:
%Here, we verify that the solution equals the solution obtained through
%linear programming.
% numRow=6;
% numCol=8;
% C=2*rand(numRow,numCol)-1;
% c=C(:);
% numEls=length(c);
% 
% %Equality constraints specifying that each constrained column is assigned 
% %assigned to exactly one row.
% AEqC=[zeros(7,6),kron(eye(7),ones(1,6))];
% bEqC=ones(7,1);
% 
% %Equality constraints specifying that each constrained row is assigned to
% %exactly one column.
% AEqR=[[zeros(5,1),eye(5)],kron(ones(1,7),[zeros(5,1),eye(5)])];
% bEqR=ones(5,1);
% 
% A=[AEqC;AEqR];
% b=[bEqC;bEqR];
% 
% %The only inequality constraint that is needed (and is not needed in a
% %normal 2D assignment problem formulated as a linear programming
% %solution) is the extra contraint on the C(1,1) term. We have to
% %explicitly say that
% %x(1)<=1.
% AIneq=[1,zeros(1,numEls-1)];
% bIneq=1;
% 
% %Solving the problem using linear programming.
% [gainLP,xLP]=linProgRevisedSimplex(A,b,AIneq,bIneq,c);
% 
% %Solving the problem using 
% [tuplesTotal,gain]=assign2DFull(C,false);
% 
% numTuples=size(tuplesTotal,2);
% x=zeros(numRow,numCol);
% for curTuple=1:numTuples
%     x(tuplesTotal(1,curTuple),tuplesTotal(2,curTuple))=x(tuplesTotal(1,curTuple),tuplesTotal(2,curTuple))+1; 
% end
% 
% %One finds that both methods produce the same assignment:
% all(x(:)==xLP(:))
% %Also, it can be verified that gainLP equals gain within certain finite
% %precision bounds. Thus, we have shown that the linear programing solution
% %is equivalent to the generalized 2D assignment solution.
%
%EXAMPLE 2:
%Here we demonstrate that the dual variables satisfy the complementary
%slackness condition for transformed versions of C This examples does this
%for algorithm 0. Example 3 shows how the constant offset needed for
%algorithm 4 differs.  If minimization is performed, one subtracts
%minC=min(0,min(C(:))) from C to get the transformed problem for
%complementary slackness; when performing maximization, subtract
%maxC=max(0,max(C(:))). The complementary slackness condition could be
%forced to hold on C by adding minC or maxC to v. The complementary
%slackness condition is described, for example, in [1].
% numRow=50;
% numCol=30;
% C=randn(numRow,numCol);
% minC=min(0,min(C(:)));
% maximize=false;
% algorithm=0;
% [tuples,gain,u,v]=assign2DFull(C,maximize,algorithm);
% slackVal=0;
% for k=1:size(tuples,2)
%     %If assigned to a constrained value.
%     if(tuples(1,k)~=1&&tuples(2,k)~=1)
%         slackVal=slackVal+abs((C(tuples(1,k),tuples(2,k))-minC)*(C(tuples(1,k),tuples(2,k))-minC-v(tuples(1,k)-1)-u(tuples(2,k)-1)));
%     end
% end
% slackVal
% %One will see that slackVal is zero (within finite precision limits),
% %which is what the complementary slackness condition says.
% %Switching to maximization with the same matrix:
% maximize=true;
% maxC=max(0,max(C(:)));
% [tuples,gain,u,v]=assign2DFull(C,maximize,algorithm);
% for k=1:size(tuples,2)
%     %If assigned to a constrained value.
%     if(tuples(1,k)~=1&&tuples(2,k)~=1)
%         slackVal=slackVal+abs((C(tuples(1,k),tuples(2,k))-maxC)*(C(tuples(1,k),tuples(2,k))-maxC-v(tuples(1,k)-1)-u(tuples(2,k)-1)));
%     end
% end
% slackVal
% %slackVal is again zero, within finite precision limits.
%
%EXAMPLE 3:
%This is the same as example 2, except algorithm 1 is used. One can note
%that the constant needed for the complementary slackness condition to hold
%is different and more complicated to compute and is given by the modC term
%below.
% numRow=50;
% numCol=30;
% C=randn(numRow,numCol);
% maximize=false;
% algorithm=1;
% [tuples,gain,u,v]=assign2DFull(C,maximize,algorithm);
% 
% modC=C;
% modC(2:numRow,2:numCol)=bsxfun(@minus,modC(2:numRow,2:numCol),C(2:numRow,1));
% minModC=min(0,min(vec(modC(1:numRow,2:numCol))));
% %We verify the correctness of the complementary slackness condition.
% slackVal=0;
% for k=1:size(tuples,2)
%     %If assigned to a constrained value
%     if(tuples(1,k)~=1&&tuples(2,k)~=1)
%         slackVal=slackVal+abs((C(tuples(1,k),tuples(2,k))-minModC)*(C(tuples(1,k),tuples(2,k))-minModC-v(tuples(1,k)-1)-u(tuples(2,k)-1)));
%     end
% end
% slackVal
% 
% %One will see that slackVal is zero (within finite precision limits),
% %which is what the complementary slackness condition says.
% %Switching to maximization with the same matrix:
% maximize=true;
% [tuples,gain,u,v]=assign2DFull(C,maximize,algorithm);
% maxModC=max(0,max(vec(modC(1:numRow,2:numCol))));
% for k=1:size(tuples,2)
%     %If assigned to a constrained value.
%     if(tuples(1,k)~=1&&tuples(2,k)~=1)
%         slackVal=slackVal+abs((C(tuples(1,k),tuples(2,k))-maxModC)*(C(tuples(1,k),tuples(2,k))-maxModC-v(tuples(1,k)-1)-u(tuples(2,k)-1)));
%     end
% end
% slackVal
% %slackVal is again zero, within finite precision limits.
%
%REFERENCES:
%[1] D. F. Crouse, "On Implementing 2D Rectangular Assignment Algorithms,"
%    IEEE Transactions on Aerospace and Electronic Systems, vol. 52, no. 4,
%    pp. 1679-1696, Aug. 2016.
%
%January 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    if(nargin<3||isempty(algorithm))
       algorithm=0; 
    end

    if(nargin<2||isempty(maximize))
        maximize=false;
    end
    
    switch(algorithm)
        case 0
            [tuples,gain,u,v]=assign2DFullStd(C,maximize);
        case 1
            [tuples,gain,u,v]=assign2DFullAlt(C,maximize);
        otherwise
            error('Unknown algorithm specified.')
    end
end

function [tuples,gain,u,v]=assign2DFullStd(C,maximize)
%%ASSIGN2DFULL This functions performs the type of 2D assignment that often
%       arises as a subproblem when solving an S-dimensional assignment
%       problem for target tracking. The cost function being minimized (or
%       maximized) is
%       \sum_{i=1}^{numRow}\sum_{j=1}^{numCol}C_{i,j}*x_{i,j}
%       subject to
%       \sum_{j=1}^{numCol}x_{i,j}=1 for i=2:numRow
%       \sum_{i=1}^{numRow}x_{i,j}=1 for j=2:numCol
%       x_{i,j}=0 or 1.
%       Note that the first row and the first column can be multiply
%       assigned and every row and every column aside from the first ones
%       must be assigned to exactly one value.
%
%INPUTS: C A (numRow)X(numCol) cost matrix where the first row is the
%          cost of not assigning each column (from number 2) and the first
%          column is the cost of not assigning each row (from number 2) and
%          C(1,1) is an assignment cost that is unconstrained and only
%          assigned if it improves the cost function (given above). The
%          matrix cannot contain any NaNs and if maximizing, cannot
%          contain any Inf values and if minimizing cannot contain any -Inf
%          values.
% maximize If true, the minimization problem is transformed into a
%          maximization problem. The default if this parameter is omitted
%          or an empty matrix is passed is false.
%
%OUTPUTS: tuples A 2XnumTuples set of assignment values. This is ordered
%               [Row Index;Column Index]. the number of tuples ranges from
%               max(numRow-1,numCol-1), for the case where as many things
%               as possible are not assigned to the uncontrained values and
%               C(1,1) is not assigned, to numRow+numCol-1 for the case
%               where all rows and columns are assigned to the
%               unconstrained hypothesis and C(1,1) is provided. If the
%               problem is infeasible, this is an empty matrix.
%          gain This is the value of the cost. This is the sum of the
%               values in C corresponding to the tuples. If
%               the problem is infeasible, this is -1.
%           u,v These are the dual variables for the constrained columns
%               and for constrained rows of C. These will satisfy the
%               complementary slackness condition of the problem as
%               transformed. See the example below.
%
%The algorithm utilizes the equivalence of the minimization for the problem
%with two unconstrained indices to one with no unconstrained indices. The
%problem described above, removing the C(1,1) element, is equivalent to the
%optimization problem 
%           minimize (or maximize)
%          \sum_{i=1}^{numAug}\sum_{j=1}^{numAug}D_{i,j}*x_{i,j}
%          subject to
%          \sum_{j=1}^{numAug}x_{i,j}=1 for all i
%          \sum_{i=1}^{numAug}x_{i,j}=1 for all j
%          x_{i,j}=0 or 1.
%where the matrix D is a numAugXnumAug matrix with numAug=(numRow+numCol-2)
%that can be formed as follows:
% D=Inf(numAug,numAug);
% D(1:(numRow-1),1:(numCol-1))=C(2:end,2:end);
% D(numRow:totalAug,numCol:totalAug)=0;
% for curRow=1:(numRow-1)
%     missedCol=(numCol-1)+curRow;
%     D(curRow,missedCol)=C(curRow+1,1);
% end
% for curCol=1:(numCol-1)
%     missedRow=(numRow-1)+curCol;
%     D(missedRow,curCol)=C(1,curCol+1);
% end
%Basically, C(2:end,2:end) is an upper-left submatrix. The upper-right and
%lower-left submatrix consist of Inf terms with diagonals that are the
%missed detection costs for the rows and columns. The lower-right submatrix
%is all zeros. This function solves the problem by implementing the 2D
%Jonker-Volgenant algorithm of [1] for the augmented problem, but the
%matrix C itself is never augmented.
%
%REFERENCES:
%[1] D. F. Crouse, "On Implementing 2D Rectangular Assignment Algorithms,"
%    IEEE Transactions on Aerospace and Electronic Systems, vol. 52, no. 4,
%    pp. 1679-1696, Aug. 2016.
%
%December 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    numRow=size(C,1);
    numCol=size(C,2);
    
    numRowsTrue=numRow;
    numColsTrue=numCol;
    %The number of rows and columns of the virtually augmented matrix.
    numRow=(numRowsTrue-1)+(numColsTrue-1);
    numCol=numRow;
    
    if(isempty(C))
        tuples=[];
        gain=0;
        u=[];
        v=[];
        return
    end
    
%The cost matrix must have all non-negative elements for the assignment
%algorithm to work. This forces all of the elements to be positive. The
%delta is added back in when computing the gain in the end.
    if(maximize==true)
        hasUnconstTuple=(C(1,1)>0);
        
        zeroOffset=max(C(:));
        
        %If C is all negative, do not shift.
        if(zeroOffset<0)
            zeroOffset=0;
        end
        
        C=-C+zeroOffset;
    else
        hasUnconstTuple=(C(1,1)<0);

        zeroOffset=min(C(:));
        
        %If C is all positive, do not shift.
        if(zeroOffset>0)
            zeroOffset=0;
        end
 
        zeroOffset=-zeroOffset;
        C=C+zeroOffset;
    end
   %The virtually augmented matrix will have to be transformed in the
   %same manner as C. zeroOffset is what the zeros in the lower-right 
   %corner have become.
        
    %These store the assignment as it is made.
    col4row=zeros(numRow,1);
    row4col=zeros(numCol,1);
    u=zeros(numCol,1);%The dual variable for the columns
    v=zeros(numRow,1);%The dual variable for the rows.
    
    %Initially, none of the columns are assigned.
    for curUnassCol=1:numCol       
        %This finds the shortest augmenting path starting at k and returns
        %the last node in the path.
        [sink,pred,u,v]=ShortestPath(curUnassCol,u,v,C,col4row,row4col,zeroOffset);
        
        %If the problem is infeasible, mark it as such and return.
        if(sink==0)
            tuples=[];
            gain=-1;
            return;
        end
        
        %We have to remove node k from those that must be assigned.
        j=sink;
        while(1)
            i=pred(j);
            col4row(j)=i;
            h=row4col(i);
            row4col(i)=j;
            j=h;
           
            if(i==curUnassCol)
                break;
            end
        end
    end
    
    %Allocate the maximum possible number of tuples.
    tuples=zeros(2,numRow+1);
    numTuples=0;
    %Deal with the completely unconstrained index.
    if(hasUnconstTuple)
        numTuples=numTuples+1;
        tuples(:,numTuples)=[1;1];
    end
    
    for curRow=1:(numRowsTrue-1)
        curCol=col4row(curRow)+1;
        numTuples=numTuples+1;
        if(curCol>numColsTrue)
            tuples(:,numTuples)=[curRow+1;1];
        else
            tuples(:,numTuples)=[curRow+1;curCol];
        end
    end
    
    for curCol=1:(numColsTrue-1)
        curRow=row4col(curCol)+1;
        if(curRow>numRowsTrue)
            numTuples=numTuples+1;
            tuples(:,numTuples)=[1;curCol+1];
        end
    end

    %Shrink to fit
    tuples=tuples(:,1:numTuples);

    %Only keep the dual variables for the real columns and rows.
    u=u(1:(numColsTrue-1));
    v=v(1:(numRowsTrue-1));
    
    %Calculate the gain that should be returned.
    if(nargout>1)
        gain=0;
        for curTuple=1:numTuples
            gain=gain+C(tuples(1,curTuple),tuples(2,curTuple));
        end
        
        %Adjust the gain for the initial offset of the cost matrix.
        if(maximize)
            gain=-gain+zeroOffset*numTuples;
            u=-u;
            v=-v;
        else
            gain=gain-zeroOffset*numTuples;
        end
    end
end

function [sink, pred, u, v]=ShortestPath(curUnassCol,u,v,C,col4row,row4col,zeroOffset)
%%SHORTESTPATH This implementation fo the shortest path algorithm for the
%           Jonker-Volgenant algorithm has been modified to handle an
%           unconstrained row and an unconstrained column by acting on a
%           virtually augmented matrix. This is called as a subroutine to
%           the assign2DFullStd function.

    %This assumes that unassigned columns go from 1:numUnassigned
    numRowTrue=size(C,1);
    numColTrue=size(C,2);
    
    %The number of rows and columns in the virtually augmented matrix. The
    %columns are augmented with (numRowTrue-1) additional rows for the
    %missed deteciton hypotheses for the rows and the rows are augmented
    %with (numColTrue-1) missed detection hypotheses for the columns.
    numRow=(numRowTrue-1)+(numColTrue-1);
    numCol=numRow;
    
    pred=zeros(numCol,1);
    
    %Initially, none of the rows and columns have been scanned.
    %This will store a 1 in every column that has been scanned.
    ScannedCols=zeros(numCol,1);
    
    %This will store a 1 in every row that has been scanned.
    ScannedRow=zeros(numRow,1);
    Row2Scan=1:numRow;
    numRow2Scan=numRow;
    
    sink=0;
    delta=0;
    curCol=curUnassCol;
    shortestPathCost=ones(numRow,1)*Inf;
    
    while(sink==0)
        %Mark the current column as having been visited.
        ScannedCols(curCol)=1;
        
        %Scan all of the rows that have not already been scanned.
        minVal=Inf;
        for curRowScan=1:numRow2Scan
            curRow=Row2Scan(curRowScan);

            %If scanning past the end of the real rows.
            if(curRow>numRowTrue-1)
                if(curCol>numColTrue-1)
                    %If we have scanned into the part of the virtually
                    %augmented matrix that is all zeros.
                    reducedCost=delta+zeroOffset-u(curCol)-v(curRow);
                elseif(curRow==numRowTrue-1+curCol)
                    %If scanning a missed detection cost.
                    reducedCost=delta+C(1,curCol+1)-u(curCol)-v(curRow);
                else
                    reducedCost=Inf;
                end
            else
                %Here, we have not scanned past the end of the real rows.
                %However, curCol might be past the end of the real columns.
                if(curCol>numColTrue-1)
                    if(curCol==numColTrue-1+curRow)
                        %If scanning a missed detection cost.
                        reducedCost=delta+C(curRow+1,1)-u(curCol)-v(curRow);
                    else
                        reducedCost=Inf;
                    end
                else
                    reducedCost=delta+C(curRow+1,curCol+1)-u(curCol)-v(curRow);
                end
            end

            if(reducedCost<shortestPathCost(curRow))
                pred(curRow)=curCol;
                shortestPathCost(curRow)=reducedCost;
            end
            
            %Find the minimum unassigned row that was scanned.
            if(shortestPathCost(curRow)<minVal)
                minVal=shortestPathCost(curRow);
                closestRowScan=curRowScan;
            end
        end

        if(~isfinite(minVal))
           %If the minimum cost row is not finite, then the problem is
           %not feasible.
           sink=0;
           return;
        end
        
        closestRow=Row2Scan(closestRowScan);
        
        %Add the row to the list of scanned rows and delete it from
        %the list of rows to scan.
        ScannedRow(closestRow)=1;
        numRow2Scan=numRow2Scan-1;
        Row2Scan(closestRowScan)=[];
        
        delta=shortestPathCost(closestRow);
        
        %If we have reached an unassigned column.
        if(col4row(closestRow)==0)
            sink=closestRow;
        else
            curCol=col4row(closestRow);
        end
    end
    
    %Dual Update Step
    
    %Update the first column in the augmenting path.
    u(curUnassCol)=u(curUnassCol)+delta;
    %Update the rest of the columns in the agumenting path.
    sel=(ScannedCols~=0);
    sel(curUnassCol)=0;
    u(sel)=u(sel)+delta-shortestPathCost(row4col(sel));
    
    %Update the scanned rows in the augmenting path.
    sel=ScannedRow~=0;
    v(sel)=v(sel)-delta+shortestPathCost(sel);
end

function [tuples,gain,u,v]=assign2DFullAlt(C,maximize)
%%ASSIGN2DFULLALT This functions performs the type of 2D assignment that often
%       arises as a subproblem when solving an S-dimensional assignment
%       problem for target tracking. The cost function being minimized (or
%       maximized) is
%       \sum_{i=1}^{numRow}\sum_{j=1}^{numCol}C_{i,j}*x_{i,j}
%       subject to
%       \sum_{j=1}^{numCol}x_{i,j}=1 for i=2:numRow
%       \sum_{i=1}^{numRow}x_{i,j}=1 for j=2:numCol
%       x_{i,j}=0 or 1.
%       Note that the first row and the first column can be multiply
%       assigned and every row and every column aside from the first ones
%       must be assigned to exactly one value.
%
%INPUTS: C A (numRow)X(numCol) cost matrix where the first row is the
%          cost of not assigning each column (from number 2) and the first
%          column is the cost of not assigning each row (from number 2) and
%          C(1,1) is an assignment cost that is unconstrained and only
%          assigned if it improves the cost function (given above). The
%          matrix cannot contain any NaNs and if maximizing, cannot
%          contain any Inf values and if minimizing cannot contain any -Inf
%          values. It is not allowed that both C(2:end,1) and C(1,2:end)
%          contain non-finite terms. However, either one individually can
%          contain non-finite terms.
% maximize If true, the minimization problem is transformed into a
%          maximization problem. The default if this parameter is omitted
%          or an empty matrix is passed is false.
%
%OUTPUTS: tuples A 2XnumTuples set of assignment values. This is ordered
%               [Row Index;Column Index]. the number of tuples ranges from
%               max(numRow-1,numCol-1), for the case where as many things as
%               possible are not assigned to the uncontrained values and
%               C(1,1) is not assigned, to numRow+numCol-1 for the case
%               where all rows and columns are assigned to the
%               unconstrained hypothesis and C(1,1) is provided. If the
%               problem is infeasible, this is an empty matrix.
%          gain This is the value of the cost. This is the sum of the
%               values in C corresponding to the tuples. If
%               the problem is infeasible, this is -1.
%           u,v These are the dual variables for the constrained columns
%               and for constrained rows of C. These will satisfy the
%               complementary slackness condition. Note that if the values
%               in C(:,1) are very large, there can be a loss of precision
%               in u and v.
%
%The algorithm is based on the fact that since every unconstrained column
%must be assigned, one can subtract the cost of the unconstrained row from
%all columns and then perform regular 2D assignment using
%assign2DMissedDetect. One can subtract it, because an assignment must be
%made for each constrained column. Subtracting the unconstrained row costs
%changes the gain, but it does not change the optimal assignment of the
%contrained columns. We find the optimal column assignment with a standard
%2D assignment algorithm. The unassigned rows then just need to be added
%back in as assigned to the unconstrained hypothesis and the gain must be
%adjusted. In general, the cost matrix will not be all-positive, so a
%constant value is added to it to make all of the elements positive. The
%hypothesis C(1,1) (from the unmodified matrix) is only added in if it
%improves the cost function. If the unconstrained column contains non-
%finite terms, then the matrix is transposed and the transposed problem is
%solved, after which the results are flipped back.
%
%October 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numRow=size(C,1);
if(numRow>1&&any(~isfinite(C(2:numRow,1))))
    C=C';
    didTranspose=true;
else
    didTranspose=false;
end

numRow=size(C,1);
numConstRow=numRow-1;
numCol=size(C,2);
numConstCol=numCol-1;

if((maximize&&C(1,1)>0)||(maximize==false&&C(1,1)<0))
    hasUnconstTuple=true;
else
    hasUnconstTuple=false;
end

%The special case where a scalar is passed (a single completely
%unconstrained  value), it is only assigned if it does not worsen the cost
%compared to zero (nothing assigned).
if(numConstCol==0&&numConstRow==0)
    if(hasUnconstTuple)
        tuples=[1;1];
        gain=C(1,1);
        u=[];
        v=[];
    else
        tuples=zeros(2,0);
        gain=0;
        u=[];
        v=[];
    end
    return;
end

%Save the original cost matrix so that the gain can be computed to a higher
%precision at the end. We also use this to force the complemantary
%slackness condition to hold to a higher degree than with standard finite
%precision effects. The transformation of the problem by subtracting off
%the unconstrained column below adds inaccuracy to the gain and the dual
%variables.
COrig=C;

if(numCol==1)
    u=[];
    v=zeros(numRow-1,1);
    tuplesRet=zeros(2,0);
else
    %Subtract the costs of the unconstrained row from all of the constrained
    %rows. The cost of assigning a column to the unconstrained row thus
    %becomes 0 (even though we do not update it), and we can then solve the
    %problem using the assign2DMissedDetect algorithm.
    C(2:numRow,2:numCol)=bsxfun(@minus,C(2:numRow,2:numCol),C(2:numRow,1));

    [tuplesRet,~,u,v]=assign2DMissedDetect(C(1:numRow,2:numCol),maximize,true);

    %If the problem appears to be infeasible.
    if(isempty(tuplesRet))
        tuples=[];
        gain=-1;
        return;
    end

    % %If one took the gain from assign2DMissedDetect, one could undo the
    % %subtraction of the unassigned row costs here. However, we choose to
    % recompute the gain from scratch at the end.
    % for k=1:size(tuplesRet,2)
    %     if(tuplesRet(1,k)~=1)
    %         gain=gain+C(tuplesRet(1,k),1);
    %     end
    % end
end

%Adjust the dual variables for the transformations that were performed
%on C.
v=v+C(2:numRow,1);

%Determine the number of rows that were assigned to something other than
%the missed detection hypothesis.
numAssignedRows=sum(tuplesRet(1,:)~=1);

%The number of rows that have to be assigned to the unconstrained
%column.
numUnconstRow=numConstRow-numAssignedRows;

numTuples=numUnconstRow+numConstCol+hasUnconstTuple;

tuples=zeros(2,numTuples);
tuples(:,1:numConstCol)=tuplesRet;
%Adjust the index, because the only rows assigned are not missed
%detections.
tuples(2,1:numConstCol)=tuples(2,1:numConstCol)+1;

%The special case for the completely unconstrained element. It is assigned
%if it does not worsen the cost.
if(hasUnconstTuple)
    tuples(:,(numConstCol+1))=[1;1];
    %If the gain were taken above, we would add in the cost of the new
    %assignment.
    %gain=gain+C(1,1);
    
    curTotalTuple=numConstCol+1;
else
    curTotalTuple=numConstCol;
end

%Sort by the assigned row, so we can tell which rows are not assigned.
[~,idx]=sort(tuples(1,1:curTotalTuple),'ascend');
tuples(:,1:curTotalTuple)=tuples(:,idx);
curTotalTuple=curTotalTuple+1;

%Next, all unassigned rows are assigned to the unconstrained row element
%and are added to the end of curTotalTuple.

%Skip the missed detection rows assigned to columns.
curTuple=1;
while(curTuple<=numTuples&&tuples(1,curTuple)==1)
    curTuple=curTuple+1;
end

%Next, find which rows are missing and add them.
curIdx=2;
while(curIdx<=numRow&&curTotalTuple<=numTuples)
    if((curTuple<=numConstCol+hasUnconstTuple)&&(tuples(1,curTuple)==curIdx))
         curTuple=curTuple+1;
         curIdx=curIdx+1;
    else
         tuples(:,curTotalTuple)=[curIdx;1];
         %If the gain were taken above, we would add in the cost of the new
         %assignment.
         %gain=gain+C(curIdx,1);

         curTotalTuple=curTotalTuple+1;
         curIdx=curIdx+1;
    end
end

%Now, we compute the gain from the original cost matrix. We could also
%adjust the dual variables so that the complementary slackness condition
%strongly holds despite finite precision effects 
gain=0;
for curTuple=1:numTuples
   gain=gain+COrig(tuples(1,curTuple),tuples(2,curTuple));
   
%     if(tuplesTotal(1,curTuple)~=1&&tuplesTotal(2,curTuple)~=1)
%         slack=COrig(tuplesTotal(1,curTuple),tuplesTotal(2,curTuple))-v(tuplesTotal(1,curTuple)-1)-u(tuplesTotal(2,curTuple)-1);
%         v(tuplesTotal(1,curTuple)-1)=v(tuplesTotal(1,curTuple)-1)+slack;
%     end
end

if(didTranspose)
    tuples=flipud(tuples);
    temp=v;
    v=u;
    u=temp;
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
