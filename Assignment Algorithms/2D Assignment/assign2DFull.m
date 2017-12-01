function [tuplesTotal,gain,u,v]=assign2DFull(C,maximize)
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
%          assigned if it does not worsen the cost function (given above).
%          The matrix cannot contain any NaNs and if maximizing, cannot
%          containing any Inf values and if minimizing cannot contain any
%          -Inf values. It is not allowed that both C(2:end,1) and
%          C(1,2:end) contain non-finite terms. However, either one
%          individually can contain non-finite terms.
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
%               unconstrained hypothesis and C(1,1) is provided.
%          gain This is the value of the cost. This is the sum of the
%               values in C corresponding to the tuples.
%           u,v These are the dual variables for the constrained columns
%               and for constrained rows of C. These will satisfy the
%               complementary slackness condition. Note that if the values
%               in C(:,1) are very large, there can be a loss of precision
%               in u and v.
%
%This algorithm uses a modified form of the Jonker-Volgenant algorithm.
%First, we note that since every unconstrained columns must be assigned, we
%can subtract the cost of the unconstrained row from all columns and then
%perform regular 2D assignment using assign2DMissedDetect. We can subtract
%it, because an assignment must be made for each unconstrained column.
%Subtracting the unconstrained row costs changes the gain, but it does not
%change the optimal assignment of the contrained columns. We find the
%optimal column assignment with a standard 2D assignment algorithm. The
%unassigned rows then just need to be added back in as assigned to the
%unconstrained hypothesis and the gain must be adjusted. In general, the
%cost matrix will not be all-positive, so a constant value is added to it
%to make all of the elements positive. The hypothesis C(1,1) (from the
%unmodified matrix) is only added in if it does not worsen the cost
%function. If the unconstrained column contains non-finite terms, then the
%matrix is transoposed and the transposed problem is solved, then the
%results are flipped back.
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
% %normal 2D assignment problem formulated as a linear programming solution)
% %is the extra contraint on the C(1,1) term. We have to explicitely say that
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
%In this example, we verify that the complementary slackness condition
%holds (i.e. that the dual variables are the values one would expect them
%to be).
% numRow=8;
% numCol=6;
% C=2*rand(numRow,numCol)-1;
%
% %Solving the problem using
% [tuples,gain,u,v]=assign2DFull(C,false);
%
% slackVal=0;
% for k=1:(numCol-1)
%     %If assigned to a constrained value.
%     if(tuples(1,k)~=1&&tuples(2,k)~=1)
%         slackVal=slackVal+abs(C(tuples(1,k),tuples(2,k))*(C(tuples(1,k),tuples(2,k))-v(tuples(1,k)-1)-u(tuples(2,k)-1)));
%     end
% end
% slackVal
%slackVal will be zero within finite precision limits, demonstrating that
%the complementary slackness conditions holds for the dual variables that
%have been returned.
%
%October 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(maximize))
    maximize=false;
end

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

if((maximize&&C(1,1)>=0)||(maximize==false&&C(1,1)<=0))
    hasZeroTuple=true;
else
    hasZeroTuple=false;
end

%The special case where a scalar is passed (a single completely
%unconstrained  value), it is only assigned if it does not worsen the cost
%compared to zero (nothing assigned).
if(numConstCol==0&&numConstRow==0)
    if(hasZeroTuple)
        tuplesTotal=[1;1];
        gain=C(1,1);
        u=[];
        v=[];
    else
        tuplesTotal=zeros(2,0);
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

%Subtract the costs of the unconstrained row from all of the constrained
%rows. The cost of assigning a column to the unconstrained row thus
%becomes 0 (even though we do not update it), and we can then solve the
%problem using the assign2DMissedDetect algorithm.
C(2:numRow,2:numCol)=bsxfun(@minus,C(2:numRow,2:numCol),C(2:numRow,1));

[tuples,~,u,v]=assign2DMissedDetect(C(1:numRow,2:numCol),maximize,true);

% %If one took the gain from assign2DMissedDetect, one could undo the
% %subtraction of the unassigned row costs here. However, we choose to
% recompute the gain from scratch at the end.
% for k=1:size(tuples,2)
%     if(tuples(1,k)~=1)
%         gain=gain+C(tuples(1,k),1);
%     end
% end

%Adjust the dual variables for the transformations that were performed on
%C.
v=v+C(2:numRow,1);

%Determine the number of rows that were assigned to something other than
%the missed detection hypothesis.
numAssignedRows=sum(tuples(1,:)~=1);

%The number of rows that have to be assigned to the unconstrained
%column.
numUnconstRow=numConstRow-numAssignedRows;

numTuples=numUnconstRow+numConstCol+hasZeroTuple;

tuplesTotal=zeros(2,numTuples);
tuplesTotal(:,1:numConstCol)=tuples;
%Adjust the index, because the only rows assigned are not missed
%detections.
tuplesTotal(2,1:numConstCol)=tuplesTotal(2,1:numConstCol)+1;

%The special case for the completely unconstrained element. It is assigned
%if it does not worsen the cost.
if(hasZeroTuple)
    tuplesTotal(:,(numConstCol+1))=[1;1];
    %If the gain were taken above, we would add in the cost of the new
    %assignment.
    %gain=gain+C(1,1);
    
    curTotalTuple=numConstCol+1;
else
    curTotalTuple=numConstCol;
end

%Sort by the assigned row, so we can tell which rows are not assigned.
[~,idx]=sort(tuplesTotal(1,1:curTotalTuple),'ascend');
tuplesTotal(:,1:curTotalTuple)=tuplesTotal(:,idx);
curTotalTuple=curTotalTuple+1;

%Next, all unassigned rows are assigned to the unconstrained row element
%and are added to the end of curTotalTuple.

%Skip the missed detection rows assigned to columns.
curTuple=1;
while(curTuple<=numTuples&&tuplesTotal(1,curTuple)==1)
    curTuple=curTuple+1;
end

%Next, find which rows are missing.
curIdx=2;
while(curIdx<=numRow&&curTotalTuple<=numTuples)
    if((curTuple<=numConstCol+hasZeroTuple)&&(tuplesTotal(1,curTuple)==curIdx))
         curTuple=curTuple+1;
         curIdx=curIdx+1;
    else
         tuplesTotal(:,curTotalTuple)=[curIdx;1];
         %If the gain were taken above, we would add in the cost of the new
         %assignment.
         %gain=gain+C(curIdx,1);
         
         curTotalTuple=curTotalTuple+1;
         curIdx=curIdx+1;
    end
end

%Now, we compute the gain from the original cost matrix. We also adjust the
%dual variables so that the complementary slackness condition strongly
%holds despite finite precision effects 
gain=0;
for curTuple=1:numTuples
   gain=gain+COrig(tuplesTotal(1,curTuple),tuplesTotal(2,curTuple));
   
    if(tuplesTotal(1,curTuple)~=1&&tuplesTotal(2,curTuple)~=1)
        slack=COrig(tuplesTotal(1,curTuple),tuplesTotal(2,curTuple))-v(tuplesTotal(1,curTuple)-1)-u(tuplesTotal(2,curTuple)-1);
        v(tuplesTotal(1,curTuple)-1)=v(tuplesTotal(1,curTuple)-1)+slack;
    end
end

if(didTranspose)
    tuplesTotal=flipud(tuplesTotal);
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
