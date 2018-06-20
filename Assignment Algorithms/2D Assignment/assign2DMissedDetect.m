function [tuples,gain,u,v]=assign2DMissedDetect(C,maximize,missedDetectByRow)
%%ASSIGN2DMISSEDDETECT Solve the two-dimensional assignment problem common
%          to target tracking, where, in a tracking context, one row (or
%          column) represents missed detection costs. The problem being
%          solved can be formulated (if missedDetectByRow=true) as
%          minimize (or maximize)
%          \sum_{i=1}^{numMeas+1}\sum_{j=1}^{numTar}C_{i,j}*x_{i,j}
%          subject to
%          \sum_{i=1}^{numMeas+1}x_{i,j}=1 for all j
%          \sum_{j=1}^{numTar}x_{i,j}<=1 for i=2:numMeas
%          x_{i,j}=0 or 1. 
%          where the index i=1 has the costs of assigning a target to a
%          missed detection, so multiple targets can be assigned to a
%          missed detection hypothesis. Measurements are indexed starting
%          at 2.
%
%INPUTS: C A (numMeas+1)XnumTar (or if targetByRow=false, a
%          numTarX(numMeas+1) cost matrix) where the first row (or column
%          if targetByRow=false) is the cost of a missed detection for each
%          target and the subsequent rows (columns) are the costs of
%          assigning a measurement to a target. The matrix cannot contain
%          any NaNs and if maximizing, cannot containing any Inf values and
%          if minimizing cannot contain any -Inf values.
% maximize If true, the minimization problem is transformed into a
%          maximization problem. The default if this parameter is omitted
%          or an empty matrix is passed is false.
% missedDetectByRow If this is true, the first row of C is the set of
%          missed detection costs. Otherwise, the first column of C is the
%          set of missed detection costs. The default if this parameter is
%          omitted or an empty matrix is passed is true.
%
%OUTPUTS: tuples A 2XnumTar set of assignment values. If
%               targetsByRow=false (or is omitted), then this is ordered
%               [Measurement Index;Target Index]; Otherwise, the ordering
%               is [Target Index; Measurement Index]. If the problem is
%               infeasible, this is an empty matrix.
%          gain This is the value of the cost. This is the sum of the
%               values in C corresponding to the tuples. If the problem is
%               infeasible, this is -1.
%           u,v These are the dual variables for the constrained columns
%               and for constrained rows of C. These will satisfy the
%               complementary slackness condition of a transformed
%               problem.
%
%This function implements a modified Jonker-Volgenant algorithm as
%described in [1] and [2]. However, unlike the assign2D function, which
%implements the same algorithm for rectangular matrices, this
%implementation modifies the ShortestPath step to allow the multiple
%assignment of targets to the missed detection hypothesis. it is
%implicitely solving an assignment problem with measurement costs augmented
%by a numTarXnumTar matrix of costs of missed detection for each target on
%the diagonal and infinite costs on the off diagonals. This changes how
%indexation is performed in the ShortestPath subroutine in this function
%compared to the assign2D function.
%
%EXAMPLE:
%Here we compare how one can represent a set of measurement assignment
%costs and missed detection values compared to using the assign2D function,
%which requires that one augment the assignment matrix. We also verify that
%the complementary slackness condition of the dual for optimality. The
%complementary slackness condition could be forced to hold on C by
%adding minC=min(0,min(C(:))) (during minimization) or
%maxC=max(0,max(C(:))) (during maximization) to v. As it is, the
%complementary slackness condition holds on C-minC for minimization and
%C-maxC for maximization (mentioned in [1]) is satisfied.
% numTar=10;
% numMeas=20;
% maximize=false;
% C=2*rand(numMeas,numTar)-1;
% missedDetect=2*rand(1,numTar)-1;
% C1=[missedDetect;C];
% minC1=min(0,min(C1(:)));
% missedDetectByRow=true;
% [tuples,gain,u,v]=assign2DMissedDetect(C1,maximize,missedDetectByRow);
% 
% %We verify the correctness of the complementary slackness condition.
% slackVal=0;
% for k=1:numTar
%     %If assigned to a constrained value
%     if(tuples(1,k)~=1&&tuples(2,k)~=1)
%         slackVal=slackVal+abs((C1(tuples(1,k),tuples(2,k))-minC1)*(C1(tuples(1,k),tuples(2,k))-minC1-v(tuples(1,k)-1)-u(tuples(2,k))));
%     end
% end
% slackVal
% %One will see that the slack value is zero, indicating that the
% %complementary slackness condition is satisfied. Note that this is on the
% %problem transformed, subtracting minC1. If maximization were performed,
% %one would have to subtract the maximum entry in C1.
% 
% CD=Inf(numTar,numTar);
% CD(diagElIdx(numTar))=missedDetect;
% C2=[C;CD];
% [~,row4col,gainAlt]=assign2D(C2,maximize);
% tuplesAlt=zeros(2,numTar);
% for curCol=1:numTar
%     rowNum=row4col(curCol);
%     if(rowNum>numMeas)
%         rowNum=1;
%     else
%         rowNum=rowNum+1;
%     end
%     tuplesAlt(:,curCol)=[rowNum;curCol];
% end
% all(tuples(:)==tuplesAlt(:))
% all(gainAlt==gain)
%One will see that both of the above are 1, demonstrating that both
%approaches yield the same solution. 
%
%REFERENCES:
%[1] D. F. Crouse, "On Implementing 2D Rectangular Assignment Algorithms,"
%    IEEE Transactions on Aerospace and Electronic Systems, vol. 52, no. 4,
%    pp. 1679-1696, Aug. 2016.
%[2] D. F. Crouse, "Advances in displaying uncertain estimates of multiple
%    targets," in Proceedings of SPIE: Signal Processing, Sensor Fusion,
%    and Target Recognition XXII, vol. 8745, Baltimore, MD, Apr. 2013.
%
%October 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    if(nargin<2||isempty(maximize))
        maximize=false;
    end
    
    if(nargin<3||isempty(missedDetectByRow))
        missedDetectByRow=true; 
    end
    
    numRow=size(C,1);
    numCol=size(C,2);
    
    if(missedDetectByRow==false)
        C=C';
        temp=numRow;
        numRow=numCol;
        numCol=temp;
    end
    numRowsTrue=numRow;
    %The number of rows of the virtually augmented matrix.
    numRow=numRow+numCol-1;

%The cost matrix must have all non-negative elements for the assignment
%algorithm to work. This forces all of the elements to be positive. The
%delta is added back in when computing the gain in the end.
    if(maximize==true)
        CDelta=max(max(C));
        
        %If C is all negative, do not shift.
        if(CDelta<0)
            CDelta=0;
        end
        
        C=-C+CDelta;
    else
        CDelta=min(min(C));
        
        %If C is all positive, do not shift.
        if(CDelta>0)
            CDelta=0;
        end
        
        C=C-CDelta;
    end

    %These store the assignment as it is made.
    col4row=zeros(numRow,1);
    row4col=zeros(numCol,1);
    u=zeros(numCol,1);%The dual variable for the columns
    v=zeros(numRow,1);%The dual variable for the rows.
    
    %Initially, none of the columns are assigned.
    for curUnassCol=1:numCol       
        %This finds the shortest augmenting path starting at k and returns
        %the last node in the path.
        [sink,pred,u,v]=ShortestPath(curUnassCol,u,v,C,col4row,row4col);
        
        %If the problem is infeasible, mark it as such and return.
        if(sink==0)
            tuples=[];
            gain=-1;
            return;
        end
        
        %We have to remove node sink from those that must be assigned.
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

    %Calculate the gain that should be returned.
    if(nargout>1)
        gain=0;
        for curCol=1:numCol
            theRow=row4col(curCol);
            
            if(theRow>numRowsTrue-1)
                gain=gain+C(1,curCol);
            else
                gain=gain+C(theRow+1,curCol);
            end
        end
        
        %Adjust the gain for the initial offset of the cost matrix.
        if(maximize)
            gain=-gain+CDelta*numCol;
            u=-u;
            v=-v;
        else
            gain=gain+CDelta*numCol;
        end
    end
    
    %Mark the columns that have been assigned to missed detections
    sel=row4col>numRowsTrue-1;
    row4col(sel)=1;%Missed detections
    row4col(~sel)=row4col(~sel)+1;%Targets

    %Get rid of the dual variables related to the "unconstrained" column.
    v=v(1:(numRowsTrue-1));

    tuples=zeros(2,numCol);
    if(missedDetectByRow==false)
        for curCol=1:numCol
           tuples(:,curCol)=[curCol;row4col(curCol)]; 
        end
        temp=u;
        u=v;
        v=temp;
    else
        for curCol=1:numCol
            tuples(:,curCol)=[row4col(curCol);curCol]; 
        end
    end
end

function [sink, pred, u, v]=ShortestPath(curUnassCol,u,v,C,col4row,row4col)
%%SHORTESTPATH This implementation of the shortest path algorithm for the
%           Jonker-Volgenant algorithm has been modified to handle an
%           unconstrained row by acting on a virtually augmented matrix.

    %This assumes that unassigned columns go from 1:numUnassigned
    numRowsTrue=size(C,1);
    numCol=size(C,2);
    pred=zeros(numCol,1);
    
    %The first row is the row containing missed detection costs for all of
    %the columns and thus is not counted. We then add numCol virtual rows
    %to the end that are the missed detection costs.
    numRow=numRowsTrue-1+numCol;
    
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
            
            if(curRow>numRowsTrue-1)
                %If scanning a missed detection cost.
                if(curRow==numRowsTrue-1+curCol)
                    reducedCost=delta+C(1,curCol)-u(curCol)-v(curRow);
                    if(reducedCost<shortestPathCost(curRow))
                        pred(curRow)=curCol;
                        shortestPathCost(curRow)=reducedCost;
                    end
                else
                    reducedCost=Inf;
                end
            else
                reducedCost=delta+C(curRow+1,curCol)-u(curCol)-v(curRow);
            end

            if(reducedCost<shortestPathCost(curRow))
                pred(curRow)=curCol;
                shortestPathCost(curRow)=reducedCost;
            end

            %Find the minimum unassigned column that was
            %scanned.
            if(shortestPathCost(curRow)<minVal)
                minVal=shortestPathCost(curRow);
                closestRowScan=curRowScan;
            end
        end
                
        if(~isfinite(minVal))
           %If the minimum cost column is not finite, then the problem is
           %not feasible.
           sink=0;
           return;
        end
        
        closestRow=Row2Scan(closestRowScan);
        
        %Add the row to the list of scanned rows and delete it from the
        %list of rows to scan.
        ScannedRow(closestRow)=1;
        numRow2Scan=numRow2Scan-1;
        Row2Scan(closestRowScan)=[];
        
        delta=shortestPathCost(closestRow);
        
        %If we have reached an unassigned row.
        if(col4row(closestRow)==0)
            sink=closestRow;
        else
            curCol=col4row(closestRow);
        end
    end
    
    %Dual Update Step
    
    %Update the first row in the augmenting path.
    u(curUnassCol)=u(curUnassCol)+delta;
    %Update the rest of the rows in the agumenting path.
    sel=(ScannedCols~=0);
    sel(curUnassCol)=0;
    u(sel)=u(sel)+delta-shortestPathCost(row4col(sel));
    
    %Update the scanned rows in the augmenting path.
    sel=(ScannedRow~=0);
    v(sel)=v(sel)-delta+shortestPathCost(sel);
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
