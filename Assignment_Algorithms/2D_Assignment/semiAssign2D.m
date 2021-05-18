function [col4row,gain]=semiAssign2D(C,b,maximize)
%%SEMIASSIGN2D Solve the two-dimensional semi-assignment problem. The
%          problem being solved can be formulated as minimize (or maximize)
%          \sum_{i=1}^{numRow}\sum_{j=1}^{numCol}C_{i,j}*x_{i,j}
%          subject to
%          \sum_{j=1}^{numCol}x_{i,j}=1 for all i
%          \sum_{i=1}^{numRow}x_{i,j}=b(j) for all j
%          x_{i,j}=0 or 1.
%          The problem can only be feasible if numRow>=numCol and
%          sum(b)=numRow.
%
%INPUTS: C A numRowXnumCol cost matrix with numRow>=numCol that does not
%          contain any NaNs and where the largest finite element minus the
%          smallest element is a finite quantity (does not overflow) when
%          performing minimization and where the smallest finite element
%          minus the largest element is finite when performing
%          maximization. Forbidden assignments can be given costs of +Inf
%          for minimization and -Inf for maximization.
%        b A 1XnumCol or numColX1 vector of constraint parameters for the
%          number of assignments to allow in a column. It is required that
%          sum(b)=numRow and that all b>0.
% maximize If true, the minimization problem is transformed into a
%          maximization problem. The default if this parameter is omitted
%          or an empty matrix is passed is false.
%
%OUTPUTS: col4row A numRowX1 vector where the entry in each element is an
%                 assignment of the element in that row to a column. If the
%                 problem is infeasible, this is an empty matrix.
%            gain The sum of the values of the assigned elements in C. If
%                 the problem is infeasible, this is -1.
%
%The 2D semi-assignment problem is just a minor modification of the
%traditional 2D assignment problem. The 2D semi-assignment problem is
%equivalent to performing traditional 2D assignment on a matrix where the
%ith column is repeated b(i) times. That makes the assignment matrix
%square. Rather than just augmenting C and calling assign2D, this function
%does the same thing without augmenting C, just by keeping track of an
%extended set of dual variables for an augmented matrix. See the comments
%to the assign2D function for more on the basics of the basic algorithm,
%which is a modified Jonker-Volgenant algorithm.
%
%Note that specialized techniques for solving the semi-assignment problem
%exist in the literature, such as in [1]. However, we do not use any of
%their heuristic methods to speed up the initialization. Additionally,
%though [1] mentions a shortest augmeting path algorithm in one section
%that is conceptually similar to what is done here, implementations details
%are given in a hard-to-obtain PhD dissertation and their claimed
%theoretical computational complexity of O(numRow^2*numCol) is less than
%the O(numRow^3) that is to be expected of this function.
%
%EXAMPLE:
%Here, we generate a simple semi-assignment problem with a random cost
%matrix. We will see that the same solution is obtain as when using
%assign2D on a fully augmented matrix.
% C=randn(12,4);
% b=[1;3;2;6];
% [col4row,gain]=semiAssign2D(C,b)
% 
% CAug=[C(:,1),C(:,2),C(:,2),C(:,2),C(:,3),C(:,3),C(:,4),C(:,4),C(:,4),C(:,4),C(:,4),C(:,4)];
% mapIdxAug2Idx=[1;2;2;2;3;3;4;4;4;4;4;4];
% [col4rowAug,~,gainAug]=assign2D(CAug);
% col4rowAug=mapIdxAug2Idx(col4rowAug)
% gainAug
%One will see that col4Row and col4RowAug are the same as are the gains.
%
%REFERENCES:
%[1] J. Kennington and Z. Wang, "A shortest augmenting path algorithm for
%    the semi-assignment problem," Operations Research, vol. 40, no. 1,
%    pp. 178-187, Jan.-Feb. 1992.
%
%June 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    if(nargin<3||isempty(maximize))
        maximize=false;
    end
    
    numRow=size(C,1);
    numCol=size(C,2);
    
    if(numRow<numCol)
        error('It is required that numRow>=numCol.')
    end

    if(sum(b)~=numRow)
        error('The sum of the b values should equal numRow.') 
    end
    
    if(length(b)~=numCol)
        error('The length of the b values should be numCol.')
    end
    
    if(any(b==0))
        error('The vector b cannot contain zeros.')
    end
    
%The cost matrix must have all non-negative elements for the assignment
%algorithm to work. This forces all of the elements to be positive. The
%delta is added back in when computing the gain in the end.
    if(maximize==true)
        CDelta=max(C(:));
        
        %If C is all negative, do not shift.
        if(CDelta<0)
            CDelta=0;
        end
        
        C=-C+CDelta;
    else
        CDelta=min(C(:));
        
        %If C is all positive, do not shift.
        if(CDelta>0)
            CDelta=0;
        end
        
        C=C-CDelta;
    end

    %These store the assignment as it is made.
    col4row=zeros(numRow,1);
    %Not numCol. It is for the number of columns of the matrix augmented
    %with entries for the columns that are effectively repeated due to b.
    row4colAug=zeros(numRow,1);
    %The dual variable for the columns --including effectively repeated
    %columns.
    u=zeros(numRow,1);
    v=zeros(numRow,1);%The dual variable for the rows.
    
    %We do not augment the C matrix, but we keep track of which columns are
    %repeated.
    fullCol2Col=zeros(numRow,1);
    curIdx=1;
    for curCol=1:numCol
        span=curIdx:(curIdx+b(curCol)-1);
        fullCol2Col(span)=curCol;
        curIdx=curIdx+b(curCol);
    end

    %Initially, none of the columns are assigned. We scan across the real
    %and the virtual repeated columns.
    for curUnassCol=1:numRow       
        %This finds the shortest augmenting path starting at k and returns
        %the last node in the path.
        [sink,pred,u,v]=ShortestPathAug(curUnassCol,u,v,C,col4row,row4colAug,fullCol2Col);

        %If the problem is infeasible, mark it as such and return.
        if(sink==0)
            col4row=[];
            gain=-1;
            return;
        end
        
        %We have to remove node k from those that must be assigned.
        j=sink;
        while(1)
            i=pred(j);
            col4row(j)=i;
            h=row4colAug(i);
            row4colAug(i)=j;
            j=h;
           
            if(i==curUnassCol)
                break;
            end
        end
    end
    
    %Map the full columns back to the original columns.
    col4row=fullCol2Col(col4row);

    %Calculate the gain that should be returned.
    if(nargout>1)
        gain=0;
        
        for curRow=1:numRow
            gain=gain+C(curRow,col4row(curRow));
        end

        %Adjust the gain for the initial offset of the augmented cost
        %matrix (why we use numRow and not numCol).
        if(maximize)
            gain=-gain+CDelta*numRow;
        else
            gain=gain+CDelta*numRow;
        end
    end
end

function [sink, pred, u, v]=ShortestPathAug(curUnassCol,u,v,C,col4row,row4colAug,fullCol2Col)
%%SHORTESTPATHAUG This shortest augmenting path algorithm assumes that
%           numCol<numRow in C and that the shortest agumenting path is
%           desired over a modified C such that certain columns have been
%           duplicated so that numCol=numRow. Rather than explicitly
%           creating the augmented C, the function fullCol2Col maps the
%           indices of the full columns to columns in C. The values in pred
%           are for the augmented C matrix.
%
%June 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.

    numRow=size(C,1);
    pred=zeros(numRow,1);
    
    %Initially, none of the rows and columns have been scanned.
    %This will store a 1 in every column that has been scanned.
    ScannedCols=zeros(numRow,1);
    
    %This will store a 1 in every row that has been scanned.
    ScannedRow=zeros(numRow,1);
    Row2Scan=1:numRow;
    numRow2Scan=numRow;
    
    sink=0;
    delta=0;
    curCol=curUnassCol;
    shortestPathCost=Inf(numRow,1);
    
    while(sink==0)        
        %Mark the current column as having been visited.
        ScannedCols(curCol)=1;
        
        %Scan all of the rows that have not already been scanned.
        minVal=Inf;
        for curRowScan=1:numRow2Scan
            curRow=Row2Scan(curRowScan);
            
            curCCol=fullCol2Col(curCol);
            
            reducedCost=delta+C(curRow,curCCol)-u(curCol)-v(curRow);
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
    %Update the rest of the columns in the augmenting path.
    sel=(ScannedCols~=0);
    sel(curUnassCol)=0;
    u(sel)=u(sel)+delta-shortestPathCost(row4colAug(sel));
    
    %Update the scanned rows in the augmenting path.
    sel=ScannedRow~=0;
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
