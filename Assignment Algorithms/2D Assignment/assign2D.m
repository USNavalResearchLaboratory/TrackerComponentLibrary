function [col4row,row4col,gain,u,v]=assign2D(C,maximize)
%%ASSIGN2D Solve the two-dimensional assignment problem with a rectangular
%          cost matrix C, scanning row-wise. The problem being solved can
%          be formulated as minimize (or maximize)
%          \sum_{i=1}^{numRow}\sum_{j=1}^{numCol}C_{i,j}*x_{i,j}
%          subject to
%          \sum_{j=1}^{numCol}x_{i,j} =1 for all i
%          \sum_{i=1}^{numRow}x_{i,j}<=1 for all j
%          x_{i,j}=0 or 1.
%          Assuming that numCol>=numRow.
%
%INPUTS: C A numRowXnumCol cost matrix that does not contain any NaNs and
%          where the largest finite element minus the smallest element is a
%          finite quantity (does not overflow) when performing minimization
%          and where the smallest finite element minus the largest
%          element is finite when performing maximization. Forbidden
%          assignments can be given costs of +Inf for minimization and -Inf
%          for maximization.
% maximize If true, the minimization problem is transformed into a
%          maximization problem. The default if this parameter is omitted
%          is false.
%
%OUTPUTS: col4row A numRowX1 vector where the entry in each element is an
%                 assignment of the element in that row to a column. 0
%                 entries signify unassigned rows.
%         row4col A numColX1 vector where the entry in each element is an
%                 assignment of the element in that column to a row. 0
%                 entries signify unassigned columns.
%            gain The sum of the values of the assigned elements in C.
%               u The dual variable for the columns. Note that this is on a
%                 transformed version of C.
%               v The dual variable for the rows. Note that this is on a
%                 transformed version of C.
%
%DEPENDENCIES: None
%
%If the number of rows is <= the number of columns, then every row is
%assigned to one column; otherwise every column is assigned to one row. The
%assignment minimizes the sum of the assigned elements (the gain).
%During minimization, assignments can be forbidden by placing Inf in
%elements. During maximization, assignment can be forbidden by placing -Inf
%in elements. The cost matrix can not contain any -Inf elements during
%minimization nor any +Inf elements during maximization to try to force an
%assignment. If no complete assignment can be made with finite cost,
%then col4row and row4col are empty and gain is set to -1.
%
%Note that the dual variables produced by a shortest path assignment
%algorithm that scans by row are not interchangeable with those of a
%shortest path assignment algorithm that scans by column. Matlab stores
%matrices row-wise. Additionally, the dual variables are only valid for the
%transformed cost matrix on which optimization is actually performed, which
%is not necessarily the original cost matrix provided.
%
%When performing minimization with all positive array elements, the initial
%preprocessing step only changes the returned dual variables by offsetting
%the u values by -min(min(C)). Adding min(min(C)) to the returned u values,
%one can get what the dual variables would have been without the
%preprocessing step.
%
%The algorithm is described in detail in [1] and [2].
%
%REFERENCES:
%[1] D. F. Crouse, "On Implementing 2D Rectangular Assignment Algorithms,"
%    IEEE Transactions on Aerospace and Electronic Systems, vol. 52, no. 4,
%    pp. 1679-1696, Aug. 2016.
%[2] D. F. Crouse, "Advances in displaying uncertain estimates of multiple
%    targets," in Proceedings of SPIE: Signal Processing, Sensor Fusion,
%    and Target Recognition XXII, vol. 8745, Baltimore, MD, Apr. 2013.
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    if(nargin<2)
        maximize=false;
    end
    
    numRow=size(C,1);
    numCol=size(C,2);
    
    didFlip=false;
    if(numCol>numRow)
        C=C';
        temp=numRow;
        numRow=numCol;
        numCol=temp;
        didFlip=true;
    end
    
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
            col4row=[];
            row4col=[];
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
    
    %Calculate the gain that should be returned.
    if(nargout>2)
        gain=0;
        for curCol=1:numCol
            gain=gain+C(row4col(curCol),curCol);
        end
        
        %Adjust the gain for the initial offset of the cost matrix.
        if(maximize==true)
            gain=-gain+CDelta*numCol;
        else
            gain=gain+CDelta*numCol;
        end
    end
    
    if(didFlip==true)
        temp=row4col;
        row4col=col4row;
        col4row=temp;
        
        temp=u;
        u=v;
        v=temp;
    end
end

function [sink, pred, u, v]=ShortestPath(curUnassCol,u,v,C,col4row,row4col)
    %This assumes that unassigned columns go from 1:numUnassigned
    numRow=size(C,1);
    numCol=size(C,2);
    pred=zeros(numCol,1);
    
    %Initially, none of the rows and columns have been scanned.
    %This will store a 1 in every column that has been scanned.
    ScannedCols=zeros(numCol,1);
    
    %This will store a 1 in every row that has been scanned.
    ScannedRow=zeros(numRow,1);
    Row2Scan=1:numRow;%Columns left to scan.
    numRow2Scan=numRow;
    
    sink=0;
    delta=0;
    curCol=curUnassCol;
    shortestPathCost=ones(numRow,1)*inf;
    
    while(sink==0)        
        %Mark the current row as having been visited.
        ScannedCols(curCol)=1;
        
        %Scan all of the columns that have not already been scanned.
        minVal=inf;
        for curRowScan=1:numRow2Scan
            curRow=Row2Scan(curRowScan);
            
            reducedCost=delta+C(curRow,curCol)-u(curCol)-v(curRow);
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
        
        %Add the column to the list of scanned columns and delete it from
        %the list of columns to scan.
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
    
    %Update the scanned columns in the augmenting path.
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
