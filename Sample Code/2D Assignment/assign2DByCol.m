function [col4row, row4col, gain, u, v]=assign2DByCol(C,maximize)
%%ASSIGN2DBYCOL     Solve the two-dimensional assignment problem with a
%                   rectangular cost matrix C, scanning column-wise.
%
%INPUTS:    C           A numRowXnumCol cost matrix that does not contain
%                       any NaNs and where the largest finite element minus
%                       the smallest element is a finite quantity (does not
%                       overflow) when performing minimization and where
%                       the smallest finite element minus the largest
%                       element is finite when performing maximization. 
%                       Forbidden assignments can be given costs of +Inf
%                       for minimization and -Inf for maximization.
%           maximize    If true, the minimization problem is transformed
%                       into a maximization problem. The default if this
%                       parameter is omitted is false.
%
%OUTPUTS:   col4row     A numRowX1 vector where the entry in each element
%                       is an assignment of the element in that row to a
%                       column. 0 entries signify unassigned rows.
%           row4col     A numColX1 vector where the entry in each element
%                       is an assignment of the element in that column to a
%                       row. 0 entries signify unassigned columns.
%           gain        The sum of the values of the assigned elements in
%                       C.
%           u           The dual variable for the rows.
%           v           The dual variable for the columns.
%
%DEPENDENCIES: None
%
%If the number of rows is <= the number of columns, then every row is
%assigned to one column; otherwise every column is assigned to one row. The
%assignment minimizes the sum of the assigned elements (the gain).
%During minimization, assignments can be forbidden by placing Inf in
%elements. During maximization, assignment can be forbidden by placing -Inf
%in elements. The cost matrix can not contain any -Inf elements during
%minimization nor any +Inf elements during mazimization to try to force an
%assignment. If no complete assignment can be made with finite cost,
%then col4row and row4col are empty and gain is set to -1.
%
%The algorithm is described in detail in [1] and [2].
%
%REFERENCES:
%[1] D. F. Crouse, "On implementing 2D rectangular assignment algorithms,"
%    IEEE Transactions on Aerospace and Electronic Systems, accepted 2016.
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

    %These store the assignment as it is made.
    col4row=zeros(numRow,1);
    row4col=zeros(numCol,1);
    u=zeros(numRow,1);%The dual variable for the rows
    v=zeros(numCol,1);%The dual variable for the columns.
    
    %Initially, none of the rows are assigned.
    for curUnassRow=1:numRow        
        %This finds the shortest augmenting path starting at k and returns
        %the last node in the path.
        [sink,pred,u,v]=ShortestPath(curUnassRow,u,v,C,col4row,row4col);
        
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
            row4col(j)=i;
            h=col4row(i);
            col4row(i)=j;
            j=h;
           
            if(i==curUnassRow)
                break;
            end
        end
    end
    
    %Calculate the gain that should be returned.
    if(nargout>2)
        gain=0;
        for curRow=1:numRow
            gain=gain+C(curRow,col4row(curRow));
        end
        %Adjust the gain for the initial offset of the cost matrix.
        if(maximize==true)
            gain=-gain+CDelta*numRow;
        else
            gain=gain+CDelta*numRow;
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

function [sink, pred, u, v]=ShortestPath(curUnassRow,u,v,C,col4row,row4col)
    %This assumes that unassigned columns go from 1:numUnassigned
    numRow=size(C,1);
    numCol=size(C,2);
    pred=zeros(numCol,1);
    
    %Initially, none of the rows and columns have been scanned.
    %This will store a 1 in every row that has been scanned.
    ScannedRows=zeros(numRow,1);
    
    %This will store a 1 in every column that has been scanned.
    ScannedCol=zeros(numCol,1);
    Col2Scan=1:numCol;%Columns left to scan.
    numCol2Scan=numCol;
    
    sink=0;
    delta=0;
    curRow=curUnassRow;
    shortestPathCost=ones(numCol,1)*inf;
    
    while(sink==0)        
        %Mark the current row as having been visited.
        ScannedRows(curRow)=1;
        
        %Scan all of the columns that have not already been scanned.
        minVal=inf;
        for curColScan=1:numCol2Scan
            curCol=Col2Scan(curColScan);
            
            reducedCost=delta+C(curRow,curCol)-u(curRow)-v(curCol);
            if(reducedCost<shortestPathCost(curCol))
                pred(curCol)=curRow;
                shortestPathCost(curCol)=reducedCost;
            end
            
            %Find the minimum unassigned column that was
            %scanned.
            if(shortestPathCost(curCol)<minVal)
                minVal=shortestPathCost(curCol);
                closestColScan=curColScan;
            end
        end
                
        if(~isfinite(minVal))
           %If the minimum cost column is not finite, then the problem is
           %not feasible.
           sink=0;
           return;
        end
        
        closestCol=Col2Scan(closestColScan);
        
        %Add the column to the list of scanned columns and delete it from
        %the list of columns to scan.
        ScannedCol(closestCol)=1;
        numCol2Scan=numCol2Scan-1;
        Col2Scan(closestColScan)=[];
        
        delta=shortestPathCost(closestCol);
        
        %If we have reached an unassigned row.
        if(row4col(closestCol)==0)
            sink=closestCol;
        else
            curRow=row4col(closestCol);
        end
    end
    
    %Dual Update Step
    
    %Update the first row in the augmenting path.
    u(curUnassRow)=u(curUnassRow)+delta;
    %Update the rest of the rows in the agumenting path.
    sel=(ScannedRows~=0);
    sel(curUnassRow)=0;
    u(sel)=u(sel)+delta-shortestPathCost(col4row(sel));
    
    %Update the scanned columns in the augmenting path.
    sel=ScannedCol~=0;
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
