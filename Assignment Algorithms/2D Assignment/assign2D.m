function [col4row,row4col,gain,u,v]=assign2D(C,maximize)
%%ASSIGN2D Solve the two-dimensional assignment problem with a rectangular
%          cost matrix C, scanning row-wise. The problem being solved can
%          be formulated as minimize (or maximize)
%          \sum_{i=1}^{numRow}\sum_{j=1}^{numCol}C_{i,j}*x_{i,j}
%          subject to
%          \sum_{j=1}^{numCol}x_{i,j}<=1 for all i
%          \sum_{i=1}^{numRow}x_{i,j}=1 for all j
%          x_{i,j}=0 or 1.
%          Assuming that numCol<=numRow. If numCol>numRow, then the
%          inequality and inequality conditions are switched. A modified
%          Jonker-Volgenant algorithm is used.
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
%          or an empty matrix is passed is false.
%
%OUTPUTS: col4row A numRowX1 vector where the entry in each element is an
%                 assignment of the element in that row to a column. 0
%                 entries signify unassigned rows. If the problem is
%                 infeasible, this is an empty matrix.
%         row4col A numColX1 vector where the entry in each element is an
%                 assignment of the element in that column to a row. 0
%                 entries signify unassigned columns. If the problem is
%                 infeasible, this is an empty matrix.
%            gain The sum of the values of the assigned elements in C. If
%                 the problem is infeasible, this is -1.
%             u,v The dual variable for the columns and for the rows. See
%                 the example below demonstrating how they relate to
%                 complementary slackness of a transformed problem for
%                 minimization or maximization.
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
%EXAMPLE 1:
% C=[Inf,  2,   Inf,Inf,3;
%      7,  Inf, 23, Inf,Inf;
%     17,  24,  Inf,Inf,Inf;
%    Inf,  6,   13, 20, Inf];
% maximize=false;
% [col4row,row4col,gain,u,v]=assign2D(C,maximize)
%One will get an optimal assignment having a gain of 47.
%
%EXAMPLE 2:
%Here we demonstrate that the dual variables satisfy the complementary
%slackness condition for transformed versions of C. If minimization is
%performed, one subtracts minC=min(0,min(C(:))) from C to get the
%transformed problem for complementary slackness; when performing
%maximization, subtract maxC=max(0,max(C(:))). The complementary slackness
%condition could be forced to hold on C by adding minC or maxC to v. The
%complementary slackness condition is described, for example, in [1].
% numRows=50;
% numCols=80;
% C=1000*randn(numRows,numCols);
% maximize=false;
% [col4row,~,gain,u,v]=assign2D(C,maximize);
% minC=min(0,min(C(:)));
% slackVal=0;
% for curRow=1:numRows
%     %Note that we are using C(curRow,col4row(curRow))-minC during
%     %minimization instead of C(curRow,col4row(curRow)).
%     slackVal=slackVal+abs((C(curRow,col4row(curRow))-minC)*(C(curRow,col4row(curRow))-minC-v(curRow)-u(col4row(curRow))));
% end
% slackVal
% %One will see that slackVal is zero (within finite precision limits),
% %which is what the complementary slackness condition says.
% %Switching to maximization with the same matrix:
% maximize=true;
% [col4row,~,gain,u,v]=assign2D(C,maximize);
% slackVal=0;
% maxC=max(0,max(C(:)));
% for curRow=1:numRows
%     %Note that we are using C(curRow,col4row(curRow))+maxC during
%     %minimization instead of C(curRow,col4row(curRow)).
%     slackVal=slackVal+abs((C(curRow,col4row(curRow))-maxC)*(C(curRow,col4row(curRow))-maxC-v(curRow)-u(col4row(curRow))));
% end
% slackVal
% %slackVal is again zero, within finite precision limits.
%
%EXAMPLE 3:
%This is the example used in [3]. Here, we demonstrate how to form
%assignment tuples from col4row.
% C=[7,   51,  52,  87,  38,  60,  74,  66,   0,   20;
%    50,  12,   0,  64,   8,  53,   0,  46,  76,  42;
%    27,  77,   0,  18,  22,  48,  44,  13,   0,  57;
%    62,   0,   3,   8,   5,   6,  14,   0,  26,  39;
%     0,  97,   0,   5,  13,   0,  41,  31,  62,  48;
%    79,  68,   0,   0,  15,  12,  17,  47,  35,  43;
%    76,  99,  48,  27,  34,   0,   0,   0,  28,   0;
%     0,  20,   9,  27,  46,  15,  84,  19,   3,  24;
%    56,  10,  45,  39,   0,  93,  67,  79,  19,  38;
%    27,   0,  39,  53,  46,  24,  69,  46,  23,   1];
% [col4row, row4col, gain]=assign2D(C);
% N=size(C,1);%It is a square matrix.
% tuples=zeros(2,N);
% for curRow=1:N
%     tuples(1,curRow)=curRow;
%     tuples(2,curRow)=col4row(curRow);
% end
% 
% gain1=0;
% gain2=0;
% gain3=0;
% for k=1:N
%     gain1=gain1+C(k,col4row(k));
%     gain2=gain2+C(row4col(k),k);
%     gain3=gain3+C(tuples(1,k),tuples(2,k));
% end
% [gain,gain1,gain2,gain3]
% tuples
%One will see that all of the gains are the same (0) and the assigned 
%tuples match what is in [3]. However, the assigned tuples is NOT obtained
%by attaching col4row to row4col.
%
%REFERENCES:
%[1] D. F. Crouse, "On Implementing 2D Rectangular Assignment Algorithms,"
%    IEEE Transactions on Aerospace and Electronic Systems, vol. 52, no. 4,
%    pp. 1679-1696, Aug. 2016.
%[2] D. F. Crouse, "Advances in displaying uncertain estimates of multiple
%    targets," in Proceedings of SPIE: Signal Processing, Sensor Fusion,
%    and Target Recognition XXII, vol. 8745, Baltimore, MD, Apr. 2013.
%[3] Murty, K. G. "An algorithm for ranking all the assignments in order of
%    increasing cost," Operations Research, vol. 16, no. 3, pp. 682-687,
%    May-Jun. 1968.
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    if(nargin<2||isempty(maximize))
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
        if(maximize)
            gain=-gain+CDelta*numCol;
            u=-u;
            v=-v;
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
            
            reducedCost=delta+C(curRow,curCol)-u(curCol)-v(curRow);
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
    u(sel)=u(sel)+delta-shortestPathCost(row4col(sel));
    
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
