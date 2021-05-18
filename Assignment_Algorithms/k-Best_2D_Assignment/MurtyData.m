classdef MurtyData < handle
%%MURTYDATA A class to hold data used in implementing Murty's Algorithm
%           using an efficient algorithm that inherits the dual variables
%           from the shortest path 2D assignment algorithm at the previous
%           step. This file includes and uses the function
%           ShortestPathUpdate to update a partial assignment. Most people
%           will be more interested in directly using the function
%           kBest2DAssign than using this helper class.
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    properties
        col4rowLCFull
        row4colLCFull
        gainFull
        u
        v
    end
    
    properties (SetAccess = private)
        A%The full assignment matrix.
        %The number of rows (from index 1 to numVarRow) that are not just
        %dummy rows to make the matrix square.
        numVarRow
        %The first row in A where new constraints will be applied when
        %splitting this hypothesis.
        activeRow
        %Columns to which the active row cannot be assigned.
        forbiddenActiveCol
    end 
    
    methods
        function newData=MurtyData(A,numVarRow,activeRow,forbiddenActiveCols,col4rowInit,row4colInit,col2Scan,uInit,vInit)
            %%MURTYDATA Create a new MurtyData object.
            %
            %INPUTS: A               The cost matrix.
            %        numVarRow       The number of rows in A before the
            %                        dummy rows.
            %        unassignedRow   The row whose assignment is not fixed.
            %                        All rows after unassignedRow can be
            %                        reassigned.
            %        forbiddenColumn A column to which the unassigned row
            %                        can not be assigned
            %        row4colInit     The initial row assignment (including
            %                        the forbidden one)
            %        col4rowInit     The initial column assignment
            %                        (including the forbidden one)
            %        col2Scan        A list of the columns that can be
            %                        considered for the new assignment
            %                        problem. It is all columns that were
            %                        not assigned in a subproblem.
            %        uInit and vInit The initial values of the dual
            %                        variables.
            %
            %OUTPUTS: newData A new MurtyData object.
            %
            %October 2013 David F. Crouse, Naval Research Laboratory,
            %Washington D.C.

            numCol=size(A,2);
            
            newData.A=A;%The FULL assignment matrix.
            %The number of rows in the full A that are not just guard rows.
            newData.numVarRow=numVarRow;
            
            if(nargin==2)
                %If this is the very first node being created.
                %Find the best solution.
                [newData.col4rowLCFull,newData.row4colLCFull,newData.gainFull, newData.u, newData.v]=assign2DByCol(A);
                
                %If a feasible solution was found.
                if(newData.gainFull~=-1)
                    newData.activeRow=1;
                    newData.forbiddenActiveCol=false(numCol,1);
                    newData.forbiddenActiveCol(newData.col4rowLCFull(1))=1;
                end
            else
                %If this is a node being split.
                
                [newData.col4rowLCFull,newData.row4colLCFull,newData.gainFull, newData.u, newData.v]=ShortestPathUpdate(A,activeRow,forbiddenActiveCols,col4rowInit,row4colInit,col2Scan,uInit,vInit);
                if(newData.gainFull~=-1)
                    newData.activeRow=activeRow;
                    newData.forbiddenActiveCol=forbiddenActiveCols;
                    newData.forbiddenActiveCol(newData.col4rowLCFull(activeRow))=1;
                end
            end
        end

        function split(data,splitList)
            %%SPLIT This splits the problem into a set of subproblems and
            %       puts the subproblem solutions into the BinaryHeap
            %       splitList.
            
            numCol=size(data.A,2);
            
            %Basically, it is all columns that have not been assigned in
            %one of the subproblems.
            col2Scan=data.col4rowLCFull(data.activeRow:end);
                       
            %The loop ends at numVarRow, because we are not going to allow
            %variable hypotheses to be applied to dummy rows.
            for curRow=data.activeRow:data.numVarRow
                %forbiddenColumns contains a 1 for every column that one is
                %not allowed to assign to curRow.
                if(curRow==data.activeRow)
                    forbiddenColumns=data.forbiddenActiveCol;
                else
                    forbiddenColumns=false(numCol,1);
                    forbiddenColumns(data.col4rowLCFull(curRow))=1;
                end
                
                %Get rid of the current assignment.
                row4colInit=data.row4colLCFull;
                col4rowInit=data.col4rowLCFull;
                row4colInit(col4rowInit(curRow))=0;
                col4rowInit(curRow)=0;                
                
                splitHyp=MurtyData(data.A,...
                                   data.numVarRow,...
                                   curRow,...
                                   forbiddenColumns,...
                                   col4rowInit,...
                                   row4colInit,...
                                   col2Scan,...
                                   data.u,...
                                   data.v);

                %If a valid assignment was produced.
                if(splitHyp.gainFull~=-1)
                    splitList.insert(splitHyp,[]);
                else
                    splitHyp.delete();
                end

                %Remove the current assigned column from the list of
                %columns that can be scanned.
                sel=col2Scan==data.col4rowLCFull(curRow);
                col2Scan(sel)=[];
            end
        end
        
        
      %Comparison functions for two data objects.
        function val=lt(data1,data2)
            %%LT Less than comparison
            if(isa(data1,'MurtyData')==true&&isa(data2,'MurtyData')==true)
                val=data1.gainFull<data2.gainFull;
            elseif(isa(data1,'MurtyData')==true)
                val=data1.gainFull<data2;
            else
                val=data1<data2.gainFull;
            end
        end
      
        function val=gt(data1,data2)
            %%GT Greater than comparison
            if(isa(data1,'MurtyData')==true&&isa(data2,'MurtyData')==true)
                val=data1.gainFull>data2.gainFull;
            elseif(isa(data1,'MurtyData')==true)
                val=data1.gainFull>data2;
            else
                val=data1.gainFull>data2;
            end
        end
      
        function val=le(data1,data2)
            %%LE Less than or equal to comparison
            if(isa(data1,'MurtyData')==true&&isa(data2,'MurtyData')==true)
                val=data1.gainFull<=data2.gainFull;
            elseif(isa(data1,'MurtyData')==true)
                val=data1.gainFull<=data2;
            else
                val=data1<=data2.gainFull;
            end
        end
      
        function val=ge(data1,data2)
            %%GE Greater than or equal to comparison
            if(isa(data1,'MurtyData')==true&&isa(data2,'MurtyData')==true)
                val=data1.gainFull>=data2.gainFull;
            elseif(isa(data1,'MurtyData')==true)
                val=data1.gainFull>=data2;
            else
                val=data1>=data2.gainFull;
            end
        end
      
        function val=ne(data1,data2)
            %%NE Not equal to comparison
            if(isa(data1,'MurtyData')==true&&isa(data2,'MurtyData')==true)
                val=data1.gainFull~=data2.gainFull;
            elseif(isa(data1,'MurtyData')==true)
                val=data1.gainFull~=data2;
            else
                val=data1~=data2.gainFull;
            end
        end
      
        function val=eq(data1,data2)
            %%EQ Equal to comparison
            if(isa(data1,'MurtyData')==true&&isa(data2,'MurtyData')==true)
                val=data1.gainFull==data2.gainFull;
            elseif(isa(data1,'MurtyData')==true)
                val=data1.gainFull==data2;
            else
                val=data1==data2.gainFull;
            end
        end
        
        function disp(data)
            %%DISP Display the MurtyData object.
            disp('Data with col4rowLC: ')
            disp(data.col4rowLCFull)
            disp('and gain: ')
            disp(data.gainFull);
        end
   end % methods
end % classdef


function [col4row, row4col, gain, u, v]=ShortestPathUpdate(C,activeRow,forbiddenActiveCols,col4row,row4col,col2Scan,u,v)
%%SHORTESTPATHUPDATE Update a partial assignment with the shortest path
%                    algorithm, inheriting the dual variables.
%
%This algorithm is based on [1].
%
%REFERENCES:
%M. L. Miller, H. S. Stone, and J. Cox, Ingemar, "Optimizing Murty's ranked
%assignment method," IEEE Transactions on Aerospace and Electronic Systems,
%vol. 33, no. 3, pp. 851-862, Jul. 1997.
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
            
    numRow=size(C,1);
    numCol=size(C,2);
    numCol2Scan=length(col2Scan);
    
    %Initially, none of the rows and columns have been scanned.
    %This will contain a 1 in every row that has been scanned.
    ScannedRows=zeros(numRow,1);
    %This will contain a 1 in every column that has been scanned.
    ScannedCol=zeros(numCol,1);
    
    sink=0;
    pred=zeros(numCol,1);
    delta=0;
    curRow=activeRow;
    shortestPathCost=ones(numCol,1)*inf;
    
    while(sink==0)        
        %Mark the current row as having been visited.
        ScannedRows(curRow)=1;
        
        %Scan all of the columns that have not already been scanned.
        minVal=inf;
        for curColScan=1:numCol2Scan
            curCol=col2Scan(curColScan);
            if(curRow==activeRow&&forbiddenActiveCols(curCol)==1)
                %Do not scan forbidden columns.
                continue;
            end
            
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
           gain=-1;
           return;
        end
        
        closestCol=col2Scan(closestColScan);
        
        %Add the column to the list of scanned columns and delete it from
        %the list of columns to scan.
        ScannedCol(closestCol)=1;
        numCol2Scan=numCol2Scan-1;
        col2Scan(closestColScan)=[];
        
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
    u(activeRow)=u(activeRow)+delta;
    %Update the rest of the rows in the agumenting path.
    sel=(ScannedRows~=0);
    sel(activeRow)=0;
    u(sel)=u(sel)+delta-shortestPathCost(col4row(sel));
    
    %Update the scanned columns in the augmenting path.
    sel=ScannedCol~=0;
    v(sel)=v(sel)-delta+shortestPathCost(sel);

%Augmentation        
    %We have to remove node k from those that must be assigned.
    j=sink;
    while(1)
        i=pred(j);
        row4col(j)=i;
        h=col4row(i);
        col4row(i)=j;
        j=h;

        if(i==activeRow)
            break;
        end
    end
    
    %Calculate the gain that should be returned. Because the duality gap is
    %zero, the gain can be most easily obtained from the dual cost
    %function.    
    gain=0;
    
    for curRow=1:numRow
        gain=gain+C(curRow,col4row(curRow));
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
