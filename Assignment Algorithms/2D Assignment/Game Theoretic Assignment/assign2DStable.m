function [col4Row,row4Col]=assign2DStable(rowPrefs,colPrefs,prefType)
%%ASSIGN2DSTABLE Perform 2D assignment to find a stable matching between
%                two sets of things where the items in set1 rank their
%                preference for the items in set2 and the items in set2
%                rank their preference for the items in set1. An
%                assignment of things in set1 to things in set2 is unstable
%                if there is any pair a,b in set1 that are assigned to
%                (A,B) in set2 even though  b prefers A to B and A prefers
%                b to a. The stable assignment is not necessarily unique.
%                This function finds the stable assignment that is optimal
%                over the set given by the rows of the preference matrices
%                (unless there are more rows than columns).
%
%INPUTS: rowPrefs,colPrefs The matrices of preferences. The format of
%               these depends on  the prefType input. Meanings depending on
%               prefType are:
%               prefType=0: B oth are numSet1XnumSet2 matrices where
%                     rowPrefs(curSet1,curSet2) is the preference of the
%                     item curSet1 in set1 for being assigned to the item
%                     curSet2 in set2, and colPrefs(curSet1,curSet2) is the
%                     preference of the item curSet2 in set 2 for being
%                     assigned to the item curSet1 in set1. The values in
%                     the matrices do not need to be intergers. Values
%                     should be finite.
%               prefType=1: rowPrefs is a numSet1XnumSet2 matrix and
%                     colPrefs is a numSet2XnumSet1 matrix.
%                     rowPrefs(curSet1,i) holds the index of the item in
%                     set2 that is ith in the order preferred by item
%                     curSet1 in set 1. Similarly, rowPrefs(curSet2,i)
%                     holds the index of the item in set1 that is ith in
%                     the order preferred by item curSet2 in set2.
%        prefType A parameter indicating the format of rowPrefs and
%                 colPrefs. If omitted or an empty matrix is passed, a
%                 default value of 0 is used.
%
%OUTPUTS:col4Row A numSet1X1 vector such that col4Row(curSet1) is the index
%                of the thing in set2 to which item curSet1 in set1 is
%                assigned. 0 means it is unassigned, which only occurs if
%                there are more things in set1 than set2.
%        row4Col A numSet2X1 vector such that row4Col(curSet2) is the index
%                of the thing in set1 to which item curSet2 in set2 is
%                assigned. 0 means it is unassigned, which only occurs if
%                there are more things in set2 than set1.
%
%This problem is typically referred to the "stable marriage problem" as it
%is often posed where set1 consists of bachelors and set2 consists of
%single women who want to get married. The algorithm then produces a stable
%assignment of bachelors to women that is optimal for the bachelors, unless
%there are more bachelors than women.
%
%The algorithm implemented is that of [1] and for sets 1 and 2 both
%consisting of n items is O(n^2) complexity. The case of differing numbers
%of elements in set1 and set 2 are also mentioned.
%
%EXAMPLES:
% This example is from [1]
% rowPrefs=[1,2,3,4;
%           1,4,3,2;
%           2,1,3,4;
%           4,2,3,1];
% colPrefs=[3,3,2,3;
%           4,1,3,2;
%           2,4,4,1;
%           1,2,1,4];
% [col4Row,row4Col]=assign2DStable(rowPrefs,colPrefs)
%Here, the optimal assignment is col4Row=[1;4;1;2];
%
%This shows how the algorithm can be used with a rectangular matrix
% rowPrefs=[1,2,3,4,5;
%           1,4,3,2,5;
%           2,1,3,4,5;
%           4,2,3,1,5];
% colPrefs=[3,3,2,3,1;
%           4,1,3,2,2;
%           2,4,4,1,3;
%           1,2,1,4,4];
% [col4Row,row4Col]=assign2DStable(rowPrefs,colPrefs)
%In this instance, col4Row is the same, but now the last item in row4col is
%zero, because one column remains unassigned.
%
%The following example is taken from Chapter 1 of [2], where the other
%format for the input is used:
% rowPrefs=[4,1,2,3;
%           2,3,1,4;
%           2,4,3,1;
%           3,1,2,4];
% colPrefs=[4,1,3,2;
%           1,3,2,4;
%           1,2,3,4;
%           4,1,3,2];
% [col4Row,row4Col]=assign2DStable(rowPrefs,colPrefs,1)
%Here, the stable assignment obtained in col4Row=[4;3;2;1];
%
%An equivalent problem is when the preference matrices are modified as per
%prefType. Then, we have
% rowPrefs=[2,3,4,1;
%           3,1,2,4;
%           4,1,3,2;
%           2,3,1,4];
% colPrefs=[2,1,1,2;
%           4,3,2,4;
%           3,2,3,3;
%           1,4,4,1];
% [col4Row,row4Col]=assign2DStable(rowPrefs,colPrefs,0)
%Which has the same solution.
%
%The final example is also from [2], and is
% rowPrefs=[5,7,1,2,6,8,4,3;
%           2,3,7,5,4,1,8,6;
%           8,5,1,4,6,2,3,7;
%           3,2,7,4,1,6,8,5;
%           7,2,5,1,3,6,8,4;
%           1,6,7,5,8,4,2,3;
%           2,5,7,6,3,4,8,1;
%           3,8,4,5,7,2,6,1];
% colPrefs=[5,3,7,6,1,2,8,4;
%           8,6,3,5,7,2,1,4;
%           1,5,6,2,4,8,7,3;
%           8,7,3,2,4,1,5,6;
%           6,4,7,3,8,1,2,5;
%           2,8,5,3,4,6,7,1;
%           7,5,2,1,8,6,4,3;
%           7,4,1,5,2,3,6,8];
% [col4Row,row4Col]=assign2DStable(rowPrefs,colPrefs,1)
%Where one gets the assignment col4Row=[5;3;8;6;7;1;2;4].
%
%REFERENCES:
%[1] D. Gale and L. S. Shapley, "College admissions and the stability of
%    marriage," The American Mathematical Monthly, vol. 69, no. 1, pp.
%    9-15, Jan. 1962.
%[2] D. Gusfield and R. W. Irving, The Stable Marriage Problem: Structure
%    and Algorithms. The MIT Press, 1989.
%
%February 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numRows=size(rowPrefs,1);

if(nargin<3||isempty(prefType))
    prefType=0;
end

didFlip=false;

%If the input was given in terms of lists where the rows rank the columns
%and the columns rank the rows.
if(prefType==1)
    numCols=size(colPrefs,1);
    %If there more rows than columns, then find a column-optimal stable
    %assignment as a row-optimal one does not exist.
    if(numCols<numRows)
       temp=numCols;
       numCols=numRows;
       numRows=temp;
       
       temp=colPrefs;
       colPrefs=rowPrefs;
       rowPrefs=temp;
       
       didFlip=true;
    end

    colPrefsNew=zeros(numRows,numCols);

    for curCol=1:numCols
        [~,idx]=sort(colPrefs(curCol,:),'ascend');
        colPrefsNew(:,curCol)=idx;
    end
    colPrefs=colPrefsNew;
    colRanking4Row=rowPrefs;
else
    numCols=size(rowPrefs,2);
    %If there more rows than columns, then find a column-optimal stable
    %assignment as a row-optimal one does not exist.
    if(numCols<numRows)
       temp=numCols;
       numCols=numRows;
       numRows=temp;
        
       temp=colPrefs;
       colPrefs=rowPrefs';
       rowPrefs=temp';
       didFlip=true;
    end
    
    colRanking4Row=zeros(numRows,numCols);
    for curRow=1:numRows
        %Go through all of the columns according to the preference of the
        %row until an assignment can be made.
        [~,colRanking]=sort(rowPrefs(curRow,:),'ascend');
        colRanking4Row(curRow,:)=colRanking;
    end
end

row4Col=zeros(numCols,1);
col4Row=zeros(numRows,1);

freeRowList=1:numRows;

%No rows have been assigned.
num2Assign=numRows;
%While rows are left unassigned.
while(num2Assign>0)
    %Find a free row. It can be any one.
    freeRow=freeRowList(num2Assign);
    %freeRow now holds the index of the current free row.
    
    %Go through all of the columns according to the preference of the row
    %until an assignment can be made.
    colRanking=colRanking4Row(freeRow,:);
    for curRank=1:numCols
        curCol=colRanking(curRank);
        
        %If the column of preference is free, then make the assignment.
        if(row4Col(curCol)==0)
            row4Col(curCol)=freeRow;
            col4Row(freeRow)=curCol;
            num2Assign=num2Assign-1;
            break;
        else
            %Otherwise, see if this assignment is better than the existing
            %assignment.
            curAssignedRow=row4Col(curCol);
            
            %If the columns prefers this row to the currently assigned one.
            if(colPrefs(freeRow,curCol)<colPrefs(curAssignedRow,curCol))
                row4Col(curCol)=freeRow;
                col4Row(freeRow)=curCol;
                
                %Unassign the previously assigned row.
                col4Row(curAssignedRow)=0;
                %Put the now unassigned row in the list of free rows.
                freeRowList(num2Assign)=curAssignedRow;
                break;
            end
        end
    end
end

if(didFlip)
    temp=row4Col;
    row4Col=col4Row;
    col4Row=temp;
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
