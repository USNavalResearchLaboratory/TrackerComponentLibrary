function [col4row,row4col]=maxCardBipartMatching(A,algorithm)
%%MAXCARDBIPARTMATCHING Solve the maximum cardinality bipartite matching
%           problem. That is, find the maximum cardinality rows of A to
%           columns of A given that a connection between the two exists
%           (meaning A(row,col)~=0).
%
%INPUTS: A A numRowsXnumCols matrix indicating connections between the rows
%          and the columns. A(i,j)=0 indicates that there is no connection
%          from row i to column j; A(i,j)=1 means there is a connection.
% algorithm An optional parameter specifying the algorithm to use. Possible
%          values are:
%          0 (The default if omitted or an empty matrix is passed) Use the
%            Hopcroft-Karp algorithm, which is described briefly in Chapter
%            26.6 of [1] and in more detail in the original paper of [2].
%          1 Use the cardinality matching algorithm given as Algorithm 3.1
%            in Chapter 3.2 of [3].
%
%OUTPUTS: col4row The column to which a particular row is assigned in the 
%                 maximum cardinality bipartite matching. col4Row(curRow)=0
%                 means that the row is not assigned to any column. This is
%                 a column vector.
%         row4col The row to which a particular column is assigned in
%                 the maximum cardinality bipartite matching.
%                 row4Col(curCol)=0 means that the column is not assigned
%                 to any row. This is a column vector.
%
%This is an implementation of the Hopcroft-Karp algorithm. The algorithm
%is described briefly in Chapter 26.6 of [1]. It is originally from [2],
%where a more detailed description is used.
%
%EXAMPLE:
% A=[1, 0, 1, 0, 0, 0, 0;
%    1, 1, 0, 0, 0, 0, 0;
%    1, 0, 1, 0, 0, 0, 0;
%    0, 0, 0, 1, 1, 1, 1;
%    0, 1, 1, 0, 0, 0, 0;
%    0, 0, 1, 1, 0, 0, 0];
% [col4row,row4col]=maxCardBipartMatching(A)
% %One finds col4row=[1;2;3;5;0;4]. 
% A=[1, 1, 0, 0, 0;
%    1, 0, 0, 0, 1;
%    0, 0, 1, 1, 0;
%    1, 0, 0, 0, 1;
%    0, 1, 0, 1, 0];
% [col4row,row4col]=maxCardBipartMatching(A)
% %One finds col4row=[2;5;3;1;4]. 
% %A completely infeasible problem.
% A=zeros(5,5);
% [col4row,row4col]=maxCardBipartMatching(A)
% col4row=[0;0;0;0;0]; 
%
%REFERENCES:
%[1] T. H. Cormen, C. E. Leiserson, R. L. Rivest, and C. Stein, 
%    Introduction to Algorithms, 2nd ed. Cambridge, MA: The MIT Press,
%    2001.
%[2] J. E. Hopcroft and R. M. Karp, "An n^(5/2) algorithm for maximum
%    matchings in bipartite graphs," SIAM Journal on Computing, vol. 2, no.
%    4, pp. 225-231, Dec. 1973.
%[3] R. Burkard, M. Dell'Amico, and S. Martello, Assignment Problems.
%    Philadelphia: Society for Industrial and Applied Mathematics, 2009.
%
%November 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(algorithm))
    algorithm=0;
end

switch(algorithm)
    case 0
        [col4row,row4col]=HopcroftKarp(A);
    case 1
        [col4row,row4col]=CardinalityMatching(A);
    otherwise
        error('Unknown algorithm specified.')
end
end

function [col4row,row4col]=HopcroftKarp(A)
%%HOPCROFTKARP Solve the maximum cardinality bipartite matching problem
%              using the Hopcroft-Karp algorithm. 
%
%The algorithm is described briefly in Chapter 26.6 of [1]. It is
%originally from [2], where a more detailed description is used.
%
%REFERENCES:
%[1] T. H. Cormen, C. E. Leiserson, R. L. Rivest, and C. Stein, 
%    Introduction to Algorithms, 2nd ed. Cambridge, MA: The MIT Press,
%    2001.
%[2] J. E. Hopcroft and R. M. Karp, "An n^(5/2) algorithm for maximum
%    matchings in bipartite graphs," SIAM Journal on Computing, vol. 2, no.
%    4, pp. 225-231, Dec. 1973.
%
%November 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.

numRows=size(A,1);
numCols=size(A,2);

row4col=zeros(numCols,1);
col4row=zeros(numRows,1);
%This holds distances to nodes for the breadth first search. distVec(1) is
%a sentinal to detect fesibility. The sentinal essentially reflects a dummy
%node to which all unmatched columns are connected.
distVec=zeros(numRows+1,1);

%These are to hold a queue of nodes for the breadth-first search.
theQueue=zeros(numRows+numCols,1);
numInQueue=0;

while(1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Begin Breadth-First Search%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%First, perform a breadth-first search to find all augmenting paths of
%assigned and unassigned edges. distVec(1) is a dummy vertex that is
%attached to all rows and all columns. Thus, if we start from the dummy
%vertex and come back to it using an alternating path of unique vertices,
%then we have found an augmenting path.

%Set all unassigned rows to a distance of 0 and all assigned rows
%to a distance of Inf. Unassigned rows are pushed onto the queue for
%investigation. All augmenting paths begin from an unassigned row.
    for curRow=1:numRows
        if(col4row(curRow)==0)
            %If the row is not already assigned.
            distVec(curRow+1)=0;
            
            %Add the row to the queue.
            numInQueue=numInQueue+1;
            theQueue(numInQueue)=curRow;
        else
            distVec(curRow+1)=Inf;
        end
    end
%The sentinal node used to indicate completion.
    distVec(1)=Inf;

%Now, we go through all of the unassigned rows on the queue.
    while(numInQueue>0)
        %Pop from the top of the queue
        curRow=theQueue(numInQueue);
        numInQueue=numInQueue-1;
        
        %If this node is not the sential node and can provide a shorter
        %path to the sentinal node than has already been found (A path to
        %the sentinal node means an alternating path was found).
        if(distVec(curRow+1)<distVec(1))
            for curCol=1:numCols
                %If there is a connection from the row to the column.
                if(A(curRow,curCol)~=0)
                    %If (row4col(curCol),curCol) is an edge that has not
                    %yet been explored.
                    if(distVec(row4col(curCol)+1)==Inf)
                        distVec(row4col(curCol)+1)=distVec(curRow+1)+1;
                        
                        %Push onto the queue
                        numInQueue=numInQueue+1;
                        theQueue(numInQueue)=row4col(curCol);
                    end
                end
            end
        end
    end
    %If the breadth-first search could not find any augmenting paths, then
    %the algorithm is complete.
    if(distVec(1)==Inf)
        break;
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%End Breadth-First Search%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %This loop tries to find/ augment on a maximal set of vertex disjoint
    %shortest augmenting paths. We try to find an augmenting path starting
    %from each row. The partial depth-first search runs until it completes
    %the first augmenting path. If it finds a path, it removes all visited
    %vertices along the way. This guarantees that paths found later will be
    %disjoint.
    for curRow=1:numRows
        %If the row has not yet ben assigned
        if(col4row(curRow)==0)
            %The depth-first search tries to find (and augment on) an
            %augmenting path beginning with curRow. It then removes the
            %nodes involved in that path.
            [~,col4row,row4col,distVec]=partialDepthFirstSearch(curRow,A,col4row,row4col,distVec);
        end
    end
end
end

function [retVal,col4row,row4col,distVec]=partialDepthFirstSearch(curRow,A,col4row,row4col,distVec)

numCols=length(row4col);

%The check for curRow has to do with the recursive use of this function.
if(curRow~=0)
    for curCol=1:numCols
       %If there is a connection from the row to the column.
        if(A(curRow,curCol)~=0)
            %Follow the distances set by the breadth-first search.
            %We can only go down things that are adjacent (differ
            %by 1).
            if(distVec(row4col(curCol)+1)==distVec(curRow+1)+1)
                [recurRetVal,col4row,row4col,distVec]=partialDepthFirstSearch(row4col(curCol),A,col4row,row4col,distVec);
                
                if(recurRetVal==true)
                    row4col(curCol)=curRow;
                    col4row(curRow)=curCol;
                    retVal=true;
                    return;
                end
            end
        end
    end
    
    %If there is no augmenting path beginning with curRow.
    distVec(curRow+1)=Inf;
    retVal=false;
    return;
else%Reached a leaf node (the sential node). we have an augmenting path.
    retVal=true;
    return
end
end

function [col4row,row4col]=CardinalityMatching(A)
%%CARDINALITYMATCHING Solve the maximum cardinality bipartite matching
%              problem using the cardinality matching algorithm. 
%
%This implements algorithm 3.1 from Chapter3.2 of [1].
%
%REFERENCES:
%[1] R. Burkard, M. Dell'Amico, and S. Martello, Assignment Problems.
%    Philadelphia: Society for Industrial and Applied Mathematics, 2009.
%
%October 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.

numRows=size(A,1);
numCols=size(A,2);

row4col=zeros(numCols,1);
col4row=zeros(numRows,1);

l=zeros(numCols,1);%Allocate for the labels and mark unlabeled.
r=zeros(numRows,1);
%Unmatched row vertices.
L=1:numRows;
numL=numRows;
R=zeros(numCols,1);%Allocate space.
numR=0;

while(numL>0||numR>0)
    if(numL>0)
        x=L(numL);
        %Scan the left vertex
        numL=numL-1;
        
        for j=1:numCols
            if(A(x,j)&&l(j)==0)%If the edge exists and it is unlabeled.
                l(j)=x;

                numR=numR+1;
                R(numR)=j;
            end
        end
    else%numR>0
        x=R(numR);
        %Scan the right vertex.
        numR=numR-1;
        
        if(row4col(x)~=0)%If this edge is in the matching.
            i=row4col(x);
            r(i)=x;%Label it.
            numL=numL+1;
            L(numL)=i;
        else
            %Backtrack the labels to get the path and keep track of
            %unmatched vertices.
            while(1)
                i=l(x);
                col4row(i)=x;
                row4col(x)=i;

                x=r(i);

                if(x==0)
                    break;
                end
            end
            
            LUnmatched=(col4row==0);
            numL=sum(LUnmatched);
            idx=1:numRows;
            L(1:numL)=idx(LUnmatched);

            numR=0;
            %Cancel all labels.
            l(:)=0;
            r(:)=0;
        end
    end
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
