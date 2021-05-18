function [col4row,row4col,gain,p,pR]=assign2DAuctionSquare(C,maximize,epsilon,algorithm)
%%ASSIGN2DAUCTIONSQUARE Solve the two-dimensional assignment problem with a
%         square cost matrix C using an auction algorithm that does not use
%         epsilon scaling. The problem being solved can be formulated as
%         minimize (or maximize)
%         \sum_{i=1}^{numRow}\sum_{j=1}^{numCol}C_{i,j}*x_{i,j}
%          subject to
%         \sum_{j=1}^{numCol}x_{i,j} =1 for all i
%         \sum_{i=1}^{numRow}x_{i,j}<=1 for all j
%         x_{i,j}=0 or 1.
%         Assuming that numCol=numRow. The Jonker-Volgenant algorithm
%         implemented in assign2D is generally better than the auction
%         algorithm. Moreover, auction algorithms with epsilon scaling as
%         described in Chapter 7.1.4 of [1] can potentially improve
%         performance.
%        
%INPUTS: C An NXN cost matrix that does not contain any NaNs. Forbidden
%          assignments can be given costs of +Inf for minimization and -Inf
%          for maximization.
% maximize If true, the minimization problem is transformed into a
%          maximization problem. The default if this parameter is omitted
%          is false.
% epsilon The epsilon value in the auction algorithm if this parameter is
%         omitted or an empty matrix is passed, then the default value if
%         the difference between the two nonequal elements in C having the
%         smallest magnitudes divided by 1.01 or eps of the maximum
%         magntiude element in C, whichever is largest.
% algorithm An optional parameter selecting which type of auction algorithm
%         is used. Possible values are:
%         0 (The default if omitted or an empty matrix is passed) Use the
%           forward auction algorithm of Chapter 7.1.1 of [1] with the
%           Jacobi version of the bid phase.
%         1 Use the reverse auction algorithm of Chapter 7.2.1 of [1] with
%           the Jacobi version of the bid phase.
%         2 Use the forward-reverse auction algorithm of Chapter 7.2.1 of
%           [1] with the Jacobi version of the bid phase using N iterations
%           for each direction.
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
%               p The price variables used in the auction algorithm. When
%                 using algorithm 2, these are the forward price variables.
%              pR When using algorithm 2, these are the price values for
%                 the reverse auction algorithm. Otherwise, this is an
%                 empty matrix.
%
%The auction algorithm is derived for use with integer costs, but can be
%used with real costs if epsilon is set sufficiently small (the default
%value). However, small values of epsilon make the algorithm slower. The
%algorithm has weak polynomial complexity. For a guaranteed optimum and a
%strong polynomial complexity execution time, use the Jonker-Volgenant
%algorithm in assign2D.
%
%Infeasibility in each implementation is detected using the criterion in
%Equation 7.13 in Chapter 7.1.5.
%
%EXAMPLE:
% C=[Inf,  2,   Inf,Inf;
%      7,  Inf, 23, Inf;
%     17,  24,  Inf,Inf;
%    Inf,  6,   13, 20];
% maximize=false;
% [col4row,row4col,gain,p]=assign2DAuctionSquare(C,maximize)
%One will get an optimal assignment having a gain of 62.
%
%REFERENCES:
%[1] D. P. Bertsekas, Network Optimization: Continuous and Discrete Models.
%    Belmont, MA: Athena Scientific, 1998.
%
%October 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    if(size(C,1)~=size(C,2))
        error('The cost matrix C must be square.')
    end

    if(nargin<2||isempty(maximize))
        maximize=false;
    end

    if(nargin<3||isempty(epsilon))
    %We have to determine the maximum epsilon that guarantees convergence.
    %This is going to be less than the (nonzero) difference between the two
    %smallest (not equal) elements of the cost matrix divided by n.
        n=length(C);
        min1=min(abs(C(:)));
        maxEps=eps(max(abs(C(:))));
        min2=min(abs(C(abs(C(:))~=min1)));
        if(~isempty(min2))%If the entires weren't all the same value.
            %A sufficiently small value.
            epsilon=max(maxEps,(min2-min1)/(1.01*n));
        else
            epsilon=1e-6;
        end
    end
    
    if(nargin<4||isempty(algorithm))
        algorithm=0;
    end

    if(maximize==false)
        C=-C; 
    end
    
    switch(algorithm)
        case 0
            [col4row,row4col,gain,p]=forwardAuction(C,epsilon);
            pR=[];
        case 1
            [col4row,row4col,gain,p]=reverseAuction(C,epsilon);
            pR=[];
        case 2
            [col4row,row4col,gain,p,pR]=FRAuction(C,epsilon);
        otherwise
            error('Unknown Algorithm specified.')
    end

    if(maximize==false)
        gain=-gain;
    end
    
    if(isempty(col4row))
        gain=-1; 
    end
end

function [x,y,gain,p]=forwardAuction(A,epsilon)
%%FORWARDAUCTION This is the basic forward auction algorithm of Chapter
%                7.1.1 of [1].
%
%REFERENCES:
%[1] D. P. Bertsekas, Network Optimization: Continuous and Discrete Models.
%    Belmont, MA: Athena Scientific, 1998.
%
%October 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.
    
    n=length(A);
    %x(row) is the column index assigned to the given row (0=unassigned).
    x=zeros(n,1);
    %y(column) is the row index assigned to the given column
    %(0=unassigned).
    y=zeros(n,1);
    %The price for each column. Initially, all prices are zero.
    p=zeros(n,1);
    
    %These are used in two checks of feasibility.
    maxA=max(A(isfinite(A)));
    minA=min(A(isfinite(A)));
    
    %If the problem is not feasible, because everything has infinite cost.
    if(isempty(maxA))
        x=[];
        y=[];
        gain=-1;
        return;
    end

    while(1)
%Bidding Phase
        %The maximum bid for each column. This is used during the
        %assignment phase.
        bMax=-Inf*ones(n,1);
        rowMax=zeros(n,1);%The maximum bidder for that column.
        
        numUnAssigned=0;
        for curRow=1:n
            %If the current row is unassigned.
            if(x(curRow)==0)
                %Find the most profitable value & column.
                [~,bestCol]=max(A(curRow,:)-p');
                %Now, find the second most profitable value.
                sel=[1:(bestCol-1) (bestCol+1):n];
                WMat=A(curRow,sel)-p(sel)';
                w=max(WMat);
                wMin=min(WMat(isfinite(WMat)));

                %Check for infeasibility as in Equation 7.13 in Chapter
                %7.1.5.
                if(isempty(wMin)||wMin-minA<-(2*n-1)*(maxA-minA)-(n-1)*epsilon)
                    x=[];
                    y=[];
                    gain=-1;
                    return;
                end
                
                b=A(curRow,bestCol)-w+epsilon;%The bid.
                if(b>bMax(bestCol))
                    bMax(bestCol)=b;
                    rowMax(bestCol)=curRow;
                end
                numUnAssigned=numUnAssigned+1;
            end
        end

%Check to see if the algorithm is finished.
        if(numUnAssigned==0)
            gain=0;
            for curRow=1:n
                gain=gain+A(curRow,x(curRow));
            end

            return
        end

%Assignment phase
        for curCol=1:n
            if(rowMax(curCol)~=0)
                %If a bid was placed on this column, raise the price to the
                %highest bid.
                p(curCol)=bMax(curCol);
                %Then, break all existing assignments and add the maximum
                %bidder.
                if(y(curCol)~=0)%If this had been previously assigned.
                    x(y(curCol))=0;
                end

                y(curCol)=rowMax(curCol);
                x(rowMax(curCol))=curCol;
            end
        end
    end
end

function [x,y,gain,pf]=reverseAuction(A,epsilon)
%%REVERSEAUCTION This is the basic reverse auction algorithm of Chapter
%                7.2.1 of [1].
%
%REFERENCES:
%[1] D. P. Bertsekas, Network Optimization: Continuous and Discrete Models.
%    Belmont, MA: Athena Scientific, 1998.
%
%October 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    %This is the basic reverse auction algorithm. It is assumed that A is
    %an nXn matrix of positive, real numbers. The returned values are the
    %row assignments, the column assignments, and the gain.
    n=length(A);
    %x(row) is the column index assigned to the given row (0=unassigned).
    x=zeros(n,1);
    %y(column) is the row index assigned to the given column
    %(0=unassigned).
    y=zeros(n,1);
    %The profit for each row. Initially, all profits are zero.
    pf=zeros(n,1);
    
    %This is used in two checks of feasibility.
    maxA=max(A(isfinite(A)));
    minA=min(A(isfinite(A)));
    
    %If the problem is not feasible, because everything has infinite cost.
    if(isempty(maxA))
        x=[];
        y=[];
        gain=-1;
        return;
    end
    
    while(1)
%Bidding Phase
        %The maximum bid for each row. This is used during the assignment
        %phase.
        bMax=-Inf*ones(n,1);
        colMax=zeros(n,1);%The maximum bidder for that row.
        
        numUnAssigned=0;
        for curCol=1:n
            %If the current column is unassigned.
            if(y(curCol)==0)
                %Find the best row.
                [~,bestRow]=max(A(:,curCol)-pf);
                %Now, find the second best value.
                sel=[1:(bestRow-1) (bestRow+1):n];
                WMat=A(sel,curCol)-pf(sel);
                w=max(WMat);
                wMin=min(WMat(isfinite(WMat)));
                
                %Check for infeasibility as in Equation 7.13 in Chapter
                %7.1.5.
                if(isempty(wMin)||wMin-minA<-(2*n-1)*(maxA-minA)-(n-1)*epsilon)
                    x=[];
                    y=[];
                    gain=-1;
                    return;
                end

                b=A(bestRow,curCol)-w+epsilon;
                if(b>bMax(bestRow))
                    bMax(bestRow)=b;
                    colMax(bestRow)=curCol;
                end
                numUnAssigned=numUnAssigned+1;
            end
        end
        
%Check to see if the algorithm is finished.
        if(numUnAssigned==0)
            gain=0;
            for curRow=1:n
                gain=gain+A(curRow,x(curRow));
            end

            return
        end
        
        %Assignment phase
        for curRow=1:n
            if(colMax(curRow)~=0)
                %If a bid was placed on this column, raise the price to the
                %higest bid.
                pf(curRow)=bMax(curRow);
                %Then, break all existing assignments and add the maximum
                %bidder.
                if(x(curRow)~=0)%If this had been previously assigned.
                    y(x(curRow))=0;
                end

                x(curRow)=colMax(curRow);
                y(colMax(curRow))=curRow;
            end
        end
        
    end
end

function [x,y,gain,p,pf]=FRAuction(A,epsilon)
%%REVERSEAUCTION This is the basic forward-reverse auction algorithm of
%                Chapter 7.2.1 of [1] run with n iterations of each
%                direction for an nXn A matrix.
%
%REFERENCES:
%[1] D. P. Bertsekas, Network Optimization: Continuous and Discrete Models.
%    Belmont, MA: Athena Scientific, 1998.
%
%October 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    %This is the basic forward-reverse auction algorithm. It is assumed
    %that A is an nXn matrix of positive, real numbers. The returned values
    %are the row assignments, the column assignments, and the gain.
    n=length(A);
    %x(row) is the column index assigned to the given row (0=unassigned).
    x=zeros(n,1);
    %y(column) is the row index assigned to the given column
    %(0=unassigned).
    y=zeros(n,1);
    %The price for each column. Initially, all prices are zero.
    p=zeros(n,1);
    %The profit for each row. These will be initialized during the forward
    %part of the algorithm.
    pf=zeros(n,1);

    %This is used in two checks of feasibility.
    maxA=max(A(isfinite(A)));
    minA=min(A(isfinite(A)));
    
    %If the problem is not feasible, because everything has infinite cost.
    if(isempty(maxA))
        x=[];
        y=[];
        gain=-1;
        return;
    end
    
    while(1)
%%%%STEP 1: Run forward auction
%Bidding Phase
        %The maximum bid for each column. This is used during the
        %assignment phase.
        bMax=-Inf*ones(n,1);
        rowMax=zeros(n,1);%The maximum bidder for that column.
        numUnAssigned=0;
        for curIter=1:n
            for curRow=1:n
                %If the current row is unassigned.
                if(x(curRow)==0)
                    %Find the most profitable value & column.
                    [~,bestCol]=max(A(curRow,:)-p');
                    %Now, find the second most profitable value.
                    sel=[1:(bestCol-1) (bestCol+1):n];
                    WMat=A(curRow,sel)-p(sel)';
                    w=max(WMat);
                    wMin=min(WMat(isfinite(WMat)));

                    %Check for infeasibility as in Equation 7.13 in Chapter
                    %7.1.5.
                    if(isempty(wMin)||wMin-minA<-(2*n-1)*(maxA-minA)-(n-1)*epsilon)
                        x=[];
                        y=[];
                        gain=-1;
                        return;
                    end

                    b=A(curRow,bestCol)-w+epsilon;%The bid.
                    if(b>bMax(bestCol))
                        bMax(bestCol)=b;
                        rowMax(bestCol)=curRow;
                    end
                    numUnAssigned=numUnAssigned+1;
                end
            end

    %Check to see if the algorithm is finished.
            if(numUnAssigned==0)
                gain=0;
                for curRow=1:n
                    gain=gain+A(curRow,x(curRow));
                end

                return
            end

    %Assignment phase
            for curCol=1:n
                if(rowMax(curCol)~=0)
                    %If a bid was placed on this column, raise the price to
                    %the highest bid.
                    p(curCol)=bMax(curCol);
                    %Set the profits to be used in the backwards algorithm.
                    pf(rowMax(curCol))=A(rowMax(curCol),curCol)-p(curCol);

                    %Then, break all existing assignments and add the
                    %maximum bidder.
                    if(y(curCol)~=0)%If this had been previously assigned.
                        x(y(curCol))=0;
                    end

                    y(curCol)=rowMax(curCol);
                    x(rowMax(curCol))=curCol;
                end
            end
        end
%%%%STEP 2: Run reverse Auction
        %The maximum bid for each row. This is used during the assignment
        %phase.
        bMax=-Inf*ones(n,1);
        colMax=zeros(n,1);%The maximum bidder for that row.
        numUnAssigned=0;

        for cutIter=1:n
            for curCol=1:n
                %If the current column is unassigned.
                if(y(curCol)==0)
                    %Find the best row.
                    [~,bestRow]=max(A(:,curCol)-pf);
                    %Now, find the second best value.
                    sel=[1:(bestRow-1) (bestRow+1):n];
                    w=max(A(sel,curCol)-pf(sel));

                    b=A(bestRow,curCol)-w+epsilon;
                    if(b>bMax(bestRow))
                        bMax(bestRow)=b;
                        colMax(bestRow)=curCol;
                    end
                    numUnAssigned=numUnAssigned+1;
                end
            end

    %Check to see if the algorithm is finished.
            if(numUnAssigned==0)
                gain=0;
                for curRow=1:n
                    gain=gain+A(curRow,x(curRow));
                end

                return
            end

            %Assignment phase
            for curRow=1:n
                if(colMax(curRow)~=0)
                    %If a bid was placed on this column, raise the price to
                    %the higest bid.
                    pf(curRow)=bMax(curRow);
                    %Set the profits to be used in the backwards algorithm.
                    p(colMax(curRow))=A(curRow,colMax(curRow))-pf(curRow);
                    %Then, break all existing assignments and add the
                    %maximum bidder.
                    if(x(curRow)~=0)%If this had been previously assigned.
                        y(x(curRow))=0;
                    end

                    x(curRow)=colMax(curRow);
                    y(colMax(curRow))=curRow;
                end
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
