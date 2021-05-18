function [col4row,row4col,gain]=assign2DHungarian(C,maximize)
%%ASSIGN2DHUNGARIAN Solve the two-dimensional assignment problem with a
%          rectangular cost matrix C using the O(n^4) complexity Hungarian
%          algorithm. The problem being solved can be formulated as
%          minimize (or maximize)
%          \sum_{i=1}^{numRow}\sum_{j=1}^{numCol}C_{i,j}*x_{i,j}
%          subject to
%          \sum_{j=1}^{numCol}x_{i,j}<=1 for all i
%          \sum_{i=1}^{numRow}x_{i,j}=1 for all j
%          x_{i,j}=0 or 1.
%          Assuming that numCol<=numRow. If numCol>numRow, then the
%          inequality and inequality conditions are switched. The function
%          assign2D uses an O(n^3) shortest augmenting path modification of
%          the Hungarian algorithm and thus should be preferred over this
%          algorithm.
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
%                 the problem is infeasible, this is -4
%
%This function implements the O(n^4) complexity Hungarian algorithm, which
%is the combination of algorithms 4.2 and 4.3 of Chapter 4.2 of [1]. The
%algorithm has been slightly modified to handle rectangular matrices and to
%transform the input matrix to have all positice elements before running
%the algorithm, so that the algorithm can work on general matrices.
%Additionally, the algorithm contains a floating-point comparison to zero,
%which has been replaced by a threshold to make the algorithm more robust
%when handling non-integer costs.
%
%See the comments to assign2D for examples. The function assign2D
%
%REFERENCES:
%[1] R. Burkard, M. Dell'Amico, and S. Martello, Assignment Problems.
%    Philadelphia: Society for Industrial and Applied Mathematics, 2009.
%
%October 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(maximize))
    maximize=false;
end

numU=size(C,1);
numV=size(C,2);

didTranspose=false;
if(numU>numV)
    C=C';
    
    temp=numU;
    numU=numV;
    numV=temp;
    
    didTranspose=true;
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

UScanned=zeros(numU,1);%Allocate space
VLabeledUnscanned=zeros(numV,1);%Allocate space
pred=zeros(numV,1);%Allocate space

VLabeled=zeros(numV,1);%Allocate space.

unassignedU=1:numU;
numUnassignedU=numU;
%Quick indicator of which vertex is assigned.
assignedBool=false(numU,1);
numAssignedU=0;

%Allocate and initialize.
u=zeros(numU,1);
v=zeros(numV,1);
row4col=zeros(numV,1);
col4row=zeros(numU,1);

while(numAssignedU<numU)
    kIdx=1;
    k=unassignedU(kIdx);
    while(assignedBool(k)==false)
        %%%GET THE ALTERNATING PATH

        %Unlabeled vertices.
        VUnlabeled=1:numV;
        numUnlabeled=numV;
        numLabeledUnscanned=0;
        numScannedU=0;
        numLabeled=0;
        
        fail=false;
        sink=0;
        i=k;
        while(fail==false&&sink==0)
            numScannedU=numScannedU+1;
            UScanned(numScannedU)=i;
            jIdx=1;
            while(jIdx<=numUnlabeled)
                j=VUnlabeled(jIdx);
                cost=C(i,j)-u(i)-v(j);
                %Ad-hoc threshold for testing if zero.
                if(cost<=4*eps(C(i,j)))
                    pred(j)=i;
                    %Remove from the list of unlabeled vertices.
                    VUnlabeled(jIdx)=VUnlabeled(numUnlabeled);
                    numUnlabeled=numUnlabeled-1;

                    %Add to the list of labeled but unscanned vertices.
                    numLabeledUnscanned=numLabeledUnscanned+1;
                    VLabeledUnscanned(numLabeledUnscanned)=j;
                    
                    %Also just add to the list of labeled vertices.
                    numLabeled=numLabeled+1;
                    VLabeled(numLabeled)=j;
                else
                    jIdx=jIdx+1; 
                end
            end
            if(numLabeledUnscanned==0)
                fail=true;
            else
                %Scan the last labeled, unscanned vertex.
                j=VLabeledUnscanned(numLabeledUnscanned);
                numLabeledUnscanned=numLabeledUnscanned-1;

                if(row4col(j)==0)
                    sink=j;
                else
                    i=row4col(j);
                end
            end
        end

        %%UPDATE USING THE ALTERNATING PATH
        if(sink>0)
            %Increase the primal solution.
            assignedBool(k)=true;
            numAssignedU=numAssignedU+1;
            unassignedU(kIdx)=unassignedU(numUnassignedU);
            numUnassignedU=numUnassignedU-1;

            j=sink;
            while(1)
                i=pred(j);
                row4col(j)=i;
                h=col4row(i);
                col4row(i)=j;
                j=h;
                if(i==k)
                    break;
                end
            end
        else
            %Update the dual solution.
            deltaVal=Inf;
            for curU=1:numScannedU
                i=UScanned(curU);
                for curV=1:numUnlabeled
                    j=VUnlabeled(curV);
                    cost=C(i,j)-u(i)-v(j);

                    if(cost<deltaVal)
                        deltaVal=cost;
                    end
                end
            end
            
            %If the problem is infeasible.
            if(~isfinite(deltaVal))
                col4row=[];
                row4col=[];
                gain=-1;
                return;
            end
            
            for curU=1:numScannedU
                i=UScanned(curU);
                u(i)=u(i)+deltaVal;
            end

            for curV=1:numLabeled
                j=VLabeled(curV);
                v(j)=v(j)-deltaVal;
            end
        end
    end
end

if(nargout>2)
    gain=0;
    for curRow=1:numU
        gain=gain+C(curRow,col4row(curRow));
    end

    %Adjust the gain for the initial offset of the cost matrix.
    if(maximize)
        gain=-gain+CDelta*numU;
    else
        gain=gain+CDelta*numU;
    end
end

if(didTranspose)
    temp=col4row;
    col4row=row4col;
    row4col=temp;
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
