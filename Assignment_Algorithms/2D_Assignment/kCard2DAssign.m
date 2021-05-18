function [tuples,gain]=kCard2DAssign(C,k,maximize)
%%KCARD2DASSIGN Solve the two-dimensional k-cardinality assignment problem
%          with a rectangular cost matrix C. The problem being solved can
%          be formulated as minimize (or maximize)
%          \sum_{i=1}^{numRow}\sum_{j=1}^{numCol}C_{i,j}*x_{i,j}
%          subject to
%          \sum_{j=1}^{numCol}x_{i,j}<=1 for all i
%          \sum_{i=1}^{numRow}x_{i,j}<=1 for all j
%          \sum_{i=1}^{numRow}\sum_{j=1}^{numCol}x_{i,j}=k
%          x_{i,j}=0 or 1.
%
%INPUTS: C A numRowXnumCol cost matrix that does not contain any NaNs and
%          where the largest finite element minus the smallest element is a
%          finite quantity (does not overflow) when performing minimization
%          and where the smallest finite element minus the largest
%          element is finite when performing maximization. Forbidden
%          assignments can be given costs of +Inf for minimization and -Inf
%          for maximization.
%        k The integer number of assignments to make.
%          k<=min(numRow,numCol).
% maximize If true, the minimization problem is transformed into a
%          maximization problem. The default if this parameter is omitted
%          or an empty matrix is passed is false.
%
%OUTPUTS: tuples A 2XnumRow set of assignment values. This is ordered
%                [rowIndex;columnIndex]. If no feasible solution exists,
%                then an empty matrix will be returned. 
%          gain This is the value of the cost. This is the sum of the
%               values in C corresponding to the tuples.
%
%As noted in [1], the k-cardinality 2D assignment problem can be
%transformed into a standard rectangular 2D assignment problem. This
%function implements that transformation and calls assign2D. 
%
%EXAMPLE:
% C=[7,   51,  52,  87;
%    50,  12,   0,  64;
%    27,  77,   0,  18;
%    62,   0,   3,   8];
% k=4;
% [tuples4,gain4]=kCard2DAssign(C,4)
% [tuples3,gain3]=kCard2DAssign(C,3)
% [tuples2,gain2]=kCard2DAssign(C,2)
% [tuples1,gain1]=kCard2DAssign(C,1)
%One will see optimal gains of 25, 7, 0, and 0.
%
%REFERENCES:
%[1] A. Volgenant, "Solving the k-cardinality assignment problem by
%    transformation," European Journal of Operational Research, vol. 157,
%    no. 2, pp. 322-331, 1 Sep. 2004.
%
%June 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    if(nargin<3||isempty(maximize))
        maximize=false;
    end
    
    numRow=size(C,1);
    numCol=size(C,2);
    
    if(k==0)
       tuples=[];
       gain=0;
       return;
    end
    
    if(k<0||k~=fix(k))
        error('k is invalid.')
    end
    
    if(k>numRow||k>numCol)
        %The problem is not feasible.
        tuples=[];
        gain=-1;
        return
    end
    
    COrig=C;
    
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
    
    %The minimum value in C is now 0 or something positive. We shall shift
    %the values to guarantee that 0 is less than every value in C.
    %Find the maximum finite value in C.
    CMax=max(max(C(isfinite(C))));

    %We shall now offset everything by a small fraction of the maximum
    %value. This is so that we guarantee that values will be assigned
    %to the zero padded rows that shall be added to C.
    COffset=CMax*1e-10;
    C=C+COffset;

    %Perform 2D assignment on the augmented matrix. An algorithm that does
    %not explicitely construct that assignment matrix is also possible.
    [~,row4Col]=assign2D([C;zeros(numCol-k,numCol)],false);
    
    if(isempty(row4Col))
        %If the problem is infeasible.
        tuples=[];
        gain=-1;
        return
    end
    
    tuples=zeros(2,k);
    curTuple=1;
    gain=0;
    for curCol=1:numCol
        
        curRow=row4Col(curCol);
        
        if(curRow<=numRow)
            gain=gain+COrig(curRow,curCol);
            tuples(:,curTuple)=[curRow;curCol];
            curTuple=curTuple+1;

           if(curTuple>k)
               break;
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
