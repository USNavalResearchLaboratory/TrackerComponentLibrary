function [rankVals,sumVals]=rankOrder(xMat,isSorted,colVars)
%%RANORDER Given an xDimX1 vector of values x, return a length xDim vector
%          where each spot contains the rank (ordered position) of the
%          value. When multiple values are equal (tied ranks), rather than
%          the positions getting an arbitrary ordering of ranks, they are
%          all set to the same rank that is the average of the span of
%          ranks. For example, if the third and forth things in order are
%          tied, then both are assigned rank (3+4)/2=3.5. This type of
%          ranking is Kendall's handling, as discusse din Chapter 3 of [1].
%          It arises in computing various robust correlation coefficients.
%
%INPUTS: xMat An xDimXN matrix of N real vectors, each of which whose rank
%             order is desired. However, if colVars=true, then this is an
%             NXxDim matrix.
%    isSorted An optional parameter indicating whether the columns of xMat
%             are already sorted. The algorithm is faster if they have
%             already been sorted. The default if omitted or an empty
%             matrix is passed is false.
%     colVars If colVars=false, then xMat is xDimXN. If colVars=true, then
%             xMat=NXxDim. The default if omitted or an empty matrix is
%             passed is false.
%
%OUTPUTS: rankVals An xDimXN matrix (or NXxDim matrix if colVars=true)
%                  where the values in each row (column) are the ranks of
%                  the corresponding values in xMat.
%          sumVals Let f_i be the number of tied values in the ith group of
%                  ties for a particular x vector. Then sum_i f_i^3-f_i is
%                  a value that arises in the computation of the Spearman
%                  rank correlation, as discussed in Chapter 3.8 of [1].
%                  sumVals(k) is the value of this sum for the kth column
%                  of xMat. This is a length xDim vector that is 1XN if
%                  colVars=false and NX1 if colVars=true.
%
%EXAMPLE:
%Here, we have a vector with no ties and a vector with ties:
% xMat=[ 4,   90;
%        6, -120;
%       12,   82;
%        9,   64;
%       12,   pi];
% [rankVals,sumVals]=rankOrder(xMat)
%One will get
%rankVals=[1.0, 5.0;
%          2.0, 1.0;
%          4.5, 4.0;
%          3.0, 3.0;
%          4.5, 2.0];
%sumVals=[6;0]
% 
%REFERENCES:
%[1] M. G. Kendall, Rank Correlation Methods, 3rd ed. LiverPool: Charles
%    Birch and Sons Ltd., 1962.
%
%August 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(isSorted))
    isSorted=false;
end

if(nargin<3||isempty(colVars))
    colVars=false;
end

if(isempty(xMat))
    rankVals=[];
    return
end

if(colVars)
    xDim=size(xMat,2);
    N=size(xMat,1);

    rankVals=zeros(N,xDim);
    sumVals=zeros(N,1);
    
    for curVec=1:N
        if(isSorted==false)
            [x,idxOrd]=sort(xMat(curVec,:),'ascend');
        else
            x=xMat(curVec,:);
        end

        idx=1;
        while(idx<xDim)
            if(x(idx+1)~=x(idx))%If there is not a tie.
                rankVals(curVec,idx)=idx;
                idx=idx+1;
            else%If there is a tie.
                %Find the last tied element.
                for k=(idx+1):xDim
                    if(x(k)~=x(idx))
                       break;
                    end
                end
                %The special case where the last element in the vector equals the
                %one before it.
                if(x(k)==x(idx))
                   k=k+1; 
                end

                rank=0.5*(idx+k-1);%The average rank.
                %All tied elements get the average rank.
                rankVals(curVec,idx:(k-1))=rank;

                f=k-idx;
                sumVals(curVec)=sumVals(curVec)+f*(f^2-1);

                idx=k;
            end
        end

        %The rank of the final element, if it was not in a tie. If it was in a tie,
        %then idx will be xDim+1.
        if(idx==xDim)
            rankVals(curVec,idx)=xDim;
        end

        if(isSorted==false)
            rankVals(curVec,:)=rankVals(curVec,idxOrd);
        end
    end
else
    xDim=size(xMat,1);
    N=size(xMat,2);
    
    rankVals=zeros(xDim,N);
    sumVals=zeros(1,N);
    for curVec=1:N
        if(isSorted==false)
            [x,idxOrd]=sort(xMat(:,curVec),'ascend');
        else
            x=xMat(:,curVec);
        end

        idx=1;
        while(idx<xDim)
            if(x(idx+1)~=x(idx))%If there is not a tie.
                rankVals(idx,curVec)=idx;
                idx=idx+1;
            else%If there is a tie.
                %Find the last tied element.
                for k=(idx+1):xDim
                    if(x(k)~=x(idx))
                       break;
                    end
                end
                %The special case where the last element in the vector equals the
                %one before it.
                if(x(k)==x(idx))
                   k=k+1; 
                end

                rank=0.5*(idx+k-1);%The average rank.
                %All tied elements get the average rank.
                rankVals(idx:(k-1),curVec)=rank;

                f=k-idx;
                sumVals(curVec)=sumVals(curVec)+f*(f^2-1);

                idx=k;
            end
        end

        %The rank of the final element, if it was not in a tie. If it was in a tie,
        %then idx will be xDim+1.
        if(idx==xDim)
            rankVals(idx,curVec)=xDim;
        end

        if(isSorted==false)
            idxOrd=inversePermutation(idxOrd);
            rankVals(:,curVec)=rankVals(idxOrd,curVec);
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
