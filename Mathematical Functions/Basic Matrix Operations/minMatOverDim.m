function M=minMatOverDim(C,minIdx)
%%MINMATOVERDIM Given a multidimensional matrix C, this function returns
%           matrix M that is obtained by minimizing C over the given
%           dimension and removing the given dimension. This is equivalent
%           to using a call to min(C,[],minIdx) in Matlab. However, here,
%           we do it without using any matrix/ vector commands so that it
%           can mirror how one might implement such a function in C.
%
%INPUTS: C A real n1Xn2Xn3..XnS S-dimensional hypermatrix.
%   minIdx The real index from 1 to S over which the minimization should be
%          performed.
%
%OUTPUTS: M An n1X...n(minIdx-1)X1Xn(minIdx+1)X...nS matrix holding the
%           minimum values over the specified dimension.
%
%The algorithm consists of two types of steps. First, a step in the linear
%indexation of C to go from one value of an index in spot minIdx is
%prod(nVals(1:(minIdx-1))), where nVals=size(C). Thus, we can go from one
%element to the next for the minimization. For a 3D matrix, the step to go
%from C(i1,1,i2) to C(i1+1,1,i2) increments by 1. However, to go from 
%C(n1,1,i2) to C(1,1,i2+1) is a step of prod(nVals(1:minIdx)). The same
%type of thing applies to a matrix with more dimensions, since all
%dimensions before the one in question and after can be collapsed into 1
%dimension. Thus, this function puts the above rules together to go
%minimize the matrix.
%
%EXAMPLE:
%Here, we just show that the results are equivalent to using the min
%command with a resize.
% C=randn(13,18,11,6,9);
% minIdx=3;
% M=minMatOverDim(C,minIdx);
% MAlt=min(C,[],minIdx);
% all(M(:)==MAlt(:))
%The result is 1, indicating that the two values are equal.
%
%March 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

nVals=size(C);
S=numel(nVals);

if(minIdx<1)
    error('Invalid value of minIdx given.')
end

if(minIdx>S)
    M=C;
    return;
end

MDims=nVals;
MDims(minIdx)=1;
M=zeros(MDims);

totalNumElsM=prod(MDims);

%Note that prod([])=1
incrMinIdx=prod(nVals(1:(minIdx-1)));

incrBigStep=incrMinIdx*nVals(minIdx);

CStartIdx=0;
curMCol=0;
for curEl=0:(totalNumElsM-1)
    curIdx=CStartIdx;
    minVal=C(curIdx+1);

    curIdx=curIdx+incrMinIdx;
    for i=2:nVals(minIdx)
        curVal=C(curIdx+1);
        
        if(curVal<minVal)
            minVal=curVal;
        end

        curIdx=curIdx+incrMinIdx;
    end
    M(curEl+1)=minVal;
    
    %If a big step has to be taken.
    if(mod(curEl+1,incrMinIdx)==0)
        curMCol=curMCol+1;
        
        CStartIdx=incrBigStep*curMCol;
    else
        CStartIdx=CStartIdx+1;
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
