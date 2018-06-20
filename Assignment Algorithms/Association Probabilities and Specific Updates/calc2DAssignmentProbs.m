function beta=calc2DAssignmentProbs(A,diagAugment)
%%CALC2DASSIGNMENTPROBS Given a matrix of all-positive likelihoods or
%                likelihood ratios, determine the probability that each row
%                (target) is assigned to each column (measurement). Whereas
%                the assign2D algorithm provides the optimal assignment of
%                all rows to columns, this provides the probability of each
%                assignment. The diagAugment parameter allows for a more
%                efficient algorithm to be used in the case that for a
%                numRow X numCol matrix where numCol>=numRow, the last
%                numRow columns are a weighted identity matrix. Such a
%                matrix structure arises when performing target-measurement
%                data association with missed detection hypotheses. In such
%                an instance, the rows are targets, the first numCol-numRow
%                columns correspond to measurements and the final numRow
%                columns have missed detection hypotheses.
%
%INPUTS: A A matrix of positive likelihoods or likelihood ratios (NOT log-
%          likelihood ratios). If diagAugment=true, then A is a
%          numTarX(numMeas+numTar) matrix of all-positive likelihoods or
%          likelihood ratios for assigning the target specified by the row
%          to the measurement specified by the column. Columns > numMeas
%          hold missed-detection likelihoods/ likelihood ratios. Thus, off-
%          diagonal terms for columns > numMeas should be set to 0 and the
%          diagonal terms set to the costs of a missed detection for each
%          given target.
% diagAugment A boolean variable indicating whether the probabilties
%          should be found for a general assignment problem or assuming A
%          ends with a diagonal submatrix of missed detection
%          probabilities. The default if omitted is false (the general
%          problem). Setting diagAugment to true changes the shape of the
%          output. The default if omitted is false. See the description of
%          the output beta for how this affects the output.
%
%OUTPUTS: beta If diagAugment is omitted or false, then beta has the same
%              dimensionality as A and hold the probability of assigning
%              each row to each column in the traditional 2D assignment
%              problem (Each row must be assigned to at least one column,
%              each column can be assigned to at most one row, or vice
%              versa if there are more rows than columns. If diagAugment
%              is true, then beta is a numTar X (numMeas+1) matrix of
%              probabilities of assigning the target given by the row to
%              the measurement given by the column. The final column is a
%              set of missed detection probabilities. Thus, in this case,
%              beta has fewer columns than A, because the final numRow
%              columns are collapsed into one column on return.
%
%The notion of using matrix permanents for evaluating 2D assignment
%probabilities, with a focus on target-measurement association
%probabilities is from [1]. The generalization to missed detections is as
%simple as the generalization of the matrix permanent to a rectangular
%matrix. The generalization is discussed in [2].
%
%The concept of the target-measurement association probability in the JPDAF
%tracker is discussed in Chapter 6.2 of [3].
%
%REFERENCES:
%[1] J. K. Uhlmann, "Matrix permanent inequalities for approximating joint 
%    assignment matrices in tracking systems," Journal of the Franklin 
%    Institute, vol. 341, pp. 569-593, 2004.
%[2] D. F. Crouse and P. Willett, "Computing Target-Measurement Association
%    Probabilities Using the Matrix Permanent," IEEE Transactions on
%    Aerospace and Electronic Systems, submitted January 2016.
%[3] Y. Bar-Shalom, P. K. Willett, and X. Tian, Tracking and Data Fusion.
%    Storrs, CT: YBS Publishing, 2011.
%
%September 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2)
    diagAugment=false;
end

if(diagAugment==true)
    numTar=size(A,1);
    numCol=size(A,2);
    numMeas=size(A,2)-numTar;

    boolRowsSkip=false(numTar,1);
    boolColsSkip=false(numCol,1);

    beta=zeros(numTar,numMeas+1);
    for curTar=1:numTar
        boolRowsSkip(curTar)=true;

         for curMeas=1:numMeas
             ati=A(curTar,curMeas);
             %The A matrix removing the row and column corresponding to the
             %selected association.
             boolColsSkip(curMeas)=true;
             beta(curTar,curMeas)=ati*perm(A,boolRowsSkip,boolColsSkip);

             boolColsSkip(curMeas)=false;
         end
         %The missed detection hypothesis
         curMeas=numMeas+curTar;
         ati=A(curTar,curMeas);

         boolColsSkip(curMeas)=true;
         beta(curTar,numMeas+1)=ati*perm(A,boolRowsSkip,boolColsSkip);

         boolColsSkip(curMeas)=false;
         boolRowsSkip(curTar)=false;
    end
else
    %If we are here, then we just want general assignment probabilities,
    %not specialized to target tracking applications.
    numRow=size(A,1);
    numCol=size(A,2);
    boolRowsSkip=false(numRow,1);
    boolColsSkip=false(numCol,1);
    
    beta=zeros(numRow,numCol);
    for curRow=1:numRow
        boolRowsSkip(curRow)=true;

         for curCol=1:numCol
             ati=A(curRow,curCol);
             %The A matrix removing the row and column corresponding to the
             %selected association.
             boolColsSkip(curCol)=true;
             beta(curRow,curCol)=ati*perm(A,boolRowsSkip,boolColsSkip);

             boolColsSkip(curCol)=false;
         end
         boolRowsSkip(curRow)=false;
    end
end

%It is faster to normalize the betas this way then to compute the
%normalization constant by finding the permanent of the entire A matrix.
beta=bsxfun(@rdivide,beta,sum(beta,2));
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
