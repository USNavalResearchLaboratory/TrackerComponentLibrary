function val=calcMOSPAError(xEst,x,w)
%%CALCMOSPAERROR Given an estimate of a set of targets as well as a set of
%                discrete state hypotheses with probabilities of being
%                correct, compute the mean optimal sub-pattern assignment 
%                error of the estimate.
%
%INPUTS: xEst The xDim X numTar estimate of numTar targets (or separate
%             vector values).
%           x An xDim X numTar X numHyp hypermatrix that holds numHyp
%             hypotheses each consisting or numTar targets (or generic
%             vectors) with xDim dimensions per target (per generic
%             vector). One should not pass x values that are a combination
%             of position and velocity components, because the OSPA error
%             in that instance has no meaning. In general, the components
%             of each target vector for each hypothesis will be position
%             only.
%           w A numHyp X 1 vector of the probabilities of each of the
%             numHyp hypotheses in x. The elements must all be positive and
%             sum to one.
%
%OUTPUTS: val The scalar MOSPA error.
%
%The definition of the OSPA statistic used here assumes a known number of
%targets, p=2 and the other measure being the square of the L2 norm. In
%other words, if the step of the algorithm changing the ordering of the
%targets were omitted, the OSPA statistic is just the root mean squared
%error summed over all targets.
%
%The concept of MOSPA error is described in detail in [1].
%
%REFERENCES:
%[1] D. F. Crouse, "Advances in displaying uncertain estimates of multiple
%    targets," in Proceedings of SPIE: Signal Processing, Sensor Fusion,
%    and Target Recognition XXII, vol. 8745, Baltimore, MD, Apr. 2013.
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numTar=size(x,2);
numHyp=size(x,3);
c=zeros(numTar,numTar);%Allocate space for the cost matrix.

val=0;
for curHyp=1:numHyp
    %Fill the cost matrix
    %Rows are from xEst.
    %Columns are from x.
    for curRow=1:numTar
       for curCol=1:numTar 
           c(curRow,curCol)=sum(xEst(:,curRow).*x(:,curCol,curHyp));
       end
    end
    %This turns the maximization optimization into a minimization.
    c=-c;
        
    %Use the shortest path algorithm to find the optimal 2D assignment.
    orderList=assign2D(c);
    
    diff=xEst-x(:,orderList,curHyp);
    val=val+w(curHyp)*sum(sum(diff.*diff));
end

val=sqrt(val/numTar);
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
