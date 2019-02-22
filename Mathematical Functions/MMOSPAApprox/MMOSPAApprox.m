function [MMOSPAEst,orderList]=MMOSPAApprox(x,w,numScans)
%%MMOSPAAPPROX Find the approximate minimum mean optimal sub-pattern
%              assignment (MMOSPA) estimate from a set of weighted 
%              discrete sets of target estimates using 2D assignment in a
%              forward-backward algorithm.
%
%INPUTS: x An xDim X numTar X numHyp hypermatrix that holds numHyp
%          hypotheses each consisting or numTar targets (or generic
%          vectors) with xDim dimensions per target (per generic vector).
%          One should not pass x values that are a combination of position
%          and velocity components, because the OSPA error in that instance
%          has no meaning. In general, the components of each target vector
%          for each hypothesis will be position only.
%        w A numHyp X 1 vector of the probabilities of each of the numHyp
%          hypotheses in x. The elements must all be positive and sum to
%          one.
% numScans An optional parameter >=1 specifying how many forward scans of
%          the approximate algorithm to perform. Even though more scans can
%          improve the estimate, global convergence is not guaranteed. The
%          default is 1 if this parameter is not provided.
%
%OUTPUTS: MMOSPAEst The approximate xDim X numTar MMOSPA estimate.
%         orderList A numTarXnumHyp matrix specifying the ordering of the
%                   targets in each hypothesis that went into the
%                   approximate MMOSPA estimate.
%
%Given a set of numHyp hypotheses, the standard expected value minimizes
%the mean squared error. The standard expected value is just
%the weighted sum of the hypothese for all of the targets. In the MMOSPA
%estimate, a weighted sum is used, but the ordering of the target states
%is modified to minimize the expected value of a specific version of the
%OSPA metric. orderList holds indices of the reordered states going into
%the MMOSPA estimate such that the MMOSPA estimate uses  
%x(:,orderList(:,curHyp),curHyp) instead of x(:,:,curHyp) for the ordered
%matrix of target states in hypothesis curHyp. Finding the MMOSPA estimate
%is generally NP-hard for more than two hypotheses. Thus, this algorithm
%approximates the MMOSPA estimate.
%
%The basic algorithm uses sequential 2D assignment going forward to
%approximate the MMOSPA estimate. If desired, assignments can be
%reevaluated in additional backward-forward passes to try to obtain an
%approximation close to the true MMOSPA estimate.
%
%The algorithm as well as the concept of MOSPA error are described in
%detail in [1].
%
%REFERENCES:
%[1] D. F. Crouse, "Advances in displaying uncertain estimates of multiple
%    targets," in Proceedings of SPIE: Signal Processing, Sensor Fusion,
%    and Target Recognition XXII, vol. 8745, Baltimore, MD, Apr. 2013.
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    if(nargin<3)
        numScans=1;
    end

    %First, get the forward solution, which might be bad.
    [MMOSPAEst,orderList]=MMOSPAApproxForward(x,w);
    numScans=numScans-1;
    
    numHyp=size(x,3);
    %Do a reverse and then a forward scan numScans times.
    while(numScans>0)
        %Re-evaluate the hypotheses going backwards.
        for curHyp=(numHyp-1):-1:2
            [MMOSPAEst, orderList]=doUpdate4Col(curHyp, x, w, MMOSPAEst,orderList);
        end
        
        %Re-evaluate the hypotheses going forwards.
        for curHyp=3:1:numHyp
            [MMOSPAEst, orderList]=doUpdate4Col(curHyp, x, w, MMOSPAEst,orderList);
        end
        
        numScans=numScans-1;
    end
end

function [MMOSPAEst,orderList]=MMOSPAApproxForward(x,w)
%%MMOSPAAPPROXFORWARD Perform a single forward step of the approximate
%                     MMOSPA algorithm without having any previous
%                     estimate.
%
%INPUTS: x An xDim X numTar X numHyp hypermatrix that holds numHyp
%          hypotheses each consisting or numTar targets (or generic
%          vectors) with xDim dimensions per target (per generic vector).
%        w A numHyp X 1 vector of the probabilities of each of the numHyp
%          hypotheses in x. The elements must all be positive and sum to
%          one.
%
%OUTPUTS: MMOSPAEst The approximate xDim X numTar MMOSPA estimate.
%         orderList A numTarXnumHyp matrix specifying the ordering of the
%                   targets in each hypothesis that went into the
%                   approximate MMOSPA estimate. 
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.

    numTar=size(x,2);
    numHyp=size(x,3);
    
    orderList=zeros(numTar,numHyp);
    
    %The initial hypothesis ordering determines the ultimate target
    %ordering.
    MMOSPAEst=w(1)*x(:,:,1);
    
    orderList(:,1)=1:numTar;
    c=zeros(numTar,numTar);%Allocate space for the cost matrix.
    for curHyp=2:numHyp
        %First, we have to calculate the MMOSPA mean up to this point.
        xM=MMOSPAEst;
        
        %Now, we must fill the cost matrix.
        %Rows are from xM.
        %Columns are from x.
        for curRow=1:numTar
           for curCol=1:numTar 
               c(curRow,curCol)=sum(xM(:,curRow).*x(:,curCol,curHyp));
           end
        end
        
        %This turns the maximization optimization into a minimization.
        c=-c;
        
        %Use the shortest path algorithm to find the optimal 2D assignment.
        orderList(:,curHyp)=assign2D(c);
        
        %Update the running MMOSPA estimate
        MMOSPAEst=MMOSPAEst+w(curHyp)*x(:,orderList(:,curHyp),curHyp);
    end
end

function [MMOSPAEst, orderList]=doUpdate4Col(varHyp, x, w, MMOSPAEst,orderList)
%%DOUPDATE4COL Fixing the ordering of all hypotheses except varHyp,
%              reevaluate the ordering of the targets in hypothesis varHyp
%              to see whether the MOSPA cost function can be reduced.
%
%INPUTS: varHyp The index of the hypothesis in x where the target ordering
%               can be changed.
%            x  An xDim X numTar X numHyp hypermatrix that holds numHyp
%               hypotheses each consisting or numTar targets (or generic
%               vectors) with xDim dimensions per target (per generic
%               vector).
%             w A numHyp X 1 vector of the probabilities of each of the
%               numHyp hypotheses in x. The elements must all be positive
%               and sum to one.
%     MMOSPAEst The current approximation to the MMOSPA estimate that was
%               found using the state orderings in orderList.
%     orderList A numTarXnumHyp matrix specifying the ordering of the
%               targets in each hypothesis that went into the current
%               approximate MMOSPA estimate.
%
%OUTPUTS: MMOSPAEst The updated approximate MMOSPA estimate.
%         orderList The updated order list corresponding to the updated
%                   approximate MMOSPA estimate
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.

        numTar=size(x,2);
        
        xOptCur=(MMOSPAEst-x(:,orderList(:,varHyp),varHyp)*w(varHyp));

        %Fill the cost matrix
        %Rows are from xOptCur.
        %Columns are from x.
        c=zeros(numTar,numTar);
        for curRow=1:numTar
           for curCol=1:numTar 
               c(curRow,curCol)=sum(xOptCur(:,curRow).*x(:,curCol,varHyp));
           end
        end
        
        %This turns the maximization optimization into a minimization.
        c=-c;
              
        %Use the shortest path algorithm to find the optimal 2D assignment.
        orderList(:,varHyp)=assign2D(c);
        
        %Update the MMOSPA estimate approximation.
        MMOSPAEst=xOptCur+w(varHyp)*x(:,orderList(:,varHyp),varHyp);
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
