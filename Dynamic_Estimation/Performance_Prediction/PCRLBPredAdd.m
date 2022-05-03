function JPred=PCRLBPredAdd(JPrior,xPrior,PPrior,Q,param5,xi,w)
%PCRLBPREDADD Update the Fisher information matrix (FIM) due to propagation
%             over a discrete time step with additive process noise.
%
%INPUTS: JPrev The xDimXxDim Fisher information matrix at the previous time
%              step. At time k=0, this should normally be the zero matrix,
%              since no information is available to the tracking algorithm.
%              However, there is no point in propagating forward an all
%              zero FIM; it will remain all zero until a measurement
%              arrives.
%       xPrior The xDimX1 mean of the distribution of the true
%              (but unknown to the tracker) possible target location at the
%              previous time. It is assumed that the distribution of the
%              target location is Gaussian. If param5  is not a function
%              handle, then xPrior is not used and an empty matrix can be
%              passed.
%       PPrior The xDimXxDim covariance of the distribution of the true
%              (but unknown to the tracker) possible target location at the
%              previous time. If the target motion is deterministic but
%              unknown to the tracker, then PPrior is a matrix of zeros. If
%              param5 is not a function handle, then PPrior is not used and
%              an empty matrix can be passed.
%            Q The xDimXxDim covariance matrix of the additive process
%              noise.
%       param5 Either the fixed xDimXxDim state transition matrix F, or a
%              function handle FJacob to get the state transition Jacobian
%              such that FJacob(x)*x is a linear approximation to the
%              transition of the state x at the given time. The function
%              must take the state as its parameter.
%        xi, w If a function handle is passed for param5 and PPrior is not
%              a zero matrix, then cubature points for evaluating expected
%              values are needed. xi and w are the cubature points and
%              weights. The cubature points must be the dimensionality of
%              the state. If these parameter are omitted or empty matrices
%              are passed, then fifthOrderCubPoints(xDim) is used. It is
%              suggested that xi and w be provided (if needed) to avoid
%              needles recomputation of the cubature points.
%
%OUTPUT: JPred The Fisher information matrix propagated forward in time
%              without a measurement.
%
%A basic description of using cubature integration for evaluating the PCRLB
%is given in [1], where it is also described how to use a prior
%distribution to evaluate the PCRLB for a simulation of random tracks
%without running Monte Carlo runs. The prior distribution parameters in
%such an instance can be calculated using the function DiscPriorPModel.
%
%REFERENCES:
%[1] David F. Crouse , "Basic tracking using nonlinear 3D monostatic and
%    bistatic measurements," IEEE Aerospace and Electronic Systems 
%    Magazine, vol. 29, no. 8, Part II, pp. 4-53, Aug. 2014.
%
%October 2013 David F.Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    xDim=size(JPrior,1);

    if(isa(param5,'function_handle'))
        if(isequal(zeros(xDim,xDim),PPrior))
            F=param5(xPrior);
            constTransMat=true;
        else
            FJacob=param5;
            
            if(nargin<8||isempty(xi))
                [xi,w]=fifthOrderCubPoints(xDim);
            end
            
            constTransMat=false;
        end
    else
        F=param5;
        constTransMat=true;
    end
    
    QInv=pinv(Q);
    if(constTransMat==true)
        D12=-F'*QInv;
        D11=F'*QInv*F;
    else
        xPoints=transformCubPoints(xi,xPrior,chol(PPrior,'lower'));
        numPoints=size(xPoints,2);
        D12=0;
        D11=0;
        for curP=1:numPoints
            F=FJacob(xPoints(:,curP));
            D12=D12-w(curP)*F'*QInv;
            D11=D11+w(curP)*F'*QInv*F;
        end
    end
    JPred=QInv-D12'*pinv(JPrior+D11)*D12;
    %Ensure symmetry is preserved.
    JPred=(JPred+JPred')/2;
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
