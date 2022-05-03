function JPost=PCRLBUpdateAddNoClutter(JPred,xCur,PCur,R,PD,paramH,xi,w)
%PCRLBUPDATEADDNOCLUTTER Update the Fisher information matrix (FIM) for a
%            measurement with additive noise (in the coordinate system of
%            the measurement), a possibly non-unity detection probability,
%            but with no clutter. The FIM for the information reduction
%            factor PCRLB is used.
%
%INPUTS: JPred The xDimXxDim Fisher information matrix predicted forward to
%              the current time-step. If this is the first step with a
%              measurement, then JPred should be the xDimXxDim zero matrix.
%         xCur The xDimX1 mean of the distribution of the true (but unknown
%              to the tracker) possible target location at the current
%              time. It is assumed that the distribution of the target
%              locations is Gaussian. This is only needed if paramH is a
%              function. Otherwise, an empty matrix can be passed.
%         PCur The xDimXxDim covariance of the distribution of the true
%              (but unknown to the tracker) possible target location at the
%              current time. If the target motion is deterministic but
%              unknown to the tracker, then PCur is a matrix of zeros. If
%              an empty matrix is passed, it is assumed that PCur is all
%              zeros.
%            R The zDimXzDim covariance matrix of the additive measurement
%              noise.
%           PD The scalar detection probability of the target at the
%              current step.
%       paramH Either the fixed measurement matrix H (which is zDimXxDim),
%              or a function handle HJacob to get the measurement Jacobian,
%              which is zDimXxDim where each row corresponds to a
%              component of the measurement vector and the column is the
%              derivative of the row component with respect the part of the
%              state vector corresponding to the column. The function must
%              take the state as its parameter.
%        xi, w If a function handle is passed for paramH and PCur is not a
%              zero matrix, then cubature points for integration are
%              needed. xi and w are the cubature points and weights. The
%              cubature points must be the dimensionality of the state. If
%              these parameter are omitted or empty matrices are passed,
%              then fifthOrderCubPoints(xDim) is used. It is suggested that
%              xi and w be provided (if needed) to avoid needless
%              recomputation of the cubature points.
%
%OUTPUT: JPost The xDimXxDim Fisher information matrix after a measurement
%              update.
%
%A basic description of using cubature integration for evaluating the PCRLB
%is given in [1], where it is also described how to use a prior
%distribution to evaluate the PCRLB for a simulation of random tracks
%without running Monte Carlo runs. The prior distribution parameters in
%such an instance can be calculated using the function DiscPriorPModel.
%
%This function implements the information reduction factor (IRF) form of
%the PCRLB for possible missed detections, which is described in Section
%III-A of [2]. The bound from this method is not as tight as the
%measurement sequence conditioned approach, which is also mentioned in [2].
%
%If multiple independent measurement are available at a particular time,
%then this function can be called multiple times to sequentially update the
%Fisher information matrix.
%
%If the state prediction model has additive process noise, then the
%PCRLBPredAdd function can be used to predict the Fisher information matrix
%forward through the dynamic model.
%
%REFERENCES:
%[1] David F. Crouse , "Basic tracking using nonlinear 3D monostatic and
%    bistatic measurements," IEEE Aerospace and Electronic Systems 
%    Magazine, vol. 29, no. 8, Part II, pp. 4-53, Aug. 2014.
%[2] M. Hernandez, B. Ristic, A. Farina, and L. Timmoneri, "A
%    comparison of two Cramer-Rao bounds for nonlinear filtering with
%    Pd<1," IEEE Transactions on Signal Processing, vol. 52, no. 9, pp.
%    2361-2370, Sep. 2004.
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    xDim=size(JPred,1);

    if(isa(paramH,'function_handle'))
        if(isempty(PCur)||isequal(zeros(xDim,xDim),PCur))
            H=paramH(xCur);
            constMeasMat=true;
        else
            HJacob=paramH;
            
            if(nargin<8||isempty(xi))
                [xi,w]=fifthOrderCubPoints(xDim);
            end
            
            constMeasMat=false;
        end
    else
        H=paramH;
        constMeasMat=true;
    end
    
    RInv=pinv(R);
    if(constMeasMat)
        JPost=JPred+PD*H'*RInv*H;
    else
        xPoints=transformCubPoints(xi,xCur,chol(PCur,'lower'));
        numPoints=size(xPoints,2);
        JPost=zeros(size(JPred));
        for curP=1:numPoints
            H=HJacob(xPoints(:,curP));
            JPost=JPost+w(curP)*(H'*RInv*H);
        end
        JPost=JPred+PD*JPost;
    end

    %Ensure symmetry is preserved.
    JPost=(JPost+JPost')/2;
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
