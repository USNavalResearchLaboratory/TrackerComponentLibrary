function JPost=PCRLBUpdateAddNoClutter(JPred,xCur,PCur,R,PD,param6,xi,w)
%PCRLBUPDATEADDNOCLUTTER Update the Fisher information matrix (FIM) for a
%            measurement with additive noise (in the coordinate system of
%            the measurement), a possibly non-unity detection probability,
%            but with no clutter. The FIM for the information reduction
%            factor PCRLB is used.
%
%INPUTS: JPred The Fisher information matrix predicted forward to the
%              current time-step. If this is the first step with a
%              measurement, then JPred should be the zero matrix.
%         xCur The (column vector) mean of the distribution of the true
%              (but unknown to the tracker) possible target location at the
%              current time. It is assumed that the distribution of the
%              target locations is Gaussian.
%         PCur The covariance of the distribution of the true (but unknown
%              to the tracker) possible target location at the current
%              time. If the target motion is deterministic but unknown to
%              the tracker, then PCur is a matrix of zeros.
%            R The covariance matrix of the additive measurement noise.
%           PD The detection probability of the target at the current step.
%       param6 Either the fixed measurement matrix H (which is zDimXxDim),
%              or a function handle HJacob to get the measurement Jacobian,
%              which is zDim X xDim where each row corresponds to a
%              component of the measurement vector and the column is the
%              derivative of the row component with respect the part of the
%              state vector corresponding to the column. The function must
%              take the state as its parameter.
%        xi, w If a function handle is passed for param6 and PCur is not a
%              zero matrix, then cubature points for integration are
%              needed. xi and w are the cubature points and weights. The
%              cubature points must be the dimensionality of the state.
%
%OUTPUT: JPost The Fisher information matrix after a measurement update.
%
%A basic description of using cubature integration for evaluating the PCRLB
%is given in [1], where it is also described how to use a prior
%distribution to evaluate the PCRLB for a simulation of random tracks
%without running Monte Carlo runs. The prior distribution parameters in
%such an instance can be calculated using the function DiscPriorPModel.
%
%The general form of the PCRLB measurement update step with a non-unity
%detection probability but no clutter is from [2].
%
%If multiple independent measurement are available at a particular time,
%then this function can be called multiple times to sequentially update the
%Fisher information matrix.
%
%REFERENCES:
%[1] David F. Crouse , "Basic tracking using nonlinear 3D monostatic and
%    bistatic measurements," IEEE Aerospace and Electronic Systems 
%    Magazine, vol. 29, no. 8, Part II, pp. 4-53, Aug. 2014.
%[2] P. Stinco, M. S. Greco, F. Gini, and A. Farina, "Posterior Cramér-Rao
%    lower bounds for passive bistatic radar tracking with uncertain target
%    measurements," Signal Processing, vol. 93, no. 12, pp. 3528-3540, Dec.
%    2013.
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    xDim=size(xCur,1);

    if(isa(param6,'function_handle'))
        if(isequal(zeros(xDim,xDim),PCur))
            H=param6(xCur);
            constMeasMat=true;
        else
            HJacob=param6;
            constMeasMat=false;
        end
    else
        H=param6;
        constMeasMat=true;
    end
    
    RInv=pinv(R);
    if(constMeasMat)
        JPost=JPred+PD*H'*RInv*H;
    else
        xPoints=transformCubPoints(xi,xCur,chol(PCur,'lower'));
        numPoints=size(xPoints,2);
        JPost=JPred;
        for curP=1:numPoints
            H=HJacob(xPoints(:,curP));
            JPost=JPost+w(curP)*(H'*RInv*H);
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
