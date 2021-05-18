function Q=QVanKeuk(T,xDim,tau,sigma)
%%QVANKEUK  Get the process noise covariance matrix for a Van Keuk dynamic
%           model. This is a direct discrete-time model such that the
%           acceleration advances in each dimension over time as
%           a[k+1]=exp(-T/tau)a[k]+sigma*sqrt(1-exp(-2*T/tau))*v[k]
%           where T is the time between steps a[k] and a[k+1] are the
%           accelerations at the two steps, v[k] is a standard normal
%           random variable and tau and sigma are parameters. The model is
%           akin to a direct discrete-time version of the Singer's
%           dynamic model.
%
%INPUTS: T The time-duration of the propagation interval in seconds.
%     xDim The dimensionality of the target state.
%          xDim is ((order+1)*numDim)X1 dimensional, where numDim is the
%          number of position dimensions of space (e.g. 2D or 3D).
%      tau The maneuver autocorrelation time in seconds. Presumably,
%          this parameter can be chosen according to the same criteria as
%          in Singer's dynamic model. That is, a reasonable range for tau
%          is between 5 and 20 seconds. This parameter can be a scalar if
%          the same parameter is used for all dimensions or a numDimX1 or
%          1XnumDim vector if different values are used across the
%          dimensions.
%    sigma The acceleration bandwidth. van Keuk shows that the RMS
%          acceleration at a given time is sigma*sqrt(2*T/theta) where T is
%          the discrete-time sampling interval. Thus, the parameter sigma
%          can be determined based on the expected average acceleration
%          magnitude.
%   numDim The number of dimensions of the simulation problem. If the
%          numDim parameter is omitted, then numDim=3 (3D motion) is
%          assumed.
%
%OUTPUTS: Q The process noise covariance matrix under a van Keuk dynamic
%           model with motion in numDim dimensions where the state is
%           stacked [position;velocity;acceleration].
%
%Van Keuk's dynamic model is described in Chapter 2.2.1 of  and is
%described in slightly more detail in [2]. The model is essentially a
%direct-discrete-time version of Singer's dynamic model.
%
%When describing Singer's dynamic model, Chapter 8.2.3 of [3] states that a
%typical value of tau for a slowly turning aircraft is 20s, and is 5s for
%an evasive maneuver. Presumably the same applies when using van Keuk's
%model.
%
%Note that the process noise covariance matrix is singular and only injects
%noise into the acceleration component of the state each discrete
%time-step.
%
%The state transition matrix associated with this model can be obtained
%using the function FVanKeuk.
%
%REFERENCES:
%[1] W. Koch, Tracking and Sensor Data Fusion: Methodological Framework and
%    Selected Applications. Heidelberg, Germany: Springer, 2014.
%[2] G. van Keuk, "Zielverfolgung, basierend auf dem Kalman-Bucy-filter,"
%    Applied Informatics, vol. 14, no. 7, pp. 302-308, 1972.
%[3] Y. Bar-Shalom, X. R. Li, and T. Kirubarajan, Estimation with
%    Applications to Tracking and Navigation. New York: John Wiley and
%    Sons, Inc, 2001.
%
%April 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numDim=xDim/3;

if(isscalar(tau))
    tau=ones(numDim,1)*tau;
end

%Get the scalar Q matrices for each dimension.
Q1=zeros(3,3,numDim);
for curDim=1:numDim
    %Equation 2.6 in Koch's book.
    G=sigma*sqrt(1-exp(-2*T/tau(curDim)))*[0;0;1];
    Q1(:,:,curDim)=G*G';
end

%Now, the elements of each of the 1D matrices just need get spread across
%multiple dimensions
xDim=3*numDim;

%Allocate space
Q=zeros(xDim,xDim);
for curRow=1:3
    rIdxMin=(curRow-1)*numDim+1;
    rIdxMax=curRow*numDim;
    spanRow=rIdxMin:rIdxMax;
    for curCol=1:3
        cIdxMin=(curCol-1)*numDim+1;
        cIdxMax=curCol*numDim;
        spanCol=cIdxMin:cIdxMax;
        
        for curDim=1:numDim
            Q(spanRow(curDim),spanCol(curDim))=Q1(curRow,curCol,curDim);
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
