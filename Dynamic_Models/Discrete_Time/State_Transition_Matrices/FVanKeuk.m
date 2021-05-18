function F=FVanKeuk(T,tau,numDim)
%%FVANKEUK  Get the state transition matrix for the van Keuk dynamic model.
%           This is a direct discrete-time model such that the acceleration
%           advances in each dimension over time as
%           a[k+1]=exp(-T/tau)a[k]+sigma*sqrt(1-exp(-2*T/tau))*v[k]
%           where T is the time between steps, a[k] and a[k+1] are the
%           accelerations at the two steps, v[k] is a standard normal
%           random variable and tau and sigma are parameters. The model is
%           akin to a direct discrete-time version of the Singer's
%           dynamic model.
%
%INPUTS: T The time-duration of the propagation interval in seconds.
%      tau The maneuver autocorrelation time in seconds. Presumably, this
%          parameter can be chosen according to the same criteria as in
%          Singer's dynamic model. That is, a reasonable range for tau is
%          between 5 and 20 seconds. This parameter can be a scalar if the
%          same parameter is used for all dimensions or a numDimX1 or
%          1XnumDim vector if different values are used across the
%          dimensions.
%   numDim The number of dimensions of the simulation problem. If the
%          numDim parameter is omitted, then numDim=3 (3D motion) is
%          assumed.
%
%OUTPUTS: F The state transition matrix under van Keuk's dynamic model in
%           numDim dimensions where the state is stacked
%           [position;velocity;acceleration].
%
%Van Keuk's dynamic model is described in Chapter 2.2.1 of [1] and is
%described in slightly more detail in [2]. The model is essentially a
%direct-discrete-time version of Singer's dynamic model.
%
%When describing Singer's dynamic model, Chapter 8.2.3 of [3] states that a
%typical value of tau for a slowly turning aircraft is 20s, and is 5s for
%an evasive maneuver. Presumably the same applies when using van Keuk's
%model.
%
%The autocorrelation of the acceleration is
%=sigma^2*exp(-(t_{k_i}-t_{k_j})/tau) where t_{k_i}-t_{k_j} are the times
%of discrete-time steps k_i and k_j. sigma is the acceleration bandwidth.
%van Keuk shows that the RMS acceleration at a given time is
%sigma*sqrt(2*T/theta) where T is the discrete-time sampling interval.
%Thus, the parameter sigma can be determined based on the expected average
%acceleration magnitude.
%
%The process noise covariance matrix associated with this model can be
%obtained using the function QVanKeuk.
%
%REFERENCES:
%[1] W. Koch,Tracking and Sensor Data Fusion: Methodological Framework and
%    Selected Applications. Heidelberg, Germany: Springer, 2014.
%[2] G. van Keuk, "Zielverfolgung, basierend auf dem Kalman-Bucy-filter,"
%    Applied Informatics, vol. 14, no. 7, pp. 302-308, 1972.
%[3] Y. Bar-Shalom, X. R. Li, and T. Kirubarajan, Estimation with
%    Applications to Tracking and Navigation. New York: John Wiley and
%    Sons, Inc, 2001.
%
%April 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(numDim))
    numDim=3;
end

if(isscalar(tau))
    tau=ones(numDim,1)*tau;
end

%Get the scalar F matrices for each dimension.
F1=zeros(3,3,numDim);
for curDim=1:numDim 
    %First, create the matrix for 1D motion
    %Equation 2.5 in Koch's book.
    F1(:,:,curDim)=[1, T, T^2/2;
                    0, 1, T;
                    0, 0, exp(-T/tau(curDim))];
end

%Now, the elements of each of the 1D matrices just need get spread across
%multiple dimensions 
xDim=3*numDim;

%Allocate space
F=zeros(xDim,xDim);
for curRow=1:3
    rIdxMin=(curRow-1)*numDim+1;
    rIdxMax=curRow*numDim;
    spanRow=rIdxMin:rIdxMax;
    for curCol=1:3
        cIdxMin=(curCol-1)*numDim+1;
        cIdxMax=curCol*numDim;
        spanCol=cIdxMin:cIdxMax;
        
        for curDim=1:numDim
            F(spanRow(curDim),spanCol(curDim))=F1(curRow,curCol,curDim);
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
