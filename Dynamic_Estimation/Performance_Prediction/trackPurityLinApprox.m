function [Pc,PInv]=trackPurityLinApprox(beta,H,F,R,Q,K,PInvInit)
%%TRACKPURITYLINAPPROX Approximate the average track purity and the final
%               track accuracy of a typical track in a cluster of closely-
%               spaced targets. it is assumed that the detection
%               probability is 100% and there are no false alarms. Track
%               purity is the probability that a measurement assigned to a
%               track is actually should be assigned to that track. Average
%               track purity is the average purity over a batch. Here, the
%               batch is taken to be the average purity from the time the
%               system is switched on going for K steps. This function
%               assumes linear dynamic and measurement models.
%
%INPUTS: beta The number of targets per unit volume in the coordinate
%             system of the measurements.
%           H The zDim X xDim measurement matrix such that H*x+w is the
%             measurement, where x is the state and w is zero-mean Gaussian
%             noise with covariance matrix R.
%           F The xDim X xDim state transition matrix The state at
%             discrete-time k+1 is modeled as F times the state at time k
%             plus zero-mean Gaussian process noise with covariance matrix
%             Q.
%           R The zDim X zDim measurement covariance matrix.
%           Q The xDim X xDim process noise covariance matrix. This can not
%             be a singular matrix.
%           K The number of scans (after the initial turn-on) over which
%             one assumes the tracker is run and over which the average
%             purity should be computed.
%    PInvInit This is an initial xDimXxDim inverse covariance matrix for
%             the tracks. One might set this to the inverse covariance
%             matrix one might use for single-point track intiation. One-
%             point differencing is discussed in [3] and Chapter 3.2.2 of
%             [4]. For example, one could use a diagonal matrix with zeros
%             for the position components and a variance for the other
%             (velocity etc.) components that makes the gate fit the
%             assumed maximum speed of a target.
%
%OUTPUTS: Pc The approximate average probability over time that a
%            particular track will have false assignments. This is the same
%            for all tracks.
%       PInv The inverse of the average target covariance matrix (which
%            would provide a figure of track accuracy) at time K (the
%            initial step is time 0) taking into account false
%            associations.  
%
%This function implements the algorithm of Section II of [1], with changes.
%Equation 5 of the model omits process noise. Rather than using Equation 5,
%information filter prediction and update routines from the flow chart in
%Appendix H of [2] are used to propagate the inverse covariance matrix. The
%necessary innovation covariance needed in the equations is just part of a
%standard Kalman filter and can be discered from [2].
%
%EXAMPLE:
% T=1/2;
% F=FPolyKal(T,6,1);
% q0=1;
% Q=QPolyKal(T,6,1,q0);
% R=diag([100;100;100]).^2;
% H=[1,0,0,0,0,0;
%    0,1,0,0,0,0;
%    0,0,1,0,0,0];
% K=5;%5 steps after the initial scan.
% beta=(1/1e3)^3;%1 target every square kilometer.
% PInvInit=diag([0;0;0;1/600^2;1/600^2;1/600^2]);
% [Pc,PInv]=trackPurityLinApprox(beta,H,F,R,Q,K,PInvInit)
% %One will get around a 65.77% average track purity figure.
%
%REFEREMCES:
%[1] S. Mori, K.-C. Chang, C.-Y. Chong, and K.-P. Dunn, "Prediction of
%    track purity and track accuracy in dense target environments," IEEE
%    Transactions on Automatic Control, vol. 40, no. 5, pp. 953-959, May
%    1995.
%[2] David F. Crouse , "Basic tracking using nonlinear 3D monostatic and
%    bistatic measurements," IEEE Aerospace and Electronic Systems 
%    Magazine, vol. 29, no. 8, Part II, pp. 4-53, Aug. 2014.
%[3] Mallick, M.,La Scala, B., "Comparison of single-point and two-point
%    difference track initiation algorithms using position measurements". 
%    Acta Automatica Sinica, 2008.
%[4] Y. Bar-Shalom, P. K. Willett, and X. Tian, Tracking and Data Fusion.
%    Storrs, CT: YBS Publishing, 2011.
%
%January 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

m=size(H,1);

%We use gammaln and exp rather than just directly taking the ratio of the
%gamma functions so as the lessen possible overlow errors.
gammaRatio=exp(gammaln((m+1)/2)-gammaln(m/2+1));
%Equation 3.
Cm=2^(m-1)*pi^((m-1)/2)*gammaRatio;

%Initial information matrix is based on a measurement without any
%misassociation error.
PInv=PInvInit+H'*inv(R)*H;

Pc=0;
for k=1:K
    %Predict the state covariance matrix forward in time.
    DInv=F'+PInv/(F)*Q;
    PInv=DInv\PInv/(F);
    
    %We initially compute innovation covariance matrix Sk (its determinant
    %is sigmak^(2*m) assuming correct association.
    Qk=H*inv(PInv)*H';
    Sk=R+Qk;
    
    %Next, we use the computed value of sigmak to compute a modified R
    %matrix as in Equation 7. This R matrix is then used to compute a new
    %sigmak.
    Omega=(2/(m+2))*((m+1)*Qk-R);%Given on page 956.
    xi=Cm*beta*sqrt(det(Sk));
    phi=xi*exp(-xi);%Given on page 956.
    
    %Equation 7.
    RCur=R+Omega*phi;
    %Now, we compute a new innovation covariance
    Sk=RCur+Qk;
    %Equation 2 added to the running sum.
    Pc=Pc+exp(-Cm*beta*sqrt(det(Sk)));
    
    %Perform the measurement update
    PInv=PInv+H'*inv(RCur)*H;
end
Pc=Pc/K;
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
