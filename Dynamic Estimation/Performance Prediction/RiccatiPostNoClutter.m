function P=RiccatiPostNoClutter(H,F,R,Q,PD,maxIter)
%%RICCATIPOSTNOCLUTTER Solve the posterior Riccati equation that has been
%                      modified for a probability of detection less
%                      than or equal to 1. This provides the average 
%                      asymptotic covariance of a linear Kalman filter with
%                      linear measurements just after a (possible with
%                      probability PD) measurement update.
%
%INPUTS:    H   The zDim X xDim measurement matrix such that H*x+w is the
%               measurement, where x is the state and w is zero-mean 
%               Gaussian noise with covariance matrix R.
%           F   The xDim X xDim state transition matrix The state at
%               discrete-time k+1 is modeled as F times the state at time k
%               plus zero-mean Gaussian process noise with covariance
%               matrix Q.
%           R   The zDim X zDim measurement covariance matrix.
%           Q   The xDim X xDim process noise covariance matrix. This can
%               not be a singular matrix.
%           PD  The optional detection probability of the target at each
%               scan. If omitted, PD is assumed to be one.
%       maxIter An optional integer specifying the maximum number of
%               iterations. By default, if omitted, this is 5000.
%               Convergence is declared if all of the components of the
%               matrix are changed by less than 1 part in 10^12 on an
%               iteration.
%
%OUTPUTS:   P The average asymptotic error covariance matrix of the
%             linear Kalman filter with a probability of detection of PD
%             after a (possible) measurement update.
%
%The notion of a Riccati equation with a detection probability less than
%one is discussed in [1]. However, they do not explicitely write out the
%posterior form. The posterior recursion goes from P_{k|k} to P_{k+1|k+1}.
%This stands in contrast to the prior Riccati equation, which expresses the
%transition from P_{k|k-1} to P_{k+1|k}.
%
%The posterior Riccati equation is
%P=F*P*F'+Q-PD*(F*P*F'*H'+Q*H')*inv(H*F*P*F'*H'+H*Q*H'+R)*(F*P*F'*H'+Q*H')'
%and is derived by substituting the prediction step into the measurement
%update step of the linear Kalman filter with additive process noise and
%additive measurement noise.
%
%REFERENCES:
%[1] Y. Boers and H. Driessen, "Modified Riccati equation and its
%    application to target tracking," IEE Proceedings Radar, Sonar and
%    Navigation, vol. 153, no. 1, pp. 7-12, Feb. 2006.
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<5)
   PD=1; 
end

if(nargin<6)
   maxIter=5000; 
end

epsilon=1e-10;%The convergence criterion.

%Get an initial estimate by solving the Riccati equation with PD=1.
AD=F';
BD=F'*H';
SD=Q*H';
RD=H*Q*H'+R;
QD=Q;
PPrev=RiccatiSolveD(AD,BD,QD,RD,SD);

if(PD~=1)
    %Now, iterate the Riccati equation with PD<1 until convergence of the
    %Frobenis norm or until a maximum number of iterations has occurred.
    curIter=0;
    while(curIter<maxIter)
        P=F*PPrev*F'+Q-PD*(F*PPrev*F'*H'+Q*H')*inv(H*F*PPrev*F'*H'+H*Q*H'+R)*(F*PPrev*F'*H'+Q*H')';
        
        diff=P-PPrev;
        if(sum(sum(abs(diff)./abs(P)>epsilon))==0)
            return;
        end
        
        PPrev=P;
        curIter=curIter+1;
    end
    display('Warning: Max Iterations Reached without Convergence.')
else
    P=PPrev;
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
