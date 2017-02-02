function J=FIMPredNoClutter(H,F,R,Q,PD)
%FIMPREDNOCLUTTER Calculate the asymptotic Fisher information matrix before
%                 a measurement update using the information reduction
%                 factor technique for a linear system with a detection
%                 probability less than or equal to 1.
%
%INPUTS:    H   The zDim X xDim measurement matrix such that H*x+w is the
%               measurement, where x is the state and w is zero-mean 
%               Gaussian noise with covariance matrix R.
%           F   The xDim X xDim state transition matrix The state at
%               discrete-time k+1 is modeled as F times the state at time k
%               plus zero-mean Gaussian process noise with covariance
%               matrix Q.
%           R   The zDim X zDim measurement covariance matrix.
%           Q   The xDim X xDim process noise covariance matrix. If this is
%               singular, an iterative solution is used. Otherwise, an
%               explicit algorithm is used. If Q is singular, then the
%               Fischer information matrix must be positive definite or
%               numerical issues will arise.
%           PD  The optional detection probability of the target at each
%               scan. If omitted, PD is assumed to be one.
%
%OUTPUTS:   J   The asymptotic prior Fisher information matrix.
%
%The inverse of the Fisher information matrix for a dynamic system is the 
%posterior Cramér-Rao lower bound (PCRLB). This finds the asymptotic value
%BEFORE a measurement update by predicting forward the posterior FIM from
%FIMPostNoClutter.
%
%This algorithm just takes the posterior FIM and puts it through the
%a recursive step for the FIM without a measurement update.
%For more information on how the function works, see the function 
%FIMPostNoClutter
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

JPost=FIMPostNoClutter(H,F,R,Q,PD);

if(rcond(Q)<1e-15)
   J=inv(Q+F*inv(JPost)*F');
else
    QInv=inv(Q);
    D11=F'*QInv*F;
    D12=-F'/Q;
    J=QInv-D12'*inv(JPost+D11)*D12;
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





