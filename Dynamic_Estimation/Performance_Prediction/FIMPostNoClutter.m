function J=FIMPostNoClutter(H,F,R,Q,PD)
%FIMPOSTNOCLUTTER Calculate the asymptotic Fisher information matrix after
%                 a measurement update using the information reduction
%                 factor technique for a linear system with a detection
%                 probability less than or equal to 1.
%
%INPUTS: H The zDim X xDim measurement matrix such that H*x+w is the
%          measurement, where x is the state and w is zero-mean Gaussian
%          noise with covariance matrix R.
%        F The xDim X xDim state transition matrix The state at discrete-
%          time k+1 is modeled as F times the state at time k plus zero-
%          mean Gaussian process noise with covariance matrix Q.
%        R The zDim X zDim measurement covariance matrix. This must be
%          positive definite.
%        Q The xDim X xDim process noise covariance matrix. If this is
%          singular, an iterative solution is used. Otherwise, an explicit
%          algorithm is used. If Q is singular, then the Fischer
%          information matrix must be positive definite or numerical issues
%          will arise.
%       PD The optional detection probability of the target at each scan.
%          If omitted, PD is assumed to be one.
%
%OUTPUTS: J The asymptotic posterior Fisher information matrix.
%
%The inverse of the Fisher information matrix for a dynamic system is the 
%posterior Cramér-Rao lower bound (PCRLB). This finds the asymptotic value
%AFTER a measurement update. For the asymptotic value BEFORE  a measurement
%update, see FIMPredNoClutter.
%
%The algorithm is from [1], which simplifies the algorithm of [2] when
%using correlated measurements. The implementation here is for the
%case where there is no clutter. As given in Section 3.2 of [3], the
%information reduction factor for PD<1 without clutter just reduces to PD.
%
%This solves for the Fischer information matrix by solving the Riccati
%equation since in Equation 49 in [4]. It is shown how the basic recursion
%can be rewritten as a Riccati equation that can be solved in a standard
%manner.
%
%If PD=1, the inverse of the FIM  is the same as the output of
%RiccatiPostNoClutter. Otherwise, the inverse FIM is less than the output
%of RiccatiPostNoClutter.
%
%EXAMPLE:
% T=1;
% F=FPolyKal(T,6,1);
% q0=1;
% Q=QPolyKal(T,6,1,q0);
% R=diag([10;10;10]);
% H=[1,0,0,0,0,0;
%    0,1,0,0,0,0;
%    0,0,1,0,0,0];
% PD=0.5;
% P=FIMPostNoClutter(H,F,R,Q,PD)
%
%REFERENCES:
%[1] Y. Bar-Shalom, X. Zhang, and P. Willett, "Simplification of the
%    dynamic Cramér-Rao bound for target tracking in clutter," IEEE
%    Transactions on Aerospace and Electronic Systems, vol. 47, no. 2, pp.
%    1481-1482, Apr. 2011.
%[2] X. Zhang, P. Willett, and Y. Bar-Shalom, "Dynamic Cramér-Rao bound for
%    target tracking in clutter," IEEE Transactions on Aerospace and
%    Electronic Systems, vol. 41, no. 4, pp. 1154-1167, Oct. 2005.
%[3] P. Stinco, M. S. Greco, F. Gini,and A.Farina, "Posterior Cramér-Rao
%    lower bounds for passive bistatic radar tracking with uncertain target
%    measurements," Signal Processing, vol. 93, no. 12, pp. 3528-3540,
%    Dec. 2013.
%[4] P. Tichavský, C. H. Muravchik, and A. Nehorai, "Posterior Cramér-Rao
%    bounds for discrete-time nonlinear filtering," IEEE Transactions on
%    Signal Processing, vol. 46, no. 5, pp. 1386-1396, May 1998.
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<5||isempty(PD))
   PD=1; 
end

n=size(H,2);

%If the Q matrix is singular, then use an iterative approach, otherwise use
%the explicit solution.
if(rcond(Q)<1e-15)
    maxIter=5000;
    
    %The convergence criteria that are used if the Q matrix is singular or
    %nearly singular.
    RelTol=1e-12;
    AbsTol=1e-15;
    
    %Iterate until convergence or until a maximum
    %number of iterations has occurred.
    curIter=0;
    JPrev=eye(size(F));
    Jz=PD*H'*inv(R)*H;
    while(curIter<maxIter)
        J=inv(Q+F*inv(JPrev)*F')+Jz;
        %Ensure symmetry is preserved.
        J=(J+J')/2;
        
        diffMag=abs(J-JPrev);
        if(all((diffMag(:)<=RelTol*abs(J(:)))|(diffMag(:)<=AbsTol)))
            return;
        end
        JPrev=J;
        curIter=curIter+1;
    end
    warning('Max (5000) Iterations Reached Without Convergence.')
else
    QInv=inv(Q);
    D11=F'*QInv*F;
    D12=-F'/Q;
    D22=QInv+PD*(H'*inv(R)*H);

    J=RiccatiSolveD(D11\D12,eye(n,n),D22-D12'*inv(D11)*D12,D11);
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
