function [PPredMissed,PPredObs,gammaVal,PG]=calcMissedGateCov(PPred,Pzz,W,PD,gammaVal,PG)
%%CALCMISSEDGATECOV When tracking performing gating, it is possible that
%            a measurement from the target exists but that the measurement
%            falls outside of the gate. Thus, when performing an update in
%            a filter for the hypothesis that no target-originated
%            measurement is present, one cannot just use the predicted
%            covariance matrix. Rather, it must have a term added to
%            account for the possibility of a true measurement being
%            outside of the gate. This matters more as the gate probability
%            declines and the target detection probability increases. This
%            function is derived in [1] for the linear-Gaussian measurement
%            case. However, it is a reasonable approximation when used with
%            parameters from other filter variants that use Gaussian
%            approximations. This function also returns the predicted
%            covariance matrix conditioned on knowing that the measurement
%            assigned will be from the target at the next time (and thus
%            must gate).
%
%INPUTS: PPred The xDimXxDim predicted target state covariance matrix.
%          Pzz The zDimXzDim innovation covariance matrix.
%            W The xDimXzDim gain from the linear Kalman filter. This also
%              arises in the extended Kalman filter and the cubature Kalman
%              filter.
%           PD The target detection probability.
% gammaVal, PG These two values  contain the same information and only one
%              needs to be provided (the other can be omitted or an empty
%              matrix passed), though passing both will avoid needless
%              recomputation of one of them. gammaVal is the bound used for
%              gating. In such an instance a measurement z gates if
%              (z-zPred)'*Pzz*(z-zPred)<gammaVal where zPred is the
%              predicted measurement and Pzz is the innovation measurement
%              covariance matrix. 0<gammaVal<=Inf and 0<PG<=1.
%
%OUTPUTS: PPredMissed The xDimXxDim covariance matrix of the state
%                 prediction conditioned on there being no detection from
%                 the target in the gate. This is Equation 28 in [1].
%        PPredObs The xDimXxDim covariance matrix of  the state prediction
%                 conditioned on the fact that the measurement that will be
%                 assigned will be from the target (and fall into the
%                 measurement gate). This is Equation 29 in [1]. This term
%                 is less useful in filtering and PPredMissed, because in
%                 Equation 33 in [1], one can see that the update
%                 conditioned on a measurement being target-originated does
%                 not need this term, though, as in Equation 34, one can
%                 use the term.
%    gammaVal, PG These have the same meaning as the input values. If one
%                 omits gammaVal or PG on the input, then one can pass the
%                 returned value on another call to avoid it being
%                 recomputed.
%
%This function implements Equation 28 in [1]. The effects of PG not being 1
%are more significant as PG decreases and PD increases. The function is
%derived assuming that the measurement z is z=H*x+w, where x is the state,
%H is the measurement matrix and w is white Gaussian noise.  
%
%REFERENCES:
%[1] X. R. Li, "Tracking in clutter with strongest neighbor measurements -
%    Part I: Theoretical analysis," IEEE Transactions on Automatic Control,
%    vol. 43, no. 11, pp. 1560-1578, Nov. 1998.
%
%May 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<6||isempty(gammaVal)||isempty(PG))
    %If the user did not provide both gammaVal and PG
    n=size(Pzz,1);

    if(isempty(gammaVal))%PG is given; determine gammaVal.
        if(PG==1)
            gammaVal=Inf;
        else
            gammaVal=ChiSquareD.invCDF(PG,n);
        end
    else%gammaVal is given; determine PG
        if(gammaVal==Inf)
            PG=1;
        else    
            %Equation 18
            PG=ChiSquareD.CDF(gammaVal,n);
        end
    end
end

zDim=size(Pzz,1);

%Equation 26 in [1].
cT=gammainc(gammaVal/2,zDim/2+1)/gammainc(gammaVal/2,zDim/2);

%Equation 28 in [1].
PPredMissed=PPred+(PD*PG*(1-cT)/(1-PD*PG))*(W*Pzz*W');

%Handle possible loss of symmetry due to how Matlab evaluates the above
%equation.
PPredMissed=(PPredMissed+PPredMissed')/2;
if(nargout>1)
    %Equation 29 in [1].
    PPredObs=PPred-(1-cT)*W*Pzz*W';
    
    %Handle possible loss of symmetry due to how Matlab evaluates the above
    %equation.
    PPredObs=(PPredObs+PPredObs')/2;
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
