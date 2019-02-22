function [zPred,PzPred,otherInfo]=separatedCovMeasPred(xPred,TPred,H)
%%SEPARATEDCOVMEASPRED Perform the measurement prediction part of the
%           measurement update step of the separated covariance filter. The
%           function separatedCovUpdateWithPred can be used to complete the
%           measurement update. Separating the measurement prediction step
%           from the rest of the update step can make the creation of
%           multiple measurement association hypotheses from a single
%           target prediction more efficient. The full measurement update
%           function is separatedCovUpdate.
%
%INPUTS: xPred The xDimX1 predicted target state.
%        TPred The xDimXxDim total error matrix. This is a combination of
%              errors due to measurement noise and filter lag.
%            H The zDim X xDim measurement matrix for a linear measurement
%              model. That is z=H*x+w, where w is measurement noise having
%              covariance matrix R.
%
%OUTPUTS: zPred The zDimX1 measurement prediction from the filter.
%        PzPred The zDimXzDim covariance matrix associated with zPred.
%     otherInfo A structure containing members of intermediate results of
%               this function that can be passed to
%               separatedCovUpdateWithPred when updating with a
%               measurement.
%
%The filter is described in [1]. See the comments to  separatedCovUpdate
%for more information on the algorithm.
%
%EXAMPLE:
%With this example, we demonstrate that one gets the same result using
%separatedCovUpdate versus calling separatedCovMeasPred and then
%separatedCovUpdateWithPred. 
% xPred=[1e3;-2e3;100;200];
% TPred=[28,   3.5,    6,  8.5;
%      3.5,    23,  8.5,   11;
%        6,   8.5,   18, 13.5;
%      8.5,    11, 13.5,   13];
% LPred=[8, 0;
%       0, 8;
%       8, 0;
%       0, 8];
% z=1e3*[-5.498856156296510;
%        1.199241491470584];
% R=eye(2);
% H=[0, 4, 9, 8;
%    6, 3, 0, 6];
% c=0.99;
% %The update in one step.
% [xUpdate,LUpdate,PUpdate,innov,Pzz,W]=separatedCovUpdate(xPred,LPred,TPred,z,R,H,c);
% %The update in two steps.
% [zPred,PzPred,otherInfo]=separatedCovMeasPred(xPred,TPred,H);
% [xUpdate1,LUpdate1,PUpdate1,innov1,Pzz1,W1]=separatedCovUpdateWithPred(z,R,zPred,PzPred,otherInfo,LPred,c);
% %One will see that the one and two step updates agree.
% max(abs([xUpdate1-xUpdate;LUpdate1(:)-LUpdate(:);PUpdate1(:)-PUpdate(:);innov1(:)-innov;Pzz1(:)-Pzz(:);W1(:)-W(:)]))
%
%REFERENCES:
%[1] G. J. Portmann, J. R. Moore, and W. G. Bath, "Separated covariance
%    filtering," in Proceedings of the IEEE International Radar Conference,
%    Arlington, VA, 7-10 May 1990, pp. 456-460.
%
%June 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

PzPred=H*TPred*H';
%Ensure symmetry
PzPred=(PzPred+PzPred')/2;
zPred=H*xPred;

Pxz=TPred*H';
%Pxz is not needed for the measurement prediction, but we compute it
%here, so that it need not be recomputed again and again if
%separatedCovUpdateWithPred is called for multiple measurements.

otherInfo.H=H;
otherInfo.Pxz=Pxz;
otherInfo.TPred=TPred;
otherInfo.xPred=xPred;
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
