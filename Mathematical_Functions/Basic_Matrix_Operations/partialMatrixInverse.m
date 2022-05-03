function PA=partialMatrixInverse(A,n)
%%PARTIALMATRIXINVERSE Given an mXm possibly singular matrix A of a rank
%   >=n, where the first n columns are independent, obtain a matrix that is
%   akin to the inverse related to the first n columns. A normal matrix
%   inverse is such that A*inv(A)=eye(m,m). Here, consider A=Q*R, where Q
%   is an orthogonal matrix and R is a lower-triangular matrix (not the
%   standard QR decomposition). The inverse relation is Q*R*inv(A)=eye(m,m)
%   or R*inv(A)=Q'*eye(m,m). This function solves for the nXn upper-left
%   submatrix of inv(A).  This function is useful for taking a singular
%   Fisher informationmatrix, where the first n rows are observable, but
%   the rest are not, and obtaining a Cramer-Rao lower bound for the
%   observable components. Typically, this might be obtaining a CRLB on
%   position components before velocity is fully observable.
% 
%INPUTS: A An mXm matrix.
%        n A number n<=m that is <= the rank of A.
%
%OUTPUTS: PA The partial inverse matrix as described above. 
%
%This function solves the problem by calling linEqSolveFirstnRows with A
%and vectors [1;0;0;0;...], [0;1;0;0;...], etc. The solution coincides with
%a submatrix of the pseudoinverse.
%
%EXAMPLE:
%This example of the CRLB of a position estimate of a target given a single
%range, direction cosine and range-rate measurement is considered. We show
%that when this function is used, the position accuracy implied by the CRLB
%is not affected by the accuracy of the range rate component (which one
%would assume). However, using the ad-hoc method of taking the inverse of a
%subset of the Fisher information matrix, one gets a different (wrong)
%answer the more accurate the range rate component is. That is due to not
%properly accounting for matrix cross terms.
% %Sensor location near Maui.
% llhRxTx=[deg2rad([20.888645;-156.022474]);0];
% %Convert to Cartesian.
% xRxTx=zeros(6,1);
% %Fill in the position. The reciever will be stationary (0 velocities).
% %Take the first one to be the transmitter.
% xRxTx(1:3,:)=ellips2Cart(llhRxTx);
% 
% %Sensor orientation.
% elAboveLevel=deg2rad(15);%The radar point slightly up.
% %Rotation matrices of the radar.
% MRxTx=findRFTransParam(llhRxTx(:,1),pi/2,elAboveLevel);%Facing East
% 
% %Construct a target state vector. 
% llhTar=[deg2rad([21.14146;-155.058839]);13e3];
% tarLocCart=ellips2Cart(llhTar);
% tarHeading=deg2rad(-140);%Radians East of North.
% angUpFromLevel=0;%Level flight.
% tarSpeed=300;%m/s
% vECEFTar=geogHeading2uVec(llhTar,tarHeading,angUpFromLevel)*tarSpeed;
% xTar=[tarLocCart;vECEFTar];
% 
% %Measuremnet and noise parameters of the radar.
% useHalfRange=false;
% sigmaR=10;%Range standard deviation.
% sigmaU=0.01;
% sigmaV=0.01;
% sigmaRR=1;%Range rate standard deviation, meters/ second.
% sigmaRRMassive=1e5;%Effectively uninformative standard deviation.
% %Lower-triangular square root measurement covariance matrix (No cross
% %terms).
% SR=diag([sigmaR;sigmaU;sigmaV;sigmaRR]);
% RInv=inv(SR*SR');
% %The same thing but with the range rate standard deviation so massive it is
% %essentially uninformative.
% SRMassive=diag([sigmaR;sigmaU;sigmaV;sigmaRRMassive]);
% RInvMassive=inv(SRMassive*SRMassive');
% 
% %Gradients for the Fisher information matrices.
% H=calcRuvRRJacob(xTar,useHalfRange,xRxTx(:,1),xRxTx(:,1),MRxTx(:,:,1));
% FIM=observedFisherInfo([],RInv,[],H);
% %The position RMSE from the CRLB found correctly using this function.
% RMSE=sqrt(sum(diag(partialMatrixInverse(FIM,3))))
% %An ad-hoc (wrong) approach that one might be tempted to take:
% RMSEWrong=sqrt(sum(diag(inv(FIM(1:3,1:3)))))
% %One can see that the position RMSE is higher when done the correct way.
% %However, a single range rate on its own doesn't inform on the position
% %accuracy. If we do the same thing cranking up the standard deviation of
% %the range rate so it is uninformative, we get:
% FIM=observedFisherInfo([],RInvMassive,[],H);
% RMSEUninfRR=sqrt(sum(diag(partialMatrixInverse(FIM,3))))
% RMSEUninfRRWrong=sqrt(sum(diag(inv(FIM(1:3,1:3)))))
% %In other words, the correct solution doesn't change, but in this case the
% %wrong solution has become correct, because the cross terms that messed it
% %up are essentially gone. 
%
%February 2022 David F.Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

m=size(A,1);
PA=zeros(n,n);

if(n>m)
    error('n cannot be > the length of A.')
end

e=zeros(m,1);
for k=1:n
    e(k)=1;
    PA(:,k)=linEqSolveFirstnRows(A,e,n);
    e(k)=0;
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
