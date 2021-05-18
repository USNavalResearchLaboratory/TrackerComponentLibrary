function [zPred,PzPred,otherInfo]=sqrtDDFMeasPred(xPred,SPred,zDim,h,algorithm)
%%SQRTDDFMEASPRED Perform the measurement prediction part of the square-
%           root version of the central difference filter (CDF), or the
%           first or second order divided difference filter (DDF) with
%           additive measurement noise. The function
%           sqrtDDFUpdateWithPred can be used to complete the measurement
%           update. Separating the measurement prediction step from the
%           rest of the update step can make the creation of multiple
%           measurement association hypotheses from a single target
%           prediction more efficient. The full measurement update function
%           is sqrtDDFUpdate.
%
%INPUTS: xPred The xDimXnumComp predicted target states
%        SPred The xDimXxDimXnumComp lower-triangular predicted state
%              covariance matrices associated with the states in xPred.
%         zDim The dimensionality of the measurement.
%            h A function handle for the measurement function that takes
%              the state as its argument.
%    algorithm An optional parameter specifying which algorithm should be
%              used for the prediction. Possible values are:
%              0 Use the CDF of [3]. This is actually the same as a first-
%                order DDF of [1] and [2], but with a different step size
%                for the finite differences.
%              1 (The default if omitted or an empty matrix is passed) Use
%                the first-order divided difference filter of [1] and [2].
%              2 Use the second-order divided difference filter of [1] and
%                [2].
%
%OUTPUTS: zPred The zDimXnumComp measurement predictions from the filter.
%        PzPred The zDimXzDimXnumComp covariance matrices associated with
%               the values in zPred.
%     otherInfo A structure containing members of intermediate results of
%               this function that can be passed to sqrtDDFUpdateWithPred
%               when updating with a measurement.
%
%The filters are discussed in [1] and [2]. See the comments to
%sqrtDDFUpdate for more information.
%
%EXAMPLE:
%With this example, we demonstrate that one gets the same result using
%sqrtDDFUpdate in one step as with using sqrtDDFMeasPred followed by
%sqrtDDFUpdateWithPred.
% xPred=[1e3;-2e3;100;200];
% h=@(x)([sum(x.^2);x(1)-x(4)^(3/2)]);
% SPred=chol([28,   3.5,    6,  8.5;
%            3.5,    23,  8.5,   11;
%              6,   8.5,   18, 13.5;
%            8.5,    11, 13.5,   13],'lower');
% z=1e6*[5.050000548964568;
%       -0.001829553054023];
% zDim=size(z,1);
% SR=eye(zDim,zDim);
% algorithm=1;
% %The update in one step.
% [xUpdate,SUpdate,innov,Szz,W]=sqrtDDFUpdate(xPred,SPred,z,SR,h,algorithm);
% %The update in two steps.
% [zPred,PzPred,otherInfo]=sqrtDDFMeasPred(xPred,SPred,zDim,h,algorithm);
% [xUpdate1,SUpdate1,innov1,Szz1,W1]=sqrtDDFUpdateWithPred(z,SR,zPred,otherInfo);
% %One will see that the one and two step updates agree.
% max(abs([xUpdate1(:)-xUpdate(:);SUpdate1(:)-SUpdate(:);innov1(:)-innov(:);Szz1(:)-Szz(:);W1(:)-W(:)]))
%
%REFERENCES:
%[1] M. Nørgaard, N. K. Poulsen, and O. Ravn, "New developments in state
%    estimation for nonlinear systems," Automatica, vol. 36, no. 11, pp.
%    1627-1638, Nov. 2000.
%[2] T. S. Schei, "A finite-difference method for linearization in
%    nonlinear estimation algorithms," Automatica, vol. 33, no. 11, pp.
%    2053-2058, Nov. 1997.
%
%June 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<5||isempty(algorithm))
    algorithm=1;
end

if(algorithm==0)
    %The central difference filter is the same as the first-order DDF,
    %except deltaH is just 1 not sqrt(3).
    deltaH=1;
else
    %The term derived from the kurtosis of the assumed normalized 
    %distribution, as derived in Section 3.3 of [2].
    deltaH=sqrt(3);
end

xDim=size(xPred,1);
numComp=size(xPred,2);

zPred=zeros(zDim,numComp);
PzPred=zeros(zDim,zDim,numComp);
Pxz=zeros(xDim,zDim,numComp);
Szx=zeros(zDim,xDim,numComp);
if(algorithm==2)
    Szx2=zeros(zDim,xDim,numComp);
else
    Szx2=[];
end

for k=1:numComp
    zPred(:,k)=h(xPred(:,k));
    Sxx=SPred(:,:,k);

    %These are to store values so that they do not have to be repeatedly
    %recomputed.
    gValsPlusX=zeros(zDim,xDim);
    gValsMinusX=zeros(zDim,xDim);
    for curDim=1:xDim
        gValsPlusX(:,curDim)=h(xPred(:,k)+deltaH*Sxx(:,curDim));
        gValsMinusX(:,curDim)=h(xPred(:,k)-deltaH*Sxx(:,curDim));
    end

    %The loops fill in Szx and Szw according to Equation 54 in [1] (equation 93
    %in [2]).
    for curDim=1:xDim
        %The zeros assumes that the mean measurement noise value is zero. If it
        %is not, then the h function could just be modified to reflect that.
        Szx(:,curDim,k)=gValsPlusX(:,curDim)-gValsMinusX(:,curDim);
    end
    Szx(:,:,k)=Szx(:,:,k)/(2*deltaH);

    if(algorithm==2)%Second order DDF
        %The loop fills in Szx2 and Szw2 according to the Equation in
        %Section 4.3 of [1] (Section 4.3 of [2]).
        for curDim=1:xDim
            Szx2(:,curDim,k)=gValsPlusX(:,curDim)+gValsMinusX(:,curDim)-2*zPred(:,k);
        end
        Szx2(:,:,k)=(sqrt(deltaH^2-1)/(2*deltaH^2))*Szx2(:,:,k);
    
        %Equation 71 in [1] (Equation 110 in [2]) without the SR term,
        %squared.
        PzPred(:,:,k)=Szx(:,:,k)*Szx(:,:,k)'+Szx2(:,:,k)*Szx2(:,:,k)';
    else
        %Equation 61 in [1] (Equation 100 in [2]) without the SR term,
        %squared.
        PzPred=Szx(:,:,k)*Szx(:,:,k)';
    end

    %Equation 63 in [1] (Equation 102 in [2]).
    Pxz=SPred(:,:,k)*Szx(:,:,k)';
    %Pxz is not needed for the measurement prediction, but we compute it
    %here, so that it need not be recomputed again and again if
    %sqrtDDFUpdateWithPred is called for multiple measurements.
end

otherInfo.algorithm=algorithm;
otherInfo.SPred=SPred;
otherInfo.Pxz=Pxz;
otherInfo.Szx=Szx;
otherInfo.Szx2=Szx2;
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
