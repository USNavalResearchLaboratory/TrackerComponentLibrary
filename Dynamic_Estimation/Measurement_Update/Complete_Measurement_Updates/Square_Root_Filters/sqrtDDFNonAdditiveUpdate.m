function [xUpdate,SUpdate,innov,Szz,W]=sqrtDDFNonAdditiveUpdate(xPred,SPred,z,SR,h,algorithm,innovTrans)
%%SQRTDDFNONADDITIVEUPDATE Perform the measurement update step in the
%               square-root version of the central difference filter (CDF),
%               or the first or second order divided difference filter
%               (DDF) with non-additive measurement noise.
%
%INPUTS: xPred The xDim X 1 predicted target state.
%        SPred The xDim X xDim lower-triangular predicted state covariance
%              matrix.
%            z The zDim X 1 vector measurement.
%           SR The zDim X zDim lower-triangular square root of the
%              measurement covariance matrix in the native coordinate
%              system of the measurement.
%            h A function handle for the measurement function. This can
%              either be of the form h(xStacked), where xStacked is [x;w],
%              where x is the state and w is the measurement noise value,
%              or this can be of the form h(x,w), where the inputs are 
%              separately specified.
%    algorithm An optional parameter specifying which algorithm should be
%              used for the prediction. Possible values are:
%              0 Use the CDF of [3]. This is actually the same as a first-
%                order DDF of [1] and [2], but with a different step size
%                for the finite differences.
%              1 (The default if omitted or an empty matrix is passed) Use
%                the first-order divided difference filter of [1] and [2].
%              2 Use the second-order divided difference filter of [1] and
%                [2].
%   innovTrans An optional function handle that computes and optionally
%              transforms the value of the difference between the
%              observation and any predicted points. This is called as
%              innovTrans(a,b) and the default if omitted or an empty
%              matrix is passed is @(a,b)bsxfun(@minus,a,b). This must be
%              able to handle sets of values. For a zDimX1 measurement,
%              either of the inputs could be zDimXN in size while one of
%              the inputs could be zDimX1 in size.  This only needs to be
%              supplied when a measurement difference must be restricted
%              to a certain range. For example, the innovation between two
%              angles will be 2*pi if one angle is zero and the other
%              2*pi, even though they are the same direction. In such an
%              instance, a function handle to the
%              wrapRange(bsxfun(@minus,a,b),-pi,pi) function with the
%              appropriate parameters should be passed for innovTrans.
%
%OUTPUTS: xUpdate The xDim X 1 updated state vector.
%         SUpdate The updated xDim X xDim lower-triangular square root
%                 state covariance matrix.
%      innov, Szz The zDimX1 innovation and the zDimXzDim square root
%                 innovation covariance matrix are returned in case one
%                 wishes to analyze the consistency of the estimator or use
%                 those values in gating or likelihood evaluation.
%               W The gain used in the update. This can be useful when
%                 gating and using the function calcMissedGateCov.
%
%The divided difference filters are implemented as described in [1] and
%[2], where they are only derived for non-additive noise and the square
%root form is the most natural form to use. The central difference filter
%of [3] is the same as a first-order divided difference filter, but using a
%different distance for the finite differences.
%
%The optional parameter innovTrans is not described in references [1]-[3]
%above reference, but allow for possible modifications to the filter as
%described in [4] The parameter has been added to allow the filter to be
%used with angular quantities. For example, if the measurement consisted of
%range and angle, z=[r;theta], then
%innovTrans=@(a,b)[bsxfun(@minus,a(1,:),b(1,:));
%                  wrapRange(bsxfun(@minus,a(2,:),b(2,:)),-pi,pi)];
%should be used to approximately deal with the circular nature of the
%measurements.
%
%REFERENCES:
%[1] M. Nørgaard, N. K. Poulsen, and O. Ravn, "New developments in state
%    estimation for nonlinear systems," Automatica, vol. 36, no. 11, pp.
%    1627-1638, Nov. 2000.
%[2] M. Nørgaard, N. K. Poulsen, and O. Ravn, "Advances in derivative-free
%    state estimation for nonlinear systems (revised edition)," Department
%    of Informatics and Mathematical Modelling, Technical University of
%    Denmark, Tech. Rep. IMM-REP-1998-15, 29 Oct. 2004. [Online].
%    Available: http://www2.imm.dtu.dk/pubdb/views/publication details.php?id=2706
%[3] T. S. Schei, "A finite-difference method for linearization in
%    nonlinear estimation algorithms," Automatica, vol. 33, no. 11, pp.
%    2053-2058, Nov. 1997.
%[4] D. F. Crouse, "Cubature/ unscented/ sigma point Kalman filtering with
%    angular measurement models," in Proceedings of the 18th International
%    Conference on Information Fusion, Washington, D.C., 6-9 Jul. 2015.
%
%October 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<6||isempty(algorithm))
    algorithm=1;
end

if(nargin<7||isempty(innovTrans))
    %The function just returns the input.
    innovTrans=@(a,b)bsxfun(@minus,a,b);
end

%This lets the user provide two different parameterization types for the
%measurement function.
if(nargin(h)<2)
    hSplit=@(x,v)h([x;v]);
else
    hSplit=h;
end

zDim=size(z,1);
xDim=size(xPred,1);

g0=hSplit(xPred,zeros(xDim,1));

if(algorithm==0)
    %The central difference filter is the same as the first-order DDF,
    %except deltaH is just 1 not sqrt(3).
    deltaH=1;
else
    %The term derived from the kurtosis of the assumed normalized 
    %distribution, as derived in Section 3.3 of [2].
    deltaH=sqrt(3);
end

Sxx=SPred;

%These are to store values so that they do not have to be repeatedly
%recomputed.
gValsPlusX=zeros(zDim,xDim);
gValsMinusX=zeros(zDim,xDim);
gValsPlusV=zeros(zDim,zDim);
gValsMinusV=zeros(zDim,zDim);
for curDim=1:xDim
    gValsPlusX(:,curDim)=hSplit(xPred+deltaH*Sxx(:,curDim),zeros(xDim,1));
    gValsMinusX(:,curDim)=hSplit(xPred-deltaH*Sxx(:,curDim),zeros(xDim,1));
end

for curDim=1:zDim
    gValsPlusV(:,curDim)=hSplit(xPred,deltaH*SR(:,curDim));
    gValsMinusV(:,curDim)=hSplit(xPred,-deltaH*SR(:,curDim));
end

%The loops fill in Szx and Szw according to Equation 54 in [1] (equation 93
%in [2]).
Szx=zeros(zDim,xDim);
Szw=zeros(zDim,zDim);
for curDim=1:xDim
    %The zeros assumes that the mean measurement noise value is zero. If it
    %is not, then the h function could just be modified to reflect that.
    Szx(:,curDim)=gValsPlusX(:,curDim)-gValsMinusX(:,curDim);
end
Szx=Szx/(2*deltaH);

for curDim=1:zDim
    Szw(:,curDim)=gValsPlusV(:,curDim)-gValsMinusV(:,curDim);
end
Szw=Szw/(2*deltaH);

if(algorithm==2)%Second order DDF
    %The loop fills in Szx2 and Szw2 according to the Equation in Section
    %4.3 of [1] (Section 4.3 of [2]).
    Szx2=zeros(zDim,xDim);
    for curDim=1:xDim
        Szx2(:,curDim)=gValsPlusX(:,curDim)+gValsMinusX(:,curDim)-2*g0;
    end
    
    Szw2=zeros(zDim,zDim);
    for curDim=1:zDim
        Szw2(:,curDim)=gValsPlusV(:,curDim)+gValsMinusV(:,curDim)-2*g0;
    end
    
    Szx2=(sqrt(deltaH^2-1)/(2*deltaH^2))*Szx2;
    Szw2=(sqrt(deltaH^2-1)/(2*deltaH^2))*Szw2;
    
    %Equation 71 in [1] (Equation 110 in [2])
    Szz=tria([Szx,Szw,Szx2,Szw2]);
else
    %Equation 61 in [1] (Equation 100 in [2])
    Szz=tria([Szx,Szw]);
end

%Equation 63 in [1] (Equation 102 in [2]).
Pxz=SPred*Szx';

%Equation 64 in [1] (Equation 103 in [2])
W=(Pxz/Szz')/Szz;

%The innovation
innov=innovTrans(z,g0);

%Equation 65 in [1] (Equation 104 in [2]).
xUpdate=xPred+W*innov;

if(algorithm==2)%Second order DDF
    %Equation 115 in [2].
    SUpdate=tria([SPred-W*Szx,W*Szw,W*Szx2,W*Szw2]);
else
    %Equation 106 in [2]
    SUpdate=tria([SPred-W*Szx,W*Szw]);
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
