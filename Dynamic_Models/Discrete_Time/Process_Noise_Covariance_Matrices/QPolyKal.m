function Q=QPolyKal(T,xDim,order,q0)
%%QPOLYKAL Get the process noise covariance matrix for a discretized
%          continuous-time linear dynamic model of the given polynomial
%          order (number of derivatives of position included) and number of
%          dimensions (generally 3 for 3D motion). order=1 means constant
%          velocity(the discretized continuous white noise acceleration
%          [DCWNA] model); order=2 means constant acceleration (the
%          discretized continuous Wiener process acceleration model
%          [DCWPA]), etc. The state is ordered in terms of position,
%          velocity, acceleration, etc. The equivalent continuous-time
%          process noise model just adds the noise to the highest order
%          derivative of position.
%
%INPUTS: T The time-duration of the propagation interval.
%     xDim The dimensionality of the target state.
%          xDim is ((order+1)*numDim)X1 dimensional, where numDim is the
%          number of position dimensions of space (e.g. 2D or 3D).
%    order The order >=0 of the filter. If order=1, then it is constant
%          velocity, 2 means constant acceleration, 3 means constant jerk,
%          etc.
%       q0 The power spectral density of the noise. If a scalar is passed,
%          it is assumed to be the same for all dimensions. Otherwise, a
%          numDimX1 or 1XnumDim vector should be passed specifying the
%          power spectral density for each dimension (dimensions of space;
%          not xDim. Unless order=0, xDim>numDim). The units are
%          length^2/time^(2*order+1).
%
%OUTPUTS: Q The process noise covariance matrix under a linear dynamic
%           model of the given order with motion in numDim dimensions
%           where the state is stacked [position;velocity;acceleration;etc]
%           where the number of derivatives of position depends on the
%           order given. Order=0 means just position.
%
%Chapter 6.2.2 of [1] describes how the covariance matrix for a linear
%dynamic model is found. This function just generalized the procedure to an
%arbitrary polynomial order. Equations for the generalization are given in
%Appendix A of [2].
%
%The state for the 1D case is assumed to be
%[position;velocity;acceleration;etc] --so in the order of increasing
%derivatives. In the multidimensional case, the same order is preserved, so
%position becomes numDim-dimensional as do velocity, acceleration, etc.
%This means that the values for the 1D case just get repeated.
%
%This process noise matrix is most commonly used with the FPolyKal state
%transition matrix.
%
%REFERENCES:
%[1] Y. Bar-Shalom, X. R. Li, and T. Kirubarajan, Estimation with
%    Applications to Tracking and Navigation. New York: John Wiley and
%    Sons, Inc, 2001.
%[2] D. F. Crouse, "Cubature/ unscented/ sigma point Kalman filtering with
%    angular measurement models," in Proceedings of the 18th International
%    Conference on Information Fusion, Washington, D.C., 6-9 Jul. 2015.
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(~isscalar(xDim)&&isvector(xDim))
    %Backwards compatibility for old version.
    xDim = length(xDim);
end
numDim=xDim/(order+1);

if(isscalar(q0))
    q0=ones(numDim,1)*q0;
end

if(length(q0)~=numDim)
    error('The size of q0 disagrees with the number of dimensions impliced by order and xDim.')
end

%First, create the matrices for 1D motion
numEl=order+1;
sel=(numEl-1):-1:0;
[colIdx,rowIdx]=meshgrid(sel,sel);

Q1=T.^(colIdx+rowIdx+1)./(factorial(colIdx).*factorial(rowIdx).*(colIdx+rowIdx+1));

%Now, the elements just get spread across identity matrices that are
%numDim dimensional to form a process noise covariance matrix of the
%desired dimensionality. This is done using a Kronecker product.
Q=kron(Q1,diag(q0));
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
