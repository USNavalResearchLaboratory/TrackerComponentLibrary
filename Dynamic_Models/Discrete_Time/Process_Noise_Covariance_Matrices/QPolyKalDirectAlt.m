function Q=QPolyKalDirectAlt(T,xDim,order,sigmaV2)
%%QPOLYKALDIRECTALT Get the process noise covariance matrix for a direct
%           discrete linear dynamic model of the given polynomial
%           order (number of derivatives of position included) and number
%           of dimensions (generally 3 for 3D motion). order=1 means
%           a nearly constant velocity model. order=2 means the
%           discrete Wiener process acceleration (DWPA) model, etc. The
%           state is ordered in terms of position, velocity, acceleration,
%           etc. Unlike the discretized continuous-time model in the
%           function QPolyKal, the direct-discrete implementation produces
%           a singular matrix and inconsistencies arise when predicting
%           over a time period of 2*T versus sequentially predicting over
%           two time perods of length T. In this function, the process
%           noise is modeled added to the same moment as the highest order
%           of the state. This is contrasted with the function
%           QPolyKalDirectDisc where the process noise is one moment higher
%           than the highest order of the state.
%
%INPUTS: T The time-duration of the propagation interval.
%     xDim The dimensionality of the target state.
%          xDim is ((order+1)*numDim)X1 dimensional, where numDim is the
%          number of position dimensions of space (e.g. 2D or 3D).
%    order The order >=0 of the filter. If order=1, then it is constant
%          velocity, 2 means constant acceleration, 3 means constant jerk,
%          etc.
%  sigmaV2 The variance driving the process noise. This has units of
%          distance^2/time^(2*order). sqrt(sigmaV2) is proportional to the
%          standard deviation of the noise affecting the highest-order
%          moment over a sampling period of length T. If a scalar is
%          passed, it is assumed to be the same for all dimensions.
%          Otherwise, a numDimX1 or 1XnumDim vector should be passed
%          specifying the value for each dimension.
%
%OUTPUTS: Q The process noise covariance matrix under the direct discrete
%           linear dynamic model of the given order with motion in numDim
%           dimensions where the state is stacked
%           [position;velocity;acceleration;etc] where the number of
%           derivatives of position depends on the order given. order=0
%           means just position.
%
%Chapter 1.5.6 of [1] presents the discrete Wiener process acceleration for
%order 2. The logic behind model is also explained in Chapter 6.3.3 of [2].
%Extending from how the standard direct discrete model of order 1 was made.
%The generalization to an arbitrary order has the matrix for 1D motion
%being sigmaV2*G*G' where the ith element in the column vector G is
%T^(order-i+1)/(order-i+1)! where i counts from 1 to order+1.
%
%Note that order=2 produces the discrete Wiener process acceleration model
%not the discrete white noise jerk model.
%
%This process noise matrix is most commonly used with the FPolyKal state
%transition matrix.
%
%REFERENCES:
%[1] Y. Bar-Shalom, P. K. Willett, and X. Tian, Tracking and Data Fusion.
%    Storrs, CT: YBS Publishing, 2011.
%[2] Y. Bar-Shalom, X. R. Li, and T. Kirubarajan, Estimation with
%    Applications to Tracking and Navigation. New York: John Wiley and
%    Sons, Inc, 2001.
%
%April 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numDim=xDim/(order+1);

if(isscalar(sigmaV2))
    sigmaV2=ones(numDim,1)*sigmaV2;
end

numEl=order+1;
i=(1:numEl)';
G=T.^(order-i+1)./factorial(order-i+1);
%These are the elements of all of the 1D submatrices, except they have to
%be multiplied by the appropriate sigmaV.
QBase=G*G';

%Now, the elements just get spread across identity matrices that are
%numDim dimensional to form a process noise covariance matrix of the
%desired dimensionality. This is done using a Kronecker product.
Q=kron(QBase,diag(sigmaV2));

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
