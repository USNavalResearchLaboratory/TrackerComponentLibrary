function xProp=sampleDiscProp(xPrev,F,SQ)
%%SAMPLEDISCPROP Assuming a linear dynamic model of
%                x_{k+1}=F*x_k+v
%                where F is a matrix and v is zero-mean Gaussian noise with
%                covariance matrix Q=SQ*SQ', or a nonlinear dynamic model
%                of
%                x_{k+1}=F(x_k)+v
%                where F is a function, this function simulate a random
%                step to propagate the trajectory.
%
%INPUTS: xPrev The xDimXN deterministic set of N target states at the
%              previous time to propagate.
%            F This is either the xDimXxDim state transition matrix or a
%              function handle to a state transition function that takes x
%              as an argument.
%           SQ The lower-triangular square root of the process noise
%              covariance matrix.
%
%OUTPUTS: xProp The xDimXN set of target states in xPrev propagated one
%               step.
%
%This function just draws a random sample fo the normal distribution to
%simulate the dynamic model.
%
%March 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

xDim=size(xPrev,1);
N=size(xPrev,2);

if(isa(F,'function_handle'))
    xProp=zeros(xDim,N);
    for k=1:N
        xProp(:,k)=F(xPrev(:,k));
    end

    xProp=xProp+randn(xDim,N);
else
    xProp=F*xPrev+SQ*randn(xDim,N);
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
