function [xi,w,CvMDistMin,exitCode,sInit]=GaussianLCDSamples(numDim,numSamples,forceCovMatch,LBFGSOptions,numIntOptions,bMax,sInit)
%%GAUSSIANLCDSAMPLES Generate cubature points for a multidimensional N(0,I)
%              Gaussian distribution using an optimization over a
%              Modified Cramer-von Mises (CvM) distance between localized
%              cumulative distributions (LCDs). As opposed to having a
%              fixed polynomial order, the cubature points are generated as
%              a best fit (based on the CvM distance) for a fixed number of
%              points. Thus, when used for cubature integration, the points
%              might not be exact for any order, but should be good over a
%              wide range of orders. The approach, as applied to track
%              filtering, is described in [2]. The symmetric algorithm of
%              [1] is used here. If an initial point pattern (sInit) is not
%              given, then running this function multiple times with the
%              same inputs will give different results.
%
%INPUTS: numDim An integer specifying the dimensionality of the cubature
%               points to be generated.
%    numSamples The number of points to use to approximate a normal zero
%               mean Gaussian distribution with the identity matrix as its
%               covariance. This must be >=2*numDim for the covariance
%               matrix of the resulting samples to be non-singular.
% forceCovMatch An optional argument specifying that the covariance matrix
%               should be forced to match that of a normal (0,I)
%               distribution. The basic algorithm best fits the
%               distribution according to certain measures that usually
%               underestimate the diagonal elements of the covariance
%               matrix. The default if this is omitted is true, unless
%               numSamples<2*numDim, in which case it is false, because
%               such a match cannot be achieved with symmetric samples.
%  LBFGSOptions A structure containing parameters for the limited-memory
%               Broyden-Fletcher-Goldfarb-Shanno algorithm, which is used
%               as a subroutine for optimization. Each member of the
%               structure is named after an option in the input to the
%               quasiNewtonLBFGS function. Possible members are numCorr,
%               epsilon, deltaTestDist, delta, lineSearchParams,
%               maxIterations, progressCallback. Note that the C and
%               l1NormRange members within lineSearchParams are not used.
%               The meanings of the parameters iaredescribed in the
%               comments to the quasiNewtonLBFGS function. If this
%               structure is omitted or an empty matrix is passed, then
%               the default values will be used. If any of the members of
%               the structure are omitted or empty matrices are passed,
%               then the default values will be used.
% numIntOptions A structure whose parameters specify the accuracy of the
%               scalar numerical integration used in this function. The
%               possible members are:
%               AbsTol The absolute error tolerance. Between a value x and
%                      its estimate xHat, this is abs(x-xHat).
%               RelTol The relative error tolerance. Between a value x and
%                      its esitimate xHat, this is
%                      abs(x-xHat)/min(abs(x),abs(xHat))
%               The default values are both 1e-14. If this parameter is
%               omitted or an empty matrix is passed, the default values
%               are used. If either member of the structure is omitted or
%               is assigned an empty matrix, the default value is used.
%          bMax The upper bound of the integrals, defined as in [1].
%               Higher numbers lead to more precise results; numbers too
%               high will run into finite precision problems with respect
%               to the integral. For numDim<=1000, the authors in [1],
%               suggest using 70. If omitted or an empty matrix is passed,
%               the default value of 70 is used.
%         sInit An optional numDim X fix(numSamples/2) set of points that
%               is used as the seed for the algorithm. If omitted, a
%               random seed is used.
%
%OUTPUTS: xi The numDim X numSamples set of cubature points for the normal
%            0-I distribution. If the exit code is not positive, then an
%            empty matrix is returned.
%          w The numDimX1 set of weights associated with the cubature
%            points. These are all just 1/numSamples. If the exit code is
%            not positive, then an empty matrix is returned.
% CvMDistMin The minimum Cramér-von Mises distance between the LCD of the
%            points and a standard multivariate gaussian distribution.
%   exitCode A number indicating the termination state of the algorithm.
%            A zero or positive number indicates success. Possible values
%            are the same as those of the quasiNewtonLBFGS function (see
%            comments to that function), with one addition:
%            -2 The points found  formed a singular covariance matrix and
%            forceCovMatch=true. If the sample covariance matrix is
%            singular, then it is not possible to force it to match the
%            true covariance matrix when using symmetric samples.
%      sInit The random seed used to seed the algorithm.
%
%Comments in the code relating pages, sections and theorems refer to those
%in [1]. The definition of the exponential integral used in [1] is the
%Cauchy principle value integral. However, in Appendix II of [3], one can
%see that the defition omits the case where the argument is zero, but this
%case arises a lot, due to the double sums. The fix added for this instance
%comes by observing the derivation of De3 in Appendices B and C, where it
%can be seen that is s_i=s_j then there is a 0*Inf term, which should just
%equal zero. Thus, the correction here is to "redefine" the exponential
%integral to return 0 when its argument is zero. that avoids the
%introduction of NaN terms to the problem.
%
%Example:
% sInit =[-0.2962    0.3036   -0.2738   -0.0276    1.9153    0.7165;
%          0.5643   -0.7903    0.2719    0.9239    0.1568    1.3373;
%          1.5826    0.8034    1.4896   -0.3213   -0.3005    2.1257;
%          2.7292   -1.3199    1.4371    0.6611   -0.5000    0.0540];
% [xi,w,CvMDistMin,exitCode,sInit]=GaussianLCDSamples(4,13,true,[],[],[],sInit)
%Using sInit as the random seed for the initialization, one obtains points
%with a minimum Cramér-von Mises distance of CvMDistMin=0.0183. The first
%point is at the origin, because an odd number of points was requested.
%
%The function is generally too slow to use to find approximations of
%Gaussian distributions in most real-time systems. However, one can
%precompute points of the desired dimensionalities and create a library
%that can then be called in lieue of using this function.
%
%REFERENCES:
%[1] J. Steinbring, M. Pander, and U. D. Hanebeck, "The smart sampling
%    Kalman filter with symmetric samples," arXiv, 10 Jun. 2015. [Online].
%    Available: http://arxiv.org/abs/1506.03254
%[2] J. Steinbring and U. D. Hanebeck, "LRKF revisited: The smart
%    sampling Kalman filter (S2KF)," Journal of Advances in Information
%    Fusion, vol. 9, no. 2, pp. 106-123, Dec. 2014.
%
%August 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%Symmetric sampling is being used.
numHalfSamples=fix(numSamples/2);
isEven=mod(numSamples,2)==0;
sDims=[numDim,numHalfSamples];

if(nargin<7||isempty(sInit))
    %If an initial set of samples is not provided.
    sInit=randn(sDims);
end

if(nargin<6||isempty(bMax))
    bMax=70;
end

%Default tolerances for numerical integration.
AbsTol=1e-14;
RelTol=1e-14;

if(nargin>=5&&~isempty(numIntOptions))
    if(isfield(numIntOptions,'AbsTol'))
        AbsTol=numIntOptions.AbsTol;
    end
    
    if(isfield(numIntOptions,'RelTol'))
        RelTol=numIntOptions.RelTol;
    end
end
    
%The default parameters for the LBGFS function: use empty matrices to
%indicate the defaults.
numCorr=[];
epsilon=[];
deltaTestDist=[];
delta=[];
lineSearchParams=[];
maxIterations=[];

if(nargin>=4&&~isempty(LBFGSOptions))
    if(isfield(LBFGSOptions,'numCorr'))
        numCorr=LBFGSOptions.numCorr;
    end
    if(isfield(LBFGSOptions,'epsilon'))
        epsilon=LBFGSOptions.epsilon;
    end
    if(isfield(LBFGSOptions,'deltaTestDist'))
        deltaTestDist=LBFGSOptions.deltaTestDist;
    end
    if(isfield(LBFGSOptions,'delta'))
        delta=LBFGSOptions.delta;
    end
    if(isfield(LBFGSOptions,'lineSearchParams'))
        lineSearchParams=LBFGSOptions.lineSearchParams;
        
        if(isfield(lineSearchParams,'C')&&~isempty(lineSearchParams.C)&&lineSearchParams.C~=0)
           error('The C option in the line search parameters cannot be nonzero.')
        end
    end
    if(isfield(LBFGSOptions,'maxIterations'))
        maxIterations=LBFGSOptions.maxIterations;
    end
end

if(nargin<3||isempty(forceCovMatch))
   if(numSamples>=2*numDim)
       forceCovMatch=true;
   else
       forceCovMatch=false;
   end
end

%Create a handle for the cost function and its gradient.
f=@(s)deal(modCvMDist(s,bMax,isEven,AbsTol,RelTol,sDims),modCvMDistGrad(s,bMax,isEven,AbsTol,RelTol,sDims));

[sMin,CvMDistMin,exitCode]=quasiNewtonLBFGS(f,sInit(:),numCorr,epsilon,deltaTestDist,delta,lineSearchParams,maxIterations);
%The vectors were stacked. This reshapes them.
sMin=reshape(sMin,sDims);

%Add in the constant terms in the CVM distance for the return value.
CvMDistMin=CvMDistMin+computeD1(bMax,numDim,AbsTol,RelTol);

if(~isEven)
   CvMDistMin=CvMDistMin-2*computeDo2ContTerm(numDim,numHalfSamples,bMax,AbsTol,RelTol);
end

%Terminate if an error occurred.
if(exitCode<0)
    xi=[];
    w=[];
    return;
end

%Mirror the points to get all of the samples and add the zero point, if an
%odd number of samples is used.
if(isEven)
    xi=[sMin,-sMin]; 
else
    xi=[zeros(numDim,1),sMin,-sMin];
end

%The samples all have the same uniform weight. For simplicity, we will use
%a scalar for the computations here, but a vector will be returned.
w=1/numSamples;

%Force the covariance matrix of the samples to be the identity matrix,
%under the assumption that the mean of the samples is zero, which is true,
%because of the symmetry constraint.
if(forceCovMatch)
    R=w*(xi*xi');%The current covariance matrix.

    if(rank(R)~=numDim)%The covariance matrix should not be singular.
        xi=[];
        w=[];
        exitCode=-2;
        return;
    end

    %These two lines force the second moment to be the identity matrix.
    SRInv=chol(inv(R),'lower');
    xi=SRInv'*xi;
end

%Duplicate the w vector to be returned.
w=w*ones(numSamples,1);
end

function D=modCvMDist(s,bMax,isEven,AbsTol,RelTol,sDims)
%Compute the modified Cramér-von Mises distance between the localized
%cumulative distribution (LCD) of a Gaussian and the LCD of a Dirac
%mixture as discussed in Section 3.2. The input s is stacked and must be
%reshaped. This modified version of the cost function omits the constant D1
%term and the constant part of the D2 term for an odd number of samples.

s=reshape(s,sDims);

%The sample independent part, D1 (pg. 9), is the same regardless of whether
%the number of dimensions is even or odd. It also does not depend on the
%samples. Thus, it is omitted from the cost function. The same goes for the
%sample-independent part of the D2 term for an odd number of samples.

%If an even number of dimensions is requested
if(isEven)
    D2=computeDe2(s,bMax,AbsTol,RelTol);
    D3=computeDe3(s,bMax);
else
    
    D2=computeDo2Simp(s,bMax,AbsTol,RelTol);
    D3=computeDo3(s,bMax);
end

%The cost without the D1 term.
D=-2*D2+D3;
end

function DGrad=modCvMDistGrad(s,bMax,isEven,AbsTol,RelTol,sDims)
%Compute the gradient of the modified CFM distance between the LCD and the
%Dirac mixture. This is the derivative of D from modCvMDist with respect to
%every single point in s. The output, DGrad, has the same dimensionality as
%s. Each entry is the derivative with respect to the corresponding entry in
%s. The input s is stacked and must be reshaped.

s=reshape(s,sDims);

%If an even number of dimensions is requested
if(isEven)
    DGrad=-2*computeDe2Grad(s,bMax,AbsTol,RelTol);
    DGrad=DGrad+computeDe3Grad(s,bMax);
else
    DGrad=-2*computeDo2Grad(s,bMax,AbsTol,RelTol);
    DGrad=DGrad+computeDo3Grad(s,bMax);
end

%Stack the results.
DGrad=DGrad(:);
end

function De1=computeD1(bMax,N,AbsTol,RelTol)
%This computes the D1 parameter used on page 9 of "The Smart Sampling
%Kalman Filter with Symmetric Samples". The integral can be expressed in
%terms of an incomplete beta function or in terms of a hypergeometric
%function, but there does not appear to be any good way to evaluate the
%integral except using numerical integration.

%The function over which integration is performed.
f=@(b)(b.*(b.^2./(1+b.^2)).^(N/2));

%Numerically evaluate the integral.
De1=integral(f,0,bMax,'AbsTol',AbsTol,'RelTol',RelTol);
end


function De2=computeDe2(s,bMax,AbsTol,RelTol)
%Compute the D2 term for an even dimensionality using the formula on page 9
%with numerical integration.

N=size(s,1);%Dimensionality of the samples.
L=size(s,2);%Number of Dirac samples.

%Numerically evaluate the integral.
De2=integral(@(b)intFun(b),0,bMax,'AbsTol',AbsTol,'RelTol',RelTol);

%The function over which integration is performed.
    function vals=intFun(bVec)
        vals=zeros(size(bVec));
        numB=length(bVec);
        
        for curB=1:numB
            b=bVec(curB);
            vals(curB)=(b/L)*((2*b^2/(1+2*b^2))^(N/2))*sum(exp(-(1/2)*sum(s.*s,1)/(1+2*b^2)));
        end
    end
end


function De3=computeDe3(s,bMax)
%Compute the D3 term for an even dimensionality using the formula in
%Theorem 3.1.

L=size(s,2);%Number of Dirac samples.

De3=0;

for i=1:L
    for j=1:L
        diffSquared=norm(s(:,i)-s(:,j))^2;
        sumSquared=norm(s(:,i)+s(:,j))^2;
        
        arg1=-(1/2)*diffSquared/(2*bMax^2);
        arg2=-(1/2)*sumSquared/(2*bMax^2);
        
        term1=(bMax^2/2)*(exp(arg1)+exp(arg2));
        term2=(1/8)*(diffSquared*Ei(arg1)+sumSquared*Ei(arg2));

        De3=De3+term1+term2;
    end
end

De3=(2/(2*L)^2)*De3;
end

function Do2=computeDo2Simp(s,bMax,AbsTol,RelTol)
%Compute the D2 term for an odd dimensionality using the formula on page
%10. The constant term involving numerical integration is omitted.

N=size(s,1);%Dimensionality of the samples.
L=size(s,2);%Number of Dirac samples.

%Evaluate the integral to 14 places and add the even term.
Do2=((2*L)/(2*L+1))*computeDe2(s,bMax,AbsTol,RelTol);
end

function val=computeDo2ContTerm(N,L,bMax,AbsTol,RelTol)
%Compute the constant part of the D2 term for an odd dimensionality using
%the formula on page 10 with numerical integration.
    val=integral(@(b)intFun(b),0,bMax,'AbsTol',AbsTol,'RelTol',RelTol);
    
%The function over which integration is performed.
    function vals=intFun(bVec)
        vals=zeros(size(bVec));
        numB=length(bVec);
        
        for curB=1:numB
            b=bVec(curB);
            vals(curB)=(b/(2*L+1))*(2*b^2/(1+2*b^2))^(N/2);
        end
    end
end

function Do3=computeDo3(s,bMax)
%Compute the D3 term for an odd dimensionality using the formula in
%Theorem 3.2.

L=size(s,2);%Number of Dirac samples.

term1=((2*L)^2/(2*L+1)^2)*computeDe3(s,bMax);
term2=bMax^2/(2*(2*L+1)^2);

term3=0;
for i=1:L
    siSquared=norm(s(:,i))^2;
    arg1=-(1/2)*siSquared/(2*bMax^2);
    
    term3=term3+(bMax^2/2)*exp(arg1)+(1/8)*siSquared*Ei(arg1);
end
term3=term3*(4/(2*L+1)^2);

Do3=term1+term2+term3;
end

function De2Grad=computeDe2Grad(s,bMax,AbsTol,RelTol)
%This computes the gradient of the D2 term for every element of every Dirac
%delta point using numerical integration for an even dimensionality. The
%formula at the bottom of page 10 is used.

N=size(s,1);%Dimensionality of the samples.
L=size(s,2);%Number of Dirac samples.

%Allocate space.
De2Grad=zeros(N,L);
for i=1:L
    f=@(b)intFun(b,i);
    %The integral term in the same for all dimensions of this point.
    intVal=integral(f,0,bMax,'AbsTol',AbsTol,'RelTol',RelTol);
    
    De2Grad(:,i)=-(s(:,i)/(2*L))*intVal;
end

%The function over which integration is performed.
    function vals=intFun(bVec,i)
        vals=zeros(size(bVec));
        numB=length(bVec);
        
        for curB=1:numB
            b=bVec(curB);
            vals(curB)=((2*b)/(1+2*b^2))*(2*b^2/(1+2*b^2))^(N/2)*exp(-(1/2)*dot(s(:,i),s(:,i))/(1+2*b^2));
        end
    end
end

function De3Grad=computeDe3Grad(s,bMax)
%This computes the gradient of the D3 term for every element of every Dirac
%delta point using numerical integration for an even dimensionality. The
%formula in Theorem 3.3 is used.

N=size(s,1);%Dimensionality of the samples.
L=size(s,2);%Number of Dirac samples.

%Allocate space.
De3Grad=zeros(N,L);

for i=1:L
    for j=1:L
        diffVal=s(:,i)-s(:,j);
        sumVal=s(:,i)+s(:,j);
        
        De3Grad(:,i)=De3Grad(:,i)+diffVal*Ei(-(1/2)*dot(diffVal,diffVal)/(2*bMax^2))+sumVal*Ei(-(1/2)*dot(sumVal,sumVal)/(2*bMax^2));
    end
end

De3Grad=(1/(2*L)^2)*De3Grad;
end

function Do2Grad=computeDo2Grad(s,bMax,AbsTol,RelTol)
%This computes the gradient of the D2 term for every element of every Dirac
%delta point using numerical integration for an odd dimensionality. The
%formula on page 11 is used.

L=size(s,2);%Number of Dirac samples.

Do2Grad=(2*L/(2*L+1))*computeDe2Grad(s,bMax,AbsTol,RelTol);
end

function Do3Grad=computeDo3Grad(s,bMax)
%This computes the gradient of the D3 term for every element of every Dirac
%delta point using numerical integration for an odd dimensionality. The
%formula in Theorem 3.4 is used.

L=size(s,2);%Number of Dirac samples.

Do3Grad=((2*L)^2/(2*L+1)^2)*computeDe3Grad(s,bMax)+bsxfun(@times,s/(2*L+1)^2,Ei(-(1/2)*sum(s.*s,1)/(2*bMax^2)));
end

function val=Ei(x)
%The Cauchy principle value exponential integral for positive values. It
%has been redefined at the zero point to deal with the i=j instances in the
%double sums.
    val=-real(expint(-x));
    val(~isfinite(val))=0;
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
