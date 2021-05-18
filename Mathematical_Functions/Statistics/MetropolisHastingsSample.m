function [xSamples,acceptRate]=MetropolisHastingsSample(xInit,numSamples,targetPDF,transPDF,transPDFSamp,numBurnIn,isSymmetric,takeEvery)
%%METROPOLISHASTINGSSAMPLE Generate samples using the Metropolis-Hastings
%               algorithm. Suppose that one wishes to determine moments or
%               other properties of a normalized target PDF given
%               targetPDF, which is the not necessarily normalized form of
%               the PDF. The Metropolis Hastings algorithm forms a Markov
%               chain of states using a transition PDF transPDF(a,b)
%               that gives the likelihood of going from state a to state b
%               such that the limiting distribution of the states is
%               targetPDF (states are column vectors). This means that
%               one can find the mean, covariance, etc. by using a large
%               number of samples from the chain. subsequent values in the
%               chain are, however, correlated, so if this algorithm is
%               used to generate "random" values of a distribution, samples
%               of the chain should be taken sufficiently far apart that
%               the autocorrelation is negligible. On the other hand, if
%               given a symmetric transition PDF, the chain always moves to
%               points of increasing likelihood in 
%
%INPUTS: xInit The xDimX1 initial value to use for the Markov chain. If
%              this is in an area of high likelihood of the target PDF,
%              then a smaller burn-in period can be used.
%   numSamples The number of samples (in a Markov chain) to generate.
%    targetPDF A function handle to the (not necessarily normalized) PDF
%              from which samples are desired. targetPDF(x) returns the
%              (nor necessarily normalized) likelihood of the point x.
%     transPDF A function handle for the transition PDF (transition
%              kernel). This has the form transPDF(x,y) and is the
%              likelihood of transitioning to state y given that one is in
%              state x. This must integrate to 1 over y. Another notable
%              requirement is that this is irreducible over time. That
%              means, that over time, one must be able to go from any state
%              to any other state in the target PDF, though that need not
%              be possible in any single step. Other conditions are
%              aperiodicity (it does not just repeat an deterministic cycle)
%              and positive recurrence. If isSymmetric=true, then an empty
%              matrix can be passed in place of this function, though
%              transPDFSamp is still required.
% transPDFSamp A function handle such that transPDFSamp(x) samples
%              transPDF(x,y) to get y.
%    numBurnIn The number of steps of the Metropolis-Hastings algorithm to
%              run before samples start being saved to be returned. If this
%              parameter is omitted or an empty matrix is passed, the
%              default of 0 is used (no burn-in at all).
%  isSymmetric This is true if transPDF(x,y)=transPDF(y,x). In such an
%              instance, an empty matrix can be passed for transPDF and a
%              more efficient algorithm can be used. If this parameter is
%              omitted or an empty matrix is passed, isSymmetric=false is
%              used.
%    takeEvery The samples produced by this function are correlated. If one
%              wishes to obtain uncorrelated samples of targetPDF, then
%              takeEvery should be set to a value at least equal to twice
%              the correlation time of the full (takeEvery=1) sequence). On
%              the other hand, if one just wishes to compute expected values
%              using the results, one should use takeEvery=1 as larger
%              values just increase the variance of the result (for the
%              same length of Markov chain generated to get the values).
%              takeEvery=n means that only every nth sample in the chain is
%              returned. If omitted or an empty matrix is passed,
%              takeEvery=1 is used.
%
%OUTPUTS: xSamples An xDimXnumSamples set of (possibly correlated) samples
%                  of targetPDF.
%       acceptRate The fraction of the total number of samples that were
%                  accepted in the acceptance-rejection step of the
%                  Metropolis-Hastings algorithm. This can be used when
%                  comparing transition densities as one generally wants as
%                  high an acceptance rate as possible.
%
%A commonly-used transition PDF is to let transPDF(x,y) be a (multivariate)
%Gaussian distribution with mean x. The covariance matrix should be
%somewhere on the order of the covariance matrix of targetPDF. Note that
%this transition PDF is symmetric.
%
%The basic Metropolis-Hastings algorithm is described in [1], Chapter 10.2
%of [2] and Chapter 3 of [3]. In Chapter 3 of [3], a more detailed analysis
%of the error performance of the algorithm is given. In Chapter 3.3.1 of
%[3], it is noted that a good simulation would tend to run fro many times
%the correlation time of the output.
%
%A use of the output of this function is based on the ergodic theorem that
%for a length-M sequence from the Metropolis-Hastings algorithm
%limit_{M->Inf} 1/M*sum_{m=1}^Mg(x_i)->integral_Xg(x)*targetPDF(x) dx
%Thus, one can use a sufficient number of terms of the output to perform
%numeric integration.
%
%EXAMPLE 1:
%One of the simplest common examples is sampling a normal distribution
%using a uniform transition density and the Metropolis Hastings algorithm.
%For simplicity, let us use the standard normal distribution. The uniform
%transition distribution must be able to cover all space (be irreducible).
%Thus, out transition density is actually a random walk -we sample from
%+/-delta around x. Thus, the algorithm is
% targetPDF=@(x)GaussianD.PDF(x);
% delta=1;%Arbitrary delta
% transPDF=@(x,y)UniformD.PDF(y,[x-delta;x+delta]);
% transPDFSamp=@(x)UniformD.rand(1,[x-delta;x+delta]);
% numSamples=25000;
% xInit=1;%An arbitrary initial value.
% %transPDF is actually not needed, because transPDF is actualyl symmetric
% %(once substituted). 
% xSamples=MetropolisHastingsSample(xInit,numSamples,targetPDF,[],transPDFSamp,[],true);
% %We will plot a histogram of the samples produced by the algorithm to see
% %if they look Gaussian
% figure()
% hold on
% h=histogram(xSamples);
% %Also plot the standard normal PDF for comparison. The scaling of the PDF
% %is to make it in line with the expected bin frequencies of random
% %samples.
% points=linspace(-4,4,1000);
% plot(points,h.BinWidth*numSamples*GaussianD.PDF(points),'-r','linewidth',2)
% %Also the variance is close to 1, as one would expect.
% var(xSamples)
%Thus, one can see the effectiveness of the algorithm. Note that there are
%much better ways to generate normal random variables. 
%
%EXAMPLE 2:
%Here, we do the same thing with a 2D Gaussian distribution and a uniform
%transition density as in Section 7.1 of [1]. In this instance, the mean mu
%and covariance matrix Sigma are not standard.
% mu=[1;2];
% Sigma=[1,  0.9;
%        0.9,1];
% delta=[0.75;1];
% targetPDF=@(x)GaussianD.PDF(x,mu,Sigma);
% transPDF=@(x,y)UniformD.PDF(y(1),[x-delta,x+delta]')
% transPDFSamp=@(x)UniformD.rand(1,[x-delta,x+delta]');
% numSamples=25000;
% xInit=[1;1];
% %transPDF is actually not needed, because transPDF is actualyl symmetric
% %(once substituted). 
% xSamples=MetropolisHastingsSample(xInit,numSamples,targetPDF,[],transPDFSamp,[],true);
% %We will plot the samples in one figure to show how well the
% %Metropolis-Hastings algorithm performs. We will then plot samples
% %generated the standard way using a Cholesky decomposition, and one can
% %compare the results.
% figure()
% scatter(xSamples(1,:),xSamples(2,:),'.b')
% %The standard method of generating normal random variables.
% cholSamples=chol(Sigma,'lower')*randn(2,numSamples);
% figure()
% scatter(cholSamples(1,:),cholSamples(2,:),'.r')
%
%EXAMPLE 3:
%In this example, we will use the Metropolis-Hastings algorithm to help
%with some Bayesian estimation. Suppose that we know the shape parameter of
%a gamma distribution k and we do not know the scale parameter theta.
%However, we have a Raylight prior on the scale parameter. We make a single
%observation obs and we would like to find the expected value of the shape
%parameter. This involves a dificult integral over a joint PDF that is not
%easy to normalize.
% lambda=0;%central gamma distribution.
% k=2;%Known shape parameter of the gamma distribution.
% obs=16;%Single observation of the gamma distribution.
% sigma=1;%Parameter of the Rayliegh prior on theta (scale parameter).
% targetPDF=@(theta)GammaD.PDF(obs,k,theta,lambda)*RayleighD.PDF(theta,sigma);
% %The transition density will be the density of the function abs(x+z) that
% %where z is a standard normal random variable. Because this transition
% %density is symmetric, we do not actually have to compute the density,
% %which would be difficult.
% transPDFSamp=@(x)abs(x+randn(1));
% numSamples=25000;
% xInit=1;%Arbitrary initial value.
% %In this instance, we use a "burn-in" length equal to the number of
% %samples that will be tkane for the mean.
% [xSamples,acceptRate]=MetropolisHastingsSample(xInit,numSamples,targetPDF,[],transPDFSamp,numSamples,true);
% thetaEst=mean(xSamples)
%One will find that the solution is in the neighborhood of a value that can
%be obtained by realizing that the solution is a ratio of two MeijerG
%functions. The solution obtained in such a manner using extended precision
%arithmetic is 
%2.5217563698462965861981671454276
%
%REFERENCES:
%[1] S. Chib and E. Greenberg, "Understanding the Metropolis-Hastings
%    algorithm," The American Statistician, vol. 49, no. 4, pp. 327-335,
%    Nov. 1995.
%[2] S. M. Ross, Simulation, 4th ed. Amsterdam: Academic Press, 2006.
%[3] Newman, M. E. J.; Barkema, G. T. Monte Carlo Methods in Statistical
%    Physics. Oxford: Clarendon Press, 2001.
%
%December 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<6||isempty(numBurnIn))
    numBurnIn=0;
end

if(nargin<7||isempty(isSymmetric))
    isSymmetric=false;
end

if(nargin<8||isempty(takeEvery))
    takeEvery=1;
end

xDim=size(xInit,1);

x=xInit;
piX=targetPDF(x);
numAccepted=0;

if(isSymmetric)%The case of a symmetric transition PDF
    %First, go through the burn-in stage.
    for curPoint=1:numBurnIn
        [x,piX,didAccept]=takeOneMHStepSymmetric(x,piX,targetPDF,transPDFSamp);
        numAccepted=numAccepted+didAccept;
    end

    %Allocate space for the samples.
    xSamples=zeros(xDim,numSamples);

    %Next, start collecting values until we have enough.
    numCollected=0;
    curStep=0;
    while(numCollected<numSamples)
        [x,piX,didAccept]=takeOneMHStepSymmetric(x,piX,targetPDF,transPDFSamp);

        if(mod(curStep,takeEvery)==0)
            numCollected=numCollected+1;
            xSamples(:,numCollected)=x;
        end
        numAccepted=numAccepted+didAccept;
        curStep=curStep+1;
    end
else%The case for a non-symmetric transition PDF.
    %First, go through the burn-in stage.
    for curPoint=1:numBurnIn
        [x,piX,didAccept]=takeOneMHStep(x,piX,targetPDF,transPDF,transPDFSamp);
        numAccepted=numAccepted+didAccept;
    end

    %Allocate space for the samples.
    xSamples=zeros(xDim,numSamples);

    %Next, start collecting values until we have enough.
    numCollected=0;
    curStep=0;
    while(numCollected<numSamples)
        [x,piX,didAccept]=takeOneMHStep(x,piX,targetPDF,transPDF,transPDFSamp);

        if(mod(curStep,takeEvery)==0)
            numCollected=numCollected+1;
            xSamples(:,numCollected)=x;
        end
        numAccepted=numAccepted+didAccept;
        curStep=curStep+1;
    end
end

%Fraction of total samples accepted
acceptRate=numAccepted/(numBurnIn+curStep);

end

function [x,piX,didAccept]=takeOneMHStep(xPrev,piPrev,targetPDF,transPDF,transPDFSamp)
%%TAKEONEMHSTEP Take one step of the Metropolis-Hastings algorithm shown in
%               Section 4 of [1].

    x=xPrev;
    y=transPDFSamp(x);
    u=rand(1);

    piY=targetPDF(y);
    piX=piPrev;
    qYX=transPDF(y,x);
    qXY=transPDF(x,y);

    alpha=min([(piY*qYX)/(piX*qXY);1]);
    
    %If the sample should be accepted.
    if(u<=alpha)
        x=y;
        piX=piY;
        didAccept=true;
    else
        didAccept=false;
    end
end

function [x,piX,didAccept]=takeOneMHStepSymmetric(xPrev,piPrev,targetPDF,transPDFSamp)
%%TAKEONEMHSTEPSYMMETRIC Take one step of the Metropolis-Hastings algorithm
%               shown in Section 4 of [1] assuming that the transition pDF
%               is symmetric.

    x=xPrev;
    y=transPDFSamp(x);
    u=rand(1);

    piY=targetPDF(y);
    piX=piPrev;

    alpha=min([piY/piX;1]);

    %If the sample should be accepted.
    if(u<=alpha)
        x=y;
        piX=piY;
        didAccept=true;
    else
        didAccept=false;
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
