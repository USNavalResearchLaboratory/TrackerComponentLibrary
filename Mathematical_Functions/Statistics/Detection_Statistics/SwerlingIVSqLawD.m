classdef SwerlingIVSqLawD
%%SWERLINGIVSQLAWD This class implements various functions and detections
%                statistics related to the Swerling IV target model when
%                used with a square law detector. All methods are static
%                (the class holds no information and does not need to be
%                instantiated). All values are with respect to a normalized
%                noise model, of which two normalization are available.
%
%Implemented methods are: mean, var, PDF, CDF, PD4Threshold,
%                         avgSNR4PDThresh, PD4PFA, rand
%
%Swerling models are given in [1] and in Chapter 11 of [2]. The Swerling IV
%model assumes that the observed power SNR of the target varies in terms of
%a central gamma distribution with shape parameter k=2 and scale parameter
%theta=avgSNR/2, that in a pulse train for detection the target SNR
%fluctuates from pulse to pulse, and that a square law detector is used to
%incoherently integrate multiple pulses The difference from a Swerling III
%model is the target amplitude fluctuations from pulse to pulse.
%
%The square law detector is explained in the comments to the NonFlucSqLawD
%class. The square law detector for N samples is
% y=sum_{i=1}^N abs(r_i)^2
%where r_i is the ith complex sampled value (including a signal and
%noise as described below) and y can be compared to a threshold. Two
%normalization are available. These are
%Definition 1:
%A is the SNR and the signal and noise model for a signle pulse is
% zSignal=sqrt(A)*exp(1j*2*pi*rand(1)) 
% zNoise=(randn(1)+1j*randn(1))/sqrt(2)
% r=zSigma+zNoise.
%so the noise is circularly symmetric complex Gaussian with variance 1.
%Definition 0:
%The model for a single pulse is now
% zSignal=sqrt(2*A)*exp(1j*2*pi*rand(1))
% zNoise=randn(1)+1j*randn(1)
% % r=zSigma+zNoise.
%So the amplitude squared of the noise is now 2, not 1.
%
%REFERENCES:
%[1] P. Swerling, "Probability of detection for fluctuating targets," The
%    RAND Corporation, Santa Monica, CA, Tech. Rep. RM-1217, 1954.
%[2] J. V. Di Franco and W. L. Rubin, Radar Detection. SciTech Publishing
%    Inc., Rayliegh, NC: 2004.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

methods(Static)
function val=mean(avgSNR,N,ampDef)
%%MEAN Obtain the mean of the distribution of the detection power
%      (normalized to a noise variance of 1) under a Swerling IV model in a
%      square-law detector.
%
%INPUTS:avgSNR The average power signal to noise ratio of the target.
%        N The number of pulses that are to be incoherently added for
%          detection (in a square-law detector). In a Swerling IV model, a
%          different realization of the target power is used across each
%          pulse and the noise corrupting the pulses varies. If this
%          parameter is omitted or an empty matrix is passed, N=1 is used.
%   ampDef This specified normalization (see help SwerlingIVSqLawD).
%          Possible values are:
%          0 The expected value of the squared magnitude of the noise is 2.
%          1 (The default if omitted or an empty matrix is passed) The
%            expected value of the squared magnitude of the noise is 1.
%            (A more common definition). 
%
%OUTPUTS: val The value of the mean for the given parameters.
%  
%The mean is given in Equation 11.5-20 of [1] where equation 11.4-8 relates
%the parameterization used there to an average squared amplitude. To
%relate the average squared amplitude ABar^2 to avgSNR as used in [2], we
%use avgSNR=ABar^2/2.
%
%EXAMPLE:
%Here, we validate the mean by generating random samples and comparing the
%computed mean to the sample mean.
% avgSNR=2;
% N=4;
% ampDef=1;
% numSamples=1e5;
% meanVal=SwerlingIVSqLawD.mean(avgSNR,N,ampDef)
% meanSampVal=mean(SwerlingIVSqLawD.rand([numSamples,1],avgSNR,N,ampDef))
%One will see that both mean values are about 12.
%
%REFERENCES:
%[1] J. V. Di Franco and W. L. Rubin, Radar Detection. SciTech Publishing
%    Inc., Rayliegh, NC: 2004.
%[2] P. Swerling, "Probability of detection for fluctuating targets," The
%    RAND Corporation, Santa Monica, CA, Tech. Rep. RM-1217, 1954.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
    
    if(nargin<3||isempty(ampDef))
        ampDef=1;
    end

    if(nargin<2||isempty(N))
        N=1;
    end
    
    if(ampDef==0)
        val=2*N*(1+avgSNR);
    else
        val=N*(1+avgSNR);
    end
end

function val=var(avgSNR,N,ampDef)
%%VAR Obtain the variance of the distribution of the detection power
%     (normalized to a noise variance of 1) under a Swerling IV model in a
%     square-law detector.
%
%INPUTS:avgSNR The average power signal to noise ratio of the target.
%        N The number of pulses that are to be incoherently added for
%          detection (in a square-law detector). In a Swerling IV model, a
%          different realization of the target power is used across each
%          pulse and the noise corrupting the pulses varies. If this
%          parameter is omitted or an empty matrix is passed, N=1 is used.
%   ampDef This specified normalization (see help SwerlingIVSqLawD).
%          Possible values are:
%          0 The expected value of the squared magnitude of the noise is 2.
%          1 (The default if omitted or an empty matrix is passed) The
%            expected value of the squared magnitude of the noise is 1.
%            (A more common definition). 
%
%OUTPUTS: val The value of the variance for the given parameters.
%
%The variance is given in Equation 11.5-22b of [1] where equation 11.4-8
%relates the parameterization used there to an average squared amplitude.
%To relate the average squared amplitude ABar^2 to avgSNR as used in [2],
%we use avgSNR=ABar^2/2. However, the variance from that expression has to
%be multiplied by 4, because as seen in Equation 10.4-4 of [1], the
%definition in is in terms of a scaled version of the square law detector.'
%
%EXAMPLE:
%Here, we validate the variance by generating random samples and comparing
%the computed variance to the sample variance.
% avgSNR=2;
% N=4;
% ampDef=0;
% numSamples=1e5;
% varVal=SwerlingIVSqLawD.var(avgSNR,N,ampDef)
% varSampVal=var(SwerlingIVSqLawD.rand([numSamples,1],avgSNR,N,ampDef))
%One will see that both variance values are about 112.
%
%REFERENCES:
%[1] J. V. Di Franco and W. L. Rubin, Radar Detection. SciTech Publishing
%    Inc., Rayliegh, NC: 2004.
%[2] P. Swerling, "Probability of detection for fluctuating targets," The
%    RAND Corporation, Santa Monica, CA, Tech. Rep. RM-1217, 1954.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.

    if(nargin<3||isempty(ampDef))
        ampDef=1;
    end

    if(nargin<2||isempty(N))
        N=1;
    end
    
    if(ampDef==0)
        val=4*N*(1+avgSNR.*(2+avgSNR/2));
    else
        val=N*(1+avgSNR.*(2+avgSNR/2));
    end
end
    
function val=PDF(x,avgSNR,N,algorithm,ampDef)
%%PDF Evaluate the scalar probability density function (PDF) of the
%     distribution of the detection power (normalized to a noise variance
%     of 1) under a Swerling IV model in a square-law detector.
%
%INPUTS: x The point or points at which the PDF should be evaluated. v>=0
%          for nonzero PDF values.
%   avgSNR The average power signal to noise ratio of the target.
%        N The number of pulses that are to be incoherently added for
%          detection (in a square-law detector). In a Swerling IV model, a
%          different realization of the target power is used across each
%          pulse and the noise corrupting the pulses varies. If this
%          parameter is omitted or an empty matrix is passed, N=1 is used.
% algorithm This selects the algorithm to use. Possible values are:
%          0 (The default if omitted or an empty matrix is passed). Invert
%            the characteristic function in [3] using a partial fraction
%            expansion and an inverse Laplace transform. This is only
%            implemented here for the case where all pulses have the same
%            SNR.
%          1 Use an algorthm that makes use of the hypergeometric1F1
%            function. This method is slow and can be inaccurate if the
%            hypergeometric1F1 does not converge.
%   ampDef This specified normalization (see help SwerlingIVSqLawD).
%          Possible values are:
%          0 The expected value of the squared magnitude of the noise is 2.
%          1 (The default if omitted or an empty matrix is passed) The
%            expected value of the squared magnitude of the noise is 1.
%            (A more common definition). 
%
%OUTPUTS: val The values of the PDF evaluated at the given points.
%
%For algorithm 1, the PDF is given in Equation 11.5-14 where equation
%11.4-8 relates the parameterization used there to an average squared
%amplitude. To relate the average squared amplitude ABar^2 to avgSNR as
%used in [2], we use avgSNR=ABar^2/2. However the PDF is derived is based
%on a scaled square-law detector, as in Equation 10.4-4 in [1]. Thus, in
%implementing the PDF, the scaling is removed.
%
%EXAMPLE:
%Here, we validate the PDF by generating random samples and comparing the
%PDF plot with a histogram of the random samples.
% avgSNR=2;
% N=4;
% ampDef=0;
% numSamples=1e4;
% 
% figure(1)
% clf
% histogram(SwerlingIVSqLawD.rand([numSamples,1],avgSNR,N,ampDef),'Normalization','pdf')
% hold on
% numPoints=1000;
% x=linspace(0,80,numPoints);
% vals=SwerlingIVSqLawD.PDF(x,avgSNR,N,[],ampDef);
% plot(x,vals,'linewidth',2)
% h1=xlabel('x');
% h2=ylabel('PDF(x)');
% set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')
%One will see that the histogram matches well with the plot.
%
%REFERENCES:
%[1] J. V. Di Franco and W. L. Rubin, Radar Detection. SciTech Publishing
%    Inc., Rayliegh, NC: 2004.
%[2] P. Swerling, "Probability of detection for fluctuating targets," The
%    RAND Corporation, Santa Monica, CA, Tech. Rep. RM-1217, 1954.
%[3] R. Kassab, T. Boutin, and C. Adnet, "Probability of detection for
%    Swerling model fluctuating targets with a square-law detector and
%    different signal to noise ratios," in 22nd International Microwave and
%    Radar Conference, Poznan, Poland, 14-17 May 2018, pp. 131-132.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.

    if(nargin<5||isempty(ampDef))
        ampDef=1;
    end

    if(nargin<4||isempty(algorithm))
        algorithm=0;
    end

    if(nargin<3||isempty(N))
        N=1;
    end

    %Initial change of variable to undo the scaling from Equation 10.4-4
    %in [1].
    if(ampDef==0)
        x=x/2;
    end

    xBar2=avgSNR/2;
    val=zeros(size(x));
    numEls=numel(x);

    switch(algorithm)
        case 0%Invert the characteristic function.
            numReps=2*N;%All denominators are squared.
            bInv=1./(1+avgSNR/2);
            scaleFactor=bInv.^(numReps);
            numPoly=1;
            for k=1:N
                numPoly=conv(numPoly,[1,1]);
            end
            [coeffs,poles]=partialFracKnownPoleDenom(-bInv,numReps,numPoly);
            coeffs=coeffs*scaleFactor;
            
            for k=1:numReps
                curPole=-poles(k);

                if(k==1)
                    val=val+coeffs(k)*exp(-curPole*x);
                else
                    val=val+coeffs(k)*(x.^(k-1)/factorial(k-1)).*exp(-curPole*x);
                end
            end

            val(x<0)=0;
            %If finite precision errors push values below zero.
            val(val<0)=0;
        case 1%Use the solution in terms of hypergeometric1F1.
            for curEl=1:numEls
                v=x(curEl);
        
                if(v>0)
                    val(curEl)=exp((N-1)*log(v)-v/(1+xBar2)-2*N*log(1+xBar2)-gammaln(N))*hypergeometric1F1(-N,N,-xBar2/(1+xBar2)*v);
                elseif(v==0&&N==1)
                    val(curEl)=1/((1+xBar2)^(2*N)*factorial(N-1));
                end
            end
        otherwise
            error('Unknown Algorithm Chosen.')
    end
    %Adjust for the change of variable to undo the scaling from Equation
    %10.4-4 in [1].
    if(ampDef==0)
        val=val/2;
    end
end

function val=CDF(v,avgSNR,N,ampDef)
%%CDF Evaluate the scalar cumulative distribution function (CDF) of the
%     distribution of the detection power (normalized to a noise variance
%     of 1) under a Swerling IV model in a square-law detector.
%
%INPUTS: v The point or points at which the CDF should be evaluated. v>0
%          for nonzero CDF values.
%   avgSNR The average power signal to noise ratio of the target.
%        N The number of pulses that are to be incoherently added for
%          detection (in a square-law detector). In a Swerling IV model, a
%          different realization of the target power is used across each
%          pulse and the noise corrupting the pulses varies. If this
%          parameter is omitted or an empty matrix is passed, N=1 is used.
%   ampDef This specified normalization (see help SwerlingIVSqLawD).
%          Possible values are:
%          0 The expected value of the squared magnitude of the noise is 2.
%          1 (The default if omitted or an empty matrix is passed) The
%            expected value of the squared magnitude of the noise is 1.
%            (A more common definition). 
%
%OUTPUTS: val The values of the CDF evaluated at the given points.
%
%The derivation of the detection probability in [1] is as 1-the value of
%the CDF. Thus, this function just evaluated 
%1-SwerlingIVSqLawD.PD4Threshold(avgSNR,v,N,ampDef).
%
%EXAMPLE:
%Here, the CDF returned by this function is plotted along with an empirical
%CDF constructed from random samples of this distribution. They agree well.
% avgSNR=5;
% N=4;
% ampDef=1;
% numSamples=1e5;
% samples=SwerlingIVSqLawD.rand([numSamples,1],avgSNR,N,ampDef);
% numPts=1000;
% x=linspace(0,80,numPts);
% CDFEmp=EmpiricalD.CDF(x,samples);
% CDF=SwerlingIVSqLawD.CDF(x,avgSNR,N,ampDef);
% figure(1)
% clf
% hold on
% plot(x,CDFEmp,'linewidth',4);
% plot(x,CDF,'linewidth',2);
% legend('Empirical CDF','Theoretical CDF')
%
%REFERENCES:
%[1] P. Swerling, "Probability of detection for fluctuating targets," The
%    RAND Corporation, Santa Monica, CA, Tech. Rep. RM-1217, 1954.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.

    if(nargin<4||isempty(ampDef))
        ampDef=1;
    end

    if(nargin<3||isempty(N))
        N=1; 
    end

    val=1-SwerlingIVSqLawD.PD4Threshold(avgSNR,v,N,ampDef);
    val(v<=0)=0;
end

function PD=PD4Threshold(avgSNR,thresh,N,ampDef)
%%PD4THRESHOLD Determine the detection probability of a Swerling IV target
%           given its signal to noise ratio and the value of the
%           detection threshold, assuming that a square-law detector is
%           used. The noise in each sample is assumed to be Gaussian with
%           variance 1 (normalized noise power).
%
%INPUTS: avgSNR A vector or matrix of average power signal to noise ratios
%          of the target at which one wishes to evaluate the detection
%          probability.
%   thresh The scalar normalized detection threshold to use. This is the
%          threshold to use if the noise variance is 1.
%        N The number of pulses that are to be incoherently added for
%          detection (in a square-law detector). In a Swerling IV model, a
%          different realization of the target power is used across each
%          pulse and the noise corrupting the pulses varies. If this
%          parameter is omitted or an empty matrix is passed, N=1 is used.
%   ampDef This specified normalization (see help SwerlingIVSqLawD).
%          Possible values are:
%          0 The expected value of the squared magnitude of the noise is 2.
%          1 (The default if omitted or an empty matrix is passed) The
%            expected value of the squared magnitude of the noise is 1.
%            (A more common definition). 
%
%OUTPUTS: PD The detection probability of the target.
%
%This function implements the sum in Equations 11.5-19 in [1]. However, the
%expression used is derived based on a scaled square-law detector, as in
%Equation 10.4-4. Thus, the threshold in this function is scaled
%appropriately depending on ampDef. There can be a loss of precision for
%0<avgSNR<eps(), though the solution is correct at 0.
%
%EXAMPLE:
%Here, we validate the results by comparing the PD from this function to
%the PD computed using random samples.
% avgSNR=8;
% thresh=30;
% N=4;
% ampDef=0;
% PD=SwerlingIVSqLawD.PD4Threshold(avgSNR,thresh,N,ampDef)
% numSamples=1e5;
% PDSamp=mean(SwerlingIVSqLawD.rand([numSamples,1],avgSNR,N,ampDef)>=thresh)
%One will see that both PD values are about 0.254.
%
%REFERENCES:
%[1] J. V. Di Franco and W. L. Rubin, Radar Detection. SciTech Publishing
%    Inc., Rayliegh, NC: 2004.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.

    if(nargin<4||isempty(ampDef))
        ampDef=1;
    end

    if(nargin<3||isempty(N))
        N=1; 
    end
    
    if(avgSNR==0)
        PD=NonFlucSqLawD.PD4Threshold(0,thresh,N,ampDef);
        return;
    end

    if(ampDef==0)
        thresh=thresh/2;
    end

    R=2*avgSNR;
    
    logR4=log(R/4);
    logCoeff=gammaln(N+1)-N*log(1+R/4);
    sumVal=0;
    for k=0:N
        curTerm=exp(logCoeff+k*logR4+log(gammainc(thresh/(1+R/4),k+N))-gammaln(k+1)-gammaln(N-k+1));

        sumVal=sumVal+curTerm;
    end

    PD=1-sumVal;
    if(~isfinite(PD)&&avgSNR<eps())
        if(ampDef==0)
            thresh=thresh*2;
        end

        PD=NonFlucSqLawD.PD4Threshold(0,thresh,N,ampDef);
        return;
    end
end

function [avgSNR,exitCode]=avgSNR4PDThresh(PD,thresh,N,ampDef,convergParams)
%%AVGSNR4PDTHRESH Given a detection probabilty and a normalized detection
%       threshold (normalized in terms of the receiver noise having a unit
%       covariance), determine the average power signal to noise ratio
%       (SNR) of the target needed under a Swerling IV model.
%
%INPUTS: PD The detection probability of the target , 0<=PD<1.
%   thresh The scalar normalized detection threshold to use. This is the
%          threshold to use if the noise variance is 1.
%        N The number of pulses that are to be incoherently added for
%          detection (in a square-law detector). In a Swerling I model, a
%          single realization of the target power is used across all
%          pulses, but the noise corrupting the pulses varies. If this
%          parameter is omitted or an empty matrix is passed, N=1 is used.
%   ampDef This specified normalization (see help SwerlingISqLawD).
%          Possible values are:
%          0 The expected value of the squared magnitude of the noise is 2.
%          1 (The default if omitted or an empty matrix is passed) The
%            expected value of the squared magnitude of the noise is 1.
%            (A more common definition). 
% convergParams An optional structure that is only used if N>1. It holds
%          parameters that define how the algorithm converges. Possible
%          members are:
%          'XTol' and 'maxIter' These have the same name as the inputs in
%           bisectionRootFind. See the comments to that function for
%           details. Default values are eps() and 100.
%          'maxIterSearch' This function uses a crude initial guess for an
%           upper bound on the solution and keeps doubling it until the PD
%           found is too large. This is the maximum number of doublings
%           that can be performed before an error is returned. The initial
%           guess is always avgSNR=100 and the default maximum number of
%           doublings, is 50.
%
%OUTPUTS: avgSNR The average SNR to give the specified PD for the specified
%                threshold. If no solution is possible (the PD is lower
%                than the PD for noise-only), then an empty matrix is
%                returne.
%       exitCode If a solution exists, then this is the exit code returned
%                by the bisectionRootFind function. Otherwise, this is -1.
%                See the comments to bisectionRootFind for more details.
%
%This function usesBisectionRootFind to invert PD minus what is equivalent
%to the output of the PD4Threshold method.
%
%EXAMPLE:
%Here, we show that the results are consistent with the PD as drawn from
%random samples.
% thresh=5;
% ampDef=1;
% PD=0.75;
% N=6;
% avgSNR=SwerlingIVSqLawD.avgSNR4PDThresh(PD,thresh,N,ampDef);
% PDBack=SwerlingIVSqLawD.PD4Threshold(avgSNR,thresh,N,ampDef)
% numSamples=1e5;
% PDSamp=mean(SwerlingIVSqLawD.rand([numSamples,1],avgSNR,N,ampDef)>=thresh)
%
%January 2024 David F. Crouse, Naval Research Laboratory, Washington D.C.

maxIterSearch=50;
maxIter=100;
XTol=eps();

if(nargin>4&&~isempty(convergParams))
    if(isfield(convergParams,'maxIterSearch'))
        maxIterSearch=convergParams.maxIterSearch;
    end
    if(isfield(convergParams,'maxIter'))
        maxIter=convergParams.maxIter;
    end
    if(isfield(convergParams,'XTol'))
        XTol=convergParams.XTol;
    end
end

if(nargin<7||isempty(maxIterSearch))
    maxIterSearch=50;
end

if(nargin<6||isempty(maxIter))
    maxIter=100;
end

if(nargin<5||isempty(XTol))
    XTol=eps();
end

if(nargin<4||isempty(ampDef))
    ampDef=1;
end

if(nargin<3||isempty(N))
    N=1;
end

if(NonFlucSqLawD.PD4Threshold(0,thresh,N,ampDef)>PD)
    %If the PD is below the PD obtained with a 0 SNR target, then no
    %solution is possible.
    avgSNR=[];
    exitCode=-1;
    return
end

%We need to find an upper bound. We start with an estimate of 100 and then
%keep doubling it until we have found an upper bound.
avgSNRUpper=100;
PDComp=SwerlingIVSqLawD.PD4Threshold(avgSNRUpper,thresh,N,ampDef);
curIter=0;
while(PDComp<PD)
    curIter=curIter+1;
    if(curIter>maxIterSearch)
        error('Unable to bracket a solution.')
    end
    avgSNRUpper=2*avgSNRUpper;
    PDComp=SwerlingIVSqLawD.PD4Threshold(avgSNRUpper,thresh,N,ampDef);
end

f=@(avgSNR)(PD-SwerlingIVSqLawD.PD4Threshold(avgSNR,thresh,N,ampDef));
[avgSNR,~,exitCode]=bisectionRootFind(f,[0;avgSNRUpper],XTol,maxIter);

end

function PD=PD4PFA(avgSNR,PFA,N,ampDef)
%%PD4PFA Determine the detection probability of a Swerling IV target given
%        its signal to noise ratio and the probability of false alarm
%        implied by a particular unspecified detection threshold, assuming
%        that a square law detector is used.
%
%INPUTS: avgSNR A vector or matrix of average power signal to noise ratios
%          of the target at which one wishes to evaluate the detection
%          probability.
%      PFA The probability of false alarm, 0<=PFA<1.
%        N The number of pulses that are to be incoherently added for
%          detection (in a square-law detector). In a Swerling IV model, a
%          different realization of the target power is used across each
%          pulse and the noise corrupting the pulses varies. If this
%          parameter is omitted or an empty matrix is passed, N=1 is used.
%   ampDef This specified normalization (see help SwerlingIVSqLawD).
%          Possible values are:
%          0 The expected value of the squared magnitude of the noise is 2.
%          1 (The default if omitted or an empty matrix is passed) The
%            expected value of the squared magnitude of the noise is 1.
%            (A more common definition). 
%
%OUTPUTS: PD The detection probability of the target, 0<=PD<=1.
%
%This function just calls PFA2SquareLawThreshold to convert the false alarm
%rate into a normalized detection threshold after which it calls
%SwerlingIVSqLawD.PD4Threshold.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.

    if(nargin<4||isempty(ampDef))
        ampDef=1;
    end

    if(nargin<3||isempty(N))
       N=1; 
    end

    thresh=PFA2SquareLawThreshold(PFA,N,ampDef);
    PD=SwerlingIVSqLawD.PD4Threshold(avgSNR,thresh,N,ampDef);
end

function vals=rand(NDims,avgSNR,N,ampDef)
%%RAND Generate random variables representing a target whose sampled signal
%      to average noise power ratio is generated according to a Swerling IV
%      model detected by a square-law detector.
%
%INPUTS: NDims If NDims is a scalar, then rand returns an NDimsXNDims
%          matrix of random variables. If NDims=[M,N1] is a two-element row
%          vector, then rand returns an MXN1 matrix of random variables.
%   avgSNR The average power signal to noise ratio of the target.
%        N The number of pulses that are to be incoherently added for
%          detection (in a square-law detector). In a Swerling IV model, a
%          different realization of the target power is used across each
%          pulse and the noise corrupting the pulses varies. If this
%          parameter is omitted or an empty matrix is passed, N=1 is used.
%   ampDef This specified normalization (see help SwerlingIVSqLawD).
%          Possible values are:
%          0 The expected value of the squared magnitude of the noise is 2.
%          1 (The default if omitted or an empty matrix is passed) The
%            expected value of the squared magnitude of the noise is 1.
%            (A more common definition). 
%
%OUTPUTS: vals A matrix whose dimensions are determined by N of the
%              generated random variables.
%
%Following the model of [1], a sample from the central gamma distribution
%is performed to obtain the complex signal to noise ratio of every single
%one of the N pulses. Then, the square law detector output conditioned on
%the signal to noise ratio is generated directly using NonFlucSqLawD.rand
%for each pulse and the values are added.
%
%REFERENCES:
%[1] P. Swerling, "Probability of detection for fluctuating targets," The
%    RAND Corporation, Santa Monica, CA, Tech. Rep. RM-1217, 1954.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.    
    
    if(nargin<4||isempty(ampDef))
        ampDef=1;
    end

    if(nargin<3||isempty(N))
        N=1; 
    end

    if(isscalar(NDims))
        dims=[NDims, NDims];
    else
        dims=NDims;
    end

    vals=zeros(dims);
    
    numEls=numel(vals);
    for curEl=1:numEls
        for curPulse=1:N        
            %Equation I.2 for  the input signal-to-noise power ratio. It is
            %different for each sample in the pulse train in the Swerling
            %IV model.
            
            curSNR=GammaD.rand(1,2,avgSNR/2);
            vals(curEl)=vals(curEl)+NonFlucSqLawD.rand(1,curSNR,1,ampDef);
        end
    end  
end
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
