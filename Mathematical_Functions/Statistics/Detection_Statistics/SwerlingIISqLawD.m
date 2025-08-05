classdef SwerlingIISqLawD
%%SWERLINGIISQLAWD This class implements various functions and detections
%                statistics related to the Swerling II target model when
%                used with a square law detector. All methods are static
%                (the class holds no information and does not need to be
%                instantiated). All values are with respect to a normalized
%                noise model, of which two normalizations are available. 
%
%Implemented methods are: mean, var, PDF, CDF, PD4Threshold, PD4PFA,
%                         avgSNR4PDThresh, avgSNR4PDPFA, thresh4AvgSNRPD,
%                         rand
%
%Swerling models are given in [1] and in Chapter 11 of [2]. The Swerling II
%model assumes that the observed power SNR of the target varies in terms of
%an exponential distribution with rate parameter (1/avgSNR), that in a
%pulse train for detection the target SNR does fluctuate from pulse to
%pulse, and that a square law detector is used to incoherently integrate
%multiple pulses. The difference from a Swerling I model is the target
%amplitude fluctuates from pulse to pulse.
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
%The comments to NonFlucSqLawD describe the normalization definitions in
%more detail.
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
%      (normalized to a noise variance of 1) under a Swerling II model in a
%      square-law detector.
%
%INPUTS:avgSNR The average power signal to noise ratio of the target.
%        N The number of pulses that are to be incoherently added for
%          detection (in a square-law detector). In a Swerling II model, a
%          different realization of the target power is used across each
%          pulse and the noise corrupting the pulses varies. If this
%          parameter is omitted or an empty matrix is passed, N=1 is used.
%   ampDef This specified normalization (see help SwerlingIISqLawD).
%          Possible values are:
%          0 The expected value of the squared magnitude of the noise is 2.
%          1 (The default if omitted or an empty matrix is passed) The
%            expected value of the squared magnitude of the noise is 1.
%            (A more common definition). 
%
%OUTPUTS: val The value of the mean for the given parameters.
%
%The expression for the mean implemented here is E{x}, where the expected
%value is taken over the PDF in Equation III.10 of [1]. However, the first
%moment from that expression has to be multiplied by 2, because as seen in
%Equation 10.4-4 of [2], the definition in [2], which is the same as in
%[1], is in terms of a scaled version of the square law detector.
%
%EXAMPLE:
%Here, we validate the mean by generating random samples and comparing the
%computed mean to the sample mean.
% avgSNR=2;
% N=4;
% ampDef=0;
% numSamples=1e5;
% meanVal=SwerlingIISqLawD.mean(avgSNR,N,ampDef)
% meanSampVal=mean(SwerlingIISqLawD.rand([numSamples,1],avgSNR,N,ampDef))
%One will see that both mean values are about 24.
%
%REFERENCES:
%[1] P. Swerling, "Probability of detection for fluctuating targets," The
%    RAND Corporation, Santa Monica, CA, Tech. Rep. RM-1217, 1954.
%[2] J. V. Di Franco and W. L. Rubin, Radar Detection. SciTech Publishing
%    Inc., Rayliegh, NC: 2004.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.

    if(nargin<3||isempty(ampDef))
        ampDef=1;
    end

    if(nargin<2||isempty(N))
        N=1;
    end

    if(ampDef==0)
        val=2*N.*(1+avgSNR);
    else
        val=N.*(1+avgSNR);
    end
end

function val=var(avgSNR,N,ampDef)
%%VAR Obtain the variance of the distribution of the detection power
%     (normalized to a noise variance of 1) under a Swerling II model in a
%     square-law detector.
%
%INPUTS:avgSNR The average power signal to noise ratio of the target.
%        N The number of pulses that are to be incoherently added for
%          detection (in a square-law detector). In a Swerling II model, a
%          different realization of the target power is used across each
%          pulse and the noise corrupting the pulses varies. If this
%          parameter is omitted or an empty matrix is passed, N=1 is used.
%   ampDef This specified normalization (see help SwerlingIISqLawD).
%          Possible values are:
%          0 The expected value of the squared magnitude of the noise is 2.
%          1 (The default if omitted or an empty matrix is passed) The
%            expected value of the squared magnitude of the noise is 1.
%            (A more common definition). 
%
%OUTPUTS: val The value of the variance for the given parameters.
%  
%To derive the formula implemented here, the first and second noncentral
%moments were obtained by finding E{x^2} and E{x} integrating over the PDF
%in Equation III.10 in [1].  One then uses the identity that
%variance=E{x^2}-E{x}^2
%However, the variance from that expression has to be multiplied by 4,
%because as seen in Equation 10.4-4 of [2], the definition in [2], which is
%the same as the definition in [1], is in terms of a scaled version of the
%square law detector (for ampDef=0).
%
%EXAMPLE:
%Here, we validate the variance by generating random samples and comparing
%the computed variance to the sample variance.
% avgSNR=2;
% N=4;
% ampDef=0;
% numSamples=1e5;
% meanVal=SwerlingIISqLawD.var(avgSNR,N,ampDef)
% meanSampVal=var(SwerlingIISqLawD.rand([numSamples,1],avgSNR,N,ampDef))
%One will see that both variance values are about 144.
%
%REFERENCES:
%[1] P. Swerling, "Probability of detection for fluctuating targets," The
%    RAND Corporation, Santa Monica, CA, Tech. Rep. RM-1217, 1954.
%[2] J. V. Di Franco and W. L. Rubin, Radar Detection. SciTech Publishing
%    Inc., Rayliegh, NC: 2004.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
    
    if(nargin<3||isempty(ampDef))
        ampDef=1;
    end

    if(nargin<2||isempty(N))
        N=1;
    end
    
    if(ampDef==0)
        val=4*N*(1+avgSNR)^2;
    else
        val=N*(1+avgSNR)^2;
    end
end

function val=PDF(v,avgSNR,N,ampDef)
%%PDF Evaluate the scalar probability density function (PDF) of the
%     distribution of the detection power (normalized to a noise variance
%     of 1) under a Swerling II model in a square-law detector.
%
%INPUTS: v The point or points at which the PDF should be evaluated. v>0
%          for nonzero PDF values.
%   avgSNR The average power signal to noise ratio of the target.
%        N The number of pulses that are to be incoherently added for
%          detection (in a square-law detector). In a Swerling II model, a
%          different realization of the target power is used across each
%          pulse and the noise corrupting the pulses varies. If this
%          parameter is omitted or an empty matrix is passed, N=1 is used.
%   ampDef This specified normalization (see help SwerlingIISqLawD).
%          Possible values are:
%          0 The expected value of the squared magnitude of the noise is 2.
%          1 (The default if omitted or an empty matrix is passed) The
%            expected value of the squared magnitude of the noise is 1.
%            (A more common definition). 
%
%OUTPUTS: val The values of the PDF evaluated at the given points.
%
%The PDF is taken from Equation III.10 of [1].
%
%EXAMPLE:
%Here, we validate the PDF by generating random samples and comparing the
%PDF plot with a histogram of the random samples.
% avgSNR=2;
% N=4;
% numSamples=1e5;
% ampDef=0;
% 
% figure(1)
% clf
% histogram(SwerlingIISqLawD.rand([numSamples,1],avgSNR,N,ampDef),'Normalization','pdf')
% hold on
% numPoints=1000;
% x=linspace(0,80,numPoints);
% vals=SwerlingIISqLawD.PDF(x,avgSNR,N,ampDef);
% plot(x,vals,'linewidth',2)
% h1=xlabel('x');
% h2=ylabel('PDF(x)');
% set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')
%One will see that the histogram matches well with the plot.
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

    if(ampDef==0)
        v=v/2;
    end

    val=exp((N-1)*log(v)-N*log(1+avgSNR)-gammaln(N)-v/(1+avgSNR));
    
    if(ampDef==0)
        val=val/2;
    end

    val(v<0)=0;
end

function val=CDF(v,avgSNR,N,ampDef)
%%CDF Evaluate the scalar cumulative distribution function (CDF) of the
%     distribution of the detection power (normalized to a noise variance
%     of 1) under a Swerling II model in a square-law detector.
%
%INPUTS: v The point or points at which the CDF should be evaluated. v>0
%          for nonzero CDF values.
%   avgSNR The average power signal to noise ratio of the target.
%        N The number of pulses that are to be incoherently added for
%          detection (in a square-law detector). In a Swerling II model, a
%          different realization of the target power is used across each
%          pulse and the noise corrupting the pulses varies. If this
%          parameter is omitted or an empty matrix is passed, N=1 is used.
%   ampDef This specified normalization (see help SwerlingIISqLawD).
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
%1-SwerlingIISqLawD.PD4Threshold(avgSNR,v,N,ampDef).
%
%EXAMPLE:
%Here, the CDF returned by this function is plotted along with an empirical
%CDF constructed from random samples of this distribution. They agree well.
% avgSNR=5;
% N=4;
% ampDef=1;
% numSamples=1e5;
% samples=SwerlingIISqLawD.rand([numSamples,1],avgSNR,N,ampDef);
% numPts=1000;
% x=linspace(0,80,numPts);
% CDFEmp=EmpiricalD.CDF(x,samples);
% CDF=SwerlingIISqLawD.CDF(x,avgSNR,N,ampDef);
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

    val=1-SwerlingIISqLawD.PD4Threshold(avgSNR,v,N,ampDef);
    val(v<=0)=0;
end

function PD=PD4Threshold(avgSNR,thresh,N,ampDef)
%%PD4THRESHOLD Determine the detection probability of a Swerling II target
%           given its signal to noise ratio and the value of the
%           detection threshold, assuming that a square-law detector is
%           used.
%
%INPUTS: avgSNR A vector or matrix of average power signal to noise ratios
%          of the target at which one wishes to evaluate the detection
%          probability.
%   thresh The scalar normalized detection threshold to use. This is the
%          threshold to use if the noise variance is 1.
%        N The number of pulses that are to be incoherently added for
%          detection (in a square-law detector). In a Swerling II model, a
%          different realization of the target power is used across each
%          pulse and the noise corrupting the pulses varies. If this
%          parameter is omitted or an empty matrix is passed, N=1 is used.
%   ampDef This specified normalization (see help SwerlingIISqLawD).
%          Possible values are:
%          0 The expected value of the squared magnitude of the noise is 2.
%          1 (The default if omitted or an empty matrix is passed) The
%            expected value of the squared magnitude of the noise is 1.
%            (A more common definition). 
%
%OUTPUTS: PD The detection probability of the target.
%
%This function implements Equation II.10 in [1]. However, the expression
%used is derived based on a scaled square-law detector, as in Equation
%10.4-4 in [2], which uses the same definition as in [1]. Thus, the
%threshold in this function is scaled appropriately.
%
%EXAMPLE:
%Here, we validate the results by comparing the PD from this function to
%the PD computed using random samples.
% avgSNR=2;
% thresh=30;
% N=4;
% ampDef=0;
% PD=SwerlingIISqLawD.PD4Threshold(avgSNR,thresh,N,ampDef)
% numSamples=1e5;
% PDSamp=mean(SwerlingIISqLawD.rand([numSamples,1],avgSNR,N,ampDef)>=thresh)
%One will see that both PD values are about 0.265.
%
%REFERENCES:
%[1] P. Swerling, "Probability of detection for fluctuating targets," The
%    RAND Corporation, Santa Monica, CA, Tech. Rep. RM-1217, 1954.
%[2] J. V. Di Franco and W. L. Rubin, Radar Detection. SciTech Publishing
%    Inc., Rayliegh, NC: 2004.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.

    if(nargin<4||isempty(ampDef))
        ampDef=1;
    end

    if(nargin<3||isempty(N))
        N=1; 
    end

    if(ampDef==0)
        thresh=thresh/2;
    end
    
    PD=1-PearsonsGammaInc(thresh./((1+avgSNR)*sqrt(N)),N-1);
end

function PD=PD4PFA(avgSNR,PFA,N,ampDef)
%%PD4PFA Determine the detection probability of a Swerling II target given
%        its signal to noise ratio and the probability of false alarm
%        implied by a particular unspecified detection threshold, assuming
%        that a square law detector is used.
%
%INPUTS: avgSNR A vector or matrix of average power signal to noise ratios
%          of the target at which one wishes to evaluate the detection
%          probability.
%      PFA The probability of false alarm, 0<=PFA<1.
%        N The number of pulses that are to be incoherently added for
%          detection (in a square-law detector). In a Swerling II model, a
%          different realization of the target power is used across each
%          pulse and the noise corrupting the pulses varies. If this
%          parameter is omitted or an empty matrix is passed, N=1 is used.
%   ampDef This specified normalization (see help SwerlingIISqLawD).
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
%SwerlingIISqLawD.PD4Threshold.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.

    if(nargin<4||isempty(ampDef))
        ampDef=1;
    end

    if(nargin<3||isempty(N))
        N=1; 
    end

    thresh=PFA2SquareLawThreshold(PFA,N,ampDef);
    PD=SwerlingIISqLawD.PD4Threshold(avgSNR,thresh,N,ampDef);
end

function avgSNR=avgSNR4PDThresh(PD,thresh,N,ampDef)
%%AVGSNR4PDTHRESH Given a detection probabilty and a normalized detection
%       threshold (normalized in terms of the receiver noise having a unit
%       covariance), determine the average power signal to noise ratio
%       (SNR) of the target needed under a Swerling II model for one or
%       more pulses in a square-law detector.
%
%INPUTS: PD The detection probability of the target , 0<=PD<1.
%    thresh The scalar normalized detection threshold to use. This is the
%           threshold to use if the noise variance is 1.
%         N The number of pulses that are to be incoherently added for
%           detection (in a square-law detector). In a Swerling II model, a
%           different realization of the target power is used across each
%           pulse and the noise corrupting the pulses varies. If this
%           parameter is omitted or an empty matrix is passed, N=1 is used.
%    ampDef This specified normalization (see help SwerlingIISqLawD).
%           Possible values are:
%          0 The expected value of the squared magnitude of the noise is 2.
%          1 (The default if omitted or an empty matrix is passed) The
%            expected value of the squared magnitude of the noise is 1.
%            (A more common definition). 
%
%OUTPUTS: avgSNR The average power signal to noise ratio needed to achieve
%                the desired probability of detection at the given
%                threshold. If this is an empty matrix, then for a given
%                PD, the threshold is so low that due to the high false
%                alarm rate, the desired PD is impossibly low.
%
%Equation II.10 in [1] gives an expression for PD in terms of thresh and
%the average SNR. Here, we have simply inverted the expression.
%
%However, the expressions used are derived based on a scaled square-law
%detector, as in Equation 10.4-4 in [2], which is the same definition used
%in [1]. Thus, the threshold in this function is scaled appropriately.
%
%EXAMPLE:
%Here, we show that the results are consistent.
% N=4;
% PD=0.6;
% thresh=1e5;
% ampDef=1;
% avgSNR=SwerlingIISqLawD.avgSNR4PDThresh(PD,thresh,N,ampDef);
% PDBack=SwerlingIISqLawD.PD4Threshold(avgSNR,thresh,N,ampDef)
% numSamples=1e5;
% PDSamp=mean(SwerlingIISqLawD.rand([numSamples,1],avgSNR,N,ampDef)>=thresh)
%One will see that PDBack is the same as PD and about the same as PDSamp.
%
%REFERENCES:
%[1] P. Swerling, "Probability of detection for fluctuating targets," The
%    RAND Corporation, Santa Monica, CA, Tech. Rep. RM-1217, 1954.
%[2] J. V. Di Franco and W. L. Rubin, Radar Detection. SciTech Publishing
%    Inc., Rayliegh, NC: 2004.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C. 
    
    if(nargin<4||isempty(ampDef))
        ampDef=1;
    end

    if(nargin<3||isempty(N))
       N=1; 
    end
    
    if(ampDef==0)
        %The division by 2 deals with the scaling in Equation 10.4-4.
        avgSNR=(thresh/2)/gammaincinv(1-PD,N)-1;
    else
        avgSNR=thresh/gammaincinv(1-PD,N)-1;
    end

    if(avgSNR<0)
        %No solution exists.
        avgSNR=[];
    end
end

function avgSNR=avgSNR4PDPFA(PD,PFA,N,ampDef)
%%AVGSNR4PDTHRESH Given a detection probabilty and the probability of false
%       alarm, determine the average power signal to noise ratio (SNR) of
%       the target needed under a  Swerling II model for one or more pulses
%       in a square-law detector.
%
%INPUTS: PD The detection probability of the target , 0<=PD<1.
%       PFA The probability of false alarm, 0<=PFA<1.
%         N The number of pulses that are to be incoherently added for
%           detection (in a square-law detector). In a Swerling II model, a
%           different realization of the target power is used across each
%           pulse and the noise corrupting the pulses varies. If this
%           parameter is omitted or an empty matrix is passed, N=1 is used.
%    ampDef This specified normalization (see help SwerlingIISqLawD).
%           Possible values are:
%          0 The expected value of the squared magnitude of the noise is 2.
%          1 (The default if omitted or an empty matrix is passed) The
%            expected value of the squared magnitude of the noise is 1.
%            (A more common definition). 
%
%OUTPUTS: avgSNR The average power signal to noise ratio needed to achieve
%                the desired probability of detection at the given
%                threshold. If this is an empty matrix, then for a given
%                PD, the threshold is so low that due to the high false
%                alarm rate, the desired PD is impossibly low.
%
%Equation II.10 in [1] gives an expression for PD in terms of a normalized
%threshold and the average SNR. Here, we call the function
%PFA2SquareLawThreshold and insert the result into
%SwerlingIISqLawD.avgSNR4PDThresh, which inverts the expression in [1].
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
    
    thresh=PFA2SquareLawThreshold(PFA,N,ampDef);
    avgSNR=SwerlingIISqLawD.avgSNR4PDThresh(PD,thresh,N,ampDef);
end

function thresh=thresh4AvgSNRPD(PD,avgSNR,N,ampDef)
%%THRESH4AVGSNRPD Given a detection probability and the average power
%           signal to noise ratio (SNR), determine the threshold needed
%           under a Swerling II model for one or more pulses in a
%           square-law detector.
%
%INPUTS: PD The detection probability of the target , 0<=PD<1.  
%    avgSNR A vector or matrix of average power signal to noise ratios
%           of the target at which one wishes to evaluate the threshold.
%         N The number of pulses that are to be incoherently added for
%           detection (in a square-law detector). In a Swerling II model, a
%           different realization of the target power is used across each
%           pulse and the noise corrupting the pulses varies. If this
%           parameter is omitted or an empty matrix is passed, N=1 is used.
%    ampDef This specified normalization (see help SwerlingIISqLawD).
%           Possible values are:
%          0 The expected value of the squared magnitude of the noise is 2.
%          1 (The default if omitted or an empty matrix is passed) The
%            expected value of the squared magnitude of the noise is 1.
%            (A more common definition). 
%
%OUTPUTS: thresh The normnalized detection threshold(s) (assuming the noise
%                variance is 1).
%
%Equation II.10 in [1] gives an expression for PD in terms of a normalized
%threshold and the average SNR. Here, we simply invert the equation.
%
%However, the expressions used are derived based on a scaled square-law
%detector, as in Equation 10.4-4 in [2], which is the same defintiion used
%in [1]. Thus, the threshold in this function is scaled appropriately.
%
%EXAMPLE:
%Here, we show that the results are consistent.
% N=4;
% PD=0.6;
% avgSNR=2;
% ampDef=1;
% thresh=SwerlingIISqLawD.thresh4AvgSNRPD(PD,avgSNR,N,ampDef);
% PDBack=SwerlingIISqLawD.PD4Threshold(avgSNR,thresh,N,ampDef)
% numSamples=1e5;
% PDSamp=mean(SwerlingIISqLawD.rand([numSamples,1],avgSNR,N,ampDef)>=thresh)
%One will see that PDBack=0.6 and PDSamp is about the same.
%
%REFERENCES:
%[1] P. Swerling, "Probability of detection for fluctuating targets," The
%    RAND Corporation, Santa Monica, CA, Tech. Rep. RM-1217, 1954.
%[2] J. V. Di Franco and W. L. Rubin, Radar Detection. SciTech Publishing
%    Inc., Rayliegh, NC: 2004.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.  
    
    if(nargin<4||isempty(ampDef))
        ampDef=1;
    end

    if(nargin<3||isempty(N))
        N=1; 
    end
    
    if(ampDef==0)
        thresh=2*gammaincinv(1-PD,N).*(avgSNR+1);
    else
        thresh=gammaincinv(1-PD,N).*(avgSNR+1);
    end
end

function vals=rand(NDims,avgSNR,N,ampDef)
%%RAND Generate random variables representing a target whose sampled signal
%      to average noise power ratio is generated according to a Swerling II
%      model detected by a square-law detector.
%
%INPUTS: NDims If NDims is a scalar, then rand returns an NDimsXNDims
%          matrix of random variables. If NDims=[M,N1] is a two-element row
%          vector, then rand returns an MXN1 matrix of random variables.
%   avgSNR The average power signal to noise ratio of the target.
%        N The number of pulses that are to be incoherently added for
%          detection (in a square-law detector). In a Swerling II model, a
%          different realization of the target power is used across each
%          pulse and the noise corrupting the pulses varies. If this
%          parameter is omitted or an empty matrix is passed, N=1 is used.
%   ampDef This specified normalization (see help SwerlingIISqLawD).
%          Possible values are:
%          0 The expected value of the squared magnitude of the noise is 2.
%          1 (The default if omitted or an empty matrix is passed) The
%            expected value of the squared magnitude of the noise is 1.
%            (A more common definition). 
%
%OUTPUTS: vals A matrix whose dimensions are determined by N of the
%              generated random variables.
%
%Following the model of [1], a sample from the exponential distribution is
%performed to obtain the complex signal to noise ratio of every single one
%of the N pulses. Then, the square law detector output conditioned on the
%signal to noise ratio is generated directly using NonFlucSqLawD.rand for
%each pulse and the values are added.
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
            %Equation I.1 for  the input signal-to-noise power ratio. It is
            %different for each sample in the pulse train in the Swerling
            %II model.
            curSNR=ExponentialD.rand(1,1/avgSNR);
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
