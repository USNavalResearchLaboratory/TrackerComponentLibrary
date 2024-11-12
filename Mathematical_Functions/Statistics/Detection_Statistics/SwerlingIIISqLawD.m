classdef SwerlingIIISqLawD
%%SWERLINGIIISQLAWD This class implements various functions and detections
%                statistics related to the Swerling III target model when
%                used with a square law detector. All methods are static
%                (the class holds no information and does not need to be
%                instantiated). All values are with respect to a normalized
%                noise model, of which two normalization are available.
%
%Implemented methods are: mean, var, PDF, CDF, PD4Threshold, PD4PFA,
%                         thresh4AvgSNRPD, rand
%
%Swerling models are given in [1] and in Chapter 11 of [2]. The Swerling
%III model assumes that the observed power SNR of the target varies in
%terms of a central gamma distribution with shape parameter k=2 and scale
%parameter theta=avgSNR/2, that in a pulse train for detection the target
%SNR is constant from pulse to pulse, and that a square law detector is
%used to incoherently integrate multiple pulses. The difference from a
%Swerling IV model is the target amplitude does not fluctuate from pulse to
%pulse.
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
%      (normalized to a noise variance of 1) under a Swerling III model in
%      a square-law detector.
%
%INPUTS: avgSNR The average power signal to noise ratio of the target.
%        N The number of pulses that are to be incoherently added for
%          detection (in a square-law detector). In a Swerling III model, a
%          single realization of the target power is used across all
%          pulses, but the noise corrupting the pulses varies. If this
%          parameter is omitted or an empty matrix is passed, N=1 is used.
%   ampDef This specified normalization (see help SwerlingIIISqLawD).
%          Possible values are:
%          0 The expected value of the squared magnitude of the noise is 2.
%          1 (The default if omitted or an empty matrix is passed) The
%            expected value of the squared magnitude of the noise is 1.
%            (A more common definition). 
%
%OUTPUTS: val The value of the mean for the given parameters.
%
%The expression for the mean implemented here is E{x}, where the expected
%value is taken over the PDF in Equation III.16 of [1]. This can be
%evaluated in Mathematica.  However, the first moment from that expression
%has to be multiplied by 2, because as seen in Equation 10.4-4 of [2], the
%definition in [2], which is the same as in [1], is in terms of a scaled
%version of the square law detector.
%
%EXAMPLE:
%Here, we validate the mean by generating random samples and comparing the
%computed mean to the sample mean.
% avgSNR=2;
% N=4;
% ampDef=1;
% numSamples=1e5;
% meanVal=SwerlingIIISqLawD.mean(avgSNR,N,ampDef)
% meanSampVal=mean(SwerlingIIISqLawD.rand([numSamples,1],avgSNR,N,ampDef))
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
        val=2*N*(1+avgSNR);
    else
        val=N*(1+avgSNR);
    end
end

function val=var(avgSNR,N,ampDef)
%%VAR Obtain the variance of the distribution of the detection power
%     (normalized to a noise variance of 1) under a Swerling III model in a
%     square-law detector.
%
%INPUTS:avgSNR The average power signal to noise ratio of the target.
%        N The number of pulses that are to be incoherently added for
%          detection (in a square-law detector). In a Swerling III model, a
%          single realization of the target power is used across all
%          pulses, but the noise corrupting the pulses varies. If this
%          parameter is omitted or an empty matrix is passed, N=1 is used.
%   ampDef This specified normalization (see help SwerlingIIISqLawD).
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
%in Equation III.16 in [1] (simplified with Mathematica). One then uses the
%identity that
%variance=E{x^2}-E{x}^2
%However, the variance from that expression has to be multiplied by 4,
%because as seen in Equation 10.4-4 in [2], the definition in [2, which is
%the same as usedin [2], is in terms of a scaled version of the square law
%detector.
%
%EXAMPLE:
%Here, we validate the variance by generating random samples and comparing
%the computed variance to the sample variance.
% avgSNR=2;
% N=4;
% ampDef=0;
% numSamples=1e5;
% varVal=SwerlingIIISqLawD.var(avgSNR,N,ampDef)
% varSampVal=var(SwerlingIIISqLawD.rand([numSamples,1],avgSNR,N,ampDef))
%One will see that both variance values are about 208.
%
%REFERENCES:
%[1] P. Swerling, "Probability of detection for fluctuating targets," The
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
        val=4*(1/2)*N*(2+avgSNR.*(4+N*avgSNR));
    else
        val=(1/2)*N*(2+avgSNR.*(4+N*avgSNR));
    end
end

function val=PDF(x,avgSNR,N,ampDef)
%%PDF Evaluate the scalar probability density function (PDF) of the
%     distribution of the detection power (normalized to a noise variance
%     of 1) under a Swerling III model in a square-law detector.
%
%INPUTS: x The point or points at which the PDF should be evaluated. x>=0
%          for nonzero PDF values.
%   avgSNR The average power signal to noise ratio of the target.
%        N The number of pulses that are to be incoherently added for
%          detection (in a square-law detector). In a Swerling III model, a
%          single realization of the target power is used across all
%          pulses, but the noise corrupting the pulses varies. If this
%          parameter is omitted or an empty matrix is passed, N=1 is used.
%   ampDef This specified normalization (see help SwerlingIIISqLawD).
%          Possible values are:
%          0 The expected value of the squared magnitude of the noise is 2.
%          1 (The default if omitted or an empty matrix is passed) The
%            expected value of the squared magnitude of the noise is 1.
%            (A more common definition). 
%
%OUTPUTS: val The values of the PDF evaluated at the given points.
%
%The PDF is taken from Equation III.16 of [1]. However, a clearer
%derivation is given in Chapter 11 of [2].
%
%EXAMPLE:
%Here, we validate the PDF by generating random samples and comparing the
%PDF plot with a histogram of the random samples.
% avgSNR=2;
% N=4;
% numSamples=1e4;
% ampDef=0;
% 
% figure(1)
% clf
% histogram(SwerlingIIISqLawD.rand([numSamples,1],avgSNR,N,ampDef),'Normalization','pdf')
% hold on
% numPoints=1000;
% x=linspace(0,80,numPoints);
% vals=SwerlingIIISqLawD.PDF(x,avgSNR,N,ampDef);
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

    %Initial change of variable to undo the scaling from Equation 10.4-4 in
    %[2].
    if(ampDef==0)
        x=x/2;
    end

    numEls=numel(x);
    val=zeros(size(x));
    for curEl=1:numEls
        v=x(curEl);
        if(v>0)
            val(curEl)=exp((N-1)*log(v)-v-gammaln(N)-2*log(1+N*avgSNR/2))*hypergeometric1F1(2,N,v/(1+2/(N*avgSNR)));
        elseif(v==0&&N==1)
            val(curEl)=(1/(factorial(N-1)*(1+N*avgSNR/2)^2));
        end
    end
    
    %Adjust for the change of variable to undo the scaling from Equation
    %10.4-4 in [2].
    if(ampDef==0)
        val=val/2;
    end
end

function val=CDF(v,avgSNR,N,ampDef)
%%CDF Evaluate the scalar cumulative distribution function (CDF) of the
%     distribution of the detection power (normalized to a noise variance
%     of 1) under a Swerling III model in a square-law detector.
%
%INPUTS: v The point or points at which the CDF should be evaluated. v>0
%          for nonzero CDF values.
%   avgSNR The average power signal to noise ratio of the target.
%        N The number of pulses that are to be incoherently added for
%          detection (in a square-law detector). In a Swerling III model, a
%          single realization of the target power is used across all
%          pulses, but the noise corrupting the pulses varies. If this
%          parameter is omitted or an empty matrix is passed, N=1 is used.
%   ampDef This specified normalization (see help SwerlingIIISqLawD).
%          Possible values are:
%          0 The expected value of the squared magnitude of the noise is 2.
%          1 (The default if omitted or an empty matrix is passed) The
%            expected value of the squared magnitude of the noise is 1.
%            (A more common definition). 
%
%OUTPUTS: val The values of the CDF evaluated at the given points.
%
%The derivation of the detection probability in [1] is as 1-the value of
%the CDF. Thus, this function just evaluates 
%1-SwerlingIIISqLawD.PD4Threshold(avgSNR,v,N,ampDef).
%
%EXAMPLE:
%Here, the CDF returned by this function is plotted along with an empirical
%CDF constructed from random samples of this distribution. They agree well.
% avgSNR=5;
% N=4;
% ampDef=0;
% numSamples=1e5;
% samples=SwerlingIIISqLawD.rand([numSamples,1],avgSNR,N,ampDef);
% numPts=1000;
% x=linspace(0,80,numPts);
% CDFEmp=EmpiricalD.CDF(x,samples);
% CDF=SwerlingIIISqLawD.CDF(x,avgSNR,N,ampDef);
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

    val=1-SwerlingIIISqLawD.PD4Threshold(avgSNR,v,N,[],ampDef);
    val(v<=0)=0;
end

function PD=PD4Threshold(avgSNR,thresh,N,method,ampDef)
%%PD4THRESHOLD Determine the detection probability of a Swerling III target
%           given it's signal to noise ratio and the value of the
%           detection threshold, assuming that a square-law detector is
%           used.
%
%INPUTS: avgSNR A vector or matrix of average power signal to noise ratios
%          of the target at which one wishes to evaluate the detection
%          probability.
%   thresh The scalar normalized detection threshold to use. This is the
%          threshold to use if the noise variance is 1.
%        N The number of pulses that are to be incoherently added for
%          detection (in a square-law detector). In a Swerling III model, a
%          single realization of the target power is used across all
%          pulses, but the noise corrupting the pulses varies. If this
%          parameter is omitted or an empty matrix is passed, N=1 is used.
%   method This specified whether an exact method or an approximation
%          should be used. For N=1 and N=2, the exact method equals the
%          approximation. Possible values are
%          0 (The default if omitted or an empty matrix is passed) Use the
%            exact solution. This is obtained by integrating Equation
%            III.16 in [1], which is the saem as Equation 11.4-18 in p2].
%            An explicit form in terms of incomplete gamma functions was
%            obtained by evaluating the integral in Mathematica and has
%            been implemented here.
%          1 Use the approximate solution of Equation II.14 of [1]. This is
%            only valid if N*avgSNR>2 or if N=1,2. This is an exact
%            solution if N=1,2.
%   ampDef This specified normalization (see help SwerlingIIISqLawD).
%          Possible values are:
%          0 The expected value of the squared magnitude of the noise is 2.
%          1 (The default if omitted or an empty matrix is passed) The
%            expected value of the squared magnitude of the noise is 1.
%            (A more common definition). 
%
%OUTPUTS: PD The detection probability of the target.
%
%In [1] and [2], only approximate expressions for the detection probability
%are given. However, an exact value was obtained here by integrating
%Equation II.16 in [1], which is equivalent to Equation 11.4-21 in [2].
%Mathematica was used for the integration. 
%
%However, the expressions used are derived based on a scaled square-law
%detector, as in Equation 10.4-4 in [2]. Thus, the threshold in this
%function is scaled appropriately.
%
%EXAMPLE:
%Here, we validate the results by comparing the PD from this function to
%the PD computed using random samples.
% avgSNR=2;
% thresh=30;
% ampDef=0;
% N=4;
% PD=SwerlingIIISqLawD.PD4Threshold(avgSNR,thresh,N,[],ampDef)
% numSamples=1e5;
% PDSamp=mean(SwerlingIIISqLawD.rand([numSamples,1],avgSNR,N,ampDef)>=thresh)
%One will see that both PD values are about 0.272.
%
%REFERENCES:
%[1] P. Swerling, "Probability of detection for fluctuating targets," The
%    RAND Corporation, Santa Monica, CA, Tech. Rep. RM-1217, 1954.
%[2] J. V. Di Franco and W. L. Rubin, Radar Detection. SciTech Publishing
%    Inc., Rayliegh, NC: 2004.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.  

    if(nargin<5||isempty(ampDef))
        ampDef=1;
    end

    if(nargin<4||isempty(method))
       method=0; 
    end

    if(nargin<3||isempty(N))
       N=1; 
    end

    if(ampDef==0)
        thresh=thresh/2;
    end
    
    switch(method)
        case 0%The exact solution
            xBar=avgSNR;
            Y=thresh;

            % term1=exp(-((2*Y)/(2+N*xBar)))*1/(N*xBar)^(N-1)*(2+N*xBar)^(-3+N)*(8+N*(-4+xBar*(N*(-2+xBar)+2*(3+Y))));
            %Using the exponent of the sums of logarithms, helps avoid overflow
            %problems.
            term1=exp(-((2*Y)/(2+N*xBar))-(N-1)*log(N*xBar)+(N-3)*log(2+N*xBar)).*(8+N*(-4+xBar*(N*(-2+xBar)+2*(3+Y))));

            term2=gammainc(Y,N,'upper');
            %term3=-(exp(-Y)*(4+N*(xBar-2))*Y^(N-1))/((2+N*xBar)*gamma(N));
            %Using the exponent of the sums of logarithms, helps avoid overflow
            %problems.
            term3=-(4+N*(xBar-2))*exp(-Y+(N-1)*log(Y)-log((2+N*xBar))-gammaln(N));
            
            %This can be numerically unstable for very small avgSNR values.
            PD=term1.*(1-gammainc(((N*xBar*Y)/(2+N*xBar)),N-1,'upper'))+(term2+term3);
            
            if(any(PD>1|PD<0|~isfinite(PD)))%||any(2^5*eps(abs(term1))>abs(term2)))
                warning('Finite precision errors detected.')
            end
            
            %Deal with finite precision errors.
            PD=min(PD,1);
            PD=max(PD,0);
        case 1%The approximate solution.
            Nx=N*avgSNR;
            
            if(any(Nx<2))
                warning('The approximate solution is not valid for N*avgSNR<=2')
            end
            PD=(1+2/Nx)^(N-2)*(1+thresh/(1+Nx/2)-2*(N-2)/Nx)*exp(-thresh/(1+Nx/2));
        otherwise
            error('Unknown method specified.')
    end
end

function PD=PD4PFA(avgSNR,PFA,N,method,ampDef)
%%PD4PFA Determine the detection probability of a Swerling III target given
%        its signal to noise ratio and the probability of false alarm
%        implied by a particular unspecified detection threshold, assuming
%        that a square law detector is used.
%
%INPUTS: avgSNR A vector or matrix of average power signal to noise ratios
%          of the target at which one wishes to evaluate the detection
%          probability.
%      PFA The probability of false alarm, 0<=PFA<1.
%        N The number of pulses that are to be incoherently added for
%          detection (in a square-law detector). In a Swerling III model, a
%          single realization of the target power is used across all
%          pulses, but the noise corrupting the pulses varies. If this
%          parameter is omitted or an empty matrix is passed, N=1 is used.
%   method This specified whether an exact method or an approximation
%          should be used. For N=1 and N=2, the exact method equals the
%          approximation. Possible values are
%          0 (The default if omitted or an empty matrix is passed) Use the
%            exact solution. This is obtained by integrating Equation
%            III.16 in [1], which is the saem as Equation 11.4-18 in p2].
%            An explicit form in terms of incomplete gamma functions was
%            obtained by evaluating the integral in Mathematica and has
%            been implemented here.
%          1 Use the approximate solution of Equation II.14 of [1]. This is
%            only valid if N*avgSNR>2 or if N=1,2. This is an exact
%            solution if N=1,2.
%   ampDef This specified normalization (see help SwerlingIIISqLawD).
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
%SwerlingIIISqLawD.PD4Threshold.
%
%REFERENCES:
%[1] P. Swerling, "Probability of detection for fluctuating targets," The
%    RAND Corporation, Santa Monica, CA, Tech. Rep. RM-1217, 1954.
%[2] J. V. Di Franco and W. L. Rubin, Radar Detection. SciTech Publishing
%    Inc., Rayliegh, NC: 2004.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.    

    if(nargin<5||isempty(ampDef))
        ampDef=1;
    end

    if(nargin<4||isempty(method))
        method=0; 
    end

    if(nargin<3||isempty(N))
       N=1; 
    end

    thresh=PFA2SquareLawThreshold(PFA,N,ampDef);
    PD=SwerlingIIISqLawD.PD4Threshold(avgSNR,thresh,N,method,ampDef);
end

function thresh=thresh4AvgSNRPD(PD,avgSNR,N,ampDef)
%%THRESH4AVGSNRPD Given a detection probability and the average power
%           signal to noise ratio (SNR), determine the threshold needed
%           under a Swerling III model. This solution is only available for
%           a N=1 and N=2 pulses in a square law detector.
%
%INPUTS: PD The detection probability of the target , 0<=PD<1.  
%    avgSNR A vector or matrix of average power signal to noise ratios
%           of the target at which one wishes to evaluate the threshold.   
%         N The number of pulses that are to be incoherently added for
%           detection (in a square-law detector). In a Swerling III
%           model, the same realization of the target power is used
%           across each pulse and the noise corrupting the pulses
%           varies. If this parameter is omitted or an empty matrix is
%           passed, N=1 is used. This function only supports N=1 and N=2.
%   ampDef This specified normalization (see help SwerlingIIISqLawD).
%          Possible values are:
%          0 The expected value of the squared magnitude of the noise is 2.
%          1 (The default if omitted or an empty matrix is passed) The
%            expected value of the squared magnitude of the noise is 1.
%            (A more common definition). 
%
%OUTPUTS: thresh The normnalized detection threshold(s) (assuming the noise
%                variance is 1).
%
%Equation II.14 in [1] gives an expression for PD in terms of a normalized
%threshold and the average SNR for N=1,2 pulses. Here, we invert the
%equation in terms of the Lambert-W function. The -1 branch of the
%Lambert-W function is chosen as it assures positive thresholds.
%
%However, the expressions used are derived based on a scaled square-law
%detector, as in Equation 10.4-4 in [2], which is the same model as in [1].
%Thus, the threshold in this function is scaled appropriately.
%
%EXAMPLE:
%Here, we show that the results are consistent.
% N=4;
% PD=0.6;
% avgSNR=2;
% ampDef=0;
% thresh=SwerlingIIISqLawD.thresh4AvgSNRPD(PD,avgSNR,N,ampDef);
% PDBack=SwerlingIIISqLawD.PD4Threshold(avgSNR,thresh,N,[],ampDef)
% numSamples=1e5;
% PDSamp=mean(SwerlingIIISqLawD.rand([numSamples,1],avgSNR,N,ampDef)>=thresh)
%One will see that PDBack=0.5 and PDSamp is about the same.
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

    xBar=avgSNR;
    
    LambWArg=-PD*exp(-1+(2*(N-2))./(N*xBar)).*((N*xBar)./(2+N*xBar)).^(N-2);
    thresh=-((2+N*xBar)./(2*N*xBar)).*(4+N*(xBar-2)+N*xBar.*LambW(LambWArg,-1));
    
    %The multiplication by 2 deals with the scaling in Equation 10.4-4.
    if(ampDef==0)
        thresh=thresh*2;
    end
end

function vals=rand(NDims,avgSNR,N,ampDef)
%%RAND Generate random variables representing a target whose sampled signal
%      to average noise power ratio is generated according to a Swerling
%      III model detected by a square-law detector.
%
%INPUTS: NDims If NDims is a scalar, then rand returns an NDimsXNDims
%          matrix of random variables. If NDims=[M,N1] is a two-element row
%          vector, then rand returns an MXN1 matrix of random variables.
%   avgSNR The average power signal to noise ratio of the target.
%        N The number of pulses that are to be incoherently added for
%          detection (in a square-law detector). In a Swerling III model,
%          the same realization of the target power is used across each
%          pulse and the noise corrupting the pulses varies. If this
%          parameter is omitted or an empty matrix is passed, N=1 is used.
%   ampDef This specified normalization (see help SwerlingIIISqLawD).
%          Possible values are:
%          0 The expected value of the squared magnitude of the noise is 2.
%          1 (The default if omitted or an empty matrix is passed) The
%            expected value of the squared magnitude of the noise is 1.
%            (A more common definition). 
%
%OUTPUTS: vals A matrix whose dimensions are determined by N of the
%              generated random variables.
%
%Following the model of [1], a single sample from the central gamma
%distribution is performed to obtain the complex signal to noise ratio of
%the target. Then, the square law detector output conditioned on the signal
%to noise ratio is generated directly using NonFlucSqLawD.rand.
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
        %Equation I.2 for  the input signal-to-noise power ratio. It will
        %be the same for all N samples that are accumulated.
        curSNR=GammaD.rand(1,2,avgSNR/2);
        vals(curEl)=NonFlucSqLawD.rand(1,curSNR,N,ampDef);
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
