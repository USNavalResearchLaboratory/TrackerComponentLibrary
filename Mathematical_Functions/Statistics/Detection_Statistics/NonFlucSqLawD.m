classdef NonFlucSqLawD
%%NONFLUCSQLAWD This class implements various functions and statistics
%   related to a square law detector applied to a non-fluctuating target.
%   The square law detector for N complex samples is
%   y=sum_{i=1}^N abs(r_i)^2
%   where r_i is the ith complex sampled value (including a signal and
%   noise as described below) and y can be compared to a threshold. All
%   functions have been normalized to a common noise power. See below for
%   the two definitions of normalization that are used. All methods are
%   static (the class holds no information and does not need to be
%   instantiated).
%
%Implemented methods are: mean, var, PDF, CDF, avgSNR4PDThresh,
%                         PD4Threshold, PD4PFA, rand
%
%This model considers detection of a non-fluctuating target with a constant
%signal power to noise ratio via non-coherent integration of multiple
%samples in a square-law detector. The model is discussed in Chapter 10 of
%[1], but significant confusion can arise on the definition of the
%normalization that is used. Here, we go over the two common definitions of
%normalization, which are options in this class.  
%
%Definition 1:
%In this instance, we assume that things have been normalized so that for a
%single sample,
% zSignal=sqrt(A)*exp(1j*2*pi*rand(1)) %The signal has a random phase.
% zNoise=(randn(1)+1j*randn(1))/sqrt(2)
%Note that zNoise is a sample of a circular complex Gaussian distribution.
%with unit variance. So E{zNoise*conj(zNoise)}=1. Similarly,
%E{zSignal*conj(zSignal)}=A. Thus, A is the signal to noise ratio (SNR) of
%the signal. The squared measurement would be 
%y=abs(r)^2=abs(zSignal+zNoise)^2=(zSignal+zNoise)*conj(zSignal+zNoise)
%Note that the random phase of zSignal does not change the statistics of y
%and one can just replace zSignal with sqrt(A). Omitting that phase and
%rewriting zNoise=(randn(1)+1j*randn(1))/sqrt(2)=(x+1j*y)/sqrt(2) where x
%and y are standard normal random variables, y can be written
%y=(sqrt(A)+(x+1j*y)/sqrt(2))*(sqrt(A)+(x-1j*y)/sqrt(2))
% =(1/2)*((sqrt(2*A)+x)^2+y^2)
%The quantity within the parentheses is distributed noncentral chi squared
%with 2 degrees of freedom and noncentrality parameter (2*A). Thus, y is
%a transformed noncentral chi squared distribution. If X2(x,nu,lambda) is
%the PDF on a noncentral chi-squared distribution with nu degrees of
%freedom and noncentrality parameter lambda, the PDF of y, y(x) is
%2*X2(x*2,2,2*A).
%Similarly, when considering the sum of N pulses, the PDF of y is
%2*X2(x*2,2*N,N*2*A)
%Note that when considering N pulses, it does not matter whether the signal
%phase is constant or changes between pulses. A simple Monte Carlo
%simulation where one can see that a histogram of the signals mentioned
%agrees with the transformed PDF as described is:
% numSamples=1e5;
% A=3;
% N=4;
% z=sqrt(A)*exp(1j*2*pi*randn(N,numSamples))+(randn(N,numSamples)+1j*randn(N,numSamples))/sqrt(2);
% r=sum(abs(z).^2,1);
% figure(1)
% clf; hold on
% histogram(r,'Normalization','pdf')
% numPts=1000;
% x=linspace(0,50,numPts);
% PDF=ChiSquareD.PDF(x*2,2*N,N*2*A)*2;
% plot(x,PDF,'linewidth',2)
%
%Definition 0:
%The motivation for this definition of the normalization is so that the
%results can be expressed in terms of a chi-squared distribution that has
%not been transformed. Here, the signal and noise models are
% zSignal=sqrt(2*A)*exp(1j*2*pi*rand(1))
% zNoise=randn(1)+1j*randn(1)
%In this instance E{abs(zNoise)^2}=2 and E{abs(zSignal)^2}=2*A, so the
%new normalization has doubled everything, but the SNR is still just A.
%Again, the phase of zSignal does not matter, and if we rewrite zNoise as
%x+1j*y where x and y are standard normal random variables, then for a
%single pulse
%y=abs(zSignal+zNoise)^2=(sqrt(2*A)+x+1j*y)*(sqrt(2*A)+x-1j*y)
% =(sqrt(2*A)+x)^2+y^2
%which is distributed noncentral chi-squared with 2 degrees of freedom and
%noncentrality parameter 2*A. Similarly, when considering N pulses, y is
%distributed noncentral chi-squared with 2*N degrees of freedom and
%noncentrality parameter 2*N*A. A simple Monte Carlo
%simulation where one can see that a histogram of the signals mentioned
%agrees with the transformed PDF as described is:
% numSamples=1e5;
% A=3;
% N=4;
% z=sqrt(2*A)+(randn(N,numSamples)+1j*randn(N,numSamples));
% r=sum(abs(z).^2,1);
% figure(1)
% clf
% hold on
% histogram(r,'Normalization','pdf')
% numPts=1000;
% x=linspace(0,100,numPts);
% PDF=ChiSquareD.PDF(x,2*N,N*(2*A));
% plot(x,PDF,'linewidth',2)
%
%So, depending on the definition of the normalization selected, random
%values, thresholds and other parameters of this class can vary.
%
%Suppose that one is given a non-normalized sample x=zSig+zNoise where
%E{abs(zNoise)^2}=B, then to scale the sample to be consistent with
%Definition 1 of this function, one should use xScaled=x*(1/sqrt(B)).
%Similarly, for Definition 0 to hold, it should be scaled as xScaled
%=x*(sqrt(2)/sqrt(B)).
%
%Chapter 10.4 of [1] discusses the PDF of the output of the square law,
%detector, such as in Equation 10.4-13. One must pay attention to how
%things are normalized.
%
%REFERENCES:
%[1] J. V. Di Franco and W. L. Rubin, Radar Detection. SciTech Publishing
%    Inc., Rayliegh, NC: 2004.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.
    
methods(Static)

function val=mean(avgSNR,N,ampDef)
%%MEAN Obtain the mean of the distribution of the detection power
%     (normalized to a noise variance of 1) under a constant amplitude
%     target model in a square-law detector. 
%
%INPUTS: avgSNR The average power signal to noise ratio of the target.
%        N The number of pulses that are to be incoherently added for
%          detection (in a square-law detector). In a nonfluctuating model,
%          the target power is constant and the noise corrupting the pulses
%          varies. If this parameter is omitted or an empty matrix is
%          passed, N=1 is used.
%   ampDef This specified normalization (see help NonFlucSqLawD). Possible
%          values are:
%          0 The expected value of the squared magnitude of the noise is 2.
%          1 (The default if omitted or an empty matrix is passed) The
%            expected value of the squared magnitude of the noise is 1.
%            (A more common definition). 
%
%OUTPUTS: val The value of the mean for the given parameters.
%
%The PDF is given by Equation 10.4-13 in [1]. However, for ampDef=0, a
%change of variables is performed to account for the values having been
%defined in terms of a scaled version of the square law in Equation 10.4-4.
%The final result is a noncentral chi-square distribution and this function
%just calls ChiSquareD.mean with the appropriate parameters.
%
%EXAMPLE:
%Here, we validate the mean by generating random samples and comparing the
%computed mean to the sample mean. 
% avgSNR=3;
% N=4;
% numSamples=1e5;
% ampDef=1;
% meanVal=NonFlucSqLawD.mean(avgSNR,N,ampDef)
% meanSampVal=mean(NonFlucSqLawD.rand([numSamples,1],avgSNR,N,ampDef))
%One will see that both values are about 16.
%
%REFERENCES:
%[1] J. V. Di Franco and W. L. Rubin, Radar Detection. SciTech Publishing
%    Inc., Rayliegh, NC: 2004.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.

    if(nargin<3||isempty(ampDef))
        ampDef=1;
    end

    if(nargin<2||isempty(N))
        N=1;
    end
    
    %Equation 10.4-13 in [1] is a noncentral Chi-Square distribution.
    nu=2*N;
    lambda=2*N*avgSNR;
    if(ampDef==0)
        val=ChiSquareD.mean(nu,lambda);
    else
        val=ChiSquareD.mean(nu,lambda)/2;
    end
end

function avgSNR=mean2AvgSNR(meanVal,N,ampDef)
%%MEAN2AVGSNR Given the mean value the distribution of the detection power
%     (normalized to a noise variance of 1), determine the average SNR.
%     This is the inverse of NonFlucSqLawS.mean.
%
%INPUTS: meanVal The mean of the distribution (normalized to unit power
%          based on the amplitude definition in ampDef).
%        N The number of pulses that are to be incoherently added for
%          detection (in a square-law detector). In a nonfluctuating model,
%          the target power is constant and the noise corrupting the pulses
%          varies. If this parameter is omitted or an empty matrix is
%          passed, N=1 is used.
%   ampDef This specified normalization (see help NonFlucSqLawD). Possible
%          values are:
%          0 The expected value of the squared magnitude of the noise is 2.
%          1 (The default if omitted or an empty matrix is passed) The
%            expected value of the squared magnitude of the noise is 1.
%            (A more common definition). 
%
%OUTPUTS: avgSNR The average power signal to noise ratio of the target.
%
%The PDF is given by Equation 10.4-13 in [1] (noncentral chi squared). This
%just solves for the SNR from the mean of such a distribution.
%
%EXAMPLE:
%This gets the mean from an SNR and then uses this method to recover the
%average SNR
% avgSNR=6;
% N=4;
% ampDef=1;
% meanVal=NonFlucSqLawD.mean(avgSNR,N,ampDef);
% avgSNRBack=NonFlucSqLawD.mean2AvgSNR(meanVal,N,ampDef)
%
%REFERENCES:
%[1] J. V. Di Franco and W. L. Rubin, Radar Detection. SciTech Publishing
%    Inc., Rayliegh, NC: 2004.
%
%August 2024 David F. Crouse, Naval Research Laboratory, Washington D.C.
    
    if(nargin<3||isempty(ampDef))
        ampDef=1;
    end

    if(nargin<2||isempty(N))
        N=1;
    end

    nu=2*N;
    if(ampDef==0)
        avgSNR=(meanVal-nu)/(2*N);
    else
        avgSNR=(2*meanVal-nu)/(2*N);
    end
end

function val=var(avgSNR,N,ampDef)
%%VAR Obtain the variance of the distribution of the detection power
%     (normalized to a noise variance of 1) under a constant amplitude
%     target model in a square-law detector.
%
%INPUTS: avgSNR The average power signal to noise ratio of the target.
%        N The number of pulses that are to be incoherently added for
%          detection (in a square-law detector). In a nonfluctuating model,
%          the target power is constant and the noise corrupting the pulses
%          varies. If this parameter is omitted or an empty matrix is
%          passed, N=1 is used.
%   ampDef This specified normalization (see help NonFlucSqLawD). Possible
%          values are:
%          0 The expected value of the squared magnitude of the noise is 2.
%          1 (The default if omitted or an empty matrix is passed) The
%            expected value of the squared magnitude of the noise is 1.
%            (A more common definition).  
%
%OUTPUTS: val The value of the variance for the given parameters.
%
%The PDF is given by Equation 10.4-13 in [1]. However, for ampDef=0, a
%change of variables is performed to account for the values having been
%defined in terms of a scaled version of the square law in Equation 10.4-4.
%The final result is a noncentral chi-square distribution and this function
%just calls ChiSquareD.var with the appropriate parameters.
%
%EXAMPLE:
%Here, we validate the variance by generating random samples and comparing
%the computed variance to the sample variance. 
% avgSNR=2;
% N=4;
% numSamples=1e5;
% ampDef=0;
% varVal=NonFlucSqLawD.var(avgSNR,N,ampDef)
% varSampVal=var(NonFlucSqLawD.rand([numSamples,1],avgSNR,N,ampDef))
%One will see that both values are about 80.
%
%REFERENCES:
%[1] J. V. Di Franco and W. L. Rubin, Radar Detection. SciTech Publishing
%    Inc., Rayliegh, NC: 2004.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
    
    if(nargin<3||isempty(ampDef))
        ampDef=1;
    end

    if(nargin<2||isempty(N))
        N=1;
    end
    
    %Equation 10.4-13 in [1] is a noncentral Chi-Square distribution.
    nu=2*N;
    lambda=2*N*avgSNR;
    if(ampDef==0)
        val=ChiSquareD.var(nu,lambda);
    else
        val=ChiSquareD.var(nu,lambda)/4;
    end
end

function val=PDF(x,avgSNR,N,ampDef)
%%PDF Evaluate the scalar probability density function (PDF) of the
%     distribution of the detection power (normalized to a noise variance
%     of 1) under a constant-amplitude model in a square-law detector.
%
%INPUTS: x The point or points at which the PDF should be evaluated. x>=0
%          for nonzero PDF values.
%   avgSNR The average power signal to noise ratio of the target.
%        N The number of pulses that are to be incoherently added for
%          detection (in a square-law detector). In a nonfluctuating model,
%          the target power is constant and the noise corrupting the pulses
%          varies. If this parameter is omitted or an empty matrix is
%          passed, N=1 is used.
%   ampDef This specified normalization (see help NonFlucSqLawD). Possible
%          values are:
%          0 The expected value of the squared magnitude of the noise is 2.
%          1 (The default if omitted or an empty matrix is passed) The
%            expected value of the squared magnitude of the noise is 1.
%            (A more common definition). 
%
%OUTPUTS: val The values of the PDF evaluated at the given points.
%
%The PDF is given by Equation 10.4-13 in [1]. However, for ampDef=0, a
%change of variables is performed to account for the values having been
%defined in terms of a scaled version of the square law in Equation 10.4-4.
%The final result is a noncentral chi-square distribution and this function
%just calls ChiSquareD.PDF with the appropriate parameters.
%
%EXAMPLE:
%Here, we validate the PDF by generating random samples and comparing the
%PDF plot with a histogram of the random samples.
% avgSNR=3;
% N=4;
% numSamples=1e5;
% ampDef=1;
% figure(1)
% clf
% hold on
% histogram(NonFlucSqLawD.rand([numSamples,1],avgSNR,N,ampDef,0),'Normalization','pdf')
% numPoints=1000;
% x=linspace(0,60,numPoints);
% vals=NonFlucSqLawD.PDF(x,avgSNR,N,ampDef);
% plot(x,vals,'linewidth',2)
% h1=xlabel('x');
% h2=ylabel('PDF(x)');
% set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')
%One will see a good match between the histogram and the PDF.
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

    %Equation 10.4-13 in [1] is a noncentral Chi-Square distribution.
    nu=2*N;
    lambda=2*N*avgSNR;
    if(ampDef==0)
        val=ChiSquareD.PDF(x,nu,lambda);
    else
        val=ChiSquareD.PDF(x*2,nu,lambda)*2;
    end
end

function val=CDF(x,avgSNR,N,ampDef)
%%CDF Evaluate the scalar cumulative distribution function (CDF) of the
%     distribution of the detection power (normalized to a noise variance
%     of 1) under a constant-amplitude model in a square-law detector.
%
%INPUTS: v The point or points at which the CDF should be evaluated. v>0
%          for nonzero CDF values.
%   avgSNR The average power signal to noise ratio of the target.
%        N The number of pulses that are to be incoherently added for
%          detection (in a square-law detector). In a nonfluctuating model,
%          the target power is constant and the noise corrupting the pulses
%          varies. If this parameter is omitted or an empty matrix is
%          passed, N=1 is used.
%   ampDef This specified normalization (see help NonFlucSqLawD). Possible
%          values are:
%          0 The expected value of the squared magnitude of the noise is 2.
%          1 (The default if omitted or an empty matrix is passed) The
%            expected value of the squared magnitude of the noise is 1.
%            (A more common definition).  
%
%OUTPUTS: val The values of the CDF evaluated at the given points.
%
%The derivation of the detection probability is in Chapter 10 of [1].
%Modifications here are related to the normalization selected by ampDef.
%
%EXAMPLE 1:
%This example shows that the CDF returned by this function matches an
%empirical CDF created using random samples.
% avgSNR=2;
% N=4;
% ampDef=0;
% numSamples=1e5;
% samples=NonFlucSqLawD.rand([numSamples,1],avgSNR,N,ampDef);
% numPts=1000;
% x=linspace(0,60,numPts);
% CDFEmp=EmpiricalD.CDF(x,samples);
% CDF=NonFlucSqLawD.CDF(x,avgSNR,N,ampDef);
% figure(1)
% clf
% hold on
% plot(x,CDFEmp,'linewidth',4);
% plot(x,CDF,'linewidth',2);
% legend('Empirical CDF','Theoretical CDF')
%
%EXAMPLE 2:
%This example demonstrates that the CDF value is one minus the probability
%of detection returned by NonFlucSqLawD.PD4Threshold with x as the
%threshold.
% avgSNR=2;
% N=4;
% ampDef=0;
% numPts=1000;
% x=linspace(0,60,numPts);
% CDF=NonFlucSqLawD.CDF(x,avgSNR,N,ampDef);
% PD=zeros(1,numPts);
% for k=1:numPts
%     PD(k)=NonFlucSqLawD.PD4Threshold(avgSNR,x(k),N,ampDef);
% end
% figure(1)
% clf
% hold on
% plot(x,CDF,'linewidth',4);
% plot(x,1-PD,'linewidth',2);
% legend('CDF','1-PD')
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

    %Equation 10.4-13 in [1] is a noncentral Chi-Square distribution.
    nu=2*N;
    lambda=2*N*avgSNR;
    if(ampDef==0)
        val=ChiSquareD.CDF(x,nu,lambda);
    else
        val=ChiSquareD.CDF(x*2,nu,lambda);
    end
end

function [avgSNR,exitCode]=avgSNR4PDThresh(PD,thresh,N,ampDef,convergParams)
%%AVGSNR4PDTHRESH Given a detection probability and the detection
% threshold, determine the power signal to noise ratio (SNR) under a
% non-fluctuating target model.
%
%INPUTS: PD The detection probability of the target , 0<=PD<1.
%   thresh The detection threshold. thresh>0
%        N The number of pulses that are to be incoherently added for
%          detection (in a square-law detector). If this parameter is
%          omitted or an empty matrix is passed, N=1 is used.
%   ampDef This specified normalization (see help NonFlucSqLawD). Possible
%          values are:
%          0 The expected value of the squared magnitude of the noise is 2.
%          1 (The default if omitted or an empty matrix is passed) The
%            expected value of the squared magnitude of the noise is 1.
%            (A more common definition).  
% convergParams An optional structure holding parameters that define how
%          the algorithm converges. Possible members are:
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
%This detection probability is computed in terms of the MarcumQ
%function. Equations 10.4-30 in [1], expresses the detection probability
%in terms of an incomplete Toronto function. In [2], it is shown that the
%incomplete Toronto function parameterized T_B(m,(m-1)/2,r) is expressed in
%terms of the MarcumQ function as 1-MarcumQ((m+1)/2,r*sqrt(2),B*sqrt(2)).
%
%EXAMPLE:
%This example demonstrates that the avgSNR computed using this method can
%be plugged into the PD4Threshold to get the same PD back as was used in
%this function (so the results are consistent). The relative error is on
%the order of finite precision limitations.
% thresh=12;
% N=4;
% ampDef=1;
% PD=0.7;
% avgSNR=NonFlucSqLawD.avgSNR4PDThresh(PD,thresh,N,ampDef);
% PDBack=NonFlucSqLawD.PD4Threshold(avgSNR,thresh,N,ampDef);
% RelErr=(PD-PDBack)/PD
%
%REFERENCES:
%[1] J. V. Di Franco and W. L. Rubin, Radar Detection. SciTech Publishing
%    Inc., Rayliegh, NC: 2004.
%[2] P. C. Sofotasios and S. Freear, "New analytic results for the
%    incomplete Toronto function and incomplete Lipschitz-Hankel
%    integrals," in Proceedings from the SMBO/IEEE MTT-S International
%    Microwave and Optoelectronics Conference, Natal, Brazil, 29 Oct. - 1
%    Nov. 2011, pp. 44-47.
%
%December 2024 David F. Crouse, Naval Research Laboratory, Washington D.C.

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

if(ampDef==1)
    thresh=2*thresh;
end
sthresh=sqrt(thresh);

%We need to find an upper bound. We start with an estimate of 100 and then
%keep doubling it until we have found an upper bound.
avgSNRUpper=100;
lambda=N*2*avgSNRUpper;
PDComp=MarcumQ(N,sqrt(lambda),sthresh);
curIter=0;
while(PDComp<PD)
    curIter=curIter+1;
    if(curIter>maxIterSearch)
        error('Unable to bracket a solution.')
    end
    avgSNRUpper=2*avgSNRUpper;
    lambda=N*2*avgSNRUpper;
    PDComp=MarcumQ(N,sqrt(lambda),sthresh);
end

f=@(avgSNR)(PD-MarcumQ(N,sqrt(N*2*avgSNR),sthresh));
[avgSNR,~,exitCode]=bisectionRootFind(f,[0;avgSNRUpper],XTol,maxIter);

end

function PD=PD4Threshold(avgSNR,thresh,N,ampDef)
%%PD4THRESHOLD Determine the detection probability of a target with a
%           constant signal to noise ratio (SNR) given the SNR and the
%           value of the detection threshold, assuming that a square-law
%           detector is used. 
%
%INPUTS: avgSNR A vector or matrix of average power signal to noise ratios
%          of the target at which one wishes to evaluate the detection
%          probability.
%   thresh The scalar normalized detection threshold to use. This is
%          the threshold to use if the noise variance is 1.
%        N The number of pulses that are to be incoherently added for
%          detection (in a square-law detector). In a nonfluctuating model,
%          the target power is constant and the noise corrupting the pulses
%          varies. If this parameter is omitted or an empty matrix is
%          passed, N=1 is used.
%   ampDef This specified normalization (see help NonFlucSqLawD). Possible
%          values are:
%          0 The expected value of the squared magnitude of the noise is 2.
%          1 (The default if omitted or an empty matrix is passed) The
%            expected value of the squared magnitude of the noise is 1.
%            (A more common definition). 
%
%OUTPUTS: PD The detection probability of the target.
%
%This function implements the detection probability in terms of the MarcumQ
%function. Equations 10.4-30 in [1], expresses the detection probability
%in terms of an incomplete Toronto function. In [2], it is shown that the
%incomplete Toronto function parameterized T_B(m,(m-1)/2,r) is expressed in
%terms of the MarcumQ function as 1-MarcumQ((m+1)/2,r*sqrt(2),B*sqrt(2)).
%This is the form in this problem. Thus, the MarcumQ function is used to
%evaluate the detection probability.
%
%EXAMPLE:
%Here, we validate the results by comparing the PD from this function to
%the PD computed using random samples.
% avgSNR=2;
% thresh=30;
% N=4;
% ampDef=0;
% PD=NonFlucSqLawD.PD4Threshold(avgSNR,thresh,N,ampDef)
% numSamples=1e5;
% PDSamp=mean(NonFlucSqLawD.rand([numSamples,1],avgSNR,N,ampDef)>=thresh)
%One will see that PD and PDSamp are both near 0.2329.
%
%REFERENCES:
%[1] J. V. Di Franco and W. L. Rubin, Radar Detection. SciTech Publishing
%    Inc., Rayliegh, NC: 2004.
%[2] P. C. Sofotasios and S. Freear, "New analytic results for the
%    incomplete Toronto function and incomplete Lipschitz-Hankel
%    integrals," in Proceedings from the SMBO/IEEE MTT-S International
%    Microwave and Optoelectronics Conference, Natal, Brazil, 29 Oct. - 1
%    Nov. 2011, pp. 44-47.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.

    if(nargin<4||isempty(ampDef))
        ampDef=1;
    end

    if(nargin<3||isempty(N))
       N=1; 
    end

    if(thresh==0)
        PD=1;
        return;
    end

    %Equation 10.4-13 in [1] is a noncentral Chi-Square distribution. We
    %can also look at Equation 10.4-4 in [1]. The point is that we can use
    %the MarcumQ function for a solution.
    lambda=N*2*avgSNR;
    nu=2*N;
    if(ampDef==0)
        PD=MarcumQ(nu/2,sqrt(lambda),sqrt(thresh));
    else
        PD=MarcumQ(nu/2,sqrt(lambda),sqrt(2*thresh));
    end
end

function PD=PD4PFA(avgSNR,PFA,N,ampDef)
%%PD4PFA Determine the detection probability of a target with a constant
%        signal to noise ratio (SNR) given the SNR and the probability of
%        false alarm, assuming that a square-law detector is used. The
%        noise in each sample is assumed to be Gaussian with variance 1
%        (normalized noise power).
%
%INPUTS: avgSNR A vector or matrix of average power signal to noise ratios
%          of the target at which one wishes to evaluate the detection
%          probability.
%      PFA The probability of false alarm, 0<=PFA<1.
%        N The number of pulses that are to be incoherently added for
%          detection (in a square-law detector). In a nonfluctuating model,
%          the target power is constant and the noise corrupting the pulses
%          varies. If this parameter is omitted or an empty matrix is
%          passed, N=1 is used.
%   ampDef This specified normalization (see help NonFlucSqLawD). Possible
%          values are:
%          0 The expected value of the squared magnitude of the noise is 2.
%          1 (The default if omitted or an empty matrix is passed) The
%            expected value of the squared magnitude of the noise is 1.
%            (A more common definition). 
%
%OUTPUTS: PD The detection probability of the target, 0<=PD<=1.
%
%This function just calls PFA2SquareLawThreshold to convert the false alarm
%rate into a normalized detection threshold after which it calls
%NonFlucSqLawD.PD4Threshold.
%
%EXAMPLE:
%This make a plot of PFA as a function of PD.
% avgSNR=2;
% N=4;
% ampDef=1;
% numPts=100;
% PFA=linspaceNoEnd(0,1,numPts,true);
% PD=zeros(numPts,1);
% for k=1:numPts
%     PD(k)=NonFlucSqLawD.PD4PFA(avgSNR,PFA(k),N,ampDef);
% end
% figure(1)
% clf
% hold on
% plot(PFA,PD,'linewidth',2);
% xlabel('PFA')
% ylabel('PD')
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.

    if(nargin<4||isempty(ampDef))
        ampDef=1;
    end

    if(nargin<3||isempty(N))
       N=1; 
    end
    
    thresh=PFA2SquareLawThreshold(PFA,N,ampDef);
    PD=NonFlucSqLawD.PD4Threshold(avgSNR,thresh,N,ampDef);
end

function vals=rand(NDims,avgSNR,N,ampDef,algorithm)
%%RAND Generate random variables representing a non-fluctuating target
%      detected by a square-law detector.
%
%INPUTS: NDims If NDims is a scalar, then rand returns an NDimsXNDims
%          matrix of random variables. If NDims=[M,N1] is a two-element row
%          vector, then rand returns an MXN1 matrix of random variables.
%   avgSNR The average power signal to noise ratio of the target.
%        N The number of pulses that are to be incoherently added for
%          detection (in a square-law detector). In a nonfluctuating model
%          the target power is constant and the noise corrupting the pulses
%          varies. If this parameter is omitted or an empty matrix is
%          passed, N=1 is used.
%   ampDef This specified normalization (see help NonFlucSqLawD). Possible
%          values are:
%          0 The expected value of the squared magnitude of the noise is 2.
%          1 (The default if omitted or an empty matrix is passed) The
%            expected value of the squared magnitude of the noise is 1.
%            (A more common definition). 
% algorithm This input can generally be omitted. It selects between an
%          efficient algorithm based on the noncentral chi-squared
%          distribution, and an inefficient algorithm that is just from 
%          the definition of summing together a nonfluctuating signal and
%          noise for each pulse. Possible values are:
%          0 (The default if omitted or an empty matrix is passed) usean
%            algorithm based on sampling a noncentral Chi-squared
%            distribution. This is given by Equation 10.4-13 in [1]. A
%            change of variables is performed when considering the
%            alternative normalizations.
%          1 Use an algorithm based on summing up the squared magnitudes of
%            multiple samples.
%
%OUTPUTS: vals A matrix whose dimensions are determined by N of the
%              generated random variables. 
%
%REFERENCES:
%[1] J. V. Di Franco and W. L. Rubin, Radar Detection. SciTech Publishing
%    Inc., Rayliegh, NC: 2004.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.

    if(nargin<5||isempty(algorithm))
        algorithm=0;
    end

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

    if(algorithm==0)
        %Equation 10.4-13 in [1] is a noncentral Chi-Square distribution,
        %but it must be transformed according to the normalziation.
        nu=2*N;
        lambda=2*N*avgSNR;
        if(ampDef==0)
            vals=ChiSquareD.rand(dims,nu,lambda);
        else
            vals=ChiSquareD.rand(dims,nu,lambda)/2;
        end
    else
        %Generate the samples (inefficiently) by directly generating a
        %nonfluctuating singla with a random phase and summing over a bunch
        %of deterministic pulses. The point of this algorithm is really
        %just as a check that the other algorithm is correct. This approach
        %is inefficient.
        if(ampDef==0)
            vals=zeros(dims);
            numPts=prod(dims);
            for k=1:numPts
                signal=sqrt(2*avgSNR)*exp(1j*2*pi*rand());
                for curPulse=1:N
                    z=signal+randn(1)+1j*randn(1);
                    vals(k)=vals(k)+abs(z)^2;
                end
            end
        else
            vals=zeros(dims);
            numPts=prod(dims);
            for k=1:numPts
                signal=sqrt(avgSNR)*exp(1j*2*pi*rand());
                for curPulse=1:N
                    z=signal+(randn(1)+1j*randn(1))/sqrt(2);
                    vals(k)=vals(k)+abs(z)^2;
                end
            end
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
