classdef SwerlingISqLawD
%%SWERLINGISQLAWD This class implements various functions and detections
%                statistics related to the Swerling I target model when
%                used with a square law detector. All methods are static
%                (the class holds no information and does not need to be
%                instantiated). All values are with respect to a normalized
%                noise model, of which two normalizations are available. 
%
%Implemented methods are: mean, var, logPDF, PDF, CDF, PD4Threshold, PD4PFA,
%                        avgSNR4PDThresh, avgSNR4PDPFA, thresh4AvgSNRPD,
%                        randSwerlingISqLawD, SNRDeriv, SNR2ndDeriv,
%                        logPDFSNRDeriv, logPDFSNR2ndDeriv, nthMoment
%
%Swerling models are given in [1] and in Chapter 11 of [2]. The Swerling I
%model assumes that the observed power SNR of the target varies in terms of
%an exponential distribution with rate parameter (1/avgSNR), that in a
%pulse train for detection the target SNR does not fluctuate from pulse to
%pulse, and that a square law detector is used to incoherently integrate
%multiple pulses. The difference from a Swerling II model is the target
%amplitude does not fluctuate from pulse to pulse.
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
%So, a Swerling 1 model basically says that for either definition, choose
%the signal power amplitude A given an average SNR according to
%A=ExponentialD.rand(1,1/avgSNR),
%the signal is sampled once given that SNR and the selected power amplitude
%definition above (the phase of the signal doesn't matter; it goes away
%when the mangitude of r is taken). Then, the measurement is
%y=sum_{i=1}^N abs(xSig+xNoise_i)^2
%so the signal value remains constant but the noise value is different at
%each sample (It doesn't matter if the noise phase changes at each sample).
%In a Swerling II model, the noise value is different at each sample.
%
%REFERENCES:
%[1] P. Swerling, "Probability of detection for fluctuating targets," The
%    RAND Corporation, Santa Monica, CA, Tech. Rep. RM-1217, 1954.
%[2] J. V. Di Franco and W. L. Rubin, Radar Detection. SciTech Publishing
%    Inc., Rayliegh, NC: 2004.
%
%August 2024 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

methods(Static)

function val=nthMoment(n,avgSNR,N,ampDef)
%%NTHMOMENT Obtain the nth noncentral moment of the Swerling I
%           distribution. This function supports up through the 15th
%           moment.
%
%INPUTS: n The number moment desired; 0<=n<=15.
%   avgSNR The average power signal to noise ratio of the target.
%        N The number of pulses that are to be incoherently added for
%          detection (in a square-law detector). In a Swerling I model, a
%          single realization of the PDFtarget power is used across all
%          pulses, but the noise corrupting the pulses varies. If this
%          parameter is omitted or an empty matrix is passed, N=1 is used.
%   ampDef Possible values are:
%          0 The expected value of the squared magnitude of the noise is 2.
%          1 (The default if omitted or an empty matrix is passed) The
%            expected value of the squared magnitude of the noise is 1.
%            (A more common definition).
%
%OUTPUT: val The value of the nth moment of the Swerling I distribution.
%
%Noncentral can be found from the charcateristic function as described in
%[1]. The characteristic function for the Swerling I distribution is given
%in Equation III.3 of [2]. Explicit solutions for n=0 to n=15 were found
%and are implemented here.
%
%EXAMPLE:
%Here, we compute the third moment and we also estimate the third moment
%via Monte Carlo simulations. The values are reasomably close.
% avgSNR=12;
% N=4;
% n=3;
% numSamples=1e6;
% ampDef=1;
% meanVal=SwerlingISqLawD.nthMoment(n,avgSNR,N,ampDef)
% meanSampVal=mean(SwerlingISqLawD.rand([numSamples,1],avgSNR,N,ampDef).^n)
%
%REFERENCES:
%[1] Rowland, Todd and Weisstein, Eric W. "Characteristic Function." From
%    MathWorld--A Wolfram Web Resource.
%    https://mathworld.wolfram.com/CharacteristicFunction.html
%[2] P. Swerling, "Probability of detection for fluctuating targets," The
%    RAND Corporation, Santa Monica, CA, Tech. Rep. RM-1217, 1954.
%
%August 2024 David F. Crouse, Naval Research Laboratory, Washington D.C.

if(nargin<4||isempty(ampDef))
    ampDef=1;
end

if(nargin<3||isempty(N))
    N=1;
end

switch(n)
    case 0
        val=1;
    case 1
        val=(1+avgSNR)*N;
    case 2
        val=N*(1+N+2*avgSNR*(1+N+avgSNR*N));
    case 3
        val=N*(6*avgSNR^3*N^2+6*avgSNR^2*N*(2+N)+(1+N)*(2+N)+3*avgSNR*(1+N)*(2+N));
    case 4
        val=N*(24*avgSNR^4*N^3+24*avgSNR^3*N^2*(3+N)+12*avgSNR^2*N*(2+N)*(3+N)+(1+N)*(2+N)*(3+N)+4*avgSNR*(1+N)*(2+N)*(3+N));
    case 5
        val=N*(120*avgSNR^5*N^4+120*avgSNR^4*N^3*(4+N)+60*avgSNR^3*N^2*(3+N)*(4+N)+20*avgSNR^2*N*(2+N)*(3+N)*(4+N)+(1+N)*(2+N)*(3+N)*(4+N)+5*avgSNR*(1+N)*(2+N)*(3+N)*(4+N));
    case 6
        val=N*(720*avgSNR^6*N^5+720*avgSNR^5*N^4*(5+N)+360*avgSNR^4*N^3*(4+N)*(5+N)+120*avgSNR^3*N^2*(3+N)*(4+N)*(5+N)+30*avgSNR^2*N*(2+N)*(3+N)*(4+N)*(5+N)+(1+N)*(2+N)*(3+N)*(4+N)*(5+N)+6*avgSNR*(1+N)*(2+N)*(3+N)*(4+N)*(5+N));
    case 7
        val=N*(5040*avgSNR^7*N^6+5040*avgSNR^6*N^5*(6+N)+2520*avgSNR^5*N^4*(5+N)*(6+N)+840*avgSNR^4*N^3*(4+N)*(5+N)*(6+N)+210*avgSNR^3*N^2*(3+N)*(4+N)*(5+N)*(6+N)+42*avgSNR^2*N*(2+N)*(3+N)*(4+N)*(5+N)*(6+N)+(1+N)*(2+N)*(3+N)*(4+N)*(5+N)*(6+N)+7*avgSNR*(1+N)*(2+N)*(3+N)*(4+N)*(5+N)*(6+N));
    case 8
        val=N*(40320*avgSNR^8*N^7+40320*avgSNR^7*N^6*(7+N)+20160*avgSNR^6*N^5*(6+N)*(7+N)+6720*avgSNR^5*N^4*(5+N)*(6+N)*(7+N)+1680*avgSNR^4*N^3*(4+N)*(5+N)*(6+N)*(7+N)+336*avgSNR^3*N^2*(3+N)*(4+N)*(5+N)*(6+N)*(7+N)+56*avgSNR^2*N*(2+N)*(3+N)*(4+N)*(5+N)*(6+N)*(7+N)+(1+N)*(2+N)*(3+N)*(4+N)*(5+N)*(6+N)*(7+N)+8*avgSNR*(1+N)*(2+N)*(3+N)*(4+N)*(5+N)*(6+N)*(7+N));
    case 9
        val=N*(362880*avgSNR^9*N^8+362880*avgSNR^8*N^7*(8+N)+181440*avgSNR^7*N^6*(7+N)*(8+N)+60480*avgSNR^6*N^5*(6+N)*(7+N)*(8+N)+15120*avgSNR^5*N^4*(5+N)*(6+N)*(7+N)*(8+N)+3024*avgSNR^4*N^3*(4+N)*(5+N)*(6+N)*(7+N)*(8+N)+504*avgSNR^3*N^2*(3+N)*(4+N)*(5+N)*(6+N)*(7+N)*(8+N)+72*avgSNR^2*N*(2+N)*(3+N)*(4+N)*(5+N)*(6+N)*(7+N)*(8+N)+(1+N)*(2+N)*(3+N)*(4+N)*(5+N)*(6+N)*(7+N)*(8+N)+9*avgSNR*(1+N)*(2+N)*(3+N)*(4+N)*(5+N)*(6+N)*(7+N)*(8+N));
    case 10
        val=(-1+N)*N*(1+N)*(2+N)*(3+N)*(4+N)*(5+N)*(6+N)*(7+N)*(8+N)+10*(-1+N)*N*(1+N)*(2+N)*(3+N)*(4+N)*(5+N)*(6+N)*(7+N)*(1+avgSNR*N)+90*(-1+N)*N*(1+N)*(2+N)*(3+N)*(4+N)*(5+N)*(6+N)*(1+avgSNR*N)^2+720*(-1+N)*N*(1+N)*(2+N)*(3+N)*(4+N)*(5+N)*(1+avgSNR*N)^3+5040*(-1+N)*N*(1+N)*(2+N)*(3+N)*(4+N)*(1+avgSNR*N)^4+30240*(-1+N)*N*(1+N)*(2+N)*(3+N)*(1+avgSNR*N)^5+151200*(-1+N)*N*(1+N)*(2+N)*(1+avgSNR*N)^6+1814400*(-1+N)*N*(1+avgSNR*N)^8+3628800*(-1+N)*(1+avgSNR*N)^9+3628800*(1+avgSNR*N)^10+604800*N*(1+avgSNR*N)^7*(-1+N^2);
    case 11
        val=(-1+N)*N*(1+N)*(2+N)*(3+N)*(4+N)*(5+N)*(6+N)*(7+N)*(8+N)*(9+N)+11*(-1+N)*N*(1+N)*(2+N)*(3+N)*(4+N)*(5+N)*(6+N)*(7+N)*(8+N)*(1+avgSNR*N)+110*(-1+N)*N*(1+N)*(2+N)*(3+N)*(4+N)*(5+N)*(6+N)*(7+N)*(1+avgSNR*N)^2+990*(-1+N)*N*(1+N)*(2+N)*(3+N)*(4+N)*(5+N)*(6+N)*(1+avgSNR*N)^3+7920*(-1+N)*N*(1+N)*(2+N)*(3+N)*(4+N)*(5+N)*(1+avgSNR*N)^4+55440*(-1+N)*N*(1+N)*(2+N)*(3+N)*(4+N)*(1+avgSNR*N)^5+332640*(-1+N)*N*(1+N)*(2+N)*(3+N)*(1+avgSNR*N)^6+1663200*(-1+N)*N*(1+N)*(2+N)*(1+avgSNR*N)^7+19958400*(-1+N)*N*(1+avgSNR*N)^9+39916800*(-1+N)*(1+avgSNR*N)^10+39916800*(1+avgSNR*N)^11+6652800*N*(1+avgSNR*N)^8*(-1+N^2);
    case 12
        val=(-1+N)*N*(1+N)*(2+N)*(3+N)*(4+N)*(5+N)*(6+N)*(7+N)*(8+N)*(9+N)*(10+N)+12*(-1+N)*N*(1+N)*(2+N)*(3+N)*(4+N)*(5+N)*(6+N)*(7+N)*(8+N)*(9+N)*(1+avgSNR*N)+132*(-1+N)*N*(1+N)*(2+N)*(3+N)*(4+N)*(5+N)*(6+N)*(7+N)*(8+N)*(1+avgSNR*N)^2+1320*(-1+N)*N*(1+N)*(2+N)*(3+N)*(4+N)*(5+N)*(6+N)*(7+N)*(1+avgSNR*N)^3+11880*(-1+N)*N*(1+N)*(2+N)*(3+N)*(4+N)*(5+N)*(6+N)*(1+avgSNR*N)^4+95040*(-1+N)*N*(1+N)*(2+N)*(3+N)*(4+N)*(5+N)*(1+avgSNR*N)^5+665280*(-1+N)*N*(1+N)*(2+N)*(3+N)*(4+N)*(1+avgSNR*N)^6+3991680*(-1+N)*N*(1+N)*(2+N)*(3+N)*(1+avgSNR*N)^7+19958400*(-1+N)*N*(1+N)*(2+N)*(1+avgSNR*N)^8+239500800*(-1+N)*N*(1+avgSNR*N)^10+479001600*(-1+N)*(1+avgSNR*N)^11+479001600*(1+avgSNR*N)^12+79833600*N*(1+avgSNR*N)^9*(-1+N^2);
    case 13
        val=(-1+N)*N*(1+N)*(2+N)*(3+N)*(4+N)*(5+N)*(6+N)*(7+N)*(8+N)*(9+N)*(10+N)*(11+N)+13*(-1+N)*N*(1+N)*(2+N)*(3+N)*(4+N)*(5+N)*(6+N)*(7+N)*(8+N)*(9+N)*(10+N)*(1+avgSNR*N)+156*(-1+N)*N*(1+N)*(2+N)*(3+N)*(4+N)*(5+N)*(6+N)*(7+N)*(8+N)*(9+N)*(1+avgSNR*N)^2+1716*(-1+N)*N*(1+N)*(2+N)*(3+N)*(4+N)*(5+N)*(6+N)*(7+N)*(8+N)*(1+avgSNR*N)^3+17160*(-1+N)*N*(1+N)*(2+N)*(3+N)*(4+N)*(5+N)*(6+N)*(7+N)*(1+avgSNR*N)^4+154440*(-1+N)*N*(1+N)*(2+N)*(3+N)*(4+N)*(5+N)*(6+N)*(1+avgSNR*N)^5+1235520*(-1+N)*N*(1+N)*(2+N)*(3+N)*(4+N)*(5+N)*(1+avgSNR*N)^6+8648640*(-1+N)*N*(1+N)*(2+N)*(3+N)*(4+N)*(1+avgSNR*N)^7+51891840*(-1+N)*N*(1+N)*(2+N)*(3+N)*(1+avgSNR*N)^8+259459200*(-1+N)*N*(1+N)*(2+N)*(1+avgSNR*N)^9+3113510400*(-1+N)*N*(1+avgSNR*N)^11+6227020800*(-1+N)*(1+avgSNR*N)^12+6227020800*(1+avgSNR*N)^13+1037836800*N*(1+avgSNR*N)^10*(-1+N^2);
    case 14
        val=(-1+N)*N*(1+N)*(2+N)*(3+N)*(4+N)*(5+N)*(6+N)*(7+N)*(8+N)*(9+N)*(10+N)*(11+N)*(12+N)+14*(-1+N)*N*(1+N)*(2+N)*(3+N)*(4+N)*(5+N)*(6+N)*(7+N)*(8+N)*(9+N)*(10+N)*(11+N)*(1+avgSNR*N)+182*(-1+N)*N*(1+N)*(2+N)*(3+N)*(4+N)*(5+N)*(6+N)*(7+N)*(8+N)*(9+N)*(10+N)*(1+avgSNR*N)^2+2184*(-1+N)*N*(1+N)*(2+N)*(3+N)*(4+N)*(5+N)*(6+N)*(7+N)*(8+N)*(9+N)*(1+avgSNR*N)^3+24024*(-1+N)*N*(1+N)*(2+N)*(3+N)*(4+N)*(5+N)*(6+N)*(7+N)*(8+N)*(1+avgSNR*N)^4+240240*(-1+N)*N*(1+N)*(2+N)*(3+N)*(4+N)*(5+N)*(6+N)*(7+N)*(1+avgSNR*N)^5+2162160*(-1+N)*N*(1+N)*(2+N)*(3+N)*(4+N)*(5+N)*(6+N)*(1+avgSNR*N)^6+17297280*(-1+N)*N*(1+N)*(2+N)*(3+N)*(4+N)*(5+N)*(1+avgSNR*N)^7+121080960*(-1+N)*N*(1+N)*(2+N)*(3+N)*(4+N)*(1+avgSNR*N)^8+726485760*(-1+N)*N*(1+N)*(2+N)*(3+N)*(1+avgSNR*N)^9+3632428800*(-1+N)*N*(1+N)*(2+N)*(1+avgSNR*N)^10+43589145600*(-1+N)*N*(1+avgSNR*N)^12+87178291200*(-1+N)*(1+avgSNR*N)^13+87178291200*(1+avgSNR*N)^14+14529715200*N*(1+avgSNR*N)^11*(-1+N^2);
    case 15
        val=(-1+N)*N*(1+N)*(2+N)*(3+N)*(4+N)*(5+N)*(6+N)*(7+N)*(8+N)*(9+N)*(10+N)*(11+N)*(12+N)*(13+N)+15*(-1+N)*N*(1+N)*(2+N)*(3+N)*(4+N)*(5+N)*(6+N)*(7+N)*(8+N)*(9+N)*(10+N)*(11+N)*(12+N)*(1+avgSNR*N)+210*(-1+N)*N*(1+N)*(2+N)*(3+N)*(4+N)*(5+N)*(6+N)*(7+N)*(8+N)*(9+N)*(10+N)*(11+N)*(1+avgSNR*N)^2+2730*(-1+N)*N*(1+N)*(2+N)*(3+N)*(4+N)*(5+N)*(6+N)*(7+N)*(8+N)*(9+N)*(10+N)*(1+avgSNR*N)^3+32760*(-1+N)*N*(1+N)*(2+N)*(3+N)*(4+N)*(5+N)*(6+N)*(7+N)*(8+N)*(9+N)*(1+avgSNR*N)^4+360360*(-1+N)*N*(1+N)*(2+N)*(3+N)*(4+N)*(5+N)*(6+N)*(7+N)*(8+N)*(1+avgSNR*N)^5+3603600*(-1+N)*N*(1+N)*(2+N)*(3+N)*(4+N)*(5+N)*(6+N)*(7+N)*(1+avgSNR*N)^6+32432400*(-1+N)*N*(1+N)*(2+N)*(3+N)*(4+N)*(5+N)*(6+N)*(1+avgSNR*N)^7+259459200*(-1+N)*N*(1+N)*(2+N)*(3+N)*(4+N)*(5+N)*(1+avgSNR*N)^8+1816214400*(-1+N)*N*(1+N)*(2+N)*(3+N)*(4+N)*(1+avgSNR*N)^9+10897286400*(-1+N)*N*(1+N)*(2+N)*(3+N)*(1+avgSNR*N)^10+54486432000*(-1+N)*N*(1+N)*(2+N)*(1+avgSNR*N)^11+653837184000*(-1+N)*N*(1+avgSNR*N)^13+1307674368000*(-1+N)*(1+avgSNR*N)^14+1307674368000*(1+avgSNR*N)^15+217945728000*N*(1+avgSNR*N)^12*(-1+N^2);
    otherwise
         error('The moment for the specified n is not supported.')
end

if(ampDef==0)
    val=2^n*val;
end

end


function val=logPDFSNRDeriv(v,avgSNR,N,ampDef)
%%LOGPDFSNRDERIV Evaluate the first derivative of the natural logarithm
%    of the PDF of the Swerling I distribution with respect to the
%    average SNR.
%
%INPUTS: v The point or points at which the first derivative of the PDF
%          with respect to avgSNR should be evaluated.
%   avgSNR The average power signal to noise ratio of the target.
%        N The number of pulses that are to be incoherently added for
%          detection (in a square-law detector). In a Swerling I model, a
%          single realization of the PDFtarget power is used across all
%          pulses, but the noise corrupting the pulses varies. If this
%          parameter is omitted or an empty matrix is passed, N=1 is used.
%   ampDef Possible values are:
%          0 The expected value of the squared magnitude of the noise is 2.
%          1 (The default if omitted or an empty matrix is passed) The
%            expected value of the squared magnitude of the noise is 1.
%            (A more common definition).
%
%OUTPUTS: val The first derivative of the logarithm of the PDF with respect
%             to avgSNR. This has the same dimensions as x.
%
%The PDF for a single pulse is given in Equation III.4 of [1] and the PDF
%for multiple pulses is Equation III.5 of [1] (Note the definition of the
%incomplete gamma function there differs from gammainc in Matlab). This
%function just implements the analytic first derivative with respect to
%avgSNR of the logarithm of the PDF.
%
%EXAMPLE:
%Here, the output of this function is plotted along with numeric
%differentiation for a variety of v values and a fixed avgSNR. One can see
%a good level of agreement.
% avgSNR=3;
% N=4;
% ampDef=1;
% 
% numPts=500;
% x=linspace(0,80,numPts);
% valNumDiff=zeros(numPts,1);
% valAnalytic=zeros(numPts,1);
% for k=1:numPts
%     f=@(SNR)SwerlingISqLawD.logPDF(x(k),SNR,N,ampDef);
%     valNumDiff(k)=numDiff(avgSNR,f,1);
%     valAnalytic(k)=SwerlingISqLawD.logPDFSNRDeriv(x(k),avgSNR,N,ampDef);
% end
% figure(1)
% clf
% hold on
% plot(x,valNumDiff,'linewidth',4)
% plot(x,valAnalytic,'linewidth',2)
% legend('Numeric Differentiation','Analytic Derivative','location','southeast')
%
%REFERENCES:
%[1] P. Swerling, "Probability of detection for fluctuating targets," The
%    RAND Corporation, Santa Monica, CA, Tech. Rep. RM-1217, 1954.
%
%August 2024 David F. Crouse, Naval Research Laboratory, Washington D.C.

if(nargin<4||isempty(ampDef))
    ampDef=1;
end

if(nargin<3||isempty(N))
    N=1;
end

if(ampDef==0)
    v=v/2;
end

if(N==1)
    val=(v-1-avgSNR)/(1+avgSNR)^2;
else
    x=(N*v*avgSNR)/((1+N*avgSNR));
    a=N-1;
    gammaInc=gammainc(x,a);
    gammaIncDeriv=(exp(-x).*x.^(a-1))/gamma(a);

    term2=N*v.*gammaIncDeriv./gammaInc;

    %Deal with x=0 by substituting the limit.
    term2(gammaInc==0)=(N-1)*(1+N*avgSNR)/avgSNR;
    
    val=((1-N*(1-v*avgSNR+N*avgSNR*(1+avgSNR)))/avgSNR+term2)/(1+N*avgSNR)^2;
end

val(v<0)=0;
end

function val=logPDFSNR2ndDeriv(v,avgSNR,N,ampDef)
%%LOGPDFSNR2NDDERIV Evaluate the second derivative of the natural logarithm
%    of the PDF of the Swerling I distribution with respect to the average
%    SNR.
%
%INPUTS: v The point or points at which the first derivative of the PDF
%          with respect to avgSNR should be evaluated.
%   avgSNR The average power signal to noise ratio of the target.
%        N The number of pulses that are to be incoherently added for
%          detection (in a square-law detector). In a Swerling I model, a
%          single realization of the PDFtarget power is used across all
%          pulses, but the noise corrupting the pulses varies. If this
%          parameter is omitted or an empty matrix is passed, N=1 is used.
%   ampDef Possible values are:
%          0 The expected value of the squared magnitude of the noise is 2.
%          1 (The default if omitted or an empty matrix is passed) The
%            expected value of the squared magnitude of the noise is 1.
%            (A more common definition).
%
%OUTPUTS: val The second derivative of the logarithm of the PDF with
%             respect to avgSNR. This has the same dimensions as x.
%
%The PDF for a single pulse is given in Equation III.4 of [1] and the PDF
%for multiple pulses is Equation III.5 of [1] (Note the definition of the
%incomplete gamma function there differs from gammainc in Matlab). This
%function just implements the analytic second derivative with respect to
%avgSNR of the logarithm of the PDF.
%
%EXAMPLE:
%Here, the output of this function is plotted along with numeric
%differentiation for a variety of v values and a fixed avgSNR. One can see
%a good level of agreement.
% avgSNR=3;
% N=4;
% ampDef=0;
% numPts=500;
% x=linspace(0,80,numPts);
% valNumDiff=zeros(numPts,1);
% valAnalytic=zeros(numPts,1);
% for k=1:numPts
%     f=@(SNR)SwerlingISqLawD.logPDFSNRDeriv(x(k),SNR,N,ampDef);
%     valNumDiff(k)=numDiff(avgSNR,f,1);
%     valAnalytic(k)=SwerlingISqLawD.logPDFSNR2ndDeriv(x(k),avgSNR,N,ampDef);
% end
% figure(1)
% clf
% hold on
% plot(x,valNumDiff,'linewidth',4)
% plot(x,valAnalytic,'linewidth',2)
% legend('Numeric Differentiation','Analytic Derivative','location','southeast')
%
%REFERENCES:
%[1] P. Swerling, "Probability of detection for fluctuating targets," The
%    RAND Corporation, Santa Monica, CA, Tech. Rep. RM-1217, 1954.
%
%August 2024 David F. Crouse, Naval Research Laboratory, Washington D.C.

if(nargin<4||isempty(ampDef))
    ampDef=1;
end

if(nargin<3||isempty(N))
    N=1;
end

if(ampDef==0)
    v=v/2;
end

if(N==1)
    val=(1-2*v+avgSNR)/(1+avgSNR)^3;
else
    x=(N*v*avgSNR)/((1+N*avgSNR));
    a=N-1;
    gammaInc=gammainc(x,a);
    gammaIncDeriv=(exp(-x).*x.^(a-1))/gamma(a);
    if(N==2)
        gammaInc2Deriv=-exp(-x).*(x-1);
    else
        gammaInc2Deriv=(exp(-x).*(a-x-1).*x.^(a-2))/gamma(a);
    end

    val=(1/((1+N*avgSNR)^4)).*((-1+N+N*avgSNR*(-4+N*(4-2*(2+v).*avgSNR+N*avgSNR*(5-2*v*avgSNR+N*avgSNR*(2+avgSNR)))))/avgSNR^2+(N^2*v.*(-v.*gammaIncDeriv.^2+gammaInc.*(-2*(1+N*avgSNR).*gammaIncDeriv+v.*gammaInc2Deriv)))./gammaInc.^2);
end

val(v<0)=0;
end

function val=SNRDeriv(v,avgSNR,N,ampDef)
%%SNRDERIV Evaluate the first derivative of the PDF of the Swerling I
%    distribution with respect to the average SNR.
%
%INPUTS: v The point or points at which the first derivative of the PDF
%          with respect to avgSNR should be evaluated.
%   avgSNR The average power signal to noise ratio of the target.
%        N The number of pulses that are to be incoherently added for
%          detection (in a square-law detector). In a Swerling I model, a
%          single realization of the PDFtarget power is used across all
%          pulses, but the noise corrupting the pulses varies. If this
%          parameter is omitted or an empty matrix is passed, N=1 is used.
%   ampDef Possible values are:
%          0 The expected value of the squared magnitude of the noise is 2.
%          1 (The default if omitted or an empty matrix is passed) The
%            expected value of the squared magnitude of the noise is 1.
%            (A more common definition).
%
%OUTPUTS: val The first derivative of the PDF with respect to avgSNR. This
%             has the same dimensions as x.
%
%The PDF for a single pulse is given in Equation III.4 of [1] and the PDF
%for multiple pulses is Equation III.5 of [1] (Note the definition of the
%incomplete gamma function there differs from gammainc in Matlab). This
%function just implements the analytic first derivative with respect to
%avgSNR.
%
%EXAMPLE:
%Here, the output of this function is plotted along with numeric
%differentiation for a variety of v values and a fixed avgSNR. One can see
%a good level of agreement.
%they agree.
% avgSNR=3;
% N=4;
% ampDef=1;
% 
% numPts=500;
% x=linspace(0,80,numPts);
% valNumDiff=zeros(numPts,1);
% valAnalytic=zeros(numPts,1);
% for k=1:numPts
%     f=@(SNR)SwerlingISqLawD.PDF(x(k),SNR,N,ampDef);
%     valNumDiff(k)=numDiff(avgSNR,f,1);
%     valAnalytic(k)=SwerlingISqLawD.SNRDeriv(x(k),avgSNR,N,ampDef);
% end
% figure(1)
% clf
% hold on
% plot(x,valNumDiff,'linewidth',4)
% plot(x,valAnalytic,'linewidth',2)
% legend('Numeric Differentiation','Analytic Derivative','location','southeast')
%
%REFERENCES:
%[1] P. Swerling, "Probability of detection for fluctuating targets," The
%    RAND Corporation, Santa Monica, CA, Tech. Rep. RM-1217, 1954.
%
%August 2024 David F. Crouse, Naval Research Laboratory, Washington D.C.

if(nargin<4||isempty(ampDef))
    ampDef=1;
end

if(nargin<3||isempty(N))
    N=1;
end

if(ampDef==0)
    v=v/2;
end

if(N==1)
    val=(exp(-(v./(1+avgSNR))).*(v-avgSNR-1))/(1+avgSNR)^3;
else
    x=(N*v*avgSNR)/((1+N*avgSNR));
    a=N-1;
    gammaInc=gammainc(x,a);
    gammaIncDeriv=(exp(-x).*x.^(a-1))/gamma(a);

    val=-((exp(-(v/(1+N*avgSNR))).*N.*(1+1/(N*avgSNR))^N.*((-1+N-N*v*avgSNR+N^2*avgSNR*(1+avgSNR)).*gammaInc-N*v.*avgSNR.*gammaIncDeriv))/(1+N*avgSNR)^4);
end

if(ampDef==0)
    val=val/2;
end

val(v<0)=0;

end

function val=SNR2ndDeriv(v,avgSNR,N,ampDef)
%%SNR2NDDERIV Evaluate the second derivative of the PDF of the Swerling I
%    distribution with respect to the average SNR.
%
%INPUTS: v The point or points at which the second derivative of the PDF
%          with respect to avgSNR should be evaluated.
%   avgSNR The average power signal to noise ratio of the target.
%        N The number of pulses that are to be incoherently added for
%          detection (in a square-law detector). In a Swerling I model, a
%          single realization of the target power is used across all
%          pulses, but the noise corrupting the pulses varies. If this
%          parameter is omitted or an empty matrix is passed, N=1 is used.
%   ampDef Possible values are:
%          0 The expected value of the squared magnitude of the noise is 2.
%          1 (The default if omitted or an empty matrix is passed) The
%            expected value of the squared magnitude of the noise is 1.
%            (A more common definition).
%
%OUTPUTS: val The second derivative of the PDF with respect to avgSNR. Tihs
%             This has the same dimensions as x.
%
%The PDF for a single pulse is given in Equation III.4 of [1] and the PDF
%for multiple pulses is Equation III.5 of [1] (Note the definition of the
%incomplete gamma function there differs from gammainc in Matlab). This
%function just implements the analytic first derivative with respect to
%avgSNR.
%
%EXAMPLE:
%Here, the output of this function is plotted along with numeric
%differentiation for a variety of v values and a fixed avgSNR. One can see
%a good level of agreement.
% avgSNR=3;
% N=4;
% ampDef=1;
% numPts=500;
% x=linspace(0,80,numPts);
% 
% valNumDiff=zeros(numPts,1);
% valAnalytic=zeros(numPts,1);
% for k=1:numPts
%     f=@(SNR)SwerlingISqLawD.SNRDeriv(x(k),SNR,N,ampDef);
%     valNumDiff(k)=numDiff(avgSNR,f,1);
%     valAnalytic(k)=SwerlingISqLawD.SNR2ndDeriv(x(k),avgSNR,N,ampDef);
% end
% figure(1)
% clf
% hold on
% plot(x,valNumDiff,'linewidth',4)
% plot(x,valAnalytic,'linewidth',2)
% legend('Numeric Differentiation','Analytic Derivative','location','southeast')
%
%REFERENCES:
%[1] P. Swerling, "Probability of detection for fluctuating targets," The
%    RAND Corporation, Santa Monica, CA, Tech. Rep. RM-1217, 1954.
%
%August 2024 David F. Crouse, Naval Research Laboratory, Washington D.C.

if(nargin<4||isempty(ampDef))
    ampDef=1;
end

if(nargin<3||isempty(N))
    N=1;
end

if(ampDef==0)
    v=v/2;
end

if(N==1)
    val=(exp(-(v/(1+avgSNR))).*(v.^2-4*v*(1+avgSNR)+2*(1+avgSNR)^2))/(1+avgSNR)^5;
else
    x=(N*v*avgSNR)/((1+N*avgSNR));
    a=N-1;
    gammaInc=gammainc(x,a);
    gammaIncDeriv=(exp(-x).*x.^(a-1))/gamma(a);
    if(N==2)
        gammaInc2Deriv=-exp(-x).*(x-1);
    else
        gammaInc2Deriv=(exp(-x).*(a-x-1).*x.^(a-2))/gamma(a);
    end
    val=(1/(avgSNR*(1+N*avgSNR)^6))*exp(-(v/(1+N*avgSNR))).*N^2.*(1+1/(N*avgSNR))^N.*((N-1+2*(N-1)*(2+N-v).*avgSNR+N*(-6+7*N+N^2-2*(1+N)*v+v.^2)*avgSNR^2+4*N^2*(N-v).*avgSNR^3+2*N^3*avgSNR^4).*gammaInc+v.*avgSNR.*(-2*(N-1+N*avgSNR*(1+N-v+2*N*avgSNR)).*gammaIncDeriv+N*v.*avgSNR.*gammaInc2Deriv));
end

if(ampDef==0)
    val=val/2;
end

val(v<0)=0;

end

function val=mean(avgSNR,N,ampDef)
%%MEAN Obtain the mean of the distribution of the detection power
%      (normalized to a noise variance of 1) under a Swerling I model in a
%      square-law detector.
%
%INPUTS:avgSNR The average power signal to noise ratio of the target.
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
%
%OUTPUTS: val The value of the mean for the given parameters.
%
%The characteristic function of the distribution is given in Equation
%11.2-23c in [1]. It is well know that 1i^(-n)*d^n C(t)/dt^n evaluated at t=0
%gives the value of the nth moment. That is, sqrt(-1) raised to the -n times
%the nth derivative of the characteristic function evaluated at zero gives
%the nth noncentral moment. The 1st moment (the mean) was obtained in such
%a manner. However, for ampDef=0,the first moment from that expression has
%to be multiplied by 2, because as seen in Equation 10.4-4, the definition
%in [1] is in terms of a scaled version of the square law detector.
%
%EXAMPLE:
%Here, we validate the mean by generating random samples and comparing the
%computed mean to the sample mean.
% avgSNR=2;
% N=4;
% numSamples=1e5;
% ampDef=1;
% meanVal=SwerlingISqLawD.mean(avgSNR,N,ampDef)
% meanSampVal=mean(SwerlingISqLawD.rand([numSamples,1],avgSNR,N,ampDef))
%One will see that both mean values are about 12.
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
    if(ampDef==0)
        val=2*N.*(1+avgSNR);
    else
        val=N.*(1+avgSNR);
    end
end

function val=var(avgSNR,N,ampDef)
%%VAR Obtain the variance of the distribution of the detection power
%     (normalized to a noise variance of 1) under a Swerling I model in a
%     square-law detector.
%
%INPUTS:avgSNR The average power signal to noise ratio of the target.
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
%
%OUTPUTS: val The value of the variance for the given parameters.
%
%The characteristic function of the distribution is given in Equation
%11.2-23c in [1]. It is well know that 1i^(-n)*d^n C(t)/dt^n evaluated at t=0
%gives the value of the nth moment. That is, sqrt(-1) raised to the -n times
%the nth derivative of the characteristic function evaluated at zero gives
%the nth noncentral moment. The 1st and second noncentral moments were
%obtained in such a manner. One then uses the identity that
%variance=E{x^2}-E{x}^2
%to get the variance. E{x^2} is the second noncentral moment and E{x} is
%the first noncentral moment. However, the values must be scaled
%differently if ampDef=0.
%
%EXAMPLE:
%Here, we validate the variance by generating random samples and comparing
%the computed variance to the sample variance.
% avgSNR=2;
% N=4;
% numSamples=1e5;
% varVal=SwerlingISqLawD.var(avgSNR,N)
% varSampVal=var(SwerlingISqLawD.rand([numSamples,1],avgSNR,N))
%One will see that both variance values are about 336.
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
    
    if(ampDef==0)
        val=4*N*(1+avgSNR*(2+N*avgSNR));
    else
        val=N*(1+avgSNR*(2+N*avgSNR));
    end
end

function val=logPDF(v,avgSNR,N,ampDef)
%%LOGPDF Evaluate the natural logarithm of the scalar probability density
%     function (PDF) of the distribution of the detection power (normalized
%     to a noise variance of 1) under a Swerling I model in a square-law
%     detector.
%
%INPUTS: v The point or points at which the PDF should be evaluated. The
%          support of the distirbution if v>0, so for anything else -Inf is
%          returned.
%   avgSNR The average power signal to noise ratio of the target.
%        N The number of pulses that are to be incoherently added for
%          detection (in a square-law detector). In a Swerling I model, a
%          single realization of the target power is used across all
%          pulses, but the noise corrupting the pulses varies. If this
%          parameter is omitted or an empty matrix is passed, N=1 is used.
%   ampDef Possible values are:
%          0 The expected value of the squared magnitude of the noise is 2.
%          1 (The default if omitted or an empty matrix is passed) The
%            expected value of the squared magnitude of the noise is 1.
%            (A more common definition). 
%
%OUTPUTS: val The values of the natural logarithm of the PDF evaluated at
%             the given points.
%
%Evaluate the natural logarithm of the PDF for a single pulse, which is
%given in Equation III.4 of [1] and the PDF for multiple pulses is
%Equation III.5 of [1] (Note the definition of the incomplete gamma 
%unction there differs from gammainc in Matlab). However,
%a clearer derivation is given in Chapter 11 of [2].
%
%EXAMPLE:
%This just verifies that the output of this fucntion is the same as what
%one gets by taking the logarithm of the output of the PDF function. The
%two lines will overlap.
% avgSNR=3;
% N=2;
% ampDef=0;
% 
% figure(1)
% clf
% hold on
% numPoints=1000;
% x=linspace(0,80,numPoints);
% vals0=log(SwerlingISqLawD.PDF(x,avgSNR,N,ampDef));
% vals=SwerlingISqLawD.logPDF(x,avgSNR,N,ampDef);
% plot(x,vals0,'linewidth',4)
% plot(x,vals,'linewidth',2)
% h1=xlabel('x');
% h2=ylabel('PDF(x)');
% set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')
%
%REFERENCES:
%[1] P. Swerling, "Probability of detection for fluctuating targets," The
%    RAND Corporation, Santa Monica, CA, Tech. Rep. RM-1217, 1954.
%[2] J. V. Di Franco and W. L. Rubin, Radar Detection. SciTech Publishing
%    Inc., Rayliegh, NC: 2004.
%
%August 2024 David F. Crouse, Naval Research Laboratory, Washington D.C.

    if(nargin<4||isempty(ampDef))
        ampDef=1;
    end

    if(nargin<3||isempty(N))
        N=1; 
    end

    if(ampDef==0)
        %Initial change of variable to undo the scaling from Equation
        %10.4-4 in [2].
        v=v/2;
    end

    if(N==1)
        %Equation III.4 in [1].
        val=-v./(1+avgSNR)-log(1+avgSNR);
    else
        NxBar=N*avgSNR;

        %The logarithm of Equation III.5 in [1]. Note that gammainc is
        %defined differently than the incomplete gamma function used in
        %[1].
        val=(N-2)*log(1+1/NxBar)-log(NxBar)-v/(1+NxBar)+log(gammainc(v/(1+1/NxBar),N-1));
    end

    if(ampDef==0)
        %Adjust for the change of variable to undo the scaling from
        %Equation 10.4-4 in [2].
        val=val-log(2);
    end
    
    val(v<0)=-Inf;
end

function val=PDF(v,avgSNR,N,ampDef)
%%PDF Evaluate the scalar probability density function (PDF) of the
%     distribution of the detection power (normalized to a noise variance
%     of 1) under a Swerling I model in a square-law detector.
%
%INPUTS: v The point or points at which the PDF should be evaluated. v>=0
%          for nonzero PDF values.
%   avgSNR The average power signal to noise ratio of the target.
%        N The number of pulses that are to be incoherently added for
%          detection (in a square-law detector). In a Swerling I model, a
%          single realization of the target power is used across all
%          pulses, but the noise corrupting the pulses varies. If this
%          parameter is omitted or an empty matrix is passed, N=1 is used.
%   ampDef Possible values are:
%          0 The expected value of the squared magnitude of the noise is 2.
%          1 (The default if omitted or an empty matrix is passed) The
%            expected value of the squared magnitude of the noise is 1.
%            (A more common definition). 
%
%OUTPUTS: val The values of the PDF evaluated at the given points.
%
%The PDF for a single pulse is given in Equation III.4 of [1] and the PDF
%for multiple pulses is Equation III.5 of [1] (Note the definition of the
%incomplete gamma function there differs from gammainc in Matlab). However,
%a clearer derivation is given in Chapter 11 of [2]. In both instances, the
%PDF derived is based on a scaled square-law detector, as in Equation
%10.4-4 in [2].
%
%EXAMPLE:
%Here, we validate the PDF by generating random samples and comparing the
%PDF plot with a histogram of the random samples.
% avgSNR=3;
% N=4;
% numSamples=10e3;
% ampDef=1;
% 
% figure(1)
% clf
% histogram(SwerlingISqLawD.rand([numSamples,1],avgSNR,N,ampDef),'Normalization','pdf')
% hold on
% numPoints=1000;
% x=linspace(0,80,numPoints);
% vals=SwerlingISqLawD.PDF(x,avgSNR,N,ampDef);
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

    if(ampDef==0)
        %Initial change of variable to undo the scaling from Equation
        %10.4-4 in [2].
        v=v/2;
    end

    if(N==1)
        %Equation III.4 in [1].
        val=(1/(1+avgSNR)).*exp(-v./(1+avgSNR));
    else
        NxBar=N*avgSNR;

        %Equation III.5 in [1]. The logarithms and exponentiation are used
        %to help reduce overflow errors. Note that gammainc is defined
        %differently than the incomplete gamma function used in [1].
        val=exp((N-2)*log(1+1/NxBar)-log(NxBar)-v/(1+NxBar)).*gammainc(v/(1+1/NxBar),N-1);
        %That is the same as
        %val=(1+1/NxBar).^(N-2).*(1/NxBar).*gammainc(v./((1+1/NxBar)),N-1).*exp(-v/(1+NxBar));
    end

    if(ampDef==0)
        %Adjust for the change of variable to undo the scaling from
        %Equation 10.4-4 in [2].
        val=val/2;
    end
    
    val(v<0)=0;
end

function val=CDF(v,avgSNR,N,ampDef)
%%CDF Evaluate the scalar cumulative distribution function (CDF) of the
%     distribution of the detection power (normalized to a noise variance
%     of 1) under a Swerling I model in a square-law detector.
%
%INPUTS: v The point or points at which the CDF should be evaluated. v>0
%          for nonzero CDF values.
%   avgSNR The average power signal to noise ratio of the target.
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
%
%OUTPUTS: val The values of the CDF evaluated at the given points.
%
%The derivation of the detection probability in [1] is as 1-the value of
%the CDF. Thus, this function just evaluated 
%1-SwerlingISqLawD.PD4Threshold(avgSNR,v,N,[],ampDef).
%
%EXAMPLE:
%Here, the CDF returned by this function is plotted along with an empirical
%CDF constructed from random samples of this distribution. They agree well.
% avgSNR=2;
% N=4;
% ampDef=1;
% numSamples=1e5;
% samples=SwerlingISqLawD.rand([numSamples,1],avgSNR,N,ampDef);
% numPts=1000;
% x=linspace(0,80,numPts);
% CDFEmp=EmpiricalD.CDF(x,samples);
% CDF=SwerlingISqLawD.CDF(x,avgSNR,N,ampDef);
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

    val=1-SwerlingISqLawD.PD4Threshold(avgSNR,v,N,[],ampDef);
    val(v<=0)=0;
end

function PD=PD4Threshold(avgSNR,thresh,N,method,ampDef)
%%PD4THRESHOLD Determine the detection probability of a Swerling I target
%           given its signal to noise ratio and the value of the
%           detection threshold, assuming that a square-law detector is
%           used. The noise in each sample is assumed to have variance 1
%           (normalized noise power).
%
%INPUTS: avgSNR A vector or matrix of average power signal to noise ratios
%          of the target at which one wishes to evaluate the detection
%          probability.
%   thresh The scalar normalized detection threshold to use. This is the
%          threshold to use if the noise SNR is 1 (subject to the
%          definition in ampDef).
%        N The number of pulses that are to be incoherently added for
%          detection (in a square-law detector). In a Swerling I model, a
%          single realization of the target power is used across all
%          pulses, but the noise corrupting the pulses varies. If this
%          parameter is omitted or an empty matrix is passed, N=1 is used.
%   method This specified whether an exact method or an approximation
%          should be used. For N=1, the exact method equals the
%          approximation. Possible values are
%          0 (The default if omitted or an empty matrix is passed) Use the
%            exact solution of Equation II.1 of [1].
%          1 Use the approximate solution of Equation II.2 of [1]. This is
%          only valid if N*avgSNR>1.
%   ampDef This specified normalization (see help SwerlingISqLawD).
%          Possible values are:
%          0 The expected value of the squared magnitude of the noise is 2.
%          1 (The default if omitted or an empty matrix is passed) The
%            expected value of the squared magnitude of the noise is 1.
%            (A more common definition). 
%
%OUTPUTS: PD The detection probability of the target.
%
%This function implements Equation II.1 and the approximation in II.2 in
%[1]. However, the implementation of Equation II.1 does not have the same
%form as in [1]. This is because in [1], Pearson's gamma function is used,
%which is defined differently than the gammainc function in Matlab.
%The value val of Pearson's incomplete gamma function of (u,p) is
%expressed in terms of gammainc as val=gammainc(u*sqrt(p+1),(p+1));
%Consequently, a number of simplifications were performed to arrive at the
%expression coded below.
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
% N=4;
% ampDef=0;
% PD=SwerlingISqLawD.PD4Threshold(avgSNR,thresh,N,[],ampDef)
% numSamples=1e5;
% PDSamp=mean(SwerlingISqLawD.rand([numSamples,1],avgSNR,N,ampDef)>=thresh)
%One will see that both PD values are about 0.268.
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
        case 0%The exact solution in II.1
            if(N==1)
                PD=exp(-thresh./(1+avgSNR));
            else
                gammaInc1=gammainc(thresh,N-1,'upper');
                if(abs(gammaInc1-1)<=eps())
                    %This avoids finite precision errors that can make the
                    %second gammainc term a NaN.
                    PD=1;
                else
                    PD=gammaInc1+...
                       ((1+1./(N*avgSNR)).^(N-1).*gammainc(thresh./(1+1./(N*avgSNR)),N-1).*exp(-thresh./(1+N*avgSNR)));
                end
            end
        case 1%The approximate solution in II.2
            NxBar=N*avgSNR;
            if(any(NxBar<1))
                warning('The approximate solution is not valid for N*avgSNR<=1')
            end

            PD=(1+1./NxBar).^(N-1).*exp(-thresh./(1+NxBar));
        otherwise
            error('Unknown method specified.')
    end
    
    %Deal with finite precision errors that might push PD slightly below 0
    %or above 1. This is primarily an issue with method 1.
    if(any(PD>1|PD<0))
        warning('Finite precision errors detected.') 
    end
    
    PD=min(PD,1);
    PD=max(PD,0);
end

function PD=PD4PFA(avgSNR,PFA,N,method,ampDef)
%%PD4PFA Determine the detection probability of a Swerling I target given
%        its signal to noise ratio and the probability of false alarm
%        implied by a particular unspecified detection threshold, assuming
%        that a square law detector is used.
%
%INPUTS: avgSNR A vector or matrix of average power signal to noise ratios
%          of the target at which one wishes to evaluate the detection
%          probability.
%      PFA The probability of false alarm, 0<=PFA<1.
%        N The number of pulses that are to be incoherently added for
%          detection (in a square-law detector). In a Swerling I model, a
%          single realization of the target power is used across all
%          pulses, but the noise corrupting the pulses varies. If this
%          parameter is omitted or an empty matrix is passed, N=1 is used.
%   method This specified whether an exact method or an approximation
%          should be used. For N=1, the exact method equals the
%          approximation. Possible values are
%          0 (The default if omitted or an empty matrix is passed) Use the
%            exact solution of Equation II.1 of [1].
%          1 Use the approximate solution of Equation II.2 of [1]. This is
%            only valid if N*avgSNR>1.
%   ampDef This specified normalization (see help SwerlingISqLawD).
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
%SwerlingISqLawD.PD4Threshold.
%
%REFERENCES:
%[1] P. Swerling, "Probability of detection for fluctuating targets," The
%    RAND Corporation, Santa Monica, CA, Tech. Rep. RM-1217, 1954.
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
    PD=SwerlingISqLawD.PD4Threshold(avgSNR,thresh,N,method,ampDef);
end

function avgSNR=avgSNR4PDThresh(PD,thresh,ampDef)
%%AVGSNR4PDTHRESH Given a detection probabilty and a normalized detection
%       threshold (normalized in terms of the receiver noise having a unit
%       covariance), determine the average power signal to noise ratio
%       (SNR) of the target needed under a Swerling I model. This solution
%       is only available for a single pulse (N=1) in a square-law
%       detector.
%
%INPUTS: PD The detection probability of the target , 0<=PD<1.
%   thresh The scalar normalized detection threshold to use. This is the
%          threshold to use if the noise variance is 1.
%   ampDef This specified normalization (see help SwerlingISqLawD).
%          Possible values are:
%          0 The expected value of the squared magnitude of the noise is 2.
%          1 (The default if omitted or an empty matrix is passed) The
%            expected value of the squared magnitude of the noise is 1.
%            (A more common definition). 
%
%OUTPUTS: avgSNR The average power signal to noise ratio needed to achieve
%                the desired probability of detection at the given
%                threshold. If this is negative, then for a given PD, the
%                threshold is so low that due to the high false alarm
%                rate, the desired PD is impossiblely low.
%
%Equation II.1 in [1] gives an expression for PD in terms of thresh and the
%average SNR. Here, we have simply inverted the expression for when N=1.
%
%However, the expressions used are derived based on a scaled square-law
%detector, as in Equation 10.4-4 in [2], which is the same model as in [1].
%Thus, the threshold in this function is scaled appropriately according to
%ampDef.
%
%EXAMPLE:
%Here, we show that the results are consistent.
%PD=0.75;
% thresh=5;
% ampDef=0;
% avgSNR=SwerlingISqLawD.avgSNR4PDThresh(PD,thresh,ampDef);
% PDBack=SwerlingISqLawD.PD4Threshold(avgSNR,thresh,1,[],ampDef)
% numSamples=1e5;
% PDSamp=mean(SwerlingISqLawD.rand([numSamples,1],avgSNR,1,ampDef)>=thresh)
%One will see that PDBack is the same as PD and about the same as PDSamp.
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

    if(ampDef==0)
        %The division by 2 deals with the scaling in Equation 10.4-4.
        avgSNR=-(thresh/2)/log(PD)-1;
    else
        avgSNR=-thresh/log(PD)-1;
    end
end
   
function avgSNR=avgSNR4PDPFA(PD,PFA,ampDef)
%%AVGSNR4PDTHRESH Given a detection probabilty and the probability of false
%       alarm, determine the average power signal to noise ratio (SNR) of
%       the target needed under a Swerling I model. This solution is only
%       available for a single pulse (N=1) in a square-law detector.
%
%INPUTS: PD The detection probability of the target , 0<=PD<1.
%      PFA The probability of false alarm, 0<=PFA<1.
%   ampDef This specified normalization (see help SwerlingISqLawD).
%          Possible values are:
%          0 The expected value of the squared magnitude of the noise is 2.
%          1 (The default if omitted or an empty matrix is passed) The
%            expected value of the squared magnitude of the noise is 1.
%            (A more common definition). 
%
%OUTPUTS: avgSNR The average power signal to noise ratio needed to achieve
%                the desired probability of detection at the given
%                threshold. If this is negative, then for a given PD, the
%                threshold is so low that due to the high false alarm
%                rate, the desired PD is impossiblely low.
%
%Equation II.1 in [1] gives an expression for PD in terms of a normalized
%threshold and the average SNR. Here, we call the function
%PFA2SquareLawThreshold and insert the result into
%SwerlingISqLawD.avgSNR4PDThresh, which inverts the expression in [1].
%
%REFERENCES:
%[1] P. Swerling, "Probability of detection for fluctuating targets," The
%    RAND Corporation, Santa Monica, CA, Tech. Rep. RM-1217, 1954.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.

    if(nargin<3||isempty(ampDef))
        ampDef=1;
    end

    thresh=PFA2SquareLawThreshold(PFA,1,ampDef);
    avgSNR=SwerlingISqLawD.avgSNR4PDThresh(PD,thresh,ampDef);
end

function thresh=thresh4AvgSNRPD(PD,avgSNR,ampDef)
%%THRESH4AVGSNRPD Given a detection probability and the average power
%           signal to noise ratio (SNR), determine the threshold needed
%           under a Swerling I model. This solution is only available for a
%           single-pulse in a square law detector.
%
%INPUTS: PD The detection probability of the target , 0<=PD<1.  
%    avgSNR A vector or matrix of average power signal to noise ratios
%           of the target at which one wishes to evaluate the threshold.
%   ampDef This specified normalization (see help SwerlingISqLawD).
%          Possible values are:
%          0 The expected value of the squared magnitude of the noise is 2.
%          1 (The default if omitted or an empty matrix is passed) The
%            expected value of the squared magnitude of the noise is 1.
%            (A more common definition). 
%
%OUTPUTS: thresh The normnalized detection threshold(s) (assuming the noise
%                variance is 1).
%
%Equation II.1 in [1] gives an expression for PD in terms of a normalized
%threshold and the average SNR. Here, we simply invert the equation.
%
%However, the expressions used are derived based on a scaled square-law
%detector, as in Equation 10.4-4 in [2]. Thus, the threshold in this
%function is scaled appropriately.
%
%EXAMPLE:
%Here, we show that the results are consistent.
% PD=0.6;
% avgSNR=2;
% ampDef=1;
% thresh=SwerlingISqLawD.thresh4AvgSNRPD(PD,avgSNR);
% PDBack=SwerlingISqLawD.PD4Threshold(avgSNR,thresh,1,[],ampDef)
% numSamples=1e5;
% PDSamp=mean(SwerlingISqLawD.rand([numSamples,1],avgSNR,1,ampDef)>=thresh)
%One will see that PDBack=0.6 and PDSamp is about the same.
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

    thresh=-(1+avgSNR).*log(PD);
    if(ampDef==0)
        %The multiplication by 2 deals with the scaling in Equation 10.4-4.
        thresh=thresh*2;
    end
end

function vals=rand(NDims,avgSNR,N,ampDef)
%%RAND Generate random variables representing a target whose sampled signal
%      to average noise power ratio is generated according to a Swerling I
%      model detected by a square-law detector.
%
%INPUTS: NDims If NDims is a scalar, then rand returns an NDimsXNDims
%          matrix of random variables. If NDims=[M,N1] is a two-element row
%          vector, then rand returns an MXN1 matrix of random variables.
%   avgSNR The average power signal to noise ratio of the target.
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
%
%OUTPUTS: vals A matrix whose dimensions are determined by N of the
%              generated random variables.
%
%Following the model of [1], a single sample from the exponential
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
        %Equation I.1 for  the input signal-to-noise power ratio. It will
        %be the same for all N samples that are accumulated.
        curSNR=ExponentialD.rand(1,1/(avgSNR));
        vals(curEl)=NonFlucSqLawD.rand(1,curSNR,N,ampDef,0);
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
