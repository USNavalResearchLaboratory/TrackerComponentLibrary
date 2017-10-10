classdef NonFlucSqLawD
%%NONFLUCMODEL This class implements various function and detection
%              statistics for targets with constant average power signal to
%              noise ratio. All methods are static (the class holds no
%              information and does not need to be instantiated).
%
%Implemented methods are: mean, var, PDF, CDF, PD4Threshold, PD4PFA, rand
%
%This model considers detection of a non-fluctuating target with a constant
%average signal power to noise ratio via non-coherent integration of
%multiple pulses in a square-law detector. The model is discussed in
%Chapter 10 of [1].
%
%In Chapter 10.4 of [1], it can be seen in Equation 10.4-13 that the PDF of
%the output of a square-law detector is related to the noncentral
%chi-square distribution. This class essentially maps the results to those
%in ChiSquareD with the appropriate parameterization.
%
%The square law detector for N samples is
%y=sum_{i=1}^N r_i^2
%where r_i is the ith real amplitude. The model for the amplitude is
%developed in Chapter 9 of [1]. The model comes from taking a sample of a
%complex amplitude coming from a filter. The squared real amplitude is
%r^2_i=y_{I,i}^2+y_{Q,i}^2 where y_{I,i} and y_{Q,i} are the in-phase and
%quadrature components of the filter output. Defining Rp=2*avgSNR as in
%Equation 9.2-34, the value y=y_{I,i}+sqrt(-1)*y_{Q,i} is modeled as being
%distributed circularly complex Gaussian with mean sqrt(Rp) and variance
%2. (See, Equation 9.3-35a in [1]). The variance being 2 simply reflects
%having normalized the variance on the I and Q components each to 1. 
%As in Equation 10.4-13, the distribution of Y=y/2 is noncentral
%chi-squared with a change of variables. The value y is noncentral chi-
%squared with nu=2*N degrees of freedom and lambda=N*Rp as the
%noncentrality parameter.
%
%Under this model, for a given average power SNR value, one can generate a
%random sample from the distribution as
% Rp=2*avgSNR;%The peak power of the signal from the SNR, Equation 9.4-8.
% A=sqrt(Rp);%Equation 9.4-8, the amplitude of the signal
% sample=sum(abs(ComplexGaussianD.rand(N,A,2)).^2);
%
%REFERENCES:
%[1] J. V. Di Franco and W. L. Rubin, Radar Detection. SciTech Publishing
%    Inc., Rayliegh, NC: 2004.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.
    
methods(Static)

function val=mean(avgSNR,N)
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
%
%OUTPUTS: val The value of the mean for the given parameters.
%
%The PDF is given by Equation 10.4-13 in [1]. However, a change of
%variables is performed to account for the values having been defined in
%terms of a scaled version of the square law in Equation 10.4-4. The final
%result is a noncentral chi-square distribution and this function just
%calls ChiSquareD.mean with the appropriate parameters.
%
%EXAMPLE:
%Here, we validate the mean by generating random samples and comparing the
%computed mean to the sample mean. 
% avgSNR=2;
% N=4;
% numSamples=1e5;
% meanVal=NonFlucSqLawD.mean(avgSNR,N)
% meanSampVal=mean(NonFlucSqLawD.rand([numSamples,1],avgSNR,N))
%One will see that both values are about 24.
%
%REFERENCES:
%[1] J. V. Di Franco and W. L. Rubin, Radar Detection. SciTech Publishing
%    Inc., Rayliegh, NC: 2004.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.

    if(nargin<2||isempty(N))
        N=1;
    end
    
    Rp=avgSNR*2;
    
    %Equation 10.4-13 in [1] is a noncentral Chi-Square distribution.
    nu=2*N;
    lambda=N*Rp;
    val=ChiSquareD.mean(nu,lambda);
end

function val=var(avgSNR,N)
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
%
%OUTPUTS: val The value of the variance for the given parameters.
%
%The PDF is given by Equation 10.4-13 in [1]. However, a change of
%variables is performed to account for the values having been defined in
%terms of a scaled version of the square law in Equation 10.4-4. The final
%result is a noncentral chi-square distribution and this function just
%calls ChiSquareD.var with the appropriate parameters.
%
%EXAMPLE:
%Here, we validate the variance by generating random samples and comparing
%the computed variance to the sample variance. 
% avgSNR=2;
% N=4;
% numSamples=1e5;
% varVal=NonFlucSqLawD.var(avgSNR,N)
% varSampVal=var(NonFlucSqLawD.rand([numSamples,1],avgSNR,N))
%One will see that both values are about 80.
%
%REFERENCES:
%[1] J. V. Di Franco and W. L. Rubin, Radar Detection. SciTech Publishing
%    Inc., Rayliegh, NC: 2004.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
    
    if(nargin<2||isempty(N))
        N=1;
    end
    
    Rp=avgSNR*2;
    
    %Equation 10.4-13 in [1] is a noncentral Chi-Square distribution.
    nu=2*N;
    lambda=N*Rp;
    val=ChiSquareD.var(nu,lambda);
end

function val=PDF(x,avgSNR,N)
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
%
%OUTPUTS: val The values of the PDF evaluated at the given points.
%
%The PDF is given by Equation 10.4-13 in [1]. However, a change of
%variables is performed to account for the values having been defined in
%terms of a scaled version of the square law in Equation 10.4-4. The final
%result is a noncentral chi-square distribution and this function just
%calls ChiSquareD.PDF with the appropriate parameters.
%
%EXAMPLE:
%Here, we validate the PDF by generating random samples and comparing the
%PDF plot with a histogram of the random samples.
% avgSNR=2;
% N=4;
% numSamples=10e3;
% 
% figure(1)
% clf
% histogram(NonFlucSqLawD.rand([numSamples,1],avgSNR,N),'Normalization','pdf')
% hold on
% numPoints=1000;
% x=linspace(0,60,numPoints);
% vals=NonFlucSqLawD.PDF(x,avgSNR,N);
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

    if(nargin<3||isempty(N))
        N=1;
    end

    Rp=avgSNR*2;
    
    %Equation 10.4-13 in [1] is a noncentral Chi-Square distribution.
    nu=2*N;
    lambda=N*Rp;
    val=ChiSquareD.PDF(x,nu,lambda);
end

function val=CDF(v,avgSNR,N)
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
%
%OUTPUTS: val The values of the CDF evaluated at the given points.
%
%The derivation of the detection probability in Chapter 10 of [1] is as
%1-the value of the CDF. Thus, this function just evaluated 
%1-NonFlucSqLawD.PD4Threshold(avgSNR,v,N).
%
%REFERENCES:
%[1] J. V. Di Franco and W. L. Rubin, Radar Detection. SciTech Publishing
%    Inc., Rayliegh, NC: 2004.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.

    if(nargin<3||isempty(N))
        N=1; 
    end

    val=1-NonFlucSqLawD.PD4Threshold(avgSNR,v,N);
    val(v<=0)=0;
end

function PD=PD4Threshold(avgSNR,thresh,N)
%%PD4THRESHOLD Determine the detection probability of a target with a
%           constant signal to noise ratio (SNR) given the SNR and the
%           value of the detection threshold, assuming that a square-law
%           detector is used. The noise in each sample is assumed to be
%           Gaussian with variance 1 (normalized noise power).
%
%INPUTS: avgSNR A vector or matrix of average power signal to noise ratios
%               of the target at which one wishes to evaluate the detection
%               probability.
%        thresh The scalar normalized detection threshold to use. This is
%               the threshold to use if the noise variance is 1.
%             N The number of pulses that are to be incoherently added for
%               detection (in a square-law detector). In a nonfluctuating
%               model, the target power is constant and the noise
%               corrupting the pulses varies. If this parameter is omitted
%               or an empty matrix is passed, N=1 is used.
%
%OUTPUTS: PD The detection probability of the target.
%
%This function implements Equations 10.4-30 in [1], which expresses the
%result in terms of an incomplete Toronto function. In [1], it is shown
%that the incomplete Toronto function parameterized T_B(m,(m-1)/2,r) is
%expressed in terms fo the MarcumQ function as
%1-MarcumQ((m+1)/2,r*sqrt(2),B*sqrt(2)). This is the form in this problem.
%Thus, the MarcumQ function is used to evaluate the detection probability.
%
%EXAMPLE:
%Here, we validate the results by comparing the PD from this function to
%the PD computed using random samples.
% avgSNR=2;
% thresh=30;
% N=4;
% PD=NonFlucSqLawD.PD4Threshold(avgSNR,thresh,N)
% numSamples=1e5;
% PDSamp=mean(NonFlucSqLawD.rand([numSamples,1],avgSNR,N)>=thresh)
%One will see that PD and PDSamp are both near 0.2329.
%
%REFERENCES:
%[1] J. V. Di Franco and W. L. Rubin, Radar Detection. SciTech Publishing
%    Inc., Rayliegh, NC: 2004.
%[2] P. C. Sofotasios and S. Freear, "New analytic results for the
%    incomplete Toronto function and imcomplete Lipschitz-Hankel
%    integrals," in Proceedings from the SMBO/IEEE MTT-S International
%    Microwave and Optoelectronics Conference, Natal, Brazil, 29 Oct. - 1
%    Nov. 2011, pp. 44-47.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.

    if(nargin<3||isempty(N))
       N=1; 
    end

    %Account for Equation 10.4-4 in [1] using a halved version of the
    %square law.
    thresh=thresh/2;
    
    B=sqrt(thresh);
    m=2*N-1;

    PD=zeros(size(avgSNR));
    numEls=numel(PD);
    
    for curEl=1:numEls
        R=2*avgSNR(curEl);
        r=sqrt(N*R/2);

        PD(curEl)=MarcumQ((m+1)/2,r*sqrt(2),B*sqrt(2));
    end
end

function PD=PD4PFA(avgSNR,PFA,N)
%%PD4PFA Determine the detection probability of a target with a constant
%        signal to noise ratio (SNR) given the SNR and the probability of
%        false alarm, assuming that a square-law detector is used. The
%        noise in each sample is assumed to be Gaussian with variance 1
%        (normalized noise power).
%
%INPUTS: avgSNR A vector or matrix of average power signal to noise ratios
%               of the target at which one wishes to evaluate the detection
%               probability.
%           PFA The probability of false alarm, 0<=PFA<1.
%             N The number of pulses that are to be incoherently added for
%               detection (in a square-law detector). In a nonfluctuating
%               model, the target power is constant and the noise
%               corrupting the pulses varies. If this parameter is omitted
%               or an empty matrix is passed, N=1 is used.
%
%OUTPUTS: PD The detection probability of the target, 0<=PD<=1.
%
%This function just calls PFA2SquareLawThreshold to convert the false alarm
%rate into a normalized detection threshold after which it calls
%NonFlucSqLawD.PD4Threshold.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.

    if(nargin<3||isempty(N))
       N=1; 
    end
    
    thresh=PFA2SquareLawThreshold(PFA,N);
    PD=NonFlucSqLawD.PD4Threshold(avgSNR,thresh,N);
end

function vals=rand(NDims,avgSNR,N)
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
%
%OUTPUTS: vals A matrix whose dimensions are determined by N of the
%              generated random variables.
%
%Though the random variables can be generated by summing complex Gaussian
%variables, that is not efficient for large values of N. Thus, the
%Chi-squared distribution is used directly. This is given by Equation
%10.4-13 in [1]. However, a change of variables is performed to account for
%the values having been defined in terms of a scaled version of the square
%law in Equation 10.4-4.
%
%REFERENCES:
%[1] J. V. Di Franco and W. L. Rubin, Radar Detection. SciTech Publishing
%    Inc., Rayliegh, NC: 2004.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.

    if(nargin<3||isempty(N))
        N=1; 
    end

    if(isscalar(NDims))
        dims=[NDims, NDims];
    else
        dims=NDims;
    end
    Rp=avgSNR*2;
    
    %Equation 10.4-13 in [1] is a noncentral Chi-Square distribution.
    %It is used here to generate the samples.
    nu=2*N;
    lambda=N*Rp;
    vals=ChiSquareD.rand(dims,nu,lambda);
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
