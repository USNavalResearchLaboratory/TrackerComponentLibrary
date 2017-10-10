classdef SwerlingIVSqLawD
%%SWERLINGIVSQLAWD This class implements various functions and detections
%                statistics related to the Swerling IV target model when
%                used with a square law detector. All methods are static
%                (the class holds no information and does not need to be
%                instantiated).
%
%Implemented methods are: mean, var, PDF, CDF, PD4Threshold, PD4PFA, rand
%
%Swerling models are given in [1] and in Chapter 11 of [2]. The Swerling IV
%model assumes that the observed power SNR of the target varies in terms of
%a central gamma distribution with shape parameter k=2 and scale parameter
%theta=avgSNR/2, that in a pulse train for detection the target SNR
%fluctuates from pulse to pulse, and that a square law detector is used to
%incoherently integrate multiple pulses The difference from a Swerling III
%model is the target amplitude fluctuations from pulse to pulse.

%The square law detector for N samples is
%y=sum_{i=1}^N r_i^2
%where r_i is the ith real amplitude. The model for the amplitude given a
%target amplitude or power is developed in Chapter 9 of [2]. The squared
%real amplitude is r^2_i=y_{I,i}^2+y_{Q,i}^2 where y_{I,i} and y_{Q,i} are
%the in-phase and quadrature components of the filter output. Define
%Rp=2*sampPowSNR, as in Equation 9.2-34. The sampled power SNR of the
%target in the Swerling IV model is different for each pulse. The sample
%power SNR is sampPowSNR=GammaD.rand(1,2,avgSNR/2) as mentioned in [1]. The
%value y=y_{I,i}+sqrt(-1)*y_{Q,i} conditioned on the sampled power SNR is
%modeled as being distributed circularly complex Gaussian with mean
%sqrt(Rp) and variance 2. (See, Equation 9.3-35a in [2]). The variance
%being 2 simply reflects having normalized the variance on the I and Q
%components each to 1.  As in Equation 10.4-13 of [2], the distribution of
%Y=y/2 conditioned on the target SNR values is noncentral chi-squared with
%a change of variables. The value y conditioned on the target SNR values is
%noncentral chi-squared with nu=2*N degrees of freedom and lambda=N*Rp as
%the noncentrality parameter. The expressions in [1] and [2] are in terms
%of Y and not y, so a transformation has to be performed to make things in
%terms of y.
%
%Under this model, for a given average power SNR value, one can generate a
%random sample from the distribution as
% sample=0;
% for curPulse=1:N
%     sampPowSNR=GammaD.rand(1,2,avgSNR/2);
%     Rp=2*sampPowSNR;%The peak power of the signal from the SNR, Equation
%                 %9.4-8 in [2].
%     A=sqrt(Rp);%Equation 9.4-8, the amplitude of the signal
%     sample=sample+sum(abs(ComplexGaussianD.rand(1,A,2)).^2);
% end
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
function val=mean(avgSNR,N)
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
%
%OUTPUTS: val The value of the mean for the given parameters.
%  
%The mean is given in Equation 11.5-20 of [1] where equation 11.4-8 relates
%the parameterization used there to an average squared amplitude. To
%relate the average squared amplitude ABar^2 to avgSNR as used in [2], we
%use avgSNR=ABar^2/2. However, the first moment from that expression has to
%be multiplied by 2, because as seen in Equation 10.4-4 of [1], the
%definition is in terms of a scaled version of the square law detector.
%
%EXAMPLE:
%Here, we validate the mean by generating random samples and comparing the
%computed mean to the sample mean.
% avgSNR=2;
% N=4;
% numSamples=1e5;
% meanVal=SwerlingIVSqLawD.mean(avgSNR,N)
% meanSampVal=mean(SwerlingIVSqLawD.rand([numSamples,1],avgSNR,N))
%One will see that both mean values are about 24.
%
%REFERENCES:
%[1] J. V. Di Franco and W. L. Rubin, Radar Detection. SciTech Publishing
%    Inc., Rayliegh, NC: 2004.
%[2] P. Swerling, "Probability of detection for fluctuating targets," The
%    RAND Corporation, Santa Monica, CA, Tech. Rep. RM-1217, 1954.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.

    if(nargin<2||isempty(N))
        N=1;
    end
    
    val=2*N*(1+avgSNR);
end

function val=var(avgSNR,N)
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
% numSamples=1e5;
% varVal=SwerlingIVSqLawD.var(avgSNR,N)
% varSampVal=var(SwerlingIVSqLawD.rand([numSamples,1],avgSNR,N))
%One will see that both variance values are about 112.
%
%REFERENCES:
%[1] J. V. Di Franco and W. L. Rubin, Radar Detection. SciTech Publishing
%    Inc., Rayliegh, NC: 2004.
%[2] P. Swerling, "Probability of detection for fluctuating targets," The
%    RAND Corporation, Santa Monica, CA, Tech. Rep. RM-1217, 1954.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.

    if(nargin<2||isempty(N))
        N=1;
    end
    
    val=4*N*(1+avgSNR.*(2+avgSNR/2));
end
    
function val=PDF(x,avgSNR,N)
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
%
%OUTPUTS: val The values of the PDF evaluated at the given points.
%
%The PDF is given in Equation 11.5-14 where equation 11.4-8 relates the
%parameterization used there to an average squared amplitude. To relate the
%average squared amplitude ABar^2 to avgSNR as used in [2], we use
%avgSNR=ABar^2/2. However the PDF is derived is based on a scaled
%square-law detector, as in Equation 10.4-4 in [1]. Thus, in implementing
%the PDF, the scaling is removed.
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
% histogram(SwerlingIVSqLawD.rand([numSamples,1],avgSNR,N),'Normalization','pdf')
% hold on
% numPoints=1000;
% x=linspace(0,80,numPoints);
% vals=SwerlingIVSqLawD.PDF(x,avgSNR,N);
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
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.

    if(nargin<3||isempty(N))
        N=1;
    end

    %Initial change of variable to undo the scaling from Equation 10.4-4 in
    %[1].
    x=x/2;

    xBar2=avgSNR/2;
    val=zeros(size(x));
    numEls=numel(x);

    for curEl=1:numEls
        v=x(curEl);

        if(v>0)
            val(curEl)=exp((N-1)*log(v)-v/(1+xBar2)-2*N*log(1+xBar2)-gammaln(N))*hypergeometric1F1(-N,N,-xBar2/(1+xBar2)*v);
        elseif(v==0&&N==1)
            val(curEl)=1/((1+xBar2)^(2*N)*factorial(N-1));
        end
    end
    
    %Adjust for the change of variable to undo the scaling from Equation
    %10.4-4 in [1].
    val=val/2;
end

function val=CDF(v,avgSNR,N)
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
%
%OUTPUTS: val The values of the CDF evaluated at the given points.
%
%The derivation of the detection probability in [1] is as 1-the value of
%the CDF. Thus, this function just evaluated 
%1-SwerlingIVSqLawD.PD4Threshold(avgSNR,v,N).
%
%REFERENCES:
%[1] P. Swerling, "Probability of detection for fluctuating targets," The
%    RAND Corporation, Santa Monica, CA, Tech. Rep. RM-1217, 1954.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.

    if(nargin<3||isempty(N))
        N=1; 
    end

    val=1-SwerlingIVSqLawD.PD4Threshold(avgSNR,v,N);
    val(v<=0)=0;
end
    
function PD=PD4Threshold(avgSNR,thresh,N)
%%PD4THRESHOLD Determine the detection probability of a Swerling IV target
%           given its signal to noise ratio and the value of the
%           detection threshold, assuming that a square-law detector is
%           used. The noise in each sample is assumed to be Gaussian with
%           variance 1 (normalized noise power).
%
%INPUTS: avgSNR A vector or matrix of average power signal to noise ratios
%               of the target at which one wishes to evaluate the detection
%               probability.
%        thresh The scalar normalized detection threshold to use. This is
%               the threshold to use if the noise variance is 1.
%             N The number of pulses that are to be incoherently added for
%               detection (in a square-law detector). In a Swerling IV
%               model, a different realization of the target power is used
%               across each pulse and the noise corrupting the pulses
%               varies. If this parameter is omitted or an empty matrix is
%               passed, N=1 is used.
%
%OUTPUTS: PD The detection probability of the target.
%
%This function implements the sum in Equations 11.5-19 in [1]. However, the
%expression used is derived based on a scaled square-law detector, as in
%Equation 10.4-4. Thus, the threshold in this function is scaled
%appropriately.
%
%EXAMPLE:
%Here, we validate the results by comparing the PD from this function to
%the PD computed using random samples.
% avgSNR=2;
% thresh=30;
% N=4;
% PD=SwerlingIVSqLawD.PD4Threshold(avgSNR,thresh,N)
% numSamples=1e5;
% PDSamp=mean(SwerlingIVSqLawD.rand([numSamples,1],avgSNR,N)>=thresh)
%One will see that both PD values are about 0.254.
%
%REFERENCES:
%[1] J. V. Di Franco and W. L. Rubin, Radar Detection. SciTech Publishing
%    Inc., Rayliegh, NC: 2004.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.

    if(nargin<3||isempty(N))
       N=1; 
    end
    
    thresh=thresh/2;
    
    R=2*avgSNR;
    
    logR4=log(R/4);
    logCoeff=gammaln(N+1)-N*log(1+R/4);
    sumVal=0;
    for k=0:N
        curTerm=exp(logCoeff+k*logR4+log(gammainc(thresh/(1+R/4),k+N))-gammaln(k+1)-gammaln(N-k+1));

        sumVal=sumVal+curTerm;
    end

    PD=1-sumVal;
end

function PD=PD4PFA(avgSNR,PFA,N)
%%PD4PFA Determine the detection probability of a Swerling IV target given
%        its signal to noise ratio and the probability of false alarm
%        implied by a particular unspecified detection threshold, assuming
%        that a square law detector is used.
%
%INPUTS: avgSNR A vector or matrix of average power signal to noise ratios
%               of the target at which one wishes to evaluate the detection
%               probability.
%           PFA The probability of false alarm, 0<=PFA<1.
%             N The number of pulses that are to be incoherently added for
%               detection (in a square-law detector). In a Swerling IV
%               model, a different realization of the target power is used
%               across each pulse and the noise corrupting the pulses
%               varies. If this parameter is omitted or an empty matrix is
%               passed, N=1 is used.
%
%OUTPUTS: PD The detection probability of the target, 0<=PD<=1.
%
%This function just calls PFA2SquareLawThreshold to convert the false alarm
%rate into a normalized detection threshold after which it calls
%SwerlingIVSqLawD.PD4Threshold.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.

    if(nargin<3||isempty(N))
       N=1; 
    end

    thresh=PFA2SquareLawThreshold(PFA,N);
    PD=SwerlingIVSqLawD.PD4Threshold(avgSNR,thresh,N);
end

function vals=rand(NDims,avgSNR,N)
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
            vals(curEl)=vals(curEl)+NonFlucSqLawD.rand(1,curSNR,1);
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
