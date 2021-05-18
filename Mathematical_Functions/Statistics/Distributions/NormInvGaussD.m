classdef NormInvGaussD
%%NORMINVGAUSSD Function related to the multivariate normal inverse
%    Gaussian distribution. This distribution has been used in applications
%    related to synthetic aperture radar (SAR), as in [1].
%Implemented methods are: mean, cov, PDF, rand
%
%REFERENCES:
%[1] T. A. Øigård, A. Hanssen, and R. E. Hansen, "The multivariate normal
%    inverse Gaussian distribution: EM-estimation and analysis of synthetic
%    aperture sonar data," in Proceedings of the 12th European Signal
%    Processing Conference, Vienna, Austria, 6-10 Sep. 2004, pp. 1433-1436.
%
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

methods(Static)
    function val=mean(mu,Gamma,beta,alpha,delta)
    %%MEAN Obtain the mean of the normal inverse Gaussian distribution.
    %
    %INPUTS:   mu The numDimX1 additive constant in the distribution.
    %       Gamma A positive semidefinite symmetric numDimXnumDim matrix
    %             such that det(Gamma)=1.
    %        beta A numDimX1 vector affecting the mean of the conditional
    %             Gaussian term.
    %       alpha A positive scalar parameter affecting the scale of the
    %             inverse Gaussian term. It is required that
    %             alpha^2>beta'*Gamma*beta
    %       delta A positive scalar parameter affecting the mean of the
    %             inverse Gaussian term.
    %
    %OUTPUTS: val The numDImX1 mean of the distribution.
    %
    %The mean of the distribution is Equation 3 in [1].
    %
    %REFERENCES:
    %[1] T. A. Øigård and A. Hanssen, "The multivariate normal inverse
    %    Gaussian heavy-tailed distribution: Simulation and estimation,"
    %    in Proceedings of the IEEE International Conference on Acoustics,
    %    Speech, and Signal Processing, vol. 2, Orlando¡ FL, 13-17 May
    %    2002, pp. 1489-1492.
    %
    %June 2016 David F. Crouse, Naval Research Laboratory, Washington D.C. 
    
        %Equation 3
        val=mu+delta*Gamma*beta/sqrt(alpha^2-beta'*Gamma*beta);
    end
    
    function val=cov(Gamma,beta,alpha,delta)
    %%MEAN Obtain the covariance matrix of the normal inverse Gaussian
    %      distribution.
    %
    %INPUTS:Gamma A positive semidefinite symmetric numDimXnumDim matrix
    %             such that det(Gamma)=1.
    %        beta A numDimX1 vector affecting the mean of the conditional
    %             Gaussian term.
    %       alpha A positive scalar parameter affecting the scale of the
    %             inverse Gaussian term. It is required that
    %             alpha^2>beta'*Gamma*beta
    %       delta A positive scalar parameter affecting the mean of the
    %             inverse Gaussian term.
    %
    %OUTPUTS: val The numDimXnumDim covariance matrix of the distribution.
    %
    %REFERENCES:
    %[1] T. A. Øigård and A. Hanssen, "The multivariate normal inverse
    %    Gaussian heavy-tailed distribution: Simulation and estimation,"
    %    in Proceedings of the IEEE International Conference on Acoustics,
    %    Speech, and Signal Processing, vol. 2, Orlando¡ FL, 13-17 May
    %    2002, pp. 1489-1492.
    %
    %June 2016 David F. Crouse, Naval Research Laboratory, Washington D.C. 
    
        alphaTerm=alpha^2-beta'*Gamma*beta;
        %Equation 4
        val=delta*sqrt(alphaTerm)*(Gamma+(1/alphaTerm)*Gamma*(beta*beta')*Gamma');
    end
    
    function vals=PDF(x,mu,Gamma,beta,alpha,delta)
    %%PDF Evaluate the probability distribution function (PDF) of the
    %     normal inverse Gassuain distribution.
    %
    %INPUTS: x The numDimXnumPoints set of points at which the PDF should
    %          be evaluated.
    %       mu The numDimX1 additive constant in the distribution.
    %    Gamma A positive semidefinite symmetric numDimXnumDim matrix such
    %          that det(Gamma)=1.
    %     beta A numDimX1 vector affecting the mean of the conditional
    %          Gaussian term.
    %    alpha A positive scalar parameter affecting the scale of the
    %          inverse Gaussian term. It is required that
    %          alpha^2>beta'*Gamma*beta
    %    delta A positive scalar parameter affecting the mean of the
    %          inverse Gaussian term.
    %
    %OUTPUTS: vals The scalar values of the PDF. If multiple points are
    %              passed (x is a matrix), then vals is a row vector.
    %
    %The PDF is given in Equation 1 in [1].
    %   
    %REFERENCES:
    %[1] T. A. Øigård and A. Hanssen, "The multivariate normal inverse
    %    Gaussian heavy-tailed distribution: Simulation and estimation,"
    %    in Proceedings of the IEEE International Conference on Acoustics,
    %    Speech, and Signal Processing, vol. 2, Orlando¡ FL, 13-17 May
    %    2002, pp. 1489-1492.
    %
    %June 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.  
    
        numDim=size(x,1);
        N=size(x,2);
        vals=zeros(1,N);

        for curVal=1:N
            diff=x(:,curVal)-mu;

            pTilde=delta*sqrt(alpha^2-beta'*Gamma*beta)+beta'*diff;
            qTilde=sqrt(delta^2+diff'*pinv(Gamma)*diff);
            vals(curVal)=(delta/(2^((numDim-1)/2)))*(alpha/(pi*qTilde))^((numDim+1)/2)*exp(pTilde)*besselk((numDim+1)/2,alpha*qTilde);
        end
    end
    
    function X=rand(N,mu,Gamma,beta,alpha,delta)
    %%RAND Generate normal inverse Gaussian random vectors.
    %
    %INPUTS: N The number of random variables to generate.
    %       mu The numDimX1 additive constant in the distribution.
    %    Gamma A positive semidefinite symmetric numDimXnumDim matrix such
    %          that det(Gamma)=1.
    %     beta A numDimX1 vector affecting the mean of the conditional
    %          Gaussian term.
    %    alpha A positive scalar parameter affecting the scale of the
    %          inverse Gaussian term. It is required that
    %          alpha^2>beta'*Gamma*beta
    %    delta A positive scalar parameter affecting the mean of the
    %          inverse Gaussian term.
    %
    %OUTPUT: X An xDimXN matrix of random instances of the normal inverse
    %          Gaussian distribution.
    %
    %As in [1], a normal inverse Gaussian random variable is expressed as
    %X=mu+Z*Gamma*beta+sqrt(Z)*cholSemiDef(Gamma,'lower')*Y
    %where Z is an inverse Gaussian distributed random variable with mean
    %delta^2 and shape parameter alpha^2-beta'*Gamma*beta and Y is a
    %standard numDim dimensional Gaussian random variable.
    %
    %EXAMPLE 1:
    %Here, we demonstrate that a histogram of many random samples will tend
    %to agree with the PDF. We choose a 1D distribution for this example.
    % alpha=12;
    % delta=2.3;
    % Gamma=1;
    % beta=4;
    % mu=7;
    % numSamp=10000;
    % samp=NormInvGaussD.rand(numSamp,mu,Gamma,beta,alpha,delta);
    % figure(1)
    % clf
    % hold on
    % h=histogram(samp,'Normalization','pdf');
    % %scatter(samp,zeros(numSamp,1),'ob')
    % %We will plot the PDF.
    % numPoints=500;
    % x=linspace(4,12,numPoints);
    % PDFVals=NormInvGaussD.PDF(x,mu,Gamma,beta,alpha,delta);
    % plot(x,PDFVals,'-r','linewidth',2)
    %
    %EXAMPLE 2:
    %This is the same as Exmaple 1, except we choose a 2D distribution.
    % alphaVal=12;
    % delta=2.3;
    % Gamma=[2/sqrt(7),1/sqrt(7);
    %        1/sqrt(7),4/sqrt(7)];%Has determinant=1.
    % beta=[-1;2];
    % mu=[0;-1];
    % numSamp=10000;
    % samp=NormInvGaussD.rand(numSamp,mu,Gamma,beta,alphaVal,delta);
    % figure(1)
    % clf
    % hold on
    % h=histogram2(samp(1,:),samp(2,:),'Normalization','pdf');
    % %scatter(samp,zeros(numSamp,1),'ob')
    % %We will plot the PDF.
    % numPoints=200;
    %  x1=linspace(-2,2,numPoints);
    % x2=linspace(-3,3,numPoints);
    % [X1,X2]=meshgrid(x1,x2);
    % x=[X1(:)';X2(:)'];
    % PDFVals=NormInvGaussD.PDF(x,mu,Gamma,beta,alphaVal,delta);
    % hs=surf(X1,X2,reshape(PDFVals,numPoints,numPoints),'EdgeColor','none');
    % alpha(hs,0.8);
    % view(45,45)
    %
    %REFERENCES:
    %[1] T. A. Øigård, A. Hanssen, and R. E. Hansen, "The multivariate
    %    normal inverse Gaussian distribution: EM-estimation and analysis
    %    of synthetic aperture sonar data," in Proceedings of the XII.
    %    European Signal Processing Conference, Vienna, Austria, 6-10 Sep.
    %    2004, pp. 1433-1436.
    %
    %June 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
    
        numDim=size(mu,1);
        X=zeros(numDim,N);

        SGamma=cholSemiDef(Gamma,'lower');

        lambdaZ=delta;
        muZ=delta/sqrt(alpha^2-beta'*Gamma*beta);
        for curVal=1:N
            Z=InvGaussianD.rand(1,muZ,lambdaZ);
            Y=randn(numDim,1);

            X(:,curVal)=mu+Z*Gamma*beta+sqrt(Z)*SGamma*Y;
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
