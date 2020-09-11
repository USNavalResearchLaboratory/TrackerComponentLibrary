classdef AlphaStableD
%%ALPHASTABLED Functions to handle univariate stable distributions
%Implemented methods are: PDF,CDF,rand (all for univariate case)
%
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.


methods(Static)
    function pdf = PDF(X,alpha,beta,gam,delta)
    %%PDF Obtain values of the alpha-stable probability density
    %     function with parameters alpha, beta, gam, and delta 
    %     evaluated at points in X. The values in X will be transformed
    %     to the corresponding range values for the standard 
    %     distribution and evaluating the pdf at those points. If only 
    %     X, alpha, and some of the input parameters are given,
    %     defaults of 0, 1, and 0 are assigned to beta, gam, and delta
    %     respectively.
    %
    %     For alpha < 0.001, the gamma function call below, 
    %     gamma(1/alpha), will return inf and therefore pdf(X==zet)
    %     will equal inf.
    %
    %     Note that for the alpha = 2 case (a Gaussian distribution) 
    %     we can relate the variance (sigma) to the scaling
    %     parameter (gam) by gam = sqrt(sigma/2) as shown in example 1.
    %     In the case of alpha = 1 (Cauchy distribution), we have that
    %     the scaling parameter of both distributions are equal as
    %     shown in example 2.
    %
    %INPUTS: X: A vector of points at which the pdf should be 
    %           evaluated.
    %    alpha: The scalar characteristic parameter of the pdf.
    %     beta: The scalar skewness parameter of the pdf.
    %      gam: The scalar dispersion parameter of the pdf.
    %    delta: The scalar location parameter of the pdf (also the
    %           mode of the distribution).
    %
    %OUTPUTS: pdf: A vector of the same size as X which contains the
    %              corresponding pdf values for a S-alpha-S
    %              distribution with parameters alpha, beta, gam, and
    %              delta.
    %
    %EXAMPLE 1: Comparing Gaussian pdf to alpha=2 S-alpha-S pdf
    % x = linspace(-10,10);
    % gpdf = GaussianD.PDF(x,0,3);
    % pdf = AlphaStableD.PDF(x,2,0,sqrt(3/2),1);
    % plot(x,gpdf,x,pdf)
    %
    %EXAMPLE 2: Comparing Cauchy pdf to alpha=1 S-alpha-S pdf
    % x = linspace(-10,10);
    % cpdf = CauchyD.PDF(x,0,2);
    % pdf = AlphaStableD.PDF(x,1,0,2,1);
    % plot(x,cpdf,x,pdf)
    %
    %EXAMPLE 3: Comparing Lévy pdf to alpha=0.5, beta=1 stable pdf
    % x = linspace(-10,10,200);
    % lpdf = LevyD.PDF(x,0,2);
    % pdf = AlphaStableD.PDF(x,0.5,1,2,1);
    % plot(x,lpdf,x,pdf)
    %
    %REFERENCES:
    %[1] Guillermo Julián-Moreno, Jorge E. López De Vergara, Iván
    %    González, Luis Pedro, Javier Royuela-Del-Val, and Federico 
    %    Simmross-Wattenberg, "Fast parallel alpha-stable distribution
    %    function evaluation and parameter estimation using OpenCL in 
    %    GPGPUs," Statistics and Computing 27 (2017), no. 5, 1365-1382.
    %
    %August 2019 Codie T. Lewis, Naval Research Laboratory, Washington D.C.
        
        %Default to standard distribution
        if(nargin<2)
            error('Must specify the characteristic parameter alpha.')
        end
        if(nargin<3)
            beta = 0;
            gam = 1;
            delta = 0;
        end
        if(nargin<4)
            gam = 1;
            delta = 0;
        end
        if(nargin<5)
            delta = 0;
        end
        if(alpha>2||alpha<=0)
            error('alpha must be in the interval (0,2]')
        end
        if(beta<-1||beta>1)
            error('beta must be in the interval [-1,1]')
        end
        if(gam<=0)
            error('gam must be greater than 0')
        end
        
        if(alpha==1)
            %Reparameterize the location parameter
            delta = delta+beta*log(gam)*2*gam/pi;
            %Apply scaling by gam and shifting location to delta
            X = (X-delta)/gam;
            if(beta==0)
                pdf = 1./(gam*pi*(1+X.^2));
            else
                V = @(theta) (2/pi)*((pi/2+beta*theta)./cos(theta)).*exp((1/beta)*(pi/2+beta*theta).*tan(theta));
                g = @(theta,x) exp(-pi*x/(2*beta)).*V(theta);
                h = @(theta,x) V(theta).*exp(-g(theta,x));
                pdf = (1/(2*abs(beta)))*exp(-pi*X/(2*beta)).*arrayfun(@(x)integral1DAdaptive(@(theta)h(theta,x),[-pi/2,pi/2]),X);
            end
        else
            %Reparameterize the location parameter
            delta = delta+beta*tan(alpha*pi/2)*gam;
            %Apply scaling by gam and shifting location to delta
            X = (X-delta)/gam;
            pdf = zeros(size(X));
            zet = -beta*tan(pi*alpha/2);
            Xpos = X;
            Xpos(X<zet) = -Xpos(X<zet);
            
            theta0above = (1/alpha)*atan(beta*tan(pi*alpha/2));
            V = @(theta) (cos(alpha*theta0above))^(1/(alpha-1))*(cos(theta)./sin(alpha*(theta0above+theta))).^(alpha/(alpha-1)).*(cos(alpha*theta0above+(alpha-1)*theta)./cos(theta));
            g = @(theta,x) V(theta).*(x-zet).^(alpha/(alpha-1));
            h = @(theta,x) g(theta,x).*exp(-g(theta,x));
            pdf(X>zet) = arrayfun(@(x)integral1DAdaptive(@(theta)h(theta,x),[-theta0above,pi/2]),Xpos(X>zet));
            pdf(X>zet) = (alpha./(pi*(Xpos(X>zet)-zet)*abs(alpha-1))).*pdf(X>zet);
            
            theta0below = (1/alpha)*atan(-beta*tan(pi*alpha/2));
            V = @(theta) (cos(alpha*theta0below))^(1/(alpha-1))*(cos(theta)./sin(alpha*(theta0below+theta))).^(alpha/(alpha-1)).*(cos(alpha*theta0below+(alpha-1).*theta)./cos(theta));
            g = @(theta,x) V(theta).*(x+zet).^(alpha/(alpha-1));
            h = @(theta,x) g(theta,x).*exp(-g(theta,x));
            pdf(X<zet) = arrayfun(@(x)integral1DAdaptive(@(theta)h(theta,x),[-theta0below,pi/2]),Xpos(X<zet));
            pdf(X<zet) = (alpha./(pi*(Xpos(X<zet)+zet)*abs(alpha-1))).*pdf(X<zet);
            
            pdf(X==zet) = (1/alpha)*gamma(1/alpha)*cos(theta0above)/(pi*(1+zet^2)^(1/(2*alpha)));
            pdf = pdf/gam;
        end
    end
    
    function cdf = CDF(X,alpha,beta,gam,delta)
    %%CDF Obtain values of the stable cumulative 
    %     density function with parameters alpha, beta, gam, and delta 
    %     evaluated at points in X. The values in X will be transformed
    %     to the corresponding range values for the standard 
    %     distribution and evaluating the pdf at those points. If only 
    %     X, alpha, and some of the input parameters are given,
    %     defaults of 0, 1, and 0 are assigned to beta, gam, and delta
    %     respectively.
    %
    %     Note that for the alpha = 2 case (a Gaussian distribution) 
    %     we can relate the variance (sigma) to the scaling
    %     parameter (gam) by gam = sqrt(sigma/2) as shown in example 1.
    %     In the case of alpha = 1 (Cauchy distribution), we have that
    %     the scaling parameter of both distributions are equal as
    %     shown in example 2.
    %
    %     In the cdf formula given by [1], there is a missing factor of
    %     1/pi for the case where alpha=1 and beta>0.
    %
    %INPUTS: X A vector of points at which the cdf should be evaluated.
    %    alpha The scalar characteristic parameter of the cdf.
    %      gam The scalar dispersion parameter of the cdf.
    %    delta The scalar location parameter of the cdf (also the
    %          mode of the distribution).
    %
    %OUTPUTS: cdf A vector of the same size as X which contains the
    %             corresponding cdf values for a S-alpha-S distribution
    %             with parameters alpha, beta, gam, and delta.
    %
    %EXAMPLE 1: Comparing Gaussian pdf to alpha=2 S-alpha-S pdf
    % x = linspace(-10,10);
    % gcdf = GaussianD.CDF(x,0,3);
    % cdf = AlphaStableD.CDF(x,2,0,sqrt(3/2),1);
    % plot(x,gcdf,x,cdf)
    %
    %EXAMPLE 2: Comparing Cauchy pdf to alpha=1 S-alpha-S pdf
    % x = linspace(-10,10);
    % ccdf = CauchyD.CDF(x,0,2);
    % cdf = AlphaStableD.CDF(x,1,0,2,1);
    % plot(x,ccdf,x,cdf)
    %
    %EXAMPLE 3: Comparing Lévy pdf to alpha=0.5, beta=1 stable pdf
    % x = linspace(0,10,200);
    % lcdf = LevyD.CDF(x,0,2);
    % cdf = AlphaStableD.CDF(x,0.5,1,2,1);
    % plot(x,lcdf,x,cdf)
    %
    %REFERENCES:
    %[1] Guillermo Julián-Moreno, Jorge E. López De Vergara, Iván
    %    González, Luis Pedro, Javier Royuela-Del-Val, and Federico 
    %    Simmross-Wattenberg, "Fast parallel alpha-stable distribution
    %    function evaluation and parameter estimation using OpenCL in 
    %    GPGPUs," Statistics and Computing 27 (2017), no. 5, 1365-1382.
    %
    %August 2019 Codie T. Lewis, Naval Research Laboratory, Washington D.C.
        
        %Default to standard distribution
        if(nargin<2)
            error('Must specify the characteristic parameter alpha.')
        end
        if(nargin<3)
            beta = 0;
            gam = 1;
            delta = 0;
        end
        if(nargin<4)
            gam = 1;
            delta = 0;
        end
        if(nargin<5)
            delta = 0;
        end
        if(alpha>2||alpha<=0)
            error('alpha must be in the interval (0,2]')
        end
        if(beta<-1||beta>1)
            error('beta must be in the interval [-1,1]')
        end
        if(gam<=0)
            error('gam must be greater than 0')
        end
        
        if(alpha==1)
            %Reparameterize the location parameter
            delta = delta+beta*log(gam)*2*gam/pi;
            %Apply scaling by gam and shifting location to delta
            X = (X-delta)/gam;
            if(beta==0)
                cdf = (1/2)+(1/pi)*atan(X);
            elseif(beta>0)
                V = @(theta) (2/pi)*((pi/2+beta*theta)./cos(theta)).*exp((1/beta)*(pi/2+beta*theta).*tan(theta));
                g = @(theta,x) exp(-pi*x/(2*beta)).*V(theta);
                cdf = (1/pi)*arrayfun(@(x)integral1DAdaptive(@(theta)exp(-g(theta,x)),[-pi/2,pi/2]),X);
            else
                beta = -beta;
                V = @(theta) (2/pi)*((pi/2+beta*theta)./cos(theta)).*exp((1/beta)*(pi/2+beta*theta).*tan(theta));
                g = @(theta,x) exp(-pi*x/(2*beta)).*V(theta);
                cdf = (1/pi)*arrayfun(@(x)integral1DAdaptive(@(theta)exp(-g(theta,x)),[-pi/2,pi/2]),X);
                cdf = 1-cdf;
            end
        else
            %Reparameterize the location parameter
            delta = delta+beta*tan(alpha*pi/2)*gam;
            %Apply scaling by gam and shifting location to delta
            X = (X-delta)/gam;
            cdf = zeros(size(X));
            zet = -beta*tan(pi*alpha/2);
            
            theta0above = (1/alpha)*atan(beta*tan(pi*alpha/2));
            theta0below = (1/alpha)*atan(-beta*tan(pi*alpha/2));
            if(alpha>1)
                c1above = 1;
                c1below = 1;
            else
                c1above = (1/pi)*((pi/2)-theta0above);
                c1below = (1/pi)*((pi/2)-theta0below);
            end
            Xpos = X;
            Xpos(X<zet) = -Xpos(X<zet);
            
            %Above
            V = @(theta) (cos(alpha*theta0above))^(1/(alpha-1))*(cos(theta)./sin(alpha*(theta0above+theta))).^(alpha/(alpha-1)).*(cos(alpha*theta0above+(alpha-1)*theta)./cos(theta));
            g = @(theta,x) V(theta).*(x-zet).^(alpha/(alpha-1));
            cdf(X>zet) = arrayfun(@(x)integral1DAdaptive(@(theta)exp(-g(theta,x)),[-theta0above,pi/2]),Xpos(X>zet));
            cdf(X>zet) = c1above+(sign(1-alpha)/pi).*cdf(X>zet);
            
            %Below
            V = @(theta) (cos(alpha*theta0below))^(1/(alpha-1))*(cos(theta)./sin(alpha*(theta0below+theta))).^(alpha/(alpha-1)).*(cos(alpha*theta0below+(alpha-1).*theta)./cos(theta));
            g = @(theta,x) V(theta).*(x+zet).^(alpha/(alpha-1));
            cdf(X<zet) = arrayfun(@(x)integral1DAdaptive(@(theta)exp(-g(theta,x)),[-theta0below,pi/2]),Xpos(X<zet));
            cdf(X<zet) = c1below+(sign(1-alpha)/pi).*cdf(X<zet);
            cdf(X<zet) = 1-cdf(X<zet);
            
            cdf(X==zet) = (1/pi)*((pi/2)-theta0above);
        end
    end
    
    function val = rand(N,alpha,beta,gam,delta)
    %%RAND Obtains an N-dimensional random sample of independent
    %      stable random variables using the formula provided in [1].
    %
    %INPUTS: N: A 1X2 vector or a scalar. If a vector [m,n], the
    %           output will be an mXn matrix with mn iid random
    %           samples. If a scalar, then the returned matrix will be
    %           NXN and contain N^2 iid random samples.
    %    alpha: The scalar shape parameter of the distribution to be
    %           sampled.
    %     beta: The scalar skewness parameter of the distribution to be
    %           sampled.
    %      gam: The scalar scale parameter of the distribution to be
    %           sampled.
    %    delta: The scalar location parameter of the distribution to be
    %           sampled.
    %
    %OUTPUTS: val: A matrix of size given by input N which contains iid
    %             random samples from the distribution.
    %
    %EXAMPLE 1: Comparing a histogram of the generated values to the
    %           exact pdf.
    % values = AlphaStableD.rand([1,1e6],0.5,1);
    %
    % % We want to see the central mass of the histogram.
    % % This will introduce some nonnegligible overestimation of the
    % % pdf near the mode of the distribution. For this reason, if you
    % % try to look at a fit to the cdf, the histogram bins to the
    % % right of the mode will be disproportionately large. More care
    % % is needed for such a plot to be accurate.
    % D = 100 % Trim to (0,100)
    % trimmed = values(values<=D);
    % 
    % x = linspace(0,D,1000);
    % pdf = AlphaStableD.PDF(x,0.5,1);
    % 
    % histogram(trimmed,'Normalization','pdf')
    % hold on
    % plot(x,pdf,'LineWidth',5)
    % hold off
    %
    %REFERENCES:
    %[1] Weron, R. (1996). "On the Chambers-Mallows-Stuck method for
    %    simulating skewed stable random variables," Statistics and
    %    Probability Letters, 28(2), 165-171.
    %
    %August 2019 Codie T. Lewis, Naval Research Laboratory, Washington D.C.
        
        %Default to standard distribution
        if(nargin<2)
            error('Must specify the characteristic parameter alpha.')
        end
        if(nargin<3)
            beta = 0;
            gam = 1;
            delta = 0;
        end
        if(nargin<4)
            gam = 1;
            delta = 0;
        end
        if(nargin<5)
            delta = 0;
        end
        if(alpha>2||alpha<=0)
            error('alpha must be in the interval (0,2]')
        end
        if(beta<-1||beta>1)
            error('beta must be in the interval [-1,1]')
        end
        if(gam<=0)
            error('gam must be greater than 0')
        end
        if(alpha>2||alpha<=0)
            error('alpha must be in the interval (0,2]')
        end
        if(beta<-1||beta>1)
            error('beta must be in the interval [-1,1]')
        end
        if(gam<0)
            error('gam must be greater than 0')
        end
        
        V = pi*(rand(N)-0.5);
        W = ExponentialD.rand(N,1);
        if(alpha~=1)
            Balbe = atan(beta*tan(pi*alpha/2))/alpha;
            Salbe = (1+beta^2*tan(pi*alpha/2))^(1/(2*alpha));
            val = Salbe*(sin(alpha*(V+Balbe))./cos(V).^(1/alpha))...
                    .*(cos(V-alpha*(V+Balbe))./W).^((1-alpha)/alpha);
            val = val*gam+delta;
        else
            val = (2/pi)*(((pi/2)+beta*V).*tan(V)-beta*log(W.*cos(V)/((pi/2)+beta*V)));
            val = val*gam+(2/pi)*beta*gam*log(gam)+delta;
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