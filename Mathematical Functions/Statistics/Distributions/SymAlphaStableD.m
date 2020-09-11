classdef SymAlphaStableD
%%SYMALPHASTABLED Functions to handle univariate symmetric alpha-stable
%                 (S-alpha-S) distributions
%Implemented methods are: PDF,CDF,rand (all for univariate case)
%
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

methods(Static)
    function pdf = PDF(X,alpha,gam,delta)
    %%PDF Obtain values of the symmetric alpha-stable (S-alpha-S)
    %     probability density function with parameters alpha, gam, and
    %     delta evaluated at points in X. The values in X will be
    %     transformed to the corresponding range values for the standard
    %     S-alpha-S distribution using (X-delta)/gam and evaluated the pdf
    %     at those points. If only X, alpha, and some of the input 
    %     parameters are given, defaults of 1 and 0 are assigned to
    %     gam and delta respectively.
    %
    %     For alpha < 0.001, the gamma function call below, 
    %     gamma(1/alpha), will return inf and therefore pdf(X==0) will 
    %     equal inf.
    %
    %     Note that for the alpha = 2 case (a Gaussian distribution) 
    %     we can relate the variance (sigma) to the scaling
    %     parameter (gam) by gam = sqrt(sigma/2) as shown in example 1.
    %     In the case of alpha = 1 (Cauchy distribution), we have that
    %     the scaling parameter of both distributions are equal as
    %     shown in example 2.
    %
    %INPUTS: X A vector of points at which the pdf should be evaluated.
    %    alpha The scalar characteristic parameter of the pdf.
    %      gam The scalar dispersion parameter of the pdf.
    %    delta The scalar location parameter of the pdf (also the
    %           mode of the distribution).
    %
    %OUTPUTS: pdf A vector of the same size as X which contains the
    %             corresponding pdf values for a S-alpha-S distribution
    %             with parameters alpha, gam, and delta.
    %
    %EXAMPLE 1: Comparing Gaussian pdf to alpha=2 S-alpha-S pdf
    % x = linspace(-10,10);
    % gpdf = GaussianD.PDF(x,0,3);
    % pdf = SymAlphaStableD.PDF(x,2,sqrt(3/2),1);
    % plot(x,gpdf,x,pdf)
    %
    %EXAMPLE 2: Comparing Cauchy pdf to alpha=1 S-alpha-S pdf
    % x = linspace(-10,10);
    % cpdf = CauchyD.PDF(x,0,2);
    % pdf = SymAlphaStableD.PDF(x,1,2,1);
    % plot(x,cpdf,x,pdf)
    %
    %REFERENCES:
    %[1] Guillermo Julián-Moreno, Jorge E. López De Vergara, Iván
    %    González, Luis Pedro, Javier Royuela-Del-Val, and Federico 
    %    Simmross-Wattenberg, "Fast parallel alpha-stable distribution
    %    function evaluation and parameter estimation using OpenCL in 
    %    GPGPUs," Statistics and Computing 27 (2017), no. 5, 1365-1382.
    %
    %August 2019 Codie T. Lewis, Naval Research Laboratory, Washington D.C.
        
        %Default to standard S-alpha-S
        if(nargin<2)
            error('Must specify the characteristic parameter alpha.')
        end
        if(nargin<3)
            gam = 1;
            delta = 0;
        end
        if(nargin<4)
            delta = 0;
        end
        if(alpha>2||alpha<=0)
            error('alpha must be in the interval (0,2]')
        end
        if(gam<=0)
            error('gam must be greater than 0')
        end
        
        %Apply scaling by gam and shifting median to delta
        X = (X-delta)./gam;
        
        if(alpha==1)
            pdf = 1./(gam*pi*(1+X.^2));
        else
            pdf = zeros(size(X));
            Xpos = abs(X(X~=0));
            V = @(theta) (cos(theta)./sin(alpha*theta)).^(alpha/(alpha-1)).*(cos((alpha-1).*theta)./cos(theta));
            g = @(theta,x) V(theta).*x.^(alpha/(alpha-1));
            h = @(theta,x) g(theta,x).*exp(-g(theta,x));
            
            pdf(X~=0) = arrayfun(@(x)integral1DAdaptive(@(theta)h(theta,x),[0,pi/2]),Xpos);
            pdf(X~=0) = (alpha./(pi*Xpos*abs(alpha-1))).*pdf(X~=0);
            pdf(X==0) = (1/alpha)*gamma(1/alpha)/pi;
            pdf = pdf/gam;
        end
    end
    
    function cdf = CDF(X,alpha,gam,delta)
    %%CDF Obtain values of the symmetric alpha-stable (S-alpha-S) 
    %     cumulative density function with parameters alpha, gam, and delta 
    %     evaluated at points in X. The values in X will be transformed
    %     to the corresponding range values for the standard S-alpha-S 
    %     distribution using (X-delta)/gam and evaluated the pdf at 
    %     those points. If only X, alpha, and some of the input 
    %     parameters are given, defaults of 1 and 0 are assigned to
    %     gam and delta respectively.
    %
    %     Note that for the alpha = 2 case (a Gaussian distribution) 
    %     we can relate the variance (sigma) to the scaling
    %     parameter (gam) by gam = sqrt(sigma/2) as shown in example 1.
    %     In the case of alpha = 1 (Cauchy distribution), we have that
    %     the scaling parameter of both distributions are equal as
    %     shown in example 2.
    %
    %INPUTS: X A vector of points at which the cdf should be evaluated.
    %    alpha The scalar characteristic parameter of the cdf.
    %      gam The scalar dispersion parameter of the cdf.
    %    delta The scalar location parameter of the cdf (also the
    %          mode of the distribution).
    %
    %OUTPUTS: cdf A vector of the same size as X which contains the
    %             corresponding cdf values for a S-alpha-S distribution
    %             with parameters alpha, gam, and delta.
    %
    %EXAMPLE 1: Comparing Gaussian cdf to alpha=2 S-alpha-S pdf
    % x = linspace(-10,10);
    % gcdf = GaussianD.CDF(x,0,3);
    % cdf = SymAlphaStableD.CDF(x,2,sqrt(3/2),1);
    % plot(x,gcdf,x,cdf)
    %
    %EXAMPLE 2: Comparing Cauchy cdf to alpha=1 S-alpha-S pdf
    % x = linspace(-10,10);
    % ccdf = CauchyD.CDF(x,0,2);
    % cdf = SymAlphaStableD.CDF(x,1,2,1);
    % plot(x,ccdf,x,cdf)
    %
    %REFERENCES:
    %[1] Guillermo Julián-Moreno, Jorge E. López De Vergara, Iván
    %    González, Luis Pedro, Javier Royuela-Del-Val, and Federico 
    %    Simmross-Wattenberg, "Fast parallel alpha-stable distribution
    %    function evaluation and parameter estimation using OpenCL in 
    %    GPGPUs," Statistics and Computing 27 (2017), no. 5, 1365-1382.
    %
    %August 2019 Codie T. Lewis, Naval Research Laboratory, Washington D.C.
        
        %Default to standard SÎ±S
        if(nargin<2)
            error('Must specify the characteristic parameter alpha.')
        end
        if(nargin<3)
            gam = 1;
            delta = 0;
        end
        if(nargin<4)
            delta = 0;
        end
        if(alpha>2||alpha<=0)
            error('alpha must be in the interval (0,2]')
        end
        if(gam<=0)
            error('gam must be greater than 0')
        end
        
        %Apply scaling by gam and shifting median to delta
        X = (X-delta)./gam;
        
        if(alpha==1)
            cdf = 0.5+(1/pi)*atan(X);
        else
            cdf = zeros(size(X));
            if(alpha<1)
                c1 = 1/2;
            else
                c1 = 1;
            end
            Xpos = abs(X(X~=0));
            V = @(theta) (cos(theta)./sin(alpha*theta)).^(alpha/(alpha-1)).*(cos((alpha-1).*theta)./cos(theta));
            g = @(theta,x) V(theta).*x.^(alpha/(alpha-1));
            cdf(X~=0) = (c1+(sign(1-alpha)/pi)*...
                arrayfun(@(x)integral1DAdaptive(@(theta)exp(-g(theta,x)),[0,pi/2]),Xpos));
            cdf(X==0) = 1/2;
            cdf(X<0) = 1-cdf(X<0);
        end
    end
    
    function val = rand(N,alpha)
    %%RAND Obtains an N-dimensional random sample of independent
    %      symmetric alpha-stable random variables using the formula
    %      provided in Chapter 2.8, page 28 of [1].
    %
    %INPUTS: N A 1X2 vector or a scalar. If a vector [m,n], the
    %          output will be an mXn matrix with mn iid random samples.
    %          If a scalar, then the returned matrix will be NXN and
    %          contain N^2 iid random samples.
    %    alpha The scalar shape parameter of the distribution to be
    %          sampled.
    %
    %OUTPUTS: val A matrix of size given by input N which contains iid
    %             random samples from the distribution.
    %
    %EXAMPLE 1: Comparing a histogram of the generated values to the
    %           exact pdf.
    % values = SymAlphaStableD.rand([1,1e6],1.5);
    %
    % % We want to see the central mass of the histogram.
    % % This will introduce some nonnegligible overestimation of the
    % % pdf near the mode of the distribution.
    % D = 10 % Trim to (-10,10)
    % trimmed = values(values>=-D);
    % trimmed = trimmed(trimmed<=D);
    % 
    % x = linspace(-D,D,1000);
    % pdf = SymAlphaStableD.PDF(x,1.5,1,0);
    % 
    % histogram(trimmed,'Normalization','pdf')
    % hold on
    % plot(x,pdf,'LineWidth',5)
    % hold off
    %
    %REFERENCES:
    % [1] C. L. Nikias and M. Shai, Signal Processing with Alpha-Stable
    %     Distributions and Applications. New York: John Wiley and Sons,
    %     1995.
    %
    %August 2019 Codie T. Lewis, Naval Research Laboratory, Washington D.C.
        
        if(alpha>2||alpha<=0)
            error('alpha must be in the interval (0,2]')
        end
        if(alpha~=1)
            Phi = pi*(rand(N)-0.5);
            W = ExponentialD.rand(N,1);
            epsilon = 1-alpha;
            a = tan(0.5*Phi);
            b = tan(0.5*epsilon*Phi);
            z = cos(epsilon*Phi)./(W.*cos(Phi));
            d = (z.^(epsilon/alpha)-1)./epsilon;
            val = (2*(a-b).*(1+a.*b))./((1-a.^2).*(1+b.^2)).*(1+epsilon.*d);
        else
            val = tan(pi*(rand(N)-0.5));
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
