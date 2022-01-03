function params = estAlphaStableParams(x,bins,p)
%%ESTALPHASTABLEPARAMS A function to estimate parameters for a symmetric
%                      alpha-stable distribution (S-alpha-S) which would
%                      likely  produce the samples in x. The characteristic 
%                      parameter must be strictly greater than 0 and 
%                      strictly less than 2, as the method from [1] uses
%                      the aymptotically algebraic tails of non-Gaussian 
%                      S-alpha-S in its derivations of the estimators.                  
%
%INPUTS: x A 1XN vector of random samples. If empty, the params members are
%          all set to inf and the function terminates.
%     bins The number of partitions for the sample when computing the
%          characteristic parameter. Defaults to 10 percent of the number
%          of samples if not given or is empty.
%        p Fractional moment to be used for computation of the dispersion
%          parameter. This should be a scalar in the interval (0,1/2) and
%          will be multiplied by the characteristic parameter. If not 
%          given, a value of 1/3 will be used as suggested in [1]. An
%          error will result if the given value is not in the allowed
%          interval.
%
%The distributions of these estimators all converge to normal distributions
%with means equal to the true parameters and variances governed by the
%parameters and size of the sample. See [1] for more details.
%
%Be aware that the number of bins is important for both the characteristic
%estimator and the dispersion estimator. Choosing too few bins will result
%in large variances for the dispersion estimator. Choosing too many bins
%will introduce bias into the characteristic estimator as well as the
%dispersion estimator.
%
%OUTPUTS: params A structure with members corresponding to the estimated
%                S-alpha-S parameters. If passed an empty vector, all
%                members are set to inf. The members are:
%                 loc An estimate of the location parameter.
%                char An estimate of the characteristic parameter.
%                disp An estimate of the dispersion parameter.
%
%EXAMPLE 1: Estimating parameters from generated random variables
% %Sample a standard S-alpha-S distribution (alpha=1.5, gam=1, delta=0)
% X = SymAlphaStableD.rand([1,1e6],1.5);
% estAlphaStableParams(X)
%
%REFERENCES:
%[1] G. A. Tsihrintzis, C. L. Nikias, "Fast estimation of the parameters of
%    alpha-stable impulsive interference", IEEE Trans. Signal Processing, 
%    vol. 44, pp. 1492-1503, June 1996.
%
%July 2019 Codie T. Lewis, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

xsrt = sort(x);
lx = length(xsrt);

if(lx==0)
    params.loc = inf;
    params.char = inf;
    params.disp = inf;
    return
end

if(~exist('bins','var')||isempty(bins))
    bins = floor(0.1*lx);
end

%Median
if(lx/2==floor(lx/2))
    params.loc = (xsrt(lx/2)+xsrt(lx/2+1))/2;
else
    params.loc = xsrt(ceil(lx/2));
end

%Characteristic Exponent
xc = x-params.loc;
if(mod(lx,bins)==0)
    K = lx/bins;
else
    K = (lx-mod(lx,bins))/bins;
end

upx = zeros([1,bins]);
downx = zeros([1,bins]);
for i = 1:bins
    upx(i) = max(xc((K*(i-1)+1):K*i));
    downx(i) = min(xc((K*(i-1)+1):K*i));
end
upx = sign(upx).*log(abs(upx));
downx = sign(downx).*log(abs(downx));

avgupx = sum(upx)/bins;
ups = sqrt((1/(bins-1))*sum((upx-avgupx).^2));
avgdownx = sum(downx)/bins;
downs = sqrt((1/(bins-1))*sum((downx-avgdownx).^2));

params.char = (pi/(2*sqrt(6)))*(ups.^(-1)+downs.^(-1));

%Dispersion
if(exist('p','var'))
    if(p>0&&p<1/2)
        p = params.char*p;
    else
        error('p must be in the interval (0,1/2)')
    end
else
    p = params.char/3;
end

C = (1/cos((pi/2)*p))*(gamma(1-(p/params.char))/gamma(1-p));
params.disp = ((1/lx)*sum(abs(xc).^p)/C).^(params.char/p);
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
