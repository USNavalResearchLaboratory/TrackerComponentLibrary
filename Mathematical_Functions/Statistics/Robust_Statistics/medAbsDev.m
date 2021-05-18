function MAD=medAbsDev(x,dim,scaleType)
%%MEDABSDEV Compute the median absolute deviation about the median (MAD) of
%           a set of points. This is a measure of the variability of a set
%           of data that is robust to outliers. This is the median of the
%           absolute value of the difference between the points and their
%           median. For scalar values, that is median(abs(x-median(x))).
%           This is more robust to outliers than the standard deviation.
%           For example, in [1], the issue of outlier detection given a
%           sample is discussed.
%
%INPUT: x A real vector, matrix, or hypermatrix of points. If an empty
%         matrix is passed, the NaN is returned. The dimension over which
%         the operation is taken is specified by dim.
%     dim The optional parameter specifying the dimension over which the
%         median absolute deviation is evaluated. The defaults if omitted
%         are as follows: If a vector is passed, then it is taken over the
%         vector. If a matrix is passed, then it is taken over the columns;
%         if a hypermatrix is passed, then the first non-unitary dimension
%         is chosen.
% scaleType This parameter specifies how the MAD should be scaled. Possible
%         values are:
%         0 (The default if omitted or an empty matrix is passed) Do not
%           scale the estimate.
%         1 Scale the estimate to be asymptotically unbiased as an
%           estimator for the standard deviation of a Gaussian
%           distribution. The estimator is biased for small n.
%         2 Scale the estimate to be approximately unbiased for the given
%           n. This uses the low-precision empirically derived
%           multiplication factors given in [2].
%
%OUTPUTS: MAD The value(s) of the median absolute deviation.
%
%The MAD of a random variable z is such that
%Pr(abs(z-median(z))<MAD)=1/2
%Let x be the value z-median(z). If z is a normal distribution, then
%median(z) corresponds to the mean, so we are just centering it. Thus:
%(1/2)=Pr(abs(x)<MAD)
%Or, getting rid of the absolute value:
%(1/2)=Pr(x<=-MAD)+Pr(x>MAD)
%(1/2)=Pr(x<=-MAD)+(1-Pr(x<=MAD))
%If X is symmetric about 0 (as in the scalar normal distribution with zero
%mean) then we know that 
%Pr(x<=-MAD)=1-Pr(x<=MAD)
%Pr(x<=MAD)=3/4
%If x is Gaussian with standard deviation sigma, then for y=x/sigma, where
%y is a normal (0,1) random variable:
%Pr(y<=MAD/sigma)=3/4
%Thus, MAD/sigma=GaussianD.invCDF(3/4);
%or
%sigma=MAD/GaussianD.invCDF(3/4)
%Thus, the MAD can be used as a robust estimator of the standard deviation
%of a set of samples that are assumed drawn from a normal distribution.
%Note that 1//GaussianD.invCDF(3/4) is about 1.482602218505601860547076529.
%However, it was empirically demonstrated in [2] that the estimator is
%based when considering a small number of samples. Thus, it must be scaled.
%Low precision heuristic scale factors are given in [2].
%
%EXAMPLE:
%Here, we sample a normal distribution and then use the scaled median
%absolute derivative to estimate the standard deviation.
% numSamples=1e4;
% sigma=20;
% mu=-8;
% x=mu+sigma*randn(numSamples,1);
% sigmaEst=medAbsDev(x,[],1)
%
%REFERENCES:
%[1] C. Leys, C. Ley, O. Klein, P. Bernard, and L. Licata, "Detecting
%    outliers: Do not use standard deviation around the mean, use absolute
%    deviation around the median," Journal of Experimental Social
%    Psychology, vol. 49, no. 4, pp. 764-766, Jul. 2013.
%[2] C. Croux and P. J. Rousseeuw, "Time-efficient algorithms for two
%    highly robust estimators of scale," in Computational Statistics.
%    Berlin: Springer-Verlag, Aug. 1992, vol. 1: Proceedings of the 10th
%    Symposium on Computational Statistics, pp. 411-428, conference
%    location: Neuchâtel, Switzerland.
%
%August 2018 David Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(isempty(x))
    MAD=NaN;
    return;
end

if(nargin<2||isempty(dim))
    if(isvector(x))
        x=x(:);
        dim=1;
    elseif(ismatrix(x))
        dim=1;
    else
        dim=find(size(x)>1,1);
        if(isempty(dim))
            dim=1; 
        end
    end
end

if(nargin<3||isempty(scaleType))
   scaleType=0; 
end

MAD=median(abs(bsxfun(@minus,x,median(x,dim))),dim);

c=1.482602218505601860547076529;
n=size(x,dim);
switch(scaleType)
    case 0%Do not scale
    case 1%Use only the asymptotic scaling.
        MAD=c*MAD;
    case 2%Use the heuristic scaling given in [2].
        switch(n)
            case 2
                bn=1.196;
            case 3
                bn=1.495;
            case 4
                bn=1.363;
            case 5
                bn=1.206;
            case 6
                bn=1.200;
            case 7
                bn=1.140;
            case 8
                bn=1.129;
            case 9
                bn=1.107;
            otherwise
                bn=n/(n-0.8);
        end
        
        MAD=bn*c*MAD;
    otherwise
        error('Unknown scale type specified.')
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
