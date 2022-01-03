function cc = estAlphaStableCovCoeff(x,y,algorithm,p,c1,c2,warn)
%%ESTALPHASTABLECOVCOEFF A function to numerically compute the covariation
%                        of two variables given samples from their
%                        respective symmetric alpha-stable distributions 
%                        (S-alpha-S). It is assumed that each pair
%                        (x(i),y(i)) is  sampled independently.
%
%INPUTS: x A vector of samples from the first S-alpha-S distribution.
%        y A vector of samples from the second S-alpha-S distribution.
% algorithm A scalar specifying which algorithm to use. Possible values
%          are:
%          0 Uses a fractional lower-order moment estimator with p=1 given
%            in Equation 6.2 of Chapter 6, page 73 of [1]. This method has
%            larger variance for alpha parameters near 1.
%          1 Uses a fractional lower-order moment estimator with chosen
%            moment p. This method has larger variance for alpha
%            parameters near 1. Increasing p within what is permitted by
%            the assumptions seems to reduce the sample standard deviation
%            for the estimator a little. However, if the drawn samples of
%            the S-alpha-S distribution have an estimated characteristic
%            below the value of p (which violates the assumptions placed on
%            p), then the standard deviation of the estimator blows up.
%            Therefore, it is advised to use algorithm 0 instead of this
%            one. See example 1 and change the p parameter passed to this
%            function to experiment with this.
%          2 Uses a screened ratio estimator given in Equation 6.4 of
%            Chapter 6, page 74 of [1]. This estimator is unbiased and
%            consistent.
%    c1,c2 Scalars used in the screened ratio estimator as a lower bound
%          and upper bound respectively on the set characteristic function
%          as defined in Chapter 6, on page 74 of [1].
%     warn A boolean value (or some value which is associated with a
%          boolean value). When true, warning messages will be printed to 
%          the screen whenever the estimated characteristic parameter of x
%          or y does not meet the assumptions for the estimators will be
%          printed to the screen. Note that these messages will occur more
%          frequently for samples from distributions with characteristic
%          parameters near 1 or 2. Defaults to false.
%
%OUTPUTS: cc A scalar representing the estimated covariaton coefficient.
%
%EXAMPLE 1: This is based on the example given in Chapter 6, pages 75-76 of
%           [1] with the purpose of highlighting the different performances
%           of the FLOM and SCR estimators. Histograms are included to 
%           better analyze the distribtion of the estimates.
% %Set model parameters
% alpha = 1.5;
% samples = 5000;
% a1 = -0.75;
% a2 = 0.25;
% b1 = 0.18;
% b2 = 0.78;
% 
% %Run trials
% trials = 50;
% FLOM = zeros([1,trials]);
% SCR = zeros([1,trials]);
% for i = 1:50
%     U1 = SymAlphaStableD.rand([1,samples],alpha);
%     U2 = SymAlphaStableD.rand([1,samples],alpha); 
% 
%     X = a1*U1+a2*U2;
%     Y = b1*U1+b2*U2;
%     FLOM(i) = estAlphaStableCovCoeff(X,Y,1,1);
%     SCR(i) = estAlphaStableCovCoeff(X,Y,2);
% end
% 
% %Compute statistics and exact value
% avgFLOM = mean(FLOM);
% stdFLOM = std(FLOM);
% avgSCR = mean(SCR);
% stdSCR = std(SCR);
% lambdaTrue = (a1*b1^(alpha-1)+a2*b2^(alpha-1))/...
%              (abs(b1)^alpha+abs(b2)^alpha);
% 
% %Print statistics
% fprintf("FLOM: avg = %0.5f, std = %0.5f\n",avgFLOM,stdFLOM)
% fprintf("SCR: avg = %0.5f, std = %0.5f\n",avgSCR,stdSCR)
% disp('------------------------------')
% fprintf('True \x03bb_{XY}: %0.5f\n',lambdaTrue)
% 
% %Display histograms of the data
% histogram(FLOM,'DisplayStyle','stairs')
% hold on
% histogram(SCR,'DisplayStyle','stairs')
% xline(lambdaTrue,'--k','exact','LabelVerticalAlignment','bottom',...
%                                'LabelHorizontalAlignment','center');
% hold off
%
%EXAMPLE 2: Estimating the covariation of X and Y using the covariation
%           coefficient and the fact that the covariation of a symmetric
%           alpha-stable random variable with itself is equal to the
%           dispersion parameter.
% %Set model parameters
% alpha = 1.5;
% samples = 5000;
% a1 = -0.75;
% a2 = 0.25;
% b1 = 0.18;
% b2 = 0.78;
% 
% %Run trials
% trials = 50;
% SCR = zeros([1,trials]);
% gamEst = zeros([1,trials]);
% for i = 1:50
%     U1 = SymAlphaStableD.rand([1,samples],alpha);
%     U2 = SymAlphaStableD.rand([1,samples],alpha); 
% 
%     X = a1*U1+a2*U2;
%     Y = b1*U1+b2*U2;
% 
%     SCR(i) = estAlphaStableCovCoeff(X,Y,2);
%     params = estAlphaStableParams(Y);
%     gamEst(i) = params.disp;
% end
% 
% %Compute average over the trials and the exact values
% avgSCR = mean(SCR);
% avggam = mean(gamEst);
% lambdaTrue = (a1*b1^(alpha-1)+a2*b2^(alpha-1))/...
%              (abs(b1)^alpha+abs(b2)^alpha);
% gam = abs(b1)^alpha+abs(b2)^alpha;
% 
% %Computeapproximate and exact covariation [X,Y]_α and print
% exact = lambdaTrue*gam;
% approx = avgSCR*avggam;
% 
% fprintf("The approximated covariation [X,Y]_\x03bb is: %0.5f\n",approx)
% fprintf("The exact covariation [X,Y]_\x03bb is: %0.5f\n",exact)
%
%REFERENCES:
%[1] Chrysostomos L. Nikias and Min Shao, Signal processing with 
%    alpha-stable distributions and applications, Wiley-Interscience, New 
%    York, NY, USA, 1995.
%
%July 2019 Codie T. Lewis, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin < 3||isempty(algorithm))
    algorithm = 0;
end
if(nargin<7)
    warn = false;
end
if(~exist('p','var')||isempty(p))
    p = 1;
end
if warn
    %Check that the assumptions are met.
    xparams = estAlphaStableParams(x);
    yparams = estAlphaStableParams(y);
    if(xparams.char<=1||yparams.char<=1)||(xparams.char>2||yparams.char>2)
        warning(strcat('The estimated characteristic of the given',...
                 ' samples does not fall on the interval (1,2].'))
    elseif(p>=xparams.char||p>=yparams.char)
        warning('p must be strictly less than the characteristic parameter.')
    end
end

switch(algorithm)
    case 0 % FLOM with p=1
        cc = sum(x.*sign(y))/sum(abs(y));
        
    case 1 % FLOM with given p
        cc = sum(x.*abs(y).^(p-1).*sign(y))/sum(abs(y).^p);
        
    case 2 % Screened Ratio Estimator
        if(~exist('c1','var'))
            c1 = 0.1;
        elseif(c1<=0)
            error('c1 must be strictly greater than 0')
        end
        if(~exist('c2','var'))
            c2 = inf;
        elseif(c2<=c1)
            error('c2 must be strictly greater than c1')
        end
        if(sum((abs(y)>c1).*(abs(y)<c2))==0)
            error('There are no samples between c1 and c2.')
        end
        
        cc = sum(x.*y.^(-1).*((abs(y)>c1).*(abs(y)<c2)))/...
             sum(((abs(y)>c1).*(abs(y)<c2)));
    otherwise
        error('Unknown algorithm specified.')
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
