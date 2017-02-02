function x=acceptRejectSample(N,f,h,hSamp,c)
%%ACCEPTREJECTSAMP Sample the probability distribution function f (which
%               need not be normalized) using acceptance-rejection
%               sampling. This technique lets one sample f based on another
%               PDF h having the same support that is easy to sample.
%
%INPUTS: N  A scalar, indicating the number of samples to generate.
%        f  A function handle for the (not necessarily normalized, possibly 
%           multivariate) PDF that is to be sampled. f can be discrete.
%        h  A function handle for the PDF that is easy to sample. This must
%           have the same support as f (be nonzero/zero in the same
%           places). If f is discrete, then h must be discrete.
%     hSamp A function handle to draw samples of the distribution h.
%           hSamp() should return a single random sample of h as a
%           numDimX1 vector.
%         c A positive scalar value such that f(x)/h(x)<=c for all x in the
%           support of f,h. Ideally, for the algorithm to be fast, c
%           should be as small as possible, meaning that c should be the
%           supremum of f(x)/h(x).
%
%OUTPUTS: x A numDimXN matrix of N samples of the distribution f that takes
%           a numDimX1 input vector.
%
%The algorithm is described in Chapters 4.4 and 5.2 of [1]. The constant c
%equals the expected number of iterations of the algorithm required until a
%random sample is generated.
%
%EXAMPLE: This is based on example 5e in 1.
% %Consider the non-normalized gamma(3/2,1) distribution.
% f=@(x)sqrt(x)*exp(-x);
% %Say that the sampling distribution is exponential with mean 3/2
% h=@(x)ExponentialD.PDF(x,2/3);
% hSamp=@()ExponentialD.rand(1,2/3);
% %A value of c such that f(x)/h(x)<=c for all x is
% c=3^(3/2)/sqrt(2*pi*exp(1));
% %Two random samples are
% acceptRejectSample(2,f,h,hSamp,c)
%
%REFERENCES:
%[1] S. M. Ross, Simulation, 4th ed. Amsterdam: Academic Press, 2006.
%
%December 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%Draw the first random sample.
xSamp=hSamp();

numDim=length(xSamp);

x=zeros(numDim,N);
curSamp=1;
while(1)
    while(1)
        u=rand(1);
    
        if(u<=f(xSamp)/(c*h(xSamp)))
            x(:,curSamp)=xSamp;
            break;
        else
            xSamp=hSamp();
        end
    end
    
    curSamp=curSamp+1;
    if(curSamp>N)
        break;
    end
    xSamp=hSamp();
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
