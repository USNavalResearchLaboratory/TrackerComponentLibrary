function f=evalContFrac(bFunc,aFunc,tinyVal,maxIter)
%%EVALCONTFRAC Evaluate a continued fraction of the form
%              f=b(1)+a(1)/(b(2)+a(2)/(b(3)+a(3)/(b(4)+a(4)/...
%              Given either matrices for a finite number of a and b terms
%              of function for the a and b values such that the required
%              number of terms can be requested.
%
%INPUTS: bFunc This can either be a function such that bFunc(i) provides
%              b(i) in the continued fraction, or this can be a maxNX1 or
%              1XmaxN vector of explicit b values.
%        aFunc This has to be of the same type as bFunc, either a function
%              handle or a vector to provide the a values. If a vector,
%              then it is 1X(maxN-1) or (maxN-1)X1 as it takes one fewer
%              call to the a terms than the b terms to evaluate the
%              continued fraction.
%      tinyVal An optional parameter that is less than typical values of
%              eps(b), but is not zero. This parameter is used to mitigate
%              divide-by zero errors. The default if this parameter is
%              omitted or an empty matrix is passed if 2e-100.
%      maxIter An optional parameter considering the maximum number of
%              terms to evaluate in the continued fraction. The default
%              value if omitted is 1000. This parameter is only used if
%              aFunc and bFunc are functions, not arrays. If they are
%              arrays, then the maximum number of terms used is determined
%              by the array length.
%
%OUTPUTS: f The value of the continued fraction. The smallest possible
%           value that this can take if b(0)=0 and all a and b are
%           non-negative is eps, due to how the function is initialized to
%           be numerically robust.
%
%The algorithm is the modified Lentz's method. A good description of the
%method, demonstrating its origin is given in Chapter II and Appendix C of
%[1]. The origin of the algorithm is [2]. However, the algorithm is harder
%to understand outside of the context of Bessel functions when consulting
%[2].
%
%An example of a continued fraction is an approximation to pi:
% bFunc=@(i)((i==1)*0+(i~=1)*(2*(i-1)-1));
% aFunc=@(i)(4*(i==1)+(i-1)^2*(i~=1));
% piApprox=evalContFrac(bFunc,aFunc)
%
%REFERENCES:
%[1] C.-Y. Shih, "An investigation of a stable algorithm for the evaluation
%    of fractional order Bessel functions," Master's thesis, Emporia
%    State University, Emporia, KS, Nov. 1993. [Online]. Available:
%    https://esirc.emporia.edu/bitstream/handle/123456789/1773/Shih%201993.pdf
%[2] W. J. Lentz, "Generating Bessel functions in Mie scattering
%    calculations using continued fractions," Applied optics, vol. 15, no.
%    3, pp. 668-671, Mar. 1976.
%
%October 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%The floating point precision.
epsVal=eps();

if(nargin<3)
    tinyVal=2^(-100);
end

if(isa(bFunc,'double'))
    maxN=length(bFunc);
    maxIter=maxN-1;
elseif(nargin<4)
    maxIter=1000;
end

b0=bFunc(1);
if(b0==0)
    f=tinyVal;
else
    f=b0;
end

CPrev=f;
DPrev=0;

j=1;
while(j<maxIter)
    bj=bFunc(j+1);
    aj=aFunc(j);
    
    DCur=bj+aj*DPrev;
    if(DCur==0)
        DCur=tinyVal;
    end
    DCur=1/DCur;
    DPrev=DCur;
    
    CCur=bj+aj/CPrev;
    if(CCur==0)
        CCur=tinyVal;
    end
    CPrev=CCur;
    
    Delta=CCur*DCur;
    f=f*Delta;
    
    if(all(abs(Delta-1)<epsVal))
        break;
    end
    j=j+1;
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
