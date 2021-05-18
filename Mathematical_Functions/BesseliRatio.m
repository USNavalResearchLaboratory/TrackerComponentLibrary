function x= BesseliRatio(nu,kappa,maxIter)
%%BESSELIRATIO Evaluate the ratio of modified Bessel functions of the first
%              kind of the form x=I_{nu}(kappa)/I_{nu-1}(kappa). 
%
%INPUTS: nu The positive integer (upper) order of the modified Bessel
%           function of the first kind in the ratio; nu>=1.
%     kappa The real argument of the modified Bessel function of the first
%           kind; kappa>=0.
%   maxIter An optional parameter specifying the maximum number of
%           iterations to use for computing the ratio using an iterative
%           method. If this parameter is omitted or an empty matrix is
%           passed, the default value of 2000 is used. Convergence usually
%           occurs long before the maximum number of iterations is reached.
%
%OUTPUTS: x The value of the ratio I_{nu}(kappa)/I_{nu-1}(kappa).
%
%Numerical precision limitations can make the evaluation of Bessel function
%radios difficult if one tries to explicitly evaluate the functions. Here,
%the algorithm of Perron described in [1] is implemented.
%
%REFERENCES:
%[1] W. Gautschi and J. Slavik, "On the computation of modified Bessel
%    function ratios," Mathematics of Computation, vol. 32, no. 143, pp.
%    865-875, Jul. 1978.
%
%September 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    if(nargin<3||isempty(maxIter))
        %It usually converges long before the maximum number of iterations.
        maxIterations=2000;
    end
    %tolerance for convergence. This should be suitable for double-
    %precision computations.
    tol=1e-15;
    
    k=1;
    cumProd=1;
    pCur=0.5*kappa*(nu+0.5)/((nu+kappa/2)*(nu+kappa+0.5)-0.5*kappa*(nu+0.5));
    cumProd=cumProd*pCur;
    cumSum=1+cumProd;
    while(k<maxIterations)
        pPrev=pCur;
        k=k+1;
        pCur=0.5*kappa*(nu+k-0.5)*(1+pPrev)/((nu+kappa+(k-1)/2)*(nu+kappa+k/2)-0.5*kappa*(nu+k-0.5)*(1+pPrev));
        cumProd=cumProd*pCur;
        
        cumSumNew=cumSum+cumProd;
        %If we have hit a precision limit; this usually happens before k
        %gets too big.
        if(abs(cumSumNew-cumSum)<tol)
            break;
        end
        cumSum=cumSumNew;
    end
    a0=(kappa+2*nu)/kappa;
    x=cumSum/a0;
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
