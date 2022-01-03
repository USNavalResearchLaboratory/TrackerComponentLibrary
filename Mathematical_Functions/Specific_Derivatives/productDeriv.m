function derivVal=productDeriv(f,fDeriv)
%%PRODUCTDERIV Evaluate the derivative of a product of function using the
%           chain rule. The product is prod_{i=1}^Nf(i). The derivatives
%           fDeriv of the functions and the functions themselves, f, are
%           provided.
%
%INPUTS: f An NX1 or 1XN vector such that f(i) is the function value of the
%          ith term in the product.
%   fDeriv An NX1 or 1XN vector such that fDeriv(i) is the derivative of
%          the ith term in the product.
%
%OUTPUTS: deirvVal The scalar derivative value of the function product.
%
%This computes the value of the derivative using the chain rules. That is
%that d/dx(f_1(x)*f_2(x))=d/dx(f_1(x))*f_2(x)+f_1(x)*d/dx(f_2(x)). This
%rule is applied repeatedly.
%
%EXAMPLE:
%The function under consideration is product_{i=0}^{N-1}(k-i)/(k-i+T) . We
%want to take the derivative with respect to T.
% N=10;
% i=0:(N-1);
% k=15;
% T=3;
% f=(k-i)./(k-i+T);
% %The derivative is with respect to T
% fDeriv=-((-i+k)./(-i+k+T).^2);
% derivVal=productDeriv(f,fDeriv)
% %One will get a derivative value of -0.053340750464453.
%
%February 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

N=length(f);

fProd=[cumprod(f(:),'reverse');1];

%This utilizes the product rule recursively.
derivVal=0;
for curLevel=N:-1:1
    derivVal=derivVal*f(curLevel)+fDeriv(curLevel)*fProd(curLevel+1);
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
