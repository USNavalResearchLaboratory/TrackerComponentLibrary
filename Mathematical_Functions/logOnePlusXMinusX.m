function val=logOnePlusXMinusX(x)
%%LOGONEPLUSXMINUSX Evaluate log(1+x)-x reducing the loss of precision that
%            occurs for values of x that are very close to 0 for real
%            values of x.
%
%INPUTS: x A matrix of real values>=-1.
%
%OUTPUTS: val The values log(1+x)-x.
%
%The function first computes z=log(1+x) using log1p(x) to reduce finite
%precision errors. This value is put into the expSlope2 function, which
%computes 2*(exp(x)-1-x)/x^2 without a low loss of finite precision.
%Multiplying the result by (-1/2)*z^2, one gets the expression log(1+x)-x.
%This approach ovoids a direct subtraction of x, which is where the finite
%precision problems occur.
%
%EXAMPLE:
%Here, we show the improvement that this algorithm offers to a direct
%evaluation.
% x=1e-150;
% valAccurate=logOnePlusXMinusX(x)
% valInaccurate=log1p(x)-x
%One will see that valAccurate is about -5e-301, whereas valInaccurate is
%just zero. The value -5e-301 can be verified using greatly extended
%precision arithmetic.
%
%October 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    if(any(x(:)<-1)||~isreal(x))
       error('This function only supports real values of x >=-1.');
    end
    
    val=zeros(size(x));
    
    selNeg=x<-0.99;
    xNeg=x(selNeg);
    xPos=x(~selNeg);
    if(~isempty(xNeg))
        valNeg=log1p(x)-x;
    else
        valNeg=[]; 
    end
    
    if(~isempty(xPos))
        z=log1p(x);
        e2Val=expSlope2(z);
        valPos=-(1/2)*e2Val.*z.^2;
    else
        valPos=[];
    end
    
    val(selNeg)=valNeg;
    val(~selNeg)=valPos;
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
