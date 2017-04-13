function val=erfComplex(z)
%%ERFCOMPLEX Evaluate the error function in a manner valid for complex
%            arguemnts. Matlab's built-in function erf cannot handle
%            complex arguments. The error function is defined as
%            (2/pi)*integral_0^z exp(-t^2) dt.
%
%INPUTS: z A scalar, vector, or matrix of values at which one wishes to
%          evaluate the error function.
%
%OUTPUTS: val A matrix having the same dimensions as z in which the values
%             of the scaled error function are computed for the values in
%             z.
%
%For values of abs(z)<0.08, six terms of the maclaurin series in [1] are
%used. These should be sufficient for convergence. For larger values, the
%Faddeeva function is used through the relation 
%erf(z)=1-exp(-z^2)*Faddeeva(1i*z);
%
%REFERENCES:
%[1] Weisstein, Eric W. "Erf." From MathWorld--A Wolfram Web Resource
%    http://mathworld.wolfram.com/Erf.html
%
%November 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numVals=numel(z);

zList=z;
val=zeros(size(z));

for curVal=1:numVals
    z=zList(curVal);
    y=imag(z);
    %If the argument is real, the built-in erf function is faster.
    if(y==0)
        x=real(z);
        val(curVal)=erf(x);
        continue;
    end

    if(abs(z)<0.08)
        %Use the Maclaurin series. For abs(z)<0.08, convergence occurs
        %within 6 terms.
        NMax=6;
        z2=z^2;

        %Initialize to the n=0 term
        sumVal=z;
        prodVal=z;
        for n=1:NMax
            prodVal=prodVal*z2*(1-2*n)/(n*(1+2*n));
            sumVal=sumVal+prodVal;
        end

        val(curVal)=2/sqrt(pi)*sumVal;
    else%Compute using the Faddeeva function.
        val(curVal)=1-exp(-z^2)*Faddeeva(1i*z);
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
