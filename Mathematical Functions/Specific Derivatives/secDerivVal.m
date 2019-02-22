function val=secDerivVal(z,n)
%%SECDERIVVAL Return the value of the nth derivative of sec(x) with respect
%             to x evaluated at x=z.
%
%INPUTS: z A matrix of real or complex values at which the nth derivative
%          of the secant function is desired.
%        n The number of derivatives to take. n>=0.
%
%OUTPUTS: val The value of the nth derivative of the secant function taken
%             at all of the points in z. val has the same dimensions as z.
%
%This function implements the algorithm of [1].
%
%REFERENCES:
%[1] M. E. Hoffman, "Derivative polynomials for tangent and secant," The
%    American Mathematical Monthly, vol. 102, no. 1, pp. 23-30, Jan. 1995.
%
%October 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numEls=numel(z);

val=zeros(size(z));
for curEl=1:numEls
    x=z(curEl);
    u=tan(x);

    P=zeros(n+1,1);
    Q=zeros(n+1,1);
    P(1)=u;
    Q(1)=1;
    
    for k=1:n
        %The loop implements Equation 5.
        sumValP=0;
        sumValQ=0;
        for i=0:(k-1)
            binomVal=binomial(k-1,i);
            sumValP=sumValP+binomVal*P(i+1)*P(k-1-i+1);
            sumValQ=sumValQ+binomVal*P(i+1)*Q(k-1-i+1);
        end
        if(k==1)
            sumValP=sumValP+1;
        end
        P(k+1)=sumValP;
        Q(k+1)=sumValQ;
    end
    val(curEl)=Q(end)*sec(x);
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
