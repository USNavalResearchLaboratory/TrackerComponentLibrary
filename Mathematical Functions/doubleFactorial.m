function val=doubleFactorial(x)
%%DOUBLEFACTORIAL Evaluate x!!, which is the double factorial (also known
%                 as semifactorial) of x. If x is odd, then the double
%                 factorial is defined as x*(x-2)*(x-4)*...*5*3*1. If x is
%                 even, then it is x(x-2)*(x-4)*...*6*4*2. The double
%                 factorials of 0 and -1 are defined to be 1.
%
%INPUTS: x A scalar or matrix, where the double factorial of each
%          element should be evaluated. The values must be positive
%          integers, or 0 or -1.   
%
%OUTPUTS: val The double factorial of each value of x.
%
%The gamma function/factorial representations of the double factorial for
%even and odd values are taken from [1].
%
%Values <=50 are tabulated. For other values, the functions are implemented
%using logarithmic functions so as to minimize the effects of overflow with
%intermediate results.
%
%REFERENCES:
%[1] Weisstein, Eric W. "Double Factorial." From MathWorld--A Wolfram Web
%    Resource. http://mathworld.wolfram.com/DoubleFactorial.html
%
%October 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numEl=length(x(:));

val=zeros(size(x));
for curEl=1:numEl
    xCur=x(curEl);
    
    if(xCur<=50)
        %If we can use a tabulated formula.
        if(xCur<-1)
            error('xCur must be >=-1.')
        end
        table=[1,1,1,2,3,8,15,48,105,384,945,3840,10395,46080,135135,645120,2027025,10321920,34459425,185794560,654729075,3715891200,13749310575,81749606400,316234143225,1961990553600,7905853580625,51011754393600,213458046676875,1428329123020800,6190283353629375,42849873690624000,191898783962510625,1371195958099968000,6332659870762850625,46620662575398912000,221643095476699771875,1678343852714360832000,8200794532637891559375,63777066403145711616000,319830986772877770815625,2551082656125828464640000,13113070457687988603440625,107145471557284795514880000,563862029680583509947946875,4714400748520531002654720000,25373791335626257947657609375,216862434431944426122117120000,1192568192774434123539907640625,10409396852733332453861621760000,58435841445947272053455474390625,520469842636666622693081088000000];
        
        val(curEl)=table(xCur+2);
    else
        if(xCur==0||xCur==-1)
            val(curEl)=1;
        elseif(mod(xCur,2)==0)%If it is divisible by 2.
            k=xCur/2;
            val(curEl)=exp(k*log(2)+gammaln(k+1));
        else%Assume that it is odd.
            n=(xCur+1)/2;
            val(curEl)=exp(gammaln(xCur+2)-(n*log(2)+gammaln(n+1)));
        end
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
