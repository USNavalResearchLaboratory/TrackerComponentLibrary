function retVals=Besseli0Inc(theta,kappa)
%%BESSELIINC Evaluate the incomplete modified Bessel function of the first
%            kind and zeroth order I_0(theta,kappa). This is equivalent to
%            the left tail area of the von Mises distribution. 
%
%INPUTS: theta A scalar or vector of angular parameters from
%              -pi<=theta<=pi. Values of theta outside that region will be
%              wrapped into that region.
%        kappa A concentration parameter. Normally 0<kappa<Inf.
%
%OUTPUTS: retVals The values of the incomplete modified Bessel function of
%                 the first kind and zeroth order evaluated at the given
%                 parameters (one result for each theta value). 
%
%The incomplete modified Bessel function of the first kind is defined to be
%inv(2*pi*besseli(0,kappa))*(the integral from -pi to theta of
%exp(kappa*cos(x)) dx where theta ranges from -pi to pi. Thus,
%Besseli0Inc(pi,kappa)=1 for all kappa.
%
%The algorithm is taken from [1]. The constants for 12 digits of accuracy
%are used.
%
%REFERENCES:
%[1] G. W. Hill, "ALGORITHM 518 incomplete Bessel function i0: The
%    von Mises distribution [S14]," ACM Transactions on Mathematical
%    Software, vol. 3, no. 3, pp. 279-284, Sep. 1977.
%
%August 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%Constants for 12 digits of accuracy, as given in Table I of the paper.
A1=28;
A2=0.5;
A3=100;
A4=5;
CK=50;
C1=50.1;

numVals=numel(theta);
retVals=zeros(size(theta));

for curVal=1:numVals
    thetaCur=theta(curVal);
    
    Z=kappa;
    %This part differs from the paper, because implementing it as in the paper
    %would cause problems for values of theta that are just below pi.
    %This is because (pi-eps(pi))+pi is numerically equal to 2*pi even though
    %pi-eps(pi) is less than pi.
    if(thetaCur==pi||(thetaCur>-pi&&thetaCur+pi<=2*pi))
        U=thetaCur+pi;
    else
        U=mod(thetaCur+pi,2*pi);
    end
    if(U<0)
        U=U+2*pi;
    end

    Y=U-pi;
    if(Z>CK)
        %For large kappa, compute the normal approximation and left tail
        C=24*Z;
        V=C-C1;
        R=sqrt((54/(347/V+26-C)-6+C)/12);
        Z=sin(Y*0.5)*R;
        S=Z*Z*2;                                           
        V=V-S+3;
        Y=(C-S-S-16)/3;
        Y=((S+1.75)*S+83.5)/V-Y;
        retVal=erf(Z-S/(Y*Y)*Z)*0.5+0.5;
    elseif(Z<=0)
        retVal=(U*0.5)/pi;
    else
        %For small kappa, sum the IP terms by backwards recursion.
        IP=fix(Z*A2-A3/(Z+A4)+A1);
        P=IP;
        S=sin(Y);
        C=cos(Y);
        Y=P*Y;
        SN=sin(Y);
        CN=cos(Y);
        R=0;
        V=0;
        Z=2/Z;
        for N=2:IP
            P=P-1;
            Y=SN;
            SN=SN*C-CN*S;
            CN=CN*C+Y*S;
            R=1/(P*Z+R);
            V=(SN/P+V)*R;
        end
        retVal=(U*0.5+V)/pi;
    end

    retVal=max([retVal,0]);
    retVals(curVal)=min([retVal,1]);
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
