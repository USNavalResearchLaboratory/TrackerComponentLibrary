function intVal=symIntSecondKind(x,y,z)
%%SYNINTSECONDKIND Compute the symmetric integral of the second kind. There
%        are two common forms of this integral. The first is (1/4) times
%        the integral from zero to infinity of
%        1/sqrt((t+x)*(t+y)*(t+z))*(x/(t+x)+y/(t+y)+z/(t+z))*t dt The
%        second form of the integral is (1/(4*pi)) times the double
%        integral from 0 to 2*pi and from 0 to pi of
%        sqrt(x*sin(theta)^2*cos(phi)^2+y*sin(theta)^2*sin(phi)^2+z*cos(theta)^2)*sin(theta) dTheta dPhi.
%        Here, at most one of x, y, and z can be zero. Sometimes, this
%        integral is referred to as RG.
%
%INPUTS: x,y,z The three parameters of the function. At most one of them
%              can be zero. When complex, the values are assumed to have
%              complex phase angle less than pi in magnitude.
%
%OUTPUTS: intVal The value of the symmetric integral of the second kind given x,
%         y, and z.
%   
%This function implements equation [7] in [1].
%   
%REFERENCES:
%[1] B. C. Carlson, "Numerical computation of real or complex elliptic
%    integrals," Numerical Algorithms, vol. 10, no. 1, pp. 13-26, 1995.
%
%July 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    %Permute so that z has the largest magnitude so as to lessen finite precision problems.
    absZ=abs(z);
    absX=abs(x);
    absY=abs(y);
    
    if(absX>absZ&&absX>absY)
        %If x is the biggest.
        temp=z;
        z=x;
        x=temp;
    elseif(absY>absZ&&absY>absX)
        %If y is the biggest.
        temp=z;
        z=y;
        y=temp;
    end
    
    intVal=(1/2)*(z*symIntFirstKind(x,y,z)-(1/3)*(x-z)*(y-z)*symIntThirdKindDegen(x,y,z)+sqrt(x*y/z));
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
