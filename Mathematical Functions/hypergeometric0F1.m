function S=hypergeometric0F1(a,z)
%%HYPERGEOMETRIC0F1 Evaluate the confluent hypergeometric function 0F1 for
%                   real or complex parameters. This is the confluent
%                   hypergeometric limit function 0F1(;a;z).
%
%INPUTS: a A real or complex scalar value.
%        z A real or complex scalar value.
%
%OUTPUTS: S The value 0F1(;a;z).
%
%This implements the Taylor series expansion of APpendix I of [1]. It was
%noted that the series generated at least 12 digits of accuracy for
%abs(z)<=1000. A check was added for values of a that are real and negtive
%as the solution should be infinite in that case. However, neighboring
%solutions are finite.
%
%REFERENCES:
%[1] J. Pearson, "Computation of hypergeometric functions," Master's
%    thesis, University of Oxford, Worcester College, 4 Sep. 2009.
%
%July 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%Check for the special case.
if(isreal(a)&&a<0)
   S=Inf;
   return;
end

A=1;
S=A;

epsVal=Inf;
k=0;
while(epsVal>eps(S))
    A=A*(1/(a+k))*(z/(k+1));
    S=S+A;
    
    epsVal=abs(A);
    
    k=k+1;
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
