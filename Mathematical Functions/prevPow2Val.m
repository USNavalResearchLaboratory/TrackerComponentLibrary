function vals=prevPow2Val(n)
%%PREVPOW2VAL Given a value n, return the value that is the closest power
%             of 2 that is <=abs(n).
%
%INPUTS: n A matrix or scalar. This can be any numeric type.
%
%OUTPUTS: vals A real matrix or scalar of the same type as n where each
%              element is the closest power of 2 that is <=abs(n).
%
%EXAMPLE:
% n=0.5:1:16.5;
% vals=prevPow2Val(n)
%One will get vals=[0,1,2,2,4,4,4,4,8,8,8,8,8,8,8,8,16] .
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(~isinteger(n))
    absN=abs(n);
    [F,E]=log2(absN);

    vals=zeros(size(n),class(n));
    %The conditional statement deals with E=0, which means that the next
    %lowest power of two is zero.
    sel=absN>=1;    
    vals(sel)=pow2(E(sel)-1);
    
    %Deal with non-finite terms
    sel=~isfinite(F);
    vals(sel)=F(sel);
else
    vals=zeros(size(n),class(n));
    absN=abs(n);

    x=bitshift(absN,-1);
    while(any(x(:))~=0)
        vals=vals+sign(x);
        x=bitshift(x,-1);
    end

    vals=bitshift(ones(class(vals)),vals);
    vals(n==0)=0;
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
