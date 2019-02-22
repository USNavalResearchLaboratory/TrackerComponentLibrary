function boolVal=isImag(x)
%%ISIMAG Returns true if x is a purely imaginary vector or matrix. That
%        is, it is complex and all real elements are exactly zero.
%
%INPUT: x A scalar or matrix. This must be of a numeric type; that is, not
%         a cell array, a string, etc.
%
%OUTPUTS: boolVal This is true is x is a complex matrix (isreal is false)
%                 and all real components are zero.
%
%December 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(isreal(x))
   boolVal=false;
   return;
end

if(all(real(x(:))==0))
    boolVal=true;
else
    boolVal=false;
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
