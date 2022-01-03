function H=angleBetweenVecsHessian(a,s,sameDeriv)
%%ANGLEBETWEENVECSHESSIAN Given 3X1 vectors a and s, this finds the
%           matrix of second derivatives of the angle between the vectors
%           taken either with respect to s twice or with respect to both a
%           and s.
%
%INPUTS: a A real 3X1 vector.
%        s A real 3X1 vector.
% sameDeriv If this is true, then the derivatives are taken with respect to
%           the vector s twice. If false, they are taken with respect to a
%           and s. The default if omitted or an empty matrix is passed is
%           true.
%
%OUTPUTS: H A 3X3 matrix of second derivatives. If sameDeriv=true, then
%           H(i,j) has the second derivative with respect to s(i) and s(j).
%           If sameDerive=false, then H(i,j) holds the second deriavtive
%           with respect to a(i) and s(j)
%
%This implements the derivatives of the gradient that is given in the 
%angleBetweenVecsGradient function.
%
%EXAMPLE 1:
%This is an example of differentiating with respect to s twice. The result
%from this function is compared to finite differencing. The agreement is
%within what one would expect due to finite precision limitations. 
% a=[6;-4;3];
% s=[-1;-2;12];
% f=@(s)vec(angleBetweenVecsGradient(a,s,1));
% HNumDiff=numDiff(s,f,3,2,1e-4*abs(s));
% H=angleBetweenVecsHessian(a,s,true);
% RelErr=max(abs((H(:)-HNumDiff(:))./H(:)))
%
%EXAMPLE 2:
%This is an example of differentiating with respect to a and s. The result
%from this function is compared to finite differencing. The agreement is
%within what one would expect due to finite precision limitations. 
% a=[6;-4;3];
% s=[-1;-2;12];
% f=@(s)vec(angleBetweenVecsGradient(a,s,0));
% HNumDiff=numDiff(s,f,3,2,1e-4*abs(s));
% sameDeriv=false;
% H=angleBetweenVecsHessian(a,s,sameDeriv);
% RelErr=max(abs((H(:)-HNumDiff(:))./H(:)))
%
%May 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(sameDeriv))
    sameDeriv=true; 
end

sx=s(1);
sy=s(2);
sz=s(3);
ax=a(1);
ay=a(2);
az=a(3);

Ma=[ax*(ay*sy+az*sz)-sx*(ay^2+az^2);
    ay*(ax*sx+az*sz)-sy*(ax^2+az^2);
    az*(ax*sx+ay*sy)-sz*(ax^2+ay^2)];
    
if(sameDeriv)
    denom1=norm(s)^2*norm(cross(a,s))^3;
    denom2=norm(s)^4*norm(cross(a,s));
    denom3=norm(s)^2*norm(cross(a,s));
    
    Ms=[sx*(ay*sy+az*sz)-ax*(sy^2+sz^2);
        sy*(ax*sx+az*sz)-ay*(sx^2+sz^2);
        sz*(ax*sx+ay*sy)-az*(sx^2+sy^2)];
    dMs=[ay*sy+az*sz, ay*sx-2*ax*sy,az*sx-2*ax*sz;
         ax*sy-2*ay*sx,ax*sx+az*sz,  az*sy-2*ay*sz;
         ax*sz-2*az*sx,ay*sz-2*az*sy,ax*sx+ay*sy];
     
    H=Ms*(Ma/denom1-2*s/denom2)'+dMs/denom3;
else
    %Take the derivative with respect to a and then the derivatives with
    %respect to s.
    dMa=[-ay^2-az^2, ax*ay,    ax*az;
          ax*ay,-ax^2-az^2,ay*az;
          ax*az,ay*az,     -ax^2-ay^2];
          
    MAlt=[ay*(ax*sy-ay*sx)+az*(ax*sz-az*sx);
          ax*(ay*sx-ax*sy)+az*(ay*sz-az*sy);
          ax*(az*sx-ax*sz)+ay*(az*sy-ay*sz)];
    denom1=norm(a)^2*norm(cross(a,s))^3;
    denom2=norm(a)^2*norm(cross(a,s));
    
    H=Ma*MAlt'/denom1+dMa/denom2;
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
