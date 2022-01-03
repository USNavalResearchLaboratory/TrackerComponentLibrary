function J=angleBetweenVecsGradient(a,s,whichVector)
%%ANGLEBETWEENVECSGRADIENT Given 3X1 vectors a and s, this finds the
%           gradient of the angle between the vectors with respect to the
%           elements of either a or s. This is the graident of the function
%           angBetweenVecs.
%
%INPUTS: a A real 3X1 vector.
%        s A real 3X1 vector.
% whichVector If this is 0, then the gradient is taken with respect to the
%          elements of a. If this is 1, then the gradient is taken with
%          respect to the elements of s. The default if omitted or an empty
%          matrix is passed is 1.
%
%OUTPUTS: J A 1X3 vector where J(i) is the derivative of the angle with
%           respect to the its element of the selected vector (a or s).
%
%This implements the gradient of the cross-product relation used in
%the angBetweenVecs function.
%
%EXAMPLE:
%The maximum absolute relative error between a gradient obtained with
%numerical differentiation and one obtained from this function is computed.
%The error is on the order of what one might expect due to finite precision
%limitiation.
% a=[5;-6;12];
% b=[-1;-2;3];
% 
% f=@(x)angBetweenVecs(x,b);
% JNumDiff=numDiff(a,f,1,1,1e-5);
% %Gradient with respect to the first input.
% J=angleBetweenVecsGradient(a,b,0);
% RelErr=max(abs((JNumDiff-J)./JNumDiff))
% 
% f=@(x)angBetweenVecs(a,x);
% JNumDiff=numDiff(b,f,1,1,1e-5);
% %Gradient with respect to the second input.
% J=angleBetweenVecsGradient(a,b,1);
% RelErr=max(abs((JNumDiff-J)./JNumDiff))
%
%May 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(whichVector))
    whichVector=1;
end

switch(whichVector)
    case 0%Derivatives with respect to the elements of a.
        temp=s;
        s=a;
        a=temp;
    case 1%Derivatives with respect to the elements of s.
    otherwise
        error('Unknown vector selected.')
end

ax=a(1);
ay=a(2);
az=a(3);
sx=s(1);
sy=s(2);
sz=s(3);

denom=norm(s)^2*norm(cross(a,s));
J=[sx*(ay*sy+az*sz)-ax*(sy^2+sz^2),...
   sy*(ax*sx+az*sz)-ay*(sx^2+sz^2),...
   sz*(ax*sx+ay*sy)-az*(sx^2+sy^2)]/denom;

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
