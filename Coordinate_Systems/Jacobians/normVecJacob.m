function J=normVecJacob(x)
%%NORVECJACOB Consider a unit vector created as u=x/norm(x). This function
%         returns the partial derivatives of u with respect to the elements
%         of x. This partial derivative 
%
%INPUTS: x A real NX1 vector.
%
%OUTPUTS: J The NXN partial derivative matrix The columns are the element
%           of x with respect to which the partial derivative is taken and
%           the rows are the elements of u.
%
%EXAMPLE:
%Here, we verify that the output of this function agrees with the numerical
%differentiation result.
% x=[4;-18;12];
% u=@(x)x/norm(x);
% JNumDiff=numDiff(x,u,3);
% J=normVecJacob(x);
% RelErr=max(abs((J(:)-JNumDiff(:))./J(:)))
%The relative error will be about 3.8488e-10, indicating a good level of
%agreement between the numerical differentiation result and the analytic
%result.
%
%November 2020 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

xDim=length(x);
normX=norm(x);

J=-x*x'/normX^3+eye(xDim,xDim)/normX;

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
