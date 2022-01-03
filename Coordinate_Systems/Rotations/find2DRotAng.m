function thetaRot=find2DRotAng(v1,v2,areNormalized)
%%FIND2DROTANG Find the angle of rotation to rotate the 2D vector v1
%              counterclockwise until it aligns with the vector v2. Such a
%              rotation could be performed with the matrix from the
%              rotMat2D function.
%
%INPUTS: v1 A 2XN set of starting vectors.
%        v2 A 2XN set of ending vectors.
% areNormalized If all of the vectors in v1 and v2 are normalized, then
%           this is true. The default if omitted or an empty matrix is
%           passed is false.
%
%OUTPUTS: thetaRot A 1XN set of counterclockwise rotation angles to rotate
%                  the vectors in v1 to those in v2.
%
%A formula using atan2 with normalized vectors is used.
%
%EXAMPLE:
%Here, a random initial vector is chosen. It is then rotated all around the
%circle. The difference between the rotation angle returned here and the
%original rotation angle is plotted. The difference is on th eorder of
%finite precision errors.
% v1=randn(2,1);
% v1=v1/norm(v1);
% 
% numPts=10;
% theta=linspace(0,2*pi,numPts);
% v2=zeros(2,numPts);
% for k=1:numPts
%     v2(:,k)=rotMat2D(theta(k))*v1;
% end
% 
% thetaBack=find2DRotAng(v1,v2);
% figure(1)
% clf
% hold on
% %Without wrapping, the difference would be 0 and then switch to 2*pi
% %because thetaBack is between -pi and pi and theta is between 0 and 2*pi.
% plot(wrapRange(theta-thetaBack,-pi,pi),'-k','linewidth',2)
%
%August 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(areNormalized))
    areNormalized=false;
end

if(areNormalized==false)
    v1=bsxfun(@rdivide,v1,sqrt(sum(v1.*v1,1)));
    v2=bsxfun(@rdivide,v2,sqrt(sum(v2.*v2,1)));
end

thetaRot=atan2(v1(1,:).*v2(2,:)-v1(2,:).*v2(1,:),v1(1,:).*v2(1,:)+v1(2,:).*v2(2,:));

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
