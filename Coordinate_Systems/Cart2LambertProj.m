function xy=Cart2LambertProj(xyz)
%%CART2LAMBERT Map a 3D point on a sphere to a 2D Lambert projection, as
%           described in [1]. This can map all points except the antipodes. 
%
%INPUTS: xyz A 3XnumPts set of 3D points.
%
%OUTPUTS: xy The 2XnumPts set of points converted from xyz.
%
%This implements equation 6 in [1].
%
%EXAMPLE:
%Here, we show that Cart2LambertProj and LambertProj2Cart are inverses of
%each other (not counting the z=r polar points).
% xyzPts=genUniformGridOnSphere(10);
% xyPts=Cart2LambertProj(xyzPts);
% xyzPtsNew=LambertProj2Cart(xyPts,1);
% absErr=max(abs(xyzPtsNew(:)-xyzPts(:)))
%The absolute error indicates more than 15 digits of accuracy
%
%REFERENCES:
%[1] D. Rosca "New uniform grids on the sphere," Astronomy & Astrophysics,
%    vol. 520, no. A63, Sep.-Oct. 2010.
%
%November 2020 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numPts=size(xyz,2);

xy=zeros(2,numPts);
for curPt=1:numPts
    x=xyz(1,curPt);
    y=xyz(2,curPt);
    z=xyz(3,curPt);
    
    if(x==0&&y==0)
        %The z=r case.
        xy(:,curPt)=[sign(x)*Inf;sign(y)*Inf];
    else
        r=sqrt(x^2+y^2+z^2);
        const=sqrt(2*r/(r-z));
        xy(:,curPt)=[x*const;y*const];
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
