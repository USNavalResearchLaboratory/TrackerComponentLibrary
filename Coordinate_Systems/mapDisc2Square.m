function ab=mapDisc2Square(ab)
%%MAPDISC2SQUARE Given Cartesian points, map them from a space on a 2D
%           disc onto a square space using the transformation of [1]. This
%           is a bijective area-preserving transformation. This function
%           is the inverse of mapSquare2Disc.
%
%INPUTS: ab A 2XnumPts list of real 2D points to transform.
%
%OUTPUTS: ab The 2XnumPts set of points transformed to the square.
%
%This function implements the area-preserving transformation of [1].
%
%EXAMPLE:
%This shows how a grid of points mapped to a disc gets mapped right back by
%this function --the grid is plotted as blue points and the
%transformed-inverse transformed points in d. The red and blue overlap. In
%second plot, just the transformation of this function is plotted, so that
%one can see that it really does change things.
% numPts=50;
% halfLine=10;
% x=linspace(-halfLine,halfLine,numPts);
% 
% pts=zeros(2,numPts,(2*halfLine+1)*2);
% curPt=1;
% for curX=-halfLine:halfLine
%     pts(:,:,curPt)=[curX*ones(1,numPts);x];
%     curPt=curPt+1;
% end
% for curY=-halfLine:halfLine
%     pts(:,:,curPt)=[x;curY*ones(1,numPts)];
%     curPt=curPt+1;
% end
% 
% figure(1)
% clf
% hold on
% scatter(pts(1,:),pts(2,:),'.b')
% pts=mapDisc2Square(mapSquare2Disc(pts(:,:)));
% scatter(pts(1,:),pts(2,:),'or')
% 
% figure(2)
% clf
% pts=mapDisc2Square(pts(:,:));
% scatter(pts(1,:),pts(2,:),'.c')
%
%REFERENCES:
%[1] D. Rosca "New uniform grids on the sphere," Astronomy & Astrophysics,
%    vol. 520, no. A63, Sep.-Oct. 2010.
%
%November 2020 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

N=size(ab,2);

for curPt=1:N
    a=ab(1,curPt);
    b=ab(2,curPt);

    if(a==0&&b==0)
        ab(:,curPt)=[0;0];
    elseif(abs(b)<=abs(a))
        %One cannot use a four-quadrant inverse tangent.
        ab(:,curPt)=sign(a)*sqrt(a^2+b^2)*[sqrt(pi)/2;(2/sqrt(pi))*atan(b/a)];
    else
        %One cannot use a four-quadrant inverse tangent.
        ab(:,curPt)=sign(b)*sqrt(a^2+b^2)*[(2/sqrt(pi))*atan(a/b);sqrt(pi)/2];
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
