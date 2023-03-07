function ab=mapSquare2Disc(ab)
%%MAPSQUARE2DISC Given points in a Cartesian space, map the points onto a
%           disc using the transformation of [1]. This is a bijective area-
%           preserving transformation. Points in a square having side-
%           length 2*L is mapped to a circle of radius 4*L^2/pi (The areas
%           are equal).
%
%INPUTS: ab A 2XnumPts list of real 2D points to transform.
%
%OUTPUTS: ab The 2XnumPts set of points transformed to the disc.
%
%This function implements the area-preserving transformation of [1].
%
%EXMAPLE:
%This plots a bunch of points forming gridlines in blue. Then, it maps them
%to the disc and plots the resulting figure in red above the blue plot.
% numPts=500;
% halfLine=10;
% x=linspace(-halfLine,halfLine,numPts);
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
% pts=mapSquare2Disc(pts(:,:));
% scatter(pts(1,:),pts(2,:),'.r')
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
        %Origin special case
        ab(:,curPt)=[0;0];
    elseif(abs(b)<=abs(a))
        %Equation 2
        ab(:,curPt)=(2*a/sqrt(pi))*[cos(b*pi/(4*a));sin(b*pi/(4*a))];
    else
        %Equation 3
        ab(:,curPt)=(2*b/sqrt(pi))*[sin(a*pi/(4*b));cos(a*pi/(4*b))];
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
