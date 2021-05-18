function xyz=genEqAreaGridOnSphere(n,algorithm)
%%GENEQAREAGRIDONSPHERE Generate a grid of Cartesian points on the
%                        surface of a unit sphere based on equal-area
%                        projections.
%
%INPUTS: n This determines the number of points to generate. For
%          algorithm=0, exectly n*n points will be generate.d For n=1,
%          approximately n*n points will be generated,
% algorithm This selects the algorithm to use. Possible values are:
%          0 Use Lambert's azimuthal projection for the entire sphere. This
%            leads to points being increasingly distorted near the North
%            pole.
%          1 (The default if omitted or an empty matrix is passed) Use
%            Lambert's azimuthal projection for each hemisphere.
%
%OUTPUTS: xyz A 3XnumPts matrix of the generated points.
%
%The algorithms of Sections 2.2 and 2.3 of [1] are used.
%
%EXAMPLE:
% n=20;
% pts=genEqAreaGridOnSphere(n);
% size(pts)
% figure(1)
% clf
% scatter3(pts(1,:),pts(2,:),pts(3,:),'.')
%
%REFERENCES:
%[1] D. Rosca "New uniform grids on the sphere," Astronomy & Astrophysics,
%    vol. 520, no. A63, Sep.-Oct. 2010.
%
%November 2020 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(algorithm))
    algorithm=1;
end

switch(algorithm)
    case 0
        %Make sure that n is even.
        n=n+mod(n,2);

        L=sqrt(pi);
        pts=linspace(-L,L,n);
        [aList,bList]=meshgrid(pts,pts);
        numPts=n*n;
        xyz=zeros(3,numPts);
        for curPt=1:numPts
            a=aList(curPt);
            b=bList(curPt);

            if(abs(b)<=abs(a))
               const=(2*a/pi)*sqrt(pi-a^2);
               arg=b*pi/(4*a);

               xyz(:,curPt)=[const*cos(arg);const*sin(arg);2*a^2/pi-1];
            else
               const=(2*b/pi)*sqrt(pi-b^2);
               arg=a*pi/(4*b);

               xyz(:,curPt)=[const*sin(arg);const*cos(arg);2*b^2/pi-1];
            end
        end
    case 1
        L=sqrt(pi/2);
     
        pts=linspace(-L,L,n);
        [xList,yList]=meshgrid(pts,pts);
        xyList=mapSquare2Disc([xList(:).';yList(:).']);

        numPts=n*n;
        xyz=zeros(3,2*numPts);
        curXYZ=1;
        for curPt=1:numPts
            x=xyList(1,curPt);
            y=xyList(2,curPt);
            
            const=sqrt(1-(x^2+y^2)/4);
            const1=(x^2+y^2)/2;
            xyz(:,curXYZ)=[const*x;const*y;1-const1];
            %Avoid there being a duplicate row along the equator.
            if(abs(const1-1)<4*eps())
                curXYZ=curXYZ+1;
            else
                xyz(:,curXYZ+1)=[xyz(1:2,curXYZ);-1+const1];
                curXYZ=curXYZ+2;
            end
        end
        %Shrink to fit. The size can vary based on the duplciate row.
        xyz=xyz(:,1:(curXYZ-1));
    otherwise
        error('Unknown algorithm selected.')
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
