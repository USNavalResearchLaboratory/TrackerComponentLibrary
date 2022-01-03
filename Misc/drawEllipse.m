function h=drawEllipse(z,A,gammaVal,varargin)
%DRAWELLIPSE Draw 2D ellipses or 3D ellipsoids on the current plot, where a
%            point zp on the ellipse/ ellipsoid satisfies the equation
%            (zp-z)'*A*(zp-z)=gammaVal. If z is a 2D/3D location estimate 
%            and inv(A) a Gaussian covariance matrix associated with the
%            estimate, then the ellipse is a probability region, where
%            gammaVal determines what amount of probability is in the
%            region. If omitted, gammaVal=16.2235 for 2D or
%            gammaVal=18.8049 for 3D, which corresponds approximately to
%            the 99.97% confidence region, is used (The value comes from
%            respectively ChiSquareD.invCDF(0.9997,2) and
%            ChiSquareD.invCDF(0.9997,3)).
%
%INPUTS: z A 2XN or 3XN vector corresponding to the centers of the N
%          ellipses or ellipsoids that should be drawn.
%        A A 2X2XN or 3X3XN set of N positive definite matrices that
%          specify the size and shape of the ellipse or ellipsoids, where
%          a point zp is on the ith ellipse/ ellipsoid if
%          (zp-z(:,i))'*A(:,:,i)*(zp-z(:,i))=gammaVal.
%          In three dimensions, A can have one zero eigenvalue, in which
%          case the fill3 command is used to draw an ellipse in 3D.
% gammaVal An optional parameter specifying the size of the ellipse/
%          ellipsoid. If omitted or an empty matrix is passed, then
%          gammaVal=16.2235 in 2D or gammaVal=18.8049 in 3D is used. These
%          are approximately the values for a 99.97% confidence region if A
%          are inverse covariance matrices of a Gaussian distribution.
%          gammaVal must be positive.
% varargin Sets of values that should be passed to the plot function to
%          format the ellipses or that will be passed to the surf function
%          to format the ellipsoids. For example, one could call the
%          function as drawEllipse(z,A,gammaVal,'--r','linewidth',2) to
%          plotellipses as thick, red lines. Note that Matlab will not
%          always properly render dashed lines due to the number of points
%          used to plot the shape. Also, if the ellipsoid is in 3D, but is
%          actually just a 2D ellipse, then if this parameter is omitted,
%          the fill3 command will be passed the option to color the
%          ellipse black. Otherwise, if this parameter is given, the user
%          must specify the ellipse color as it is not done here.
%
%OUTPUTS: h A NX1 cell array containing the plot objects for each of the
%           ellipses. This can be useful if, for example, one wishes to
%           change the transparency of object i to 50%, one can use the
%           alpha(h{i},0.5) command.
%
%When drawing 2D ellipses, the function getEllipsePoints is used to get the
%points on the ellipse, which are subsequently plotted.
%
%To find the points on an axis-aligned ellipsoid, we will use the ellipsoid
%function that is built into Matlab. Suppose that l=[x;y;z] and
%D=diag(a,b,c). Then the ellipsoid equation is x^2*a+y^2*b+z^2c=gammaVal.
%This means that the semi-axis lengths of this centered, axis-aligned
%ellipsoid are sqrt(gammaVal/a), sqrt(gammaVal/b), and sqrt(gammaVal/c).
%2501 vertices are used in the ellipsoids.
%
%July 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numEllipse=size(z,2);

%Save the value of hold on the plot so that it can be reset to its previous
%state after plotting all of the ellipses.
holdVal=ishold();

numDim=size(z,1);

if(nargout>0)
    h=cell(numEllipse,1);
end

switch(numDim)
    case 2
        if(nargin<3||isempty(gammaVal))
            gammaVal=16.2235;
        end

        numPoints=1000;
        for curEllip=1:numEllipse            
            zp=getEllipsePoints(z(:,curEllip),A(:,:,curEllip),gammaVal,numPoints,false);
            
            %Make sure that all of the ellipses are printed.
            if(curEllip~=1)
                hold on
            end
            x=plot(zp(1,:),zp(2,:),varargin{:});
            
            if(nargout>0)
               h{curEllip}=x;
            end
        end
    case 3
        if(nargin<3||isempty(gammaVal))
            gammaVal=18.8049;
        end

        for curEllip=1:numEllipse
            if(rank(A(:,:, curEllip))<3)
                numPoints=500;
                %If the ellipsoid is actually just an ellipse, then we will
                %draw the 2D ellipse in 3D.
                
                %Perform an eigenvalue decomposition of A.
                [V,D]=eig(A(:,:,curEllip));
                %We will put the zero eigenvalue at the start, if it is
                %not already there.
                [~,idx]=sort(diag(D),'ascend');
                V=V(idx,idx);
                D=D(idx,idx);
                
                a=D(2,2);
                b=D(3,3);
                
                xBound=sqrt(gammaVal/a);
                x=linspace(-xBound,xBound,numPoints);
                %The real command deals with possible finite precision
                %issues.
                y=real(sqrt((gammaVal-x.^2*a)/b));
                %The centered, axis-aligned ellipse points.
                l=[zeros(1,2*numPoints);
                   x,fliplr(x);
                   y,-fliplr(y)];
                zp=bsxfun(@plus,V*l,z(idx,curEllip));
                idxInv=inversePermutation(idx);
                %Undo the reordering.
                zp=zp(idxInv);

                %Make sure that all of the ellipses are printed.
                if(curEllip~=1)
                    hold on
                end
                %Draw the ellipse as a filled shape.
                if(isempty(varargin))
                    x=fill3(zp(1,:),zp(2,:),zp(3,:),'k');
                else
                    x=fill3(zp(1,:),zp(2,:),zp(3,:),varargin{:});
                end
                if(nargout>0)
                    h{curEllip}=x;
                end
            else
                numPoints=50;
                %Perform an eigenvalue decomposition of A.
                [V,D]=eig(A(:,:,curEllip));
                zCur=z(:,curEllip);

                [xp,yp,zp]=ellipsoid(0,0,0,sqrt(gammaVal/D(1,1)),sqrt(gammaVal/D(2,2)),sqrt(gammaVal/D(3,3)),numPoints);

                %Rotate to the correct orientation
                xyzPoints=V*[xp(:)';yp(:)';zp(:)'];
                xp=reshape(zCur(1)+xyzPoints(1,:),numPoints+1,numPoints+1);
                yp=reshape(zCur(2)+xyzPoints(2,:),numPoints+1,numPoints+1);
                zp=reshape(zCur(3)+xyzPoints(3,:),numPoints+1,numPoints+1);

                %Make sure that all of the ellipses are printed.
                if(curEllip~=1)
                    hold on
                end

                x=surf(xp,yp,zp,varargin{:});
                if(nargout>0)
                    h{curEllip}=x;
                end
            end
        end
    otherwise
        error('Invalid Dimensionality')
end
        
%Restore the hold value to its original setting.
if(~holdVal)
    hold off
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
