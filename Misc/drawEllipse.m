function drawEllipse(z,A,gamma,varargin)
%DRAWELLIPSE Draw 2D ellipses or 3D ellipsoids on the curren plot, where a
%            point zp on the ellipse/ ellipsoid satisfies the equation
%            (zp-z)'*A*(zp-z)=gamma. If z is a 2D/3D location estimate and 
%            inv(A) a Gaussian covariance matrix associated with the
%            estimate, then the ellipse is a probability region, where
%            gamma determines what amount of probability is in the region.
%            If omitted, gamma=16.2235 for 2D or gamma=18.8049 for 3D,
%            which corresponds to approximately to the 99.97% confidence
%            region, is used (The value comes from respectively
%            ChiSquareD.invCDF(0.9997,2) and ChiSquareD.invCDF(0.9997,3)).
%
%INPUTS: z  A 2XN or 3XNvector corresponding to the centers of the N
%           ellipses or ellipsoids that should be drawn.
%        A  A 2X2XN or 3X3XN set of N positive definite matrices that
%           specify the size and shape of the ellipse or ellipsoids, where
%           a point zp is on the ith ellipse/ ellipsoid if
%           (zp-z(:,i))'*A(:,:,i)*(zp-z(:,i))=gamma.
%           In three dimensions, A can have one zero eigenvalue, in which
%           case the fill3 command is used to draw an ellipse in 3D.
%     gamma An optional parameter specifying the size of the ellipse/
%           ellipsoid. If omitted or an empty matrix is passed, then
%           gamma=16.2235 in 2D or gamma=18.8049 in 3D is used. gamma must be
%           positive.
%  varargin Sets of values that should be passed to the plot function to
%           format the ellipses or that will be passed to the surf function
%           to format the ellipsoids. For example, one could call the
%           function as drawEllipse(z,A,gamma,'--r','linewidth',2) to plot
%           ellipses as thick, red lines. Note that Matlab will not always
%           properly render dashed lines due to the number of points used
%           to plot the shape. Also, if the ellipsoid is in 3D, but is
%           actually just a 2D ellipse, then if this parameter is omitted,
%           the fill3 command will be passed the option to color the
%           ellipse black. Otherwise, if this parameter is given, the user
%           must specify the ellipse color as it is not done here.
%
%OUTPUTS: None. The results are plotted.
%
%An eigendecomposition of a matrix A breaks it into parts V and D such that
%V*D*V'=A, where V is a rotation matrix and D is a diagonal matrix.
%Considering the problem here, this means that (zp-z)'*A*(zp-z)=
%(V'*zp-V'*z)'*D*(V'*zp-V'*z). If we substitute kp=V'*zp and k=V'*z, then
%the equation is (kp-k)'*D*(kp-k)=gamma, which is the equation for an
%ellipse where the axes are aligned with the coordinate axes. Thus, one can
%first find the points for an ellipse that is aligned with the coordinate
%axes (the kp points), and then rotate it back to the proper alignment (the
%zp points). However, it is simpler to first find an ellipse centered about
%zero. Thus, substitute l=kp-k and find the l points, then shift and rotate
%them to get the z-points.
%
%When considering an ellipse, suppose that l=[x;y] and D=diag(a,b). Then
%the ellipse equation is
%x^2*a+y^2*b=gamma
%Solving for y in terms of x, we have
%y=+/-sqrt((gamma-x^2*a)/b).
%The x values are limited so that the argument of the square root is
%positive. Thus the x values range from -sqrt(gamma/a) to +sqrt(gamma/a).
%Thus, a simple way to plot this centered, axis-aligned ellipse is to
%generate x-values in the valid range and then find the two sets of y
%values. This function does that and transforms the results to get the
%corresponding z values. 1000 points are used in the plot.
%
%To find the points on an axis-aligned ellipsoid, we will use the ellipsoid
%function that is built into Matlab. Suppose that l=[x;y;z] and
%D=diag(a,b,c). Then the ellipsoid equation is x^2*a+y^2*b+z^2c=gamma. This
%means that the semi-axis lengths of this centered, axis-aligned ellipsoid
%are sqrt(gamma/a), sqrt(gamma/b), and sqrt(gamma/c). 2501 vertices are
%used in the ellipsoids.
%
%July 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numEllipse=size(z,2);

%Save the value of hold on the plot so that it can be reset to its previous
%state after plotting all of the ellipses.
holdVal=ishold();

numDim=size(z,1);

switch(numDim)
    case 2
        if(nargin<3||isempty(gamma))
            gamma=16.2235;
        end

        numPoints=1000;
        
        for curEllip=1:numEllipse
            %Perform an eigenvalue decomposition of A.
            [V,D]=eig(A(:,:,curEllip));
            
            a=D(1,1);
            b=D(2,2);

            xBound=sqrt(gamma/a);
            x=linspace(-xBound,xBound,numPoints);
            %The real command deals with possible finite precision issues.
            y=real(sqrt((gamma-x.^2*a)/b));
            %The centered, axis-aligned ellipse points.
            l=[x,fliplr(x);
               y,-fliplr(y)];
            zp=bsxfun(@plus,V*l,z(:,curEllip));

            %Make sure that all of the ellipses are printed.
            if(curEllip~=1)
                hold on
            end
            plot(zp(1,:),zp(2,:),varargin{:});
        end
    case 3
        if(nargin<3||isempty(gamma))
            gamma=18.8049;
        end

        for curEllip=1:numEllipse
            if(rank(A)<3)
                numPoints=500;
                %If the ellipsoid is actually just an ellipse, then we will
                %draw the 2D ellipse in 3D.
                
                %Perform an eigenvalue decomposition of A.
                [V,D]=eig(A(:,:,curEllip));
                %We will put the zero eigenvalue at the statrt, if it is
                %not already there.
                [~,idx]=sort(diag(D),'ascend');
                V=V(idx,idx);
                D=D(idx,idx);
                
                a=D(2,2);
                b=D(3,3);
                
                xBound=sqrt(gamma/a);
                x=linspace(-xBound,xBound,numPoints);
                %The real command deals with possible finite precision
                %issues.
                y=real(sqrt((gamma-x.^2*a)/b));
                %The centered, axis-aligned ellipse points.
                l=[zeros(1,2*numPoints);
                   x,fliplr(x);
                   y,-fliplr(y)];
                zp=bsxfun(@plus,V*l,z(:,curEllip));
                %Make sure that all of the ellipses are printed.
                if(curEllip~=1)
                    hold on
                end
                %Draw the ellipse as a filled shape.
                if(isempty(varargin))
                    fill3(zp(1,:),zp(2,:),zp(3,:),'k');
                else
                    fill3(zp(1,:),zp(2,:),zp(3,:),varargin{:});
                end
            else
                numPoints=50;
                %Perform an eigenvalue decomposition of A.
                [V,D]=eig(A(:,:,curEllip));
                zCur=z(:,curEllip);

                [xp,yp,zp] = ellipsoid(0,0,0,sqrt(gamma/D(1,1)),sqrt(gamma/D(2,2)),sqrt(gamma/D(3,3)),numPoints);

                %Rotate to the correct orientation
                xyzPoints=V*[xp(:)';yp(:)';zp(:)'];
                xp=reshape(zCur(1)+xyzPoints(1,:),numPoints+1,numPoints+1);
                yp=reshape(zCur(2)+xyzPoints(2,:),numPoints+1,numPoints+1);
                zp=reshape(zCur(3)+xyzPoints(3,:),numPoints+1,numPoints+1);

                %Make sure that all of the ellipses are printed.
                if(curEllip~=1)
                    hold on
                end

                surf(xp,yp,zp,varargin{:});
            
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
