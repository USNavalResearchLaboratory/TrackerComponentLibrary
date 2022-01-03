function vals=randPointInEllipsoid(N,invS,z0,algorithm,gamma)
%%RANDPOINTINELLIPSOID Generate a random point z uniformly distributed in
%              a circle, sphere, ellipse or ellipsoid such that
%              (z-z0)'*invS*(z-z0)<=gamma
%
%INPUTS: N The number of random points in the ellipse to generate.
%     invS A numDimXnumDim positive definite matrix that specifies the
%          size and shape of the ellipse or ellipsoid, where
%          a point z is in the ellipse if (z-z0)'*invS*(z-z0)<=gamma.
%       z0 The optional numDimX1 vector defining the center of the ellipse.
%          If this parameter is omitted or an empty matrix is passed, then 
%          an all zero vector will be used.
% algorithm An optional parameter specifying the algorithm to used.
%          Possible values are:
%          0 Use the algorithm of Section 3.1.1 of [1], modified for more
%            than two dimensions.
%          1 Use algorithm B1 of Section 3.2.1 of [1].
%          2 Use the algorithm of Section 3.3.1 of [1]. This is the default
%            if this parameter is omitted or an empty matrix is passed.
%    gamma The optional threshold for determining the size of the ellipse.
%          If this parameter is omitted or an empty matrix is passed, then
%          the default of 1 will be used.
%
%OUTPUTS: vals A numDimXN set of random samples within the ellipse.
%
%All of the algorithms are taken from [1].
%
%EXAMPLE 1:
%Here, we use the 2D ellipse of Section 4 of [1] with the similar
%parameters. The random points generated from each method are displayed.
% S=[10,6;
%     6,10];
% invS=inv(S);
% z0=[100;100];
% gamma=9.2103;
% N=231;
% algorithm=0;
% tic;vals0=randPointInEllipsoid(N,invS,z0,algorithm,gamma);toc
% algorithm=1;
% tic;vals1=randPointInEllipsoid(N,invS,z0,algorithm,gamma);toc
% algorithm=2;
% tic;vals2=randPointInEllipsoid(N,invS,z0,algorithm,gamma);toc
% figure(1)
% clf
% hold on
% scatter(vals0(1,:),vals0(2,:),'.r')
% scatter(vals1(1,:),vals1(2,:),'.g')
% scatter(vals2(1,:),vals2(2,:),'.b')
%
%REFERENCES:
%[1] H. Sun and M. Farooq, "Note on the generation of random points
%    uniformly distributed in hyper-ellipsoids," in Proceedings of the
%    Fifth International Conference on Information Fusion, Annapolis, MD,
%    8-11 Jul. 2002, pp. 489-496.
%
%October 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numDim=size(invS,1);

if(nargin>4&&~isempty(gamma))
    invS=invS/gamma; 
end

if(nargin<4||isempty(algorithm))
    algorithm=2;
end

if(nargin<3||isempty(z0))
    z0=zeros(numDim,1);
end

vals=zeros(numDim,N);
switch(algorithm)
    case 0%The algorithm of Section 3.1.1 of [1], modified for more than
          %two dimensions.
        S=inv(invS);
        minBounds=-sqrt(diag(S));
        maxBounds=sqrt(diag(S));
        span=maxBounds-minBounds;
        
        for curVal=1:N
            while(1)
                x=minBounds+span.*rand(numDim,1);
                if(x'*invS*x<=1)
                    vals(:,curVal)=z0+x;
                    break
                end
            end
        end
    case 1%Algorithm B1 of Section 3.2.1 of [1].
        S=inv(invS);
        for curVal=1:N
            theta=2*pi*rand(1);
            nu=rand(1);
            rho=nthroot(nu,numDim);

            phi=zeros(numDim-2,1);
            for i=1:(numDim-2)
                phi(i)=SinPowD.rand(1,numDim-i-1);
            end

            y=hypersphere2Cart([rho;phi;theta]);
            A=chol(S,'lower');

            vals(:,curVal)=z0+A*y;
        end
    case 2%The algorithm of Section 3.3.1 of [1].
        %All eigenvalues must be positive and real.
        [L,Lambda]=eig(invS);
        lambda=diag(Lambda);
        invLambdaRoots=1./sqrt(lambda);
        
        for curVal=1:N
            while(1)
                x=UniformD.rand(1,[-invLambdaRoots';invLambdaRoots']);

                if(x'*Lambda*x<=1)
                    vals(:,curVal)=z0+L*x;
                    break;
                end
            end
        end
    otherwise
        error('Unknown algorithm specified')
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
