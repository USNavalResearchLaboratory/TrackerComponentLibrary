function [xi,w]=fifthOrderSpherSurfCubPoints(numDim,algorithm)
%%FIFTHORDERSPHERSURFCUBPOINTS Generate fifth-order cubature points for
%               integration over the surface of a unit (hyper-)sphere
%               (weighting function is just 1).
%
%INPUTS: numDim An integer specifying the dimensionality of the points
%               to be generated.
%     algorithm A value indicating which algorithm should be used.
%               Possible values are
%               0 (The default if omitted or an empty matrix is passed)
%                 Use algorithm Un 5-1 from [1]. This requires 2*numDim^2
%                 points. Negative weights present for numDim>4.
%               1 Algorithm Un 5-2 of [1] requiring 2^numDim+2*numDim
%                 points.
%               2 Algorithm Un 5-3 of [1] requiring 2^(numDim+1)-2 points.
%               3 Algorithm Un 5-4 of [1] requiring numDim*2^numDim
%                 points.
%               4 The algorithm of [2] as taken from [3] requiring
%                 (numDim+1)*(numDim+2) points. Negative weights present
%                 for numDim>7.
%               5 Algorithm U3 5-1 of [1], requiring 12 points and that
%                 numDim=3.
%               6 Algorithm U3 5-2 of [1], requiring 14 points and that
%                 numDim=3.
%               7 Algorithm U3 5-3 of [1], requiring 18 points and that
%                 numDim=3.
%               8 Algorithm U3 5-4 of [1], requiring 20 points and that
%                 numDim=3.
%               9 Algorithm U3 5-5 of [1], requiring 30 points and that
%                 numDim=3.
%
%OUTPUTS: xi A numDim X numCubaturePoints matrix containing the cubature
%            points. (Each "point" is a vector)
%          w A numCubaturePoints X 1 vector of the weights associated with
%            the cubature points.
%
%REFERENCES:
%[1] A.H. Stroud, Approximate Calculation of Multiple Integrals. Cliffs,
%    NJ: Prentice-Hall, Inc., 1971.
%[2] I. P. Mysovskikh, The approximation of multiple integrals by using
%    interpolatory cubature formulae, in Quantitative Approximation, R. A.
%    DeVore and K. Scherer, eds., Academic Press, New York, 1980, pp.
%    217-243.
%[3] A. Genz and J. Monahan, "Stochastic integration rules for infinite
%    regions," SIAM Journal on Scientific Computing, vol. 19, no. 2, pp.
%    426-439, Mar. 1998.
%
%October 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(algorithm))
    algorithm=0;
end

switch(algorithm)
    case 0%Stroud, algorithm 5-1, 2*numDim^2 points
        r=1;
        s=1/sqrt(2);
        V=2*pi^(numDim/2)/gamma(numDim/2);
        B1=(4-numDim)/(2*numDim*(numDim+2))*V;
        B2=1/(numDim*(numDim+2))*V;
        
        w=[repmat(B1,[2*numDim,1]);repmat(B2,[2*(numDim^2-numDim),1])];
        
        xi=zeros(numDim,2*numDim^2);
        
        %The first half of the points
        curPoint=1;
        for curDim=1:numDim
            xi(curDim,curPoint)=r;
            xi(curDim,curPoint+1)=-r;
            curPoint=curPoint+2;
        end
        
        %The second half of the points
        for curIdx1=1:(numDim-1)
            for curIdx2=(curIdx1+1):numDim
                sel=[curIdx1;curIdx2];
                xi(sel,curPoint)=[s;s];
                xi(sel,curPoint+1)=[s;-s];
                xi(sel,curPoint+2)=[-s;s];
                xi(sel,curPoint+3)=[-s;-s];

                curPoint=curPoint+4;
            end
        end
    case 1%Stroud, algorithm 5-2, 2^numDim+2*numDim points
        r=1;
        s=1/sqrt(numDim);
        V=2*pi^(numDim/2)/gamma(numDim/2);
        B1=1/(numDim*(numDim+2))*V;
        B2=numDim/(2^numDim*(numDim+2))*V;
        
        w=[repmat(B1,[2*numDim,1]);repmat(B2,[2^numDim,1])];
        
        xi=zeros(numDim,2*numDim+2^numDim);
        %The first half of the points
        curPoint=1;
        for curDim=1:numDim
            xi(curDim,curPoint)=r;
            xi(curDim,curPoint+1)=-r;
            curPoint=curPoint+2;
        end
        
        xi(:,curPoint:end)=PMCombos(repmat(s,[numDim,1]));
    case 2%Stroud, algorithm 5-3, 2^(numDim+1)-2 points
        numPoints=2^(numDim+1)-2;
        
        w=zeros(numPoints,1);
        xi=zeros(numDim,numPoints);
        
        s=1/sqrt(numDim+2);
        V=2*pi^(numDim/2)/gamma(numDim/2);
        
        valVec=s*ones(numDim,1);
        curStart=1;
        for k=1:numDim
            r=sqrt((k+2)/(numDim+2));
            B=V*2^(k-numDim)*(numDim+2)/(numDim*(k+1)*(k+2));
            valVec(k)=r;
            
            points=PMCombos(valVec(k:end));
            num2Add=size(points,2);
            
            xi(k:end,curStart:(curStart+num2Add-1))=points;
            w(curStart:(curStart+num2Add-1))=B;
            
            curStart=curStart+num2Add;
        end  
    case 3%Stroud, algorithm 5-4, numDim*2^numDim points
        u=sqrt((numDim+2+(numDim-1)*sqrt(2*(numDim+2)))/(numDim*(numDim+2)));
        v=sqrt((numDim+2-sqrt(2*(numDim+2)))/(numDim*(numDim+2)));
        V=2*pi^(numDim/2)/gamma(numDim/2);
        
        numPoints=numDim*2^numDim;
        
        w=repmat(V/numPoints,[numPoints,1]);
        
        xi=zeros(numDim,numPoints);
        
        vec=repmat(v,[numDim,1]);
        sel=1:(2^numDim);
        for curSpot=1:numDim
            vecCur=vec;
            vecCur(curSpot)=u;
            
            xi(:,sel)=PMCombos(vecCur);
            sel=sel+2^numDim;
        end
    case 4%Mysovskikh, (numDim+1)*(numDim+2) point.
        Um=2*pi^(numDim/2)/gamma(numDim/2);
        
        numPoints=(numDim+1)*(numDim+2);
        xi=zeros(numDim,numPoints);
        w=zeros(numPoints,1);
        
        %The first sum
        w(1:(2*(numDim+1)))=Um*(7-numDim)*numDim/(2*(numDim+1)^2*(numDim+2));
        %Get the vertices of an n-simplex on the surface of the unit
        %sphere.
        v=regularNSimplexCoords(numDim);
        xi(:,1:(numDim+1))=-v;
        xi(:,(numDim+2):(2*(numDim+1)))=v;
        
        %The second sum. We have to get the edges of the m-simplex. An edge
        %exists between each vertex.
        curPoint=2*(numDim+1)+1;
        w(curPoint:end)=Um*2*(numDim-1)^2/(numDim*(numDim+1)^2*(numDim+2));
        
        for vertex1=1:numDim
           for vertex2=(vertex1+1):(numDim+1)
               y=(v(:,vertex1)+v(:,vertex2))/2;%Midpoint of the edge
               y=y/norm(y);%Project onto the unit sphere.
               
               xi(:,curPoint)=-y;
               xi(:,curPoint+1)=y;
               curPoint=curPoint+2;
           end
        end
    case 5%Algorithm U3 5-1 of [1], requiring 12 points
        if(numDim~=3)
           error('This formula requires that numDim=3') 
        end
        
        xi=regularIcosahedronCoords();
        V=2*pi^(3/2)/gamma(3/2);
        w=(V/12)*ones(12,1);
    case 6%Algorithm U3 5-2 of [1], requiring 14 points
        if(numDim~=3)
           error('This formula requires that numDim=3') 
        end
        
        r=1;
        s=1/sqrt(3);
        V=2*pi^(3/2)/gamma(3/2);
        B1=(8/120)*V;
        B2=(9/120)*V;

        xi=[fullSymPerms([r;0;0]),PMCombos([s;s;s])];
        w=[B1*ones(6,1);B2*ones(8,1)];
    case 7%Algorithm U3 5-3 of [1], requiring 18 points 
        if(numDim~=3)
           error('This formula requires that numDim=3') 
        end
        
        r=1;
        s=1/sqrt(2);
        V=2*pi^(3/2)/gamma(3/2);
        B1=(1/30)*V;
        B2=(2/30)*V;

        xi=[fullSymPerms([r;0;0]),fullSymPerms([s;s;0])];
        w=[B1*ones(6,1);B2*ones(12,1)];
    case 8%Algorithm U3 5-5 of [1], requiring 30 points
        if(numDim~=3)
           error('This formula requires that numDim=3') 
        end
        
        V=2*pi^(3/2)/gamma(3/2);
        r=1/2;
        s=(sqrt(5)+1)/4;
        t=(sqrt(5)-1)/4;
        B=V/30;
        
        xi=[fullSymPerms([1;0;0]),PMCombos([r;s;t]),PMCombos([t;r;s]),PMCombos([s;t;r])];
        w=B*ones(30,1);
    otherwise
        error('Unknown algorithm specified');
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
