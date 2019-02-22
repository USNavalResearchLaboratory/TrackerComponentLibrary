function [xi,w]=thirdOrderNDimCubPoints(numDim,algorithm)
%%THIRDORDERNDIMCUBPOINTS Generate third-order cubature points for
%               integration over a numDim-dimensional cube with bounds in
%               coordinates of (-1,1).
%
%INPUTS: numDim An integer specifying the dimensionality of the points to
%               be generated. numDim>1.
%     algorithm An optional parameter specifying the algorithm to be
%               used to generate the points. Possible values are:
%               0  (The default if omitted or an empty matrix is passed)
%                  Formula Cn 3-1 in [1], pg. 230, 2*numDim points, with
%                  the correction listed in Table I of [2].
%               1  Formula Cn 3-3 in [1], pg. 230, 2*numDim+1 points.
%               2  Formula Cn 3-4 in [1], pg. 231, 2^numDim points.
%               3  Formula Cn 3-6 in [1], pg. 231, 3^numDim points.
%               4  Formula Cn 3-6 in [1], pg. 231, 3^numDim points.
%               5  Formula C2 3-1 in [1], pg. 243, 4 points, numDim=2,
%                  with the correction listed in Table I of [2].
%               6  Formula C2 3-3 in [1], pg. 244, 9 points, numDim=2.
%               7  Formula C2 3-4 in [1], pg. 245, 9 points, numDim=2.
%               8  Formula C3 3-1 in [1], pg. 261, 6 points, numDim=3.
%               9  Formula C3 3-3 in [1], pg. 261, 9 points, numDim=3.
%               10  Formula C3 3-4 in [1], pg.261, 13 points, numDim=3.
%               11 Formula C3 3-5 in [1], pg. 262, 15 points, numDim=3.
%               12 Formula C3 3-6 in [1], pg. 262, 19 points, numDim=3.
%               13 Formula C3 3-7 in [1], pg. 263, 20 points, numDim=3.
%
%OUTPUTS: xi A numDim X numCubaturePoints matrix containing the cubature
%            points (Each "point" is a vector).
%          w A numCubaturePoints X 1 vector of the weights associated with
%            the cubature points.
%
%REFERENCES:
%[1] A.H. Stroud, Approximate Calculation of Multiple Integrals. Cliffs,
%    NJ: Prentice-Hall, Inc., 1971.
%[2] R. Cools, "An encyclopedia of cubature formulas," Journal of
%    Complexity, vol. 19, no. 3, pp. 445-453, Jun. 2003.
%
%October 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(algorithm))
    algorithm=0; 
end

switch(algorithm)
    case 0%Cn 3-1 in [1], pg. 230, 2*numDim points.
        V=2^numDim;
        
        w=(1/(2*numDim))*ones(2*numDim,1)*V;
        xi=zeros(numDim,2*numDim);
        i=1:(2*numDim);
        for k=1:fix(numDim/2)
            xi(2*k-1,i)=sqrt(2/3)*cos((2*k-1)*i*pi/numDim);
            xi(2*k,i)=sqrt(2/3)*sin((2*k-1)*i*pi/numDim);
        end
        
        if(mod(numDim,2)~=0)
            xi(numDim,i+1)=(-1).^i/sqrt(3);
        end
    case 1%Cn 3-3 in [1], pg. 230, 2*numDim+1 points.
        V=2^numDim;
        
        B0=V*(3-numDim)/3;
        B1=V/6;
        
        xi=[zeros(numDim,1),fullSymPerms([1;zeros(numDim-1,1)])];
        w=[B0;B1*ones(2*numDim,1)];
    case 2%Cn 3-4 in [1], pg. 231, 2^numDim points.
        V=2^numDim;
        
        r=sqrt(3)/3;
        w=V/2^numDim*ones(2^numDim,1);
        xi=PMCombos(r*ones(numDim,1));
    case 3%Cn 3-5 in [1], pg. 231, 2^numDim+1 points.
        V=2^numDim;
        B0=(2/3)*V;
        B1=(1/(3*2^numDim))*V;
        
        xi=[zeros(numDim,1),PMCombos(ones(numDim,1))];
        w=[B0;B1*ones(2^numDim,1)];
    case 4%Cn 3-6 in [1], pg. 231, 3^numDim points.
        r=[-1;0;1];
        A=[1/3;4/3;1/3];
        
        xi=zeros(numDim,3^numDim);
        w=zeros(3^numDim,1);
        
        dims=3*ones(numDim,1);
        for curIdx=1:(3^numDim)
            idxVec=index2NDim(dims,curIdx);
            
            xi(:,curIdx)=r(idxVec);
            w(curIdx)=prod(A(idxVec));
        end
    case 5%Formula C2 3-1 in [1], pg. 243, 4 points, numDim=2.
        if(numDim~=2)
           error('The selected algorithm requires that numDim=2') 
        end
        
        r=1/3;
        xi=PMCombos([sqrt(r);sqrt(r)]);
        
        
        V=2^numDim;
        w=(1/4)*V*ones(4,1);
    case 6%C2 3-3 in [1], pg. 244, 9 points, numDim=2.
        if(numDim~=2)
           error('The selected algorithm requires that numDim=2') 
        end
        
        V=2^numDim;
        xi=[[0;0],fullSymPerms([1;0]),PMCombos([1;1])];
        w=[4/9;1/9*ones(4,1);1/36*ones(4,1)]*V;
    case 7%C2 3-4 in [1], pg. 245, 9 points, numDim=2.
        if(numDim~=2)
           error('The selected algorithm requires that numDim=2') 
        end
        
        V=2^numDim;
        xi=[[0;0],fullSymPerms([1;0]),PMCombos([1;1])];
        w=[5/12;1/8*ones(4,1);1/48*ones(4,1)]*V;
    case 8%C3 3-1 in [1], pg. 261, 6 points, numDim=3.
        if(numDim~=3)
           error('The selected algorithm requires that numDim=3') 
        end
        
        V=2^numDim;
        xi=fullSymPerms([1;0;0]);
        w=V/6*ones(6,1);
    case 9%C3 3-3 in [1], pg. 261, 9 points, numDim=3.
        if(numDim~=3)
           error('The selected algorithm requires that numDim=3') 
        end
        
        V=2^numDim;

        xi=[zeros(3,1),PMCombos([1;1;1])];
        w=[2/3;1/24*ones(8,1)]*V;
    case 10%C3 3-4 in [1], pg.261, 13 points, numDim=3.
        if(numDim~=3)
           error('The selected algorithm requires that numDim=3') 
        end
        
        V=2^numDim;
        xi=[zeros(3,1),fullSymPerms([1;1;0])];
        w=[1/2;1/24*ones(12,1)]*V;
    case 11%C3 3-5 in [1], pg. 262, 15 points, numDim=3.
        if(numDim~=3)
           error('The selected algorithm requires that numDim=3') 
        end
        
        V=2^numDim;
        xi=[[0;0;0],fullSymPerms([1;0;0]),PMCombos([1;1;1])];
        w=[2/9;1/9*ones(6,1);1/72*ones(8,1)]*V;
    case 12%C3 3-6 in [1], pg. 262, 19 points, numDim=3.
        if(numDim~=3)
           error('The selected algorithm requires that numDim=3') 
        end
        
        V=2^numDim;
        xi=[[0;0;0],fullSymPerms([1;0;0]),fullSymPerms([1;1;0])];
        w=[1/4;1/12*ones(6,1);1/48*ones(12,1)]*V;
    case 13%C3 3-7 in [1], pg. 263, 20 points, numDim=3.
        if(numDim~=3)
           error('The selected algorithm requires that numDim=3') 
        end
        
        V=2^numDim;
        xi=[fullSymPerms([1;1;0]),PMCombos([1;1;1])];
        w=[1/6*ones(12,1);-1/8*ones(8,1)]*V;
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
