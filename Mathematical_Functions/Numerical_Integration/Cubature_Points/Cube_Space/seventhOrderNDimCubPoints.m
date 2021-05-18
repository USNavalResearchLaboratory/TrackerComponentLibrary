function [xi,w]=seventhOrderNDimCubPoints(numDim,algorithm)
%%SEVENTHORDERNDIMCUBPOINTS Generate seventh-order cubature points for
%               integration over a numDim-dimensional cube with bounds in
%               coordinates of (-1,1).
%
%INPUTS: numDim An integer specifying the dimensionality of the points to
%               be generated. Currently, only numDim=2 and numDim=3 are
%              supported.
%    algorithm An optional parameter specifying the algorithm to be
%              used to generate the points. Possible values are:
%              0 (The default if omitted or an empty matrix is passed and
%                numDim=2) Formula C2 7-1 in [1], pg. 252, 12 points,
%                numDim=2.
%              1 Formula C2 7-3 in [1], pg. 253, 13 points, numDim=2, with
%                the correction in Table I of [2].
%              2 Formula C2 7-4 in [1], pg. 255, 16 points, numDim=2.
%              3 Formula C2 7-5 in [1], pg. 255, 21 points, numDim=2.
%              4 Formula C2 7-6 in [1], pg. 255, 25 points, numDim=2.
%              5 Formula C3 7-2 in [1], pg. 265, 34 points, numDim=3.
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
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2)
   if(numDim==2)
       algorithm=0;
   else
       algorithm=5;
   end
end

switch(algorithm)
    case 0%Formula C2 7-1 in [1], pg. 252, 12 points, numDim=2.
        if(numDim~=2)
           error('The selected algorithm requires that numDim=2') 
        end
        
        V=2^numDim;
        
        r=sqrt(6/7);
        s=sqrt((114-3*sqrt(583))/287);
        t=sqrt((114+3*sqrt(583))/287);
        
        xi=[fullSymPerms([r;0]),PMCombos([s;s]),PMCombos([t;t])];
        
        B1=(49/810)*V;
        B2=((178981+2769*sqrt(583))/1888920)*V;
        B3=((178981-2769*sqrt(583))/1888920)*V;
        
        w=[B1*ones(4,1);B2*ones(4,1);B3*ones(4,1)];
    case 1%Formula C2 7-3 in [1], pg. 253, 13 points, numDim=2.
        if(numDim~=2)
           error('The selected algorithm requires that numDim=2') 
        end
        V=2^numDim;
        
        r=sqrt(12/35);
        s=sqrt((93+3*sqrt(186))/155);
        t=sqrt((93-3*sqrt(186))/155);
        
        xi=[[0;0],fullSymPerms([r;0]),fullSymPerms([s;t])];
        w=[(1/81)*V;(49/324)*V*ones(4,1);(31/648)*V*ones(8,1)];
    case 2%Formula C2 7-4 in [1], pg. 255, 16 points, numDim=2.
        if(numDim~=2)
           error('The selected algorithm requires that numDim=2') 
        end
        V=2^numDim;
        
        r=sqrt((15-2*sqrt(30))/35);
        s=sqrt((15+2*sqrt(30))/35);
        
        xi=[PMCombos([r;r]),PMCombos([s;s]),fullSymPerms([r;s])];
        
        B1=((59+6*sqrt(30))/864)*V;
        B2=((59-6*sqrt(30))/864)*V;
        B3=(49/864)*V;
        
        w=[B1*ones(4,1);B2*ones(4,1);B3*ones(8,1)];
    case 3%Formula C2 7-5 in [1], pg. 255, 21 points, numDim=2.
        if(numDim~=2)
           error('The selected algorithm requires that numDim=2') 
        end
        V=2^numDim;
        
        r=2/3;
        s=1/3;
        t=1/2;
        
        xi=[[0;0],fullSymPerms([1;0]),fullSymPerms([r;0]),fullSymPerms([s;0]),PMCombos([1;1]),PMCombos([t;t])];
        w=[(449/315)*V;
           (37/1260)*V*ones(4,1);
           (3/28)*V*ones(4,1);
           (-69/140)*V*ones(4,1);
           (7/540)*V*ones(4,1);
           (32/135)*V*ones(4,1)];
        
    case 4%Formula C2 7-6 in [1], pg. 255, 25 points, numDim=2.
        if(numDim~=2)
           error('The selected algorithm requires that numDim=2') 
        end
        V=2^numDim;
        
        r=2/3;
        s=1/3;
        
        xi=[[0;0],PMCombos([r;r]),fullSymPerms([r;0]),PMCombos([s;s]),fullSymPerms([1;s]),PMCombos([1;1])];
        
        w=[(1024/6720)*V;
           (576/6720)*V*ones(8,1);
           (-9/6720)*V*ones(4,1);
           (117/6720)*V*ones(8,1);
           (47/6720)*V*ones(4,1)];
    case 5%Formula C3 7-2 in [1], pg. 265, 34 points, numDim=3.
        if(numDim~=3)
           error('The selected algorithm requires that numDim=3') 
        end
        V=2^numDim;
        
        r=sqrt(6/7);
        s=sqrt((960-3*sqrt(28798))/2726);
        t=sqrt((960+3*sqrt(28798))/2726);
        
        xi=[fullSymPerms([r;0;0]),fullSymPerms([r;r;0]),PMCombos([s;s;s]),PMCombos([t;t;t])];

        B1=(1078/29160)*V;
        B2=(343/29160)*V;
        B3=((774*t^2-230)/(9720*(t^2-s^2)))*V;
        B4=((230-774*s^2)/(9720*(t^2-s^2)))*V;
        
        w=[B1*ones(6,1);B2*ones(12,1);B3*ones(8,1);B4*ones(8,1)];
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
