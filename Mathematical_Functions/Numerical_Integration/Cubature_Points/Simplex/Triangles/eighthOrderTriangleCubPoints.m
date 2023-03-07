function [xi,w]=eighthOrderTriangleCubPoints(algorithm)
%%EIGHTHORDERTRIANGLECUBPOINTS Obtain eighth-order cubature points for
%   integration over a triangle in 2D. The points and weights are for the
%   triangle with vertices (1,0), (0,1), (0,0), but can be transformed to
%   any triangle using transformSimplexTriPoints.
%
%INPUTS: algorithm An optional parameter selecting the algorithm for the
%                  specific point set. Possible values are:
%                  0 (The default if omitted or an empty matrix is passed)
%                    Use the algorithm of [1] (16 points).
%                  1 Use the 8th order points found using the algorithm of
%                    [2], given in the supplementary material of [2] (16
%                    points).
%
%OUTPUTS: xi A 2XnumCubPoints set of points for the standard triangle.
%          w A 1XnumCubPoints set of cubature weights. This sums to the
%            volume of the triangle (1/2).
%
%This function implements the points given in [1] and [2].
%
%EXAMPLE:
%Given the vertices of the simplex, we compare an eighth-order moment
%computed using these cubature points to one computed using
%monomialIntSimplex. The results are the same within typical finite
%precision limits.
% [xi,w]=eighthOrderTriangleCubPoints();
% alpha=[6;2];
% theMoment=findMomentFromSamp(alpha,xi,w)
% intVal=monomialIntSimplex(alpha)
%
%REFERENCES:
%[1]L. Zhang and T. Cui, "A set of symmetric quadrature rules on triangles
%   and tetrahedra," Journal of Computational Mathematics, vol. 27, no. 1,
%   pp. 89-96, Jan. 2009.
%[2] F. D. Witherden and P. E. Vincent, "On the identification of symmetric
%    quadrature rules for finite element methods," Computer and Mathematics
%    with Applications, vol. 69, no. 10, pp. 1232-1241, May 2015.
%
%October 2022 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<1||isempty(algorithm))
    algorithm=0;
end

switch(algorithm)
    case 0
        w1=0.1443156076777871682510911104890646;
        w2=0.1032173705347182502817915502921290;
        w3=0.0324584976231980803109259283417806;
        w4=0.0950916342672846247938961043885843;
        w5=0.0272303141744349942648446900739089;
        w=[w1;w2;w2;w2;w3;w3;w3;w4;w4;w4;w5;w5;w5;w5;w5;w5];
        
        p1=0.3333333333333333333333333333333333;
        p2=0.1705693077517602066222935014914645;
        p3=0.0505472283170309754584235505965989;
        p4=0.4592925882927231560288155144941693;
        p5a=0.2631128296346381134217857862846436;
        p5b=0.0083947774099576053372138345392944;
        
        xiBary=zeros(3,16);
        xiBary(:,1)=[p1;p1;p1];
        xiBary(:,2:4)=genAllMultisetPermutations([p2;p2;1-2*p2]);
        xiBary(:,5:7)=genAllMultisetPermutations([p3;p3;1-2*p3]);
        xiBary(:,8:10)=genAllMultisetPermutations([p4;p4;1-2*p4]);
        xiBary(:,11:16)=genAllMultisetPermutations([p5a;p5b;1-p5a-p5b]);
        
        %Adjust w for the area of the standard triangle.
        w=w/2;
        
        %Convert the barycentric points into normal cubature points for the
        %standard triangle given the vertices.
        vertices=[1,0,0;
                  0,1,0];
        xi=barycentricCoords2Pt(xiBary,vertices);
    case 1
        M=[-0.33333333333333333333333333333333333333,   -0.33333333333333333333333333333333333333,     0.2886312153555743365021822209781292496;
          -0.081414823414553687942368971011661355879,   -0.83717035317089262411526205797667728824,     0.1901832685345692495877922087771686332;
           -0.83717035317089262411526205797667728824,  -0.081414823414553687942368971011661355879,     0.1901832685345692495877922087771686332;
          -0.081414823414553687942368971011661355879,  -0.081414823414553687942368971011661355879,     0.1901832685345692495877922087771686332;
           -0.65886138449647958675541299701707099796,    0.31772276899295917351082599403414199593,    0.20643474106943650056358310058425806003;
            0.31772276899295917351082599403414199593,   -0.65886138449647958675541299701707099796,    0.20643474106943650056358310058425806003;
           -0.65886138449647958675541299701707099796,   -0.65886138449647958675541299701707099796,    0.20643474106943650056358310058425806003;
           -0.89890554336593804908315289880680210631,    0.79781108673187609816630579761360421262,   0.064916995246396160621851856683561193593;
            0.79781108673187609816630579761360421262,   -0.89890554336593804908315289880680210631,   0.064916995246396160621851856683561193593;
           -0.89890554336593804908315289880680210631,   -0.89890554336593804908315289880680210631,   0.064916995246396160621851856683561193593;
           -0.98321044518008478932557233092141110162,    0.45698478591080856248200075835212392604,    0.05446062834886998852968938014781784832;
            0.45698478591080856248200075835212392604,   -0.98321044518008478932557233092141110162,    0.05446062834886998852968938014781784832;
           -0.47377434073072377315642842743071282442,    0.45698478591080856248200075835212392604,    0.05446062834886998852968938014781784832;
            0.45698478591080856248200075835212392604,   -0.47377434073072377315642842743071282442,    0.05446062834886998852968938014781784832;
           -0.47377434073072377315642842743071282442,   -0.98321044518008478932557233092141110162,    0.05446062834886998852968938014781784832;
           -0.98321044518008478932557233092141110162,   -0.47377434073072377315642842743071282442,    0.05446062834886998852968938014781784832];

        w=M(:,3);
        xi=M(:,1:2)';
        %Transform the points to the standard triangle.
        v1=[-1,-1, 1;
            -1, 1,-1];
        v2=[1,0,0;
            0,1,0];
        [A,d]=affineTransBetweenTriangles(v1,v2);
        xi=bsxfun(@plus,A*xi,d);
        w=w/4;
    otherwise
        error('Unknown Algorithm Specified')
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
