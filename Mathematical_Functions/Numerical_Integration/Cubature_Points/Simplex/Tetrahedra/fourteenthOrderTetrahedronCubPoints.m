function [xi,w]=fourteenthOrderTetrahedronCubPoints()
%%FOURTEENTHORDERTETRAHEDRONCUBPOINTS Obtain fourteenth-order cubature
%   points for integration over a tetrahedron in 3D. The points and weights
%   are for the tetrahedron with vertices (1,0,0), (0,1,0), (0,0,1), and
%   (0,0,0), but can be transformed to any tetrahedron using
%   transformSimplexTetrahedronPts.
%
%INPUTS: None
%
%OUTPUTS: xi This is a 3XnumCubPoints set of points for the standard
%            tetrahedron.
%          w A 1XnumCubPoints set of cubature weights. This sums to the
%            volume of the standard tetrahedron (1/6).
%
%This function implements the points given in [1] (236 points).
%
%EXAMPLE:
%Given the vertices of the simplex, we compare an fourteenth-order moment
%computed using these cubature points to one computed using
%monomialIntSimplex. The results are the same within typical finite
%precision limits.
% [xi,w]=fourteenthOrderTetrahedronCubPoints();
% alpha=[4;6;4];
% theMoment=findMomentFromSamp(alpha,xi,w);
% intVal=monomialIntSimplex(alpha);
% RelErr=(theMoment-intVal)/intVal
%
%REFERENCES:
%[1]L. Zhang and T. Cui, "A set of symmetric quadrature rules on triangles
%   and tetrahedra," Journal of Computational Mathematics, vol. 27, no. 1,
%   pp. 89-96, Jan. 2009.
%
%October 2022 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

w1=0.0040651136652707670436208836835636;
w2=0.0022145385334455781437599569500071;
w3=0.0058134382678884505495373338821455;
w4=0.0196255433858357215975623333961715;
w5=0.0003875737905908214364538721248394;
w6=0.0116429719721770369855213401005552;
w7=0.0052890429882817131317736883052856;
w8=0.0018310854163600559376697823488069;
w9=0.0082496473772146452067449669173660;
w10=0.0030099245347082451376888748208987;
w11=0.0008047165617367534636261808760312;
w12=0.0029850412588493071187655692883922;
w13=0.0056896002418760766963361477811973;
w14=0.0041590865878545715670013980182613;
w15=0.0007282389204572724356136429745654;
w16=0.0054326500769958248216242340651926;

w=[w1*ones(4,1);
   w2*ones(4,1);
   w3*ones(4,1);
   w4*ones(4,1);
   w5*ones(4,1);
   w6*ones(12,1);
   w7*ones(12,1);
   w8*ones(12,1);
   w9*ones(12,1);
   w10*ones(24,1);
   w11*ones(24,1);
   w12*ones(24,1);
   w13*ones(24,1);
   w14*ones(24,1);
   w15*ones(24,1);
   w16*ones(24,1)];

p1=0.3272533625238485639093096692685289;
p2=0.0447613044666850808837942096478842;
p3=0.0861403311024363536537208740298857;
p4=0.2087626425004322968265357083976176;
p5=0.0141049738029209600635879152102928;
p6a=0.1021653241807768123476692526982584;
p6b=0.5739463675943338202814002893460107;
p7a=0.4075700516600107157213295651301783;
p7b=0.0922278701390201300000000000000000;
p8a=0.0156640007402803585557586709578084;
p8b=0.7012810959589440327139967673208426;
p9a=0.2254963562525029053780724154201103;
p9b=0.4769063974420887115860583354107011;
p10a=0.3905984281281458000000000000000000;
p10b=0.2013590544123922168123077327235092;
p10c=0.0161122880710300298578026931548371;
p11a=0.1061350679989021455556139029848079;
p11b=0.0327358186817269284944004077912660;
p11c=0.0035979076537271666907971523385925;
p12a=0.5636383731697743896896816630648502;
p12b=0.2302920722300657454502526874135652;
p12c=0.1907199341743551862712487790637898;
p13a=0.3676255095325860844092206775991167;
p13b=0.2078851380230044950717102125250735;
p13c=0.3312104885193449000000000000000000;
p14a=0.7192323689817295295023401840796991;
p14b=0.1763279118019329762157993033636973;
p14c=0.0207602362571310090754973440611644;
p15a=0.5278249952152987298409240075817276;
p15b=0.4372890892203418165526238760841918;
p15c=0.0092201651856641949463177554949220;
p16a=0.5483674544948190728994910505607746;
p16b=0.3447815506171641228703671870920331;
p16c=0.0867217283322215394629438740085828;

xiBary=zeros(4,236);
xiBary(:,1:4)=genAllMultisetPermutations([p1;p1;p1;1-3*p1]);
xiBary(:,5:8)=genAllMultisetPermutations([p2;p2;p2;1-3*p2]);
xiBary(:,9:12)=genAllMultisetPermutations([p3;p3;p3;1-3*p3]);
xiBary(:,13:16)=genAllMultisetPermutations([p4;p4;p4;1-3*p4]);
xiBary(:,17:20)=genAllMultisetPermutations([p5;p5;p5;1-3*p5]);

xiBary(:,21:32)=genAllMultisetPermutations([p6a;p6a;p6b;1-2*p6a-p6b]);
xiBary(:,33:44)=genAllMultisetPermutations([p7a;p7a;p7b;1-2*p7a-p7b]);
xiBary(:,45:56)=genAllMultisetPermutations([p8a;p8a;p8b;1-2*p8a-p8b]);
xiBary(:,57:68)=genAllMultisetPermutations([p9a;p9a;p9b;1-2*p9a-p9b]);

xiBary(:,69:92)=genAllMultisetPermutations([p10a;p10b;p10c;1-p10a-p10b-p10c]);
xiBary(:,93:116)=genAllMultisetPermutations([p11a;p11b;p11c;1-p11a-p11b-p11c]);
xiBary(:,117:140)=genAllMultisetPermutations([p12a;p12b;p12c;1-p12a-p12b-p12c]);
xiBary(:,141:164)=genAllMultisetPermutations([p13a;p13b;p13c;1-p13a-p13b-p13c]);
xiBary(:,165:188)=genAllMultisetPermutations([p14a;p14b;p14c;1-p14a-p14b-p14c]);
xiBary(:,189:212)=genAllMultisetPermutations([p15a;p15b;p15c;1-p15a-p15b-p15c]);
xiBary(:,213:236)=genAllMultisetPermutations([p16a;p16b;p16c;1-p16a-p16b-p16c]);

%Adjust w for the area of the standard tetrahedron.
w=w/6;

%Convert the barycentric points into normal cubature points for the
%standard tetrahedron given the vertices.
vertices=[1,0,0,0;
          0,1,0,0;
          0,0,1,0];
xi=barycentricCoords2Pt(xiBary,vertices);

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
