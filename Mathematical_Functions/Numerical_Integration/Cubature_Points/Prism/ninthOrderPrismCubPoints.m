function [xi,w]=ninthOrderPrismCubPoints()
%%NINTHORDERPRISMCUBPOINTS Generate ninth-order cubature points for
%  integration over a standard prism in 3D. This is a triangle that has
%  been extruded upward. The vertices are (1,-1,-1), (-1,-1,-1), (-1,1,-1),
%  (1,-1,1), (-1,-1,1), and (-1,1,1).
% 
%INPUTS: None
%
%OUTPUTS: xi This is a 3XnumCubPoints set of points for the standard prism.
%          w A 1XnumCubPoints set of cubature weights. This sums to the
%            volume of the standard prism (4).
%
%This function implements the points given in [1] (60 points).
%
%EXAMPLE:
%We compare a 9th-order moment computed using these cubature points
%to one computed using monomialIntPrism. The results are the same within
%typical finite precision limits.
% [xi,w]=ninthOrderPrismCubPoints();
% alpha=[3;2;4];
% theMoment=findMomentFromSamp(alpha,xi,w);
% intVal=monomialIntPrism(alpha);
% RelErr=(theMoment-intVal)/intVal
%
%REFERENCES:
%[1] F. D. Witherden and P. E. Vincent, "On the identification of symmetric
%    quadrature rules for finite element methods," Computer and Mathematics
%    with Applications, vol. 69, no. 10, pp. 1232-1241, May 2015.
%
%October 2022 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

M=[-0.33333333333333333333333333333333333333,    -0.33333333333333333333333333333333333333,                                            0,     0.19114798844620108508537046584425003736;
  -0.33333333333333333333333333333333333333,    -0.33333333333333333333333333333333333333,    -0.87583915879987761660772171906042408394,      0.1060969195893273614653564550346218398;
  -0.33333333333333333333333333333333333333,    -0.33333333333333333333333333333333333333,     0.87583915879987761660772171906042408394,      0.1060969195893273614653564550346218398;
-0.0010699296586281549894456411507705393926,    -0.99786014068274369002110871769845892121,                                            0,    0.049462405441895223670279761487537762472;
  -0.99786014068274369002110871769845892121,  -0.0010699296586281549894456411507705393926,                                            0,    0.049462405441895223670279761487537762472;
-0.0010699296586281549894456411507705393926,  -0.0010699296586281549894456411507705393926,                                            0,    0.049462405441895223670279761487537762472;
  -0.10986499991499267966311447997302837166,    -0.78027000017001464067377104005394325668,    -0.52963355520607675696300739352733925666,     0.15080793261275860374457825128943921743;
  -0.78027000017001464067377104005394325668,    -0.10986499991499267966311447997302837166,    -0.52963355520607675696300739352733925666,     0.15080793261275860374457825128943921743;
  -0.10986499991499267966311447997302837166,    -0.10986499991499267966311447997302837166,    -0.52963355520607675696300739352733925666,     0.15080793261275860374457825128943921743;
  -0.10986499991499267966311447997302837166,    -0.78027000017001464067377104005394325668,     0.52963355520607675696300739352733925666,     0.15080793261275860374457825128943921743;
  -0.78027000017001464067377104005394325668,    -0.10986499991499267966311447997302837166,     0.52963355520607675696300739352733925666,     0.15080793261275860374457825128943921743;
  -0.10986499991499267966311447997302837166,    -0.10986499991499267966311447997302837166,     0.52963355520607675696300739352733925666,     0.15080793261275860374457825128943921743;
  -0.57769113257897895683607913964943464903,     0.15538226515795791367215827929886929806,    -0.32397688216877005165250824998597314548,     0.09338246974553840864616210022748459774;
   0.15538226515795791367215827929886929806,    -0.57769113257897895683607913964943464903,    -0.32397688216877005165250824998597314548,     0.09338246974553840864616210022748459774;
  -0.57769113257897895683607913964943464903,    -0.57769113257897895683607913964943464903,    -0.32397688216877005165250824998597314548,     0.09338246974553840864616210022748459774;
  -0.57769113257897895683607913964943464903,     0.15538226515795791367215827929886929806,     0.32397688216877005165250824998597314548,     0.09338246974553840864616210022748459774;
   0.15538226515795791367215827929886929806,    -0.57769113257897895683607913964943464903,     0.32397688216877005165250824998597314548,     0.09338246974553840864616210022748459774;
  -0.57769113257897895683607913964943464903,    -0.57769113257897895683607913964943464903,     0.32397688216877005165250824998597314548,     0.09338246974553840864616210022748459774;
  -0.91231503711701435232450671034134809014,     0.82463007423402870464901342068269618029,    -0.40165166857078872104499428459378049449,     0.03838828638520642757038156568015552951;
   0.82463007423402870464901342068269618029,    -0.91231503711701435232450671034134809014,    -0.40165166857078872104499428459378049449,     0.03838828638520642757038156568015552951;
  -0.91231503711701435232450671034134809014,    -0.91231503711701435232450671034134809014,    -0.40165166857078872104499428459378049449,     0.03838828638520642757038156568015552951;
  -0.91231503711701435232450671034134809014,     0.82463007423402870464901342068269618029,     0.40165166857078872104499428459378049449,     0.03838828638520642757038156568015552951;
   0.82463007423402870464901342068269618029,    -0.91231503711701435232450671034134809014,     0.40165166857078872104499428459378049449,     0.03838828638520642757038156568015552951;
  -0.91231503711701435232450671034134809014,    -0.91231503711701435232450671034134809014,     0.40165166857078872104499428459378049449,     0.03838828638520642757038156568015552951;
  -0.03450013665826621486161695118676603904,    -0.93099972668346757027676609762646792192,     -0.9677921254642105035410819146400981344,    0.026528431202320433845464036710085263608;
  -0.93099972668346757027676609762646792192,    -0.03450013665826621486161695118676603904,     -0.9677921254642105035410819146400981344,    0.026528431202320433845464036710085263608;
  -0.03450013665826621486161695118676603904,    -0.03450013665826621486161695118676603904,     -0.9677921254642105035410819146400981344,    0.026528431202320433845464036710085263608;
  -0.03450013665826621486161695118676603904,    -0.93099972668346757027676609762646792192,      0.9677921254642105035410819146400981344,    0.026528431202320433845464036710085263608;
  -0.93099972668346757027676609762646792192,    -0.03450013665826621486161695118676603904,      0.9677921254642105035410819146400981344,    0.026528431202320433845464036710085263608;
  -0.03450013665826621486161695118676603904,    -0.03450013665826621486161695118676603904,      0.9677921254642105035410819146400981344,    0.026528431202320433845464036710085263608;
   -0.6456435117766454460710959796587760834,     0.29128702355329089214219195931755216681,    -0.89078174045697066578383528302064143989,    0.066719942573448407007300296577585702464;
   0.29128702355329089214219195931755216681,     -0.6456435117766454460710959796587760834,    -0.89078174045697066578383528302064143989,    0.066719942573448407007300296577585702464;
   -0.6456435117766454460710959796587760834,     -0.6456435117766454460710959796587760834,    -0.89078174045697066578383528302064143989,    0.066719942573448407007300296577585702464;
   -0.6456435117766454460710959796587760834,     0.29128702355329089214219195931755216681,     0.89078174045697066578383528302064143989,    0.066719942573448407007300296577585702464;
   0.29128702355329089214219195931755216681,     -0.6456435117766454460710959796587760834,     0.89078174045697066578383528302064143989,    0.066719942573448407007300296577585702464;
   -0.6456435117766454460710959796587760834,     -0.6456435117766454460710959796587760834,     0.89078174045697066578383528302064143989,    0.066719942573448407007300296577585702464;
  -0.90753763917811234690103494985844043039,     0.81507527835622469380206989971688086078,     -0.9577062038929269606887790227035291939,    0.012131173836407148944270419030858660558;
   0.81507527835622469380206989971688086078,    -0.90753763917811234690103494985844043039,     -0.9577062038929269606887790227035291939,    0.012131173836407148944270419030858660558;
  -0.90753763917811234690103494985844043039,    -0.90753763917811234690103494985844043039,     -0.9577062038929269606887790227035291939,    0.012131173836407148944270419030858660558;
  -0.90753763917811234690103494985844043039,     0.81507527835622469380206989971688086078,      0.9577062038929269606887790227035291939,    0.012131173836407148944270419030858660558;
   0.81507527835622469380206989971688086078,    -0.90753763917811234690103494985844043039,      0.9577062038929269606887790227035291939,    0.012131173836407148944270419030858660558;
  -0.90753763917811234690103494985844043039,    -0.90753763917811234690103494985844043039,      0.9577062038929269606887790227035291939,    0.012131173836407148944270419030858660558;
  -0.88126598051326876444819587362986475675,     0.44939053910599895722040631866985064464,                                            0,    0.087691161312839084087580790508946261895;
   0.44939053910599895722040631866985064464,    -0.88126598051326876444819587362986475675,                                            0,    0.087691161312839084087580790508946261895;
  -0.56812455859273019277221044503998588789,     0.44939053910599895722040631866985064464,                                            0,    0.087691161312839084087580790508946261895;
   0.44939053910599895722040631866985064464,    -0.56812455859273019277221044503998588789,                                            0,    0.087691161312839084087580790508946261895;
  -0.56812455859273019277221044503998588789,    -0.88126598051326876444819587362986475675,                                            0,    0.087691161312839084087580790508946261895;
  -0.88126598051326876444819587362986475675,    -0.56812455859273019277221044503998588789,                                            0,    0.087691161312839084087580790508946261895;
   -0.5517449084932878730050536869300298641,     0.49681091913207739922639584227322548803,    -0.69522055387163478698981196744187399375,      0.0495312141698622864915543816230467997;
   0.49681091913207739922639584227322548803,     -0.5517449084932878730050536869300298641,    -0.69522055387163478698981196744187399375,      0.0495312141698622864915543816230467997;
  -0.94506601063878952622134215534319562393,     0.49681091913207739922639584227322548803,    -0.69522055387163478698981196744187399375,      0.0495312141698622864915543816230467997;
   0.49681091913207739922639584227322548803,    -0.94506601063878952622134215534319562393,    -0.69522055387163478698981196744187399375,      0.0495312141698622864915543816230467997;
  -0.94506601063878952622134215534319562393,     -0.5517449084932878730050536869300298641,    -0.69522055387163478698981196744187399375,      0.0495312141698622864915543816230467997;
   -0.5517449084932878730050536869300298641,    -0.94506601063878952622134215534319562393,    -0.69522055387163478698981196744187399375,      0.0495312141698622864915543816230467997;
   -0.5517449084932878730050536869300298641,     0.49681091913207739922639584227322548803,     0.69522055387163478698981196744187399375,      0.0495312141698622864915543816230467997;
   0.49681091913207739922639584227322548803,     -0.5517449084932878730050536869300298641,     0.69522055387163478698981196744187399375,      0.0495312141698622864915543816230467997;
  -0.94506601063878952622134215534319562393,     0.49681091913207739922639584227322548803,     0.69522055387163478698981196744187399375,      0.0495312141698622864915543816230467997;
   0.49681091913207739922639584227322548803,    -0.94506601063878952622134215534319562393,     0.69522055387163478698981196744187399375,      0.0495312141698622864915543816230467997;
  -0.94506601063878952622134215534319562393,     -0.5517449084932878730050536869300298641,     0.69522055387163478698981196744187399375,      0.0495312141698622864915543816230467997;
   -0.5517449084932878730050536869300298641,    -0.94506601063878952622134215534319562393,     0.69522055387163478698981196744187399375,      0.0495312141698622864915543816230467997];
w=M(:,4);
xi=M(:,1:3)';

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