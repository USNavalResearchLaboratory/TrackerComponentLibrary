function [xi,w]=fourthOrderPrismCubPoints()
%%FOURTHORDERPRISMCUBPOINTS Generate fourth-order cubature points for
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
%This function implements the points given in [1] (11 points).
%
%EXAMPLE:
%We compare a 4th-order moment computed using these cubature points
%to one computed using monomialIntPrism. The results are the same within
%typical finite precision limits.
% [xi,w]=fourthOrderPrismCubPoints();
% alpha=[2;2;0];
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

M=[-0.33333333333333333333333333333333333333,   -0.33333333333333333333333333333333333333,   -0.86686197400903479841583223937151297546,    0.43164789926214201915560571539043719684;
 -0.33333333333333333333333333333333333333,   -0.33333333333333333333333333333333333333,    0.86686197400903479841583223937151297546,    0.43164789926214201915560571539043719684;
-0.062688380276009575455822705856210492672,   -0.87462323944798084908835458828757901466,                                           0,    0.54565845042191057226153838655080393041;
 -0.87462323944798084908835458828757901466,  -0.062688380276009575455822705856210492672,                                           0,    0.54565845042191057226153838655080393041;
-0.062688380276009575455822705856210492672,  -0.062688380276009575455822705856210492672,                                           0,    0.54565845042191057226153838655080393041;
 -0.79851918840217872745314255427215462609,    0.59703837680435745490628510854430925218,   -0.67563982368225971739661759355428029031,    0.24995480836833070748402890159445230251;
  0.59703837680435745490628510854430925218,   -0.79851918840217872745314255427215462609,   -0.67563982368225971739661759355428029031,    0.24995480836833070748402890159445230251;
 -0.79851918840217872745314255427215462609,   -0.79851918840217872745314255427215462609,   -0.67563982368225971739661759355428029031,    0.24995480836833070748402890159445230251;
 -0.79851918840217872745314255427215462609,    0.59703837680435745490628510854430925218,    0.67563982368225971739661759355428029031,    0.24995480836833070748402890159445230251;
  0.59703837680435745490628510854430925218,   -0.79851918840217872745314255427215462609,    0.67563982368225971739661759355428029031,    0.24995480836833070748402890159445230251;
 -0.79851918840217872745314255427215462609,   -0.79851918840217872745314255427215462609,    0.67563982368225971739661759355428029031,    0.24995480836833070748402890159445230251];

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
