function [xi,w]=fifthOrderPyramidCubPoints()
%%FIFTHORDERPYRAMIDCUBPOINTS Generate fifth-order cubature points for
%  integration over a 3-dimensional pyramid with a square base with the
%  peak vertex at (0,0,1) and the base vertices at (1,-1,-1), (-1,-1,-1),
%  (-1,1,-1), and (1,1,-1).
% 
%INPUTS: None
%
%OUTPUTS: xi This is a 3XnumCubPoints set of points for the standard
%            square pyramid.
%          w A 1XnumCubPoints set of cubature weights. This sums to the
%            volume of the standard square pyramid (8/3).
%
%This function implements the points given in [1] (15 points).
%
%EXAMPLE:
%We compare a 5th-order moment computed using these cubature points
%to one computed using monomialIntPyramid. The results are the same within
%typical finite precision limits.
% [xi,w]=fifthOrderPyramidCubPoints();
% alpha=[2;2;1];
% theMoment=findMomentFromSamp(alpha,xi,w);
% intVal=monomialIntPyramid(alpha);
% RelErr=(theMoment-intVal)/intVal
%
%REFERENCES:
%[1] F. D. Witherden and P. E. Vincent, "On the identification of symmetric
%    quadrature rules for finite element methods," Computer and Mathematics
%    with Applications, vol. 69, no. 10, pp. 1232-1241, May 2015.
%
%October 2022 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

M=[                                     0,                                          0,   0.45971576156501338586164265377920811314,   0.18249431975770692138374895897213800931;
                                        0,                                          0,  -0.39919795837246198593385139590712322914,   0.45172563864726406056400285032640105704;
                                        0,                                          0,  -0.99999998701645569241460017355590234925,   0.15654542887619877154120304977336547704;
 0.70652603154632457420722562974792066862,                                          0,                                      -0.75,   0.20384344839498724639142514342645843799;
                                        0,   0.70652603154632457420722562974792066862,                                      -0.75,   0.20384344839498724639142514342645843799;
-0.70652603154632457420722562974792066862,                                          0,                                      -0.75,   0.20384344839498724639142514342645843799;
                                        0,  -0.70652603154632457420722562974792066862,                                      -0.75,   0.20384344839498724639142514342645843799;
 0.70511712277882760181079385797948261057,   0.70511712277882760181079385797948261057,  -0.87777618587595407108464357252416911085,   0.10578907087905457654220386143818487109;
 0.70511712277882760181079385797948261057,  -0.70511712277882760181079385797948261057,  -0.87777618587595407108464357252416911085,   0.10578907087905457654220386143818487109;
-0.70511712277882760181079385797948261057,   0.70511712277882760181079385797948261057,  -0.87777618587595407108464357252416911085,   0.10578907087905457654220386143818487109;
-0.70511712277882760181079385797948261057,  -0.70511712277882760181079385797948261057,  -0.87777618587595407108464357252416911085,   0.10578907087905457654220386143818487109;
 0.43288286410354097685000790909815143591,   0.43288286410354097685000790909815143591,  -0.15279732576055038842025517341026975071,   0.15934280057233240536079894703404722173;
 0.43288286410354097685000790909815143591,  -0.43288286410354097685000790909815143591,  -0.15279732576055038842025517341026975071,   0.15934280057233240536079894703404722173;
-0.43288286410354097685000790909815143591,   0.43288286410354097685000790909815143591,  -0.15279732576055038842025517341026975071,   0.15934280057233240536079894703404722173;
-0.43288286410354097685000790909815143591,  -0.43288286410354097685000790909815143591,  -0.15279732576055038842025517341026975071,   0.15934280057233240536079894703404722173];

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
