function retVec=ICRS2J2000F(vec,TT1,TT2)
%%ICRS2J2000F Convert vectors from the International Celestial Reference
%             System (ICRS) to the J2000.0 dynamical frame. This involves
%             removing a frame bias rotation.
%
%INPUTS: vec The 3XN matrix of N vectors that are to be rotated rom the
%            ICRS to the J2000.0 dynamical frame.
%   TT1, TT2 Two parts of a Julian date given in terrestrial time (TT).
%            The units of the date are days. The full date is the sum of
%            both terms. The date is broken into two parts to provide more
%            bits of precision. It does not matter how the date is split.
%
%OUTPUTS: retVec The 3XN set of vectors rotated into the J2000.0 dynamical
%                frame.
%
%As described in [1], the ICRS is ofset from the J2000 dynamical frame by
%a bias rotation.
%
%This is a wrapper for the function iauBp06 and some matrix operations in
%the International Astronomical Union's (IAU) Standard's of Fundamental
%Astronomy library.
%
%The algorithm can be compiled for use in Matlab  using the
%CompileCLibraries function.
%
%The algorithm is run in Matlab using the command format
%retVec=ICRS2J2000F(vec,TT1,TT2);
%
%REFERENCES:
%[1] G. Petit and B. Luzum, IERS Conventions (2010), International Earth
%    Rotation and Reference Systems Service Std. 36, 2010.
%
%March 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

error('This function is only implemented as a mexed C or C++ function. Please run CompileCLibraries.m to compile the function for use.')

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
