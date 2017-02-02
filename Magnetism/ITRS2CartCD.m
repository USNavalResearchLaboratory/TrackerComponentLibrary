function [zCartCD,rotMat]=ITRS2CartCD(zCart,param1,param2,S11)
%%ITRS2CARTCD  Convert points given in the International Terrestrial
%              Reference System (ITRS), a standard Earth-Cenetred-Earth-
%              fixed (ECEF) coordinate system, in 3D Cartesian coordinates
%              into centered dipole (CD) coordinates, a type of cooridnate
%              system where the z-axis is aligned with the magnetic dipole
%              of the Earth. This coordinate system is useful when studying
%              the ionosphere. Note that the inputs and outputs of this
%              function are in Cartesian coordinates, but many papers only
%              use CD coordinates in spherical form. The function
%              Cart2Sphere could be used to perform the conversion to
%              spherical form.
%
%INPUTS: zCart  One or more points given in Cartesian coordinates in the
%               ITRS. zCart is a 3XN matrix with each column
%               having the format [x;y;z].
%param1,param2,S11 These optional parameters specify the magnetic field
%               parameters from which the direction of the CD pole is
%               derived. These aprameters are impied set base don how the
%               function is called:
%               1) ITRS2CartCD(zCart)
%               If all of the other parameters are omitted, then the CD
%               pole parameters for the latest epoch of the International
%               Geomagnetic Reference Field (IGRF) are used; the Schmidt
%               semi-normalized coefficients are obtained using the
%               function getIGRFCoeffs.
%               2) ITRS2CartCD(zCart, year)
%               The parameters of the IGRF for the specified epoch year are
%               used. The year is in the Gregorian calendar and is
%               specified as noted in the comments to the getIGRFCoeffs
%               function.
%               3) ITRS2CartCD(zCart, year, model). This is the
%               same as spher2CenteredDipole(zSpher, year), except model
%               can be 'IGRF' or 'WMM' to specify which magnetic field
%               model to use. 'WMM' refers to the World magnetic model
%               using the function getWMMCoeffs. Note than an empty year
%               matrix can be passed to get the latest epoch year of either
%               model.
%               4) ITRS2CartCD(zCart,C10,C11,S11). In this
%               instance, the first few Schmidt semi-normalized
%               coefficients for a magnetic field model are given. These
%               are the only three coefficients used to define the
%               orientation of the CD pole. If using the getWMMCoeffs or
%               getIGRFCoeffs function to get the coordinates, then from
%               the outputs of those functions (specifying in the second
%               input that they should not be fully normalized), the
%               outputs C and S of the functions are used here as
%               C10=C(1+1,0+1), C11=C(1+1,1+1), S11=C(1+1,1+1).
%
%OUTPUTS: zCartCD The points in zCart rotated into centered dipole
%               coordinates under the selected model. The ITRS and the CD
%               coordinate systems are related by a rotation.
%        rotMat The 3X3 rotation matrix such that
%               zCartCD(:,i)=rotMat*zCart(:,i).
%
%The definition of centered dipole (CD) coordinates is taken from [1],
%though the relations given in spherical coordinates are not used.
%
%REFERENCES:
%[1] A. C. Fraser-Smith, "Centered and eccentric geomagnetic dipoles and
%    their poles 1600-1985," Reviews of Geophysics, vol. 25, no. 1, pp.
%    1-16, Feb. 1987.
%
%June 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%%Determine the centered dipole location from a magnetic field model.
if(nargin==1||nargin==2||nargin==3)
    if(nargin==1||isempty(param1))
    %If only position components were given, then just use the values for
    %the reference epoch of the IGRF. The unnormalized coefficients are
    %desired.
        [C,S]=getIGRFCoeffs([],false);
    elseif(nargin==2)
    %If a reference year is given, then use that in the IGRF. 
        [C,S]=getIGRFCoeffs(param1,false);
    elseif(nargin==3)%If one wishes to select the magnetic field model.
        switch(param2)
            case 'IGRF'
                if(isempty(param1))
                    [C,S]=getIGRFCoeffs([],false);
                else%If the year is given
                    [C,S]=getIGRFCoeffs(param1,false);
                end
            case 'WMM'
                if(isempty(param1))
                    [C,S]=getWMMCoeffs([],false);
                else%If the year is given
                    [C,S]=getWMMCoeffs(param1,false);
                end
            otherwise
                error('Unknown magnetic field model selected')
        end
    end
 
    g10=C(1+1,1+0);
    g11=C(1+1,1+1);
    h11=S(1+1,1+1);
elseif(nargin==4)
    %If the required coefficients are directly provided.
    g10=param1;
    g11=param2;
    h11=S11;
end

%The paper does not explicitely state it, but from Equation 5 of [1]
%(assuming one uses a four-quadrant inverse tangent), one can infer that
%the pole being used is given by the following vector in terms of
%geomagnetic components:
uCDPole=[-g11;-h11;-g10];
B0=norm(uCDPole);%Reference magnetic field, the "reduced moment" (Tesla)
uCDPole=uCDPole/B0;%Turn the pole direction into a unit vector.

%The sphericla geometric relations in the paper contain a number of
%singularities. Thus, it is best if the rotations take place in Cartesian
%coordinates.
zSpherPol=Cart2Sphere(uCDPole);
phiN=zSpherPol(2);%Azimuth of the CD pole
thetaN=pi/2-zSpherPol(3);%Colatitude of CD pole 

%Rotation matrix from geocentric coordinates to centered-dipole
%coordinates. This rotates uCDPole to [0;0;1] --The z axis. The first
%rotation about the z axis makes the point where the pole is have zero
%longitude in CD coordinates.
rotMat=Euler2Ang2RotMat(-thetaN,-phiN,'yz');

%Rotate the points. Doing this in Cartesian space avoids the singularities
%in Equation 3 of [1].
zCartCD=rotMat*zCart;
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
