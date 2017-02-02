function zSpher=spherCD2SpherITRS(zCD,param1,param2,S11)
%%SPHERCD2SPHERITRS Convert a point given in spherical cooridnates in
%          centered dipole (CD) coordinates into spherical coordinates in
%          the International Terrestrial Reference System (ITRS), a
%          standard Earth-Cenetred-Earth-fixed (ECEF) coordinate system. CD
%          coordinates is a type of coordinate system where the z-axis is
%          aligned with the the magnetic dipole of the Earth (under a
%          dipole model). This coordinate system is useful when studying
%          the ionosphere.
%
%INPUTS: zCD One or more points given in terms of range, azimuth and
%            elevation, with the angles in radians, or in terms of just
%            azimuth and elevation if a conversion of directions is
%            desired, in CD coordinates. To convert N points, zCD is a 3XN
%            matrix with each column having the format [range;azimuth;
%            elevation] or it is a 2XN matrix with each column having
%            format [azimuth; elevation]. Azimuth is measured
%            counterclockwise from the x-axis in the x-y plane. Elevation
%            is measured up from the x-y plane (towards the z-axis).
%param1,param2,S11 These optional parameters specify the magnetic field
%               parameters from which the direction of the CD pole is
%               derived. These aprameters are impied set base don how the
%               function is called:
%               1) spherCD2SpherITRS(zCD)
%               If all of the other parameters are omitted, then the CD
%               pole parameters for the latest epoch of the International
%               Geomagnetic Reference Field (IGRF) are used; the Schmidt
%               semi-normalized coefficients are obtained using the
%               function getIGRFCoeffs.
%               2) spherCD2SpherITRS(zCD, year)
%               The parameters of the IGRF for the specified epoch year are
%               used. The year is in the Gregorian calendar and is
%               specified as noted in the comments to the getIGRFCoeffs
%               function.
%               3) spherCD2SpherITRS(zCD, year, model). This is the
%               same as centeredDipole2Sphere(zCD, year), except model
%               can be 'IGRF' or 'WMM' to specify which magnetic field
%               model to use. 'WMM' refers to the World magnetic model
%               using the function getWMMCoeffs. Note than an empty year
%               matrix can be passed to get the latest epoch year of either
%               model.
%               4) spherCD2SpherITRS(zCD,C10,C11,S11). In this
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
%OUTPUTS: zSphere The points in zCD converted into spherical coordinates 
%                 under the selected model (systemType 0 in the Cart2Sphere
%                 function). The z axis is algined with the If zCD is 3XN,
%                 then zSphere consists of range as well as azimuth and
%                 elevation. If zCD was 2XN, then zSphere just consists of
%                 azimuth and elevation. All angles are in radians.
%
%The definition of centered dipole (CD) coordinates is taken from [1].
%However, as opposed to using the spherical trignonometric equations in the
%paper, the rotation is performed in Cartesian coordinates so as to avoid
%singularities.
%
%REFERENCES:
%[1] A. C. Fraser-Smith, "Centered and eccentric geomagnetic dipoles and
%    their poles 1600-1985," Reviews of Geophysics, vol. 25, no. 1, pp.
%    1-16, Feb. 1987.
%
%June 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numDim=size(zCD,1);

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

zITRS=CartCD2ITRS(spher2Cart(zCD),g10,g11,h11);

%Convert the points back to spherical coordinates to return.
zSpher=Cart2Sphere(zITRS);

if(numDim<3)
%If range was not provided, discard the unit ranges that arose using the
%Spher2Cart conversion.
    zSpher=zSpher(2:3,:);
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
