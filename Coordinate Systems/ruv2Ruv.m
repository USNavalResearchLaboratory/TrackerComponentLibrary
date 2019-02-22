function zNew=ruv2Ruv(z,useHalfRange,zTx1,zRx1,M1,zTx2,zRx2,M2,includeW)
%%RUV2RUV Convert from bistatic r-u-v (or r-u-v-w) coordinates with respect
%         to one bistatic radar pair into those for another. This function
%         can be useful for converting bistatically received measurements
%         into the coordinate system of the transmitter so that they can be
%         gated with the transmit beam of a radar. If r-u-v coordinates are
%         used, the target is assumed to be in front of the receiver (local
%         z coordinate is positive). r-u-v coordinates consist of a
%         bistatic range and direction cosines at the receiver. The
%         "direction cosines" u and v are just the x and y coordinates of a
%         unit vector from the receiver to the target in the coordinate
%         system at the receiver. This basically assumes that the boresight
%         direction of the receiver is the z axis. Assuming the target is
%         in front of the receiver, the third unit vector coordinate is not
%         needed. However, r-u-v-w coordinates include the third component.
%
%INPUTS: z A 3XN matrix of vectors with elements [r;u;v], where r is the
%          bistatic range from the transmitter to the target to the
%          receiver, and u and v are direction cosines. Each u,v pair
%          should have a magnitude less than or equal to one. If the
%          magnitude is greater than one, then the pair is normalized
%          before conversion to avoid imaginary results. Alternatively,
%          one can pass a 4XN matrix of [r;u;v;w] vectors where [u;v;w]
%          form a full unit vector in the receiver's local 3D Cartesian
%          coordinates.
% useHalfRange A scalar boolean value or a 2X1 or 1X2 vector of boolean
%          values specifying whether the bistatic range value should be
%          divided by two. This normally comes up when operating in
%          monostatic mode, so that the range reported is a one-way range.
%          The default if an empty matrix is provided is false.
%          useHalfRange(1) applies to the first (source) bistatic channel
%          and useHalfRange(2) to the second (destination) bistatic
%          channel. If a scalar is passed, then both values are taken to
%          be the same.
%     zTx1 The 3XN [x;y;z] location vectors of the transmitters
%          originating the measurements z in global Cartesian coordinates.
%          If only a single vector is passed, then the transmitter
%          location is assumed the same for all of the target states being
%          converted. zTx1 can have more than 3 rows; additional rows are
%          ignored.
%     zRx1 The 3XN [x;y;z] location vectors of the receivers in Cartesian
%          coordinates that produced the measurements z . If only a single
%          vector is passed, then the receiver location is assumed the
%          same for all of the target states being converted. zRx1 can
%          have more than 3 rows; additional rows are ignored.
%       M1 A 3X3XN hypermatrix of the rotation matrices to go from the
%          alignment of the global coordinate system to that at the
%          receiver that produced the measurements z. The z-axis of the
%          local coordinate system of the receiver is the pointing
%          direction of the receiver. If omitted, then it is assumed that
%          the local coordinate system is aligned with the global and
%          M1=eye(3) --the identity matrix is used. If only a single 3X3
%          matrix is passed, then it is assumed to be the same for all of
%          the N conversions.
% zTx2,zRx2,M2 These are the same as aTx1,zRx1, and M1, but for the
%          coordinate system into which the bistatic r-u-v(-w)
%          measurements should be converted.
% includeW An optional boolean value indicating whether a third direction
%          cosine component should be included in the output (regardless of
%          whether it was provided in the input). The u and v direction
%          cosines are two parts of a 3D unit vector. Generally, one might
%          assume that the target is in front of the sensor, so the third
%          component would be positive and is not needed. However, the
%          third component can be included if ambiguity exists. The
%          default if this parameter is omitted or an empty matrix is
%          passed is false if the input z is 3XN and true if z is 4XN.
%
%OUTPUTS: zNew The measurements converted into the other coordinate system.
%              If z was 3XN, then this is 3XN. If z was 4XN, then this is
%              4XN, reflecting the use of r-u-v-w coordinates in place of
%              just r-u-v coordinates.
%
%The function just calls the ruv2Cart and Cart2Ruv functions with the
%appropriate parameters.
%
%October 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numDim=size(z,1);

if(isscalar(useHalfRange))
    useHalfRange=[useHalfRange;useHalfRange];
end

if(nargin<9||isempty(includeW))
    %4D means r-u-v-w coordinates, so assume that the w coordinate should
    %be included.
    includeW=(numDim==4);
end

zC=ruv2Cart(z,useHalfRange(1),zTx1,zRx1,M1);
zNew=Cart2Ruv(zC,useHalfRange(2),zTx2,zRx2,M2,includeW);

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
