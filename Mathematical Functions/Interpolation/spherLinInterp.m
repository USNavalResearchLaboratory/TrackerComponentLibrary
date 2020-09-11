function pInterp=spherLinInterp(u,p1,p2,normalize)
%%SPHERLININTERP Perform spherical linear interpolation (SLERP) between two
%                vectors at a number of points. SLERP is useful for
%                obtaining a smooth interpolation between two 3D directions
%                as unit vectors or between two orientations, represented
%                in 4D using  quaternions. The interpolation is performed
%                the short way on the sphere between the vectors. However,
%                note that with quaternions q and -q represent the same
%                orientation, thus if the last elements of the quaternions
%                have opposite signs, then the rotation will be the long
%                way. The interpolation routine can also be used with non-
%                unit vectors.
%
%INPUTS: u  A 1XN or NX1 vector of the n points at which interpolation
%           should be performed. These are values from 0 to 1 that
%           interpolate along the way from p1 to p2.
%     p1,p2 The pDimX1 starting and ending interpolation points. This are
%           typically 3D for directions and 4D for interpolating
%           orientations via quaternions. They do not have to be unit
%           vectors. pDim>=2.
% normalize An optional parameter indicating whether the output pInterp
%           should be normalized. If p1 and p2 have the same magnitude,
%           then barring finite precision problems, normalization is not
%           necessary, because the output will have the same magnitude as
%           the inputs. The default if this parameter is omitted or an
%           empty matrix is passed is false.
%
%OUTPUTS: pInterp A pDimXn matrix of the vectors interpolated at the N
%                 fractions given in u.
%
%The interpolation path is a constant-velocity torque-minimal path, making
%it an interpolation method consistent with physics when considering
%rotations. The SLERP algorithm originates in [1] in the "great arc
%in-betweening" section.
%
%Note that spherLinInterp2Pt(u,p1,p2) equals spherLinInterp2Pt(1-u,p2,p1),
%as one would expect.
%
%REFERENCES:
%[1] K. Shoemake, "Animating rotations with quaternion curves," in
%    Proceedings of the 12th Annual Conference on Computer graphics and
%    Interactive Techniques (SIGGRAPH '85), vol. 19, no. 3, San Francisco,
%    CA, 22-26 1985, pp. 245-254.
%
%August 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4||isempty(normalize))
    normalize=false;
end

if(length(p1)==3)
    theta=angBetweenVecs(p1,p2);
else%Assume that we are dealing with
    cosTheta=dot(p1,p2)/(norm(p1)*norm(p2));
    %Deal with finite precision issues pushing the magnitude over 1.
    if(abs(cosTheta)>1)
        cosTheta=sign(cosTheta);
    end
    
    theta=acos(cosTheta);
    if(~isfinite(theta))%Deal with small magnitudes.
        theta=0;
    end
end

sinTheta=sin(theta);
u=u(:)';%Make it a row vector.
if(sinTheta==0)
    %Deal with parallel vectors by inserting the limit as theta->0. From
    %l'Hôpital's rule, this is just linear interpolation.
    pInterp=bsxfun(@times,(1-u),p1)+bsxfun(@times,u,p2);
else
    pInterp=bsxfun(@times,(sin((1-u)*theta)/sinTheta),p1)+bsxfun(@times,(sin(u*theta)/sinTheta),p2);
end

if(normalize)
    mags=sqrt(sum(pInterp.*pInterp,1));
    pInterp=bsxfun(@rdivide,pInterp,mags);
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
