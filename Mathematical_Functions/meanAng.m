function angAvg=meanAng(ang,w,dim)
%%MEANANG Find the weighted or unweighted mean of a set of angles in
%         radians. This function is needed, because one cannot simply
%         average angles, because issues can arise at the boundary of
%         0/(2*pi) or -pi/pi or any other areas that are multiples of
%         2*pi apart.
%
%INPUTS: ang A vector, matrix or hypermatrix of angles, whereby the
%            average over a particular dimension (rows=1, columns=2), or
%            all elements is desired.
%          w An optional vector of weights having the same dimensionality
%            as ang, if a weighted average is desired. The weights do not
%            have to sum to one. If omitted or an empty matrix is passed,
%            all angles are treated as having the same weight.
%        dim An optional parameter specifying the dimension over which
%            averaging is performed. This can be 0 (average averything) 1,
%            average across rows, 2, average across columns, or higher
%            numbers to average across higher indices of hypermatrices.
%
%OUTPUTS: angAvg The average in radians of the angles in ang, taken
%                across dimensions dim (0 for all), weighted by
%                w, if w was provided. angAvg is from -pi to pi. If dim=0,
%                then angAvg is scalar. Otherwise, the dimensionality of 
%                angAvg is the same as ang, except dimension dim has been
%                collapsed to one.
%
%The method for performing a weighted average of points comes from the
%definition of population mean angle used in Chapter 9.2.1 of [1].
%Basically, the angles are converted into 2D unit vectors, the weighted
%average of the unit vectors is taken, and an inverse tangent is used to
%transform the 2D unit vectors back into an angle. In the event that the
%unit vectors all cancel out, for example averaging -pi/2 and pi/2, the
%result can be anywhere on the unit circle (due to numerical issues not
%completely cancelling the vectors). If the vectors cancel perfectly, then
%atan2(0,0) is the result, which currently returns 0.
%
%REFERENCES:
%[1] K. V. Mardia and P. E. Jupp, Directional Statistics. Chichester,
%    England: John Wiley & Sons, 2000.
%
%September 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3)
    dim=0;
end

if(nargin<2)
   w=[]; 
end

u=cos(ang);
v=sin(ang);

%Weight the values.
if(~isempty(w))
   u=u(:).*w(:);
   v=v(:).*w(:);
end

if(dim==0)
    u=sum(u(:));
    v=sum(v(:));
else
    u=sum(u,dim);
    v=sum(v,dim);
end

angAvg=atan2(v,u);
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
