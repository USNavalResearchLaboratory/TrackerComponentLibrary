function NsVals=reduceStdRefrac2Spher(NMeas,height,expConst,multConst,varargin)
%%REDUCESTDREFRAC2SPHER Given a measurement of the atmospheric refractivity
%                       at a particular height above sea level, determine
%                       the equivalent refractivity at sea level using  a
%                       standard exponential atmospheric model. The
%                       algorithm can handle sea surface refractivities
%                       from 200 to 450. 
%
%INPUTS: NMeas The measured atmospheric refractivity. Note that the
%              refractivity is (n-1)*1e6, where n is the index of
%              refraction.
%       height The height in meters at which the measurement of NMeas was
%              taken. One would expect this to be an orthometic height
%              (with respect to mean sea level). However, the model is low
%              fidelity, so using a height with respect to the Earth's
%              reference ellipsoid would presumably not introduce a
%              significant amount of error.
%  expConst, multConst The two optional, positive  parameters
%              parameterizing the decay constant in the model. Assume that
%              at a point, the refractivity is N. Increasing the height by
%              1km, the model is that the refractivity changes by
%              deltaN=-multConst*exp(expConst*N). deltaN cannot be
%              negative. If these parameters are omitted or empty matrices
%              are passed, then the values fitting the data in [1] are
%              used: expConst=0.005577 and multConst=7.32. Note that the
%              value ce=log(Ns/(Ns+DeltaN))/1000 is the decay constant in
%              inverse meters.
%     varargin Parameters to pass to the fminbnd function. These are
%              standard comma-separated keys and values, such as
%              'TolX',1e-8.
%
%OUTPUTS: Ns The indices of refraction reduced to sea level under a
%            standard exponential atmospheric model. There can be 1-2
%            solutions.
%
%The standard exponential atmospheric model is derived in [1]. This
%function determines the refractivity at the given altitudes for sea
%surface refractivities from 200 to 450 and subtracts the observed
%refractivity. Any zero crossings indicate potential solutions. The
%function fminbnd is then used to find each of the zeros given bounded
%regions for the zero crossings.
%
%The use of this type of very basic refraction model is discussed in [2].
%
%REFERENCES:
%[1] B. R. Bean and G. D. Thayer, CRPL Exponential Reference Atmosphere.
%     Washington, D.C.: U. S. Department of Commerce, National Bureau of
%     Standards, Oct. 1959. [Online]. Available:
%     http://digicoll.manoa.hawaii.edu/techreports/PDF/NBS4.pdf
%[2] D. F. Crouse, "Basic tracking using 3D monostatic and bistatic
%    measurements in refractive environments," IEEE Aerospace and
%    Electronic Systems Magazine, vol. 29, no. 8, Part II, pp. 54-75, Aug.
%    2014.
%
%June 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(expConst))
    expConst=0.005577; 
end

if(nargin<4||isempty(multConst))
	multConst=7.32; 
end

%Evaluate the refractivity on a grid.
numPoints=100;
Ns=linspace(200,450,numPoints);
theVals=Ns.*(Ns./(Ns-multConst*exp(expConst*Ns))).^(-height/1000)-NMeas;

%Find where the zero crossings are
diffIdx=find(diff(theVals>0));

%Each crossing gives a region to search to find the zero. Because it is
%bounded, we will use fminbnd to find the minimum of the squared value in
%the bounded region.
numNs=length(diffIdx);
NsVals=zeros(numNs,1);

costFun=@(Ns)(Ns.*(Ns./(Ns-multConst*exp(expConst*Ns))).^(-height/1000)-NMeas)^2;
for curNs=1:numNs
    NsCur=fminbnd(costFun,Ns(diffIdx(curNs)),Ns(diffIdx(curNs)+1),varargin{:});
    NsVals(curNs)=NsCur;
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
