function vol=hypersphereVolume(r,numDim)
%%HYPERSPHEREVOLUME Determine the volume of a sphere of a given radius in
%                   multiple dimensions.
%                   
%INPUTS: r The radius of the hypersphere. r>0.
%   numDim The number of dimensions of the hypersphere. numDim>=1. One
%          dimension is a line, two a circle, three a sphere, etc.
%
%OUTPUTS: vol The volume of the hypersphere.
%
%The volume of a hypersphere is given in terms of the gamma function in
%[1].
%
%EXAMPLE:
%We can verify the correctness of this function by comparing the volume
%obtained by this function to that obtained using A Monte Carlo simulation.
% r=5;
% numDim=4;
% numSamples=1e5;
% 
% numInSphere=0;
% for curSamp=1:numSamples
%     %Generate a point unifomrly in the hypercube where one side has
%     %length 2*r.
%     p=2*r*rand(numDim,1)-r;
%     
%     if(norm(p)<=r)
%         numInSphere=numInSphere+1; 
%     end
% end
% %The volume of the hypercube is (2*r)^numDim. The fraction of the samples
% %in the sphere times the volume of the cube is the approximate volume of
% %the hypersphere.
% approxVol=(numInSphere/numSamples)*(2*r)^numDim
% exactVol=hypersphereVolume(r,numDim)
%One will see that the approximate and exact volumes are close.
%
%REFERENCES:
%[1] S. Li, "Concise formulas for the area and volume of a hyperspherical
%    cap," Asian Journal of Mathematics and Statistics, vol. 4, no. 1, pp.
%    66-70, 2011.
%
%October 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    vol=pi^(numDim/2)/gamma(numDim/2+1)*r^numDim;
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