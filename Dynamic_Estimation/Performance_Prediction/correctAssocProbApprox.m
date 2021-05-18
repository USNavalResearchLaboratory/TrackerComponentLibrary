function Pc=correctAssocProbApprox(zDim,beta,detSPred)
%%CORRECTASSOCPROBAPPROX Given the average target density per unit volume
%           over a certain region and the average determiniant of the
%           innovation covariance matrix of the targets, get the
%           approximate probability of any particular track being assigned
%           to the correct measurement assuming 100% detection probability
%           and no false alarms. This can help assess the difficult of a
%           scenario.
%
%INPUTS: zDim The dimensionality of the measurements.
%        beta The number of targets per unit volume in the coordinate
%             system of the measurements. This could be Cartesian
%             coordinates, but it need not be.
%    detSPred The determinant of the average innovation covariance matrix
%             of the targets. here, it is assumed that track accuracy for
%             all targets is about the same, so the innovation covariance
%             matrix should be about the same.
%
%OUTPUTS: Pc The approximate probability that target i is assigned to the
%            correct measurement in an ideal tracking algorithm. This
%            propbability is the same for all targets since this function
%            is considering an "average" scenario.
%
%This function implements the approximation given in Equation 2 of [1].
%One application of this function could be to determine whether a scenario
%appears to be difficult for a tracking algorithm. The value detSPred could
%be set the determinant of the measurement variance if one wants to
%consider track initiation with a rapid update rate. Then, if Pc is low, it
%is clear that track initiation will be challenging.
%
%EXAMPLE:
% zDim=3;%3D measurements.
% %Ten targets every kilometer in 3D.
% beta=(10/1e3)^zDim;
% %Innovation covariance on average to 1km meters accuracy in all
% %dimensions.
% detSPred=det(diag([1e3;1e3;1e3]));
% Pc=correctAssocProbApprox(zDim,beta,detSPred)
%One gets a correct association probability of about 74.16%.
%
%REFEREMCES:
%[1] S. Mori, K.-C. Chang, C.-Y. Chong, and K.-P. Dunn, "Prediction of
%    track purity and track accuracy in dense target environments," IEEE
%    Transactions on Automatic Control, vol. 40, no. 5, pp. 953-959, May
%    1995.
%
%January 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

m=zDim;

%We use gammaln and exp rather than just directly taking the ratio of the
%gamma functions so as the lessen possible overlow errors.
gammaRatio=exp(gammaln((m+1)/2)-gammaln(m/2+1));

%Equation 3
Cm=2^(m-1)*pi^((m-1)/2)*gammaRatio;

%Equation 2
Pc=exp(-Cm*beta*sqrt(detSPred));

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
