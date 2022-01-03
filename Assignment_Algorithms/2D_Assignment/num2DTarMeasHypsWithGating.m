function numHyps=num2DTarMeasHypsWithGating(GateMat)
%%NUMTARMEASHYPSWITHGATING Determine the maximum number of target-
%           measurement association hypotheses of associating targets to 
%           measurements when all targets have a missed detection
%           hypothesis, but not all targets gate with all measurements.
%           This is the number of joint association events in a JPDAF with
%           gating. In the worst-case scenario, where everything gates with
%           everything else, the function num2DTarMeasHyps is significantly
%           faster.
%
%INPUTS: GateMat A numTarXnumMeas 0-1 matrix, where GateMat(i,j)=1 if
%                target i gates with measurement j and it is 0 if target i
%                does not gate with measurement j.
%
%OUTPUTS: numHyps The total number of possible joint association events.
%
%The permanent of an NXM matrix with M>=N sums all possible combinations of
%the product of one item selected from each of the N rows and similarly for
%an MXN matrix with M>=N, it sums all possible combinations of the
%product of one item selected from each of the N columns. If we have a
%numTarXnumMeas binary matrix  where 1 indicated that a target gates with a
%measurment and 0 that a target does not gate, then the permanent of the
%matrix is the number of possible target-measurement association hypotheses
%not considering the possibility of missed detections. To handle the
%possibility of missed detections, we augment the numTarXnumMeas matrix
%with a numTarXnumTar identity matrix.
%
%While there are indubitably optimizations that could be performed to make
%this more general function faster, it is unlikely that the computational
%complexity could become polynomial in terms of the size of the problem,
%because the computation of the permanent of a general 0-1 matrix is #-P
%complete as proven in [1] and [2]. The matrix permanent has also found use
%for the computation of target-measurement association probabilities in
%[3].
%
%EXAMPLE:
%Here, we verify that the function produces the correct solution in two
%simple scenarios. In the first scenario, if GateMat is all ones, then the
%solution must be the same as num2DTarMeasHyps. In the second scenario, if
%GateMat is an NXN identity matrix, then the solution must be 2^N, because
%the only possible hypotheses for each target are observed and not
%observed.
% numTar=4;
% numMeas=5;
% GateMat=ones(numTar,numMeas);
% %These two values should be the same
% num2DTarMeasHypsWithGating(GateMat)
% num2DTarMeasHyps(numTar,numMeas)
% %These two values should be the same.
% 2^numTar
% num2DTarMeasHypsWithGating(eye(numTar,numTar))
%
%REFERENCES:
%[1] L. G. Valiant, "The complexity of computing the permanent,"
%    Theoretical Computer Science, vol. 8, no. 2, pp. 189-201, 1979.
%[2] A. Ben-Dor and S. Halevi, "Zero-one permanent is #P-complete, a
%    simpler proof," in Proceedings of the 2nd Israel Symposium on Theory
%    and Computing Systems, Natanya, Israel, 7-9 Jun. 1993, pp. 108-117.
%[3] D. F. Crouse and P. Willett, "Computation of target-measurement
%    association probabilities using the matrix permanent," IEEE
%    Transactions on Aerospace and Electronic Systems, vol. 53, no. 2, pp.
%    698-702, Apr. 2017.
%
%July 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numTar=size(GateMat,1);
numHyps=perm([GateMat,eye(numTar,numTar)]);

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
