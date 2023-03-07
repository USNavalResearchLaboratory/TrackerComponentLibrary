function numHyps=numJPDAStarTarMeasHyps(numTar,numMeas)
%%NUMJPDASTARTARMEASHYPS Determine the maximum number of target-measurement
%             association hypotheses that are considered in the JPDA*
%             algorithm of [1]. The JPDA* algorithm throws out many of the
%             joint association events (associating targets to
%             measurements, considering missed detections), as compared to
%             the JPDAF. This is not the maximum complexity of the update
%             step, because  in the implementation of the JPDA*, there is a
%             sub-step involving 2D assignment. The complexity and the
%             number of hypotheses is discussed in [2].
%
%INPUTS: numTar The number of target present (numTar>=0).
%        numMeas The number of measurements present (numMeas>=0).
%
%OUTPUTS: numHyps The total number of hypotheses considered in the JPDA*
%                 algorithm.
%
%The number is equal to Equation 16a in [1] removing the l! term (because
%only one hypothesis for each combination of observed targets and target-
%originated measurements is considered). The result can be expressed as a
%ratio of factorials, which we evaluate as the exponential of the
%logarithms of gamma functions so as to be more robust to finite precision
%limitations.
%
%REFERENCES:
%[1] H. A. P. Blom, E. A. Bloem, and D. Musicki, "JIPDA*: Automatic
%    target tracking avoiding track coalescence," IEEE Transactions on
%    Aerospace and Electronic Systems, vol. 51, no. 2, pp. 962-974, Apr.
%    2015.
%[2] D. F. Crouse, Y. Bar-Shalom, P. Willett, and L. Svensson, "The JPDAF
%    in practical systems: Computation and snake oil," in Proceedings of
%    SPIE: Signal and Data Processing of Small Targets, vol. 7698, Orlando,
%    FL, 5 Apr. 2010.
%
%July 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%Use gammaln and exp instead of factorials to handle big numbers and then
%round so as to deal with finite precision limitations.
numHyps=round(exp(gammaln(numMeas+numTar+1)-gammaln(numMeas+1)-gammaln(numTar+1)));

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
