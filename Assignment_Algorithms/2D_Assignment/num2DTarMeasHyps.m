function numHyp=num2DTarMeasHyps(numTar,numMeas)
%%NUM2DTARMEASHYPS Determine the maximum number of target-measurement
%             association hypotheses when all targets could be associated
%             with all measurements and all targets have a missed detection
%             hypothesis. This is the maximum complexity of the update step
%             of a JPDAF. It is the case where everything gates with
%             everything else.
%
%INPUTS:  numTar The number of target present (numTar>=0).
%        numMeas The number of measurements present (numMeas>=0).
%
%OUTPUTS: numHyp The number of possible ways of assigning targets to
%                measurements taking into account that not all targets will
%                be observed.
%
%The formula for the number of hypotheses is given in [1]. This function
%avoid overflows in the numerator and denominator of each term by
%computating the next term by multipying the previous term by an
%appropriate value.
%
%REFERENCES:
%[1] D. F. Crouse, Y. Bar-Shalom, P. Willett, and L. Svensson, "The JPDAF
%    in practical systems: Computation and snake oil," in Proceedings of
%    SPIE: Signal and Data Processing of Small Targets, vol. 7698, Orlando,
%    FL, 5 Apr. 2010.
%
%October 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numHyp=1;
curTerm=1;
for l=1:min(numMeas,numTar)
    curTerm=curTerm*(numMeas-l+1)*(numTar-l+1)/l;
    numHyp=numHyp+curTerm;
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
