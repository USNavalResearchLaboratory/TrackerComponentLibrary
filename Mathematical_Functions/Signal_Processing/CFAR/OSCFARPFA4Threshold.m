function PFA=OSCFARPFA4Threshold(T,N,k)
%%OSCFARPFA4THRESHOLD Determine the probability of false alarm in uniform
%           clutter when using the order-statistics constant false alarm
%           (OS-CFAR) detector with a simple exponential noise model as used
%           in [1].
%
%INPUTS: T The positive threshold to use for the CFAR detector. A matrix
%          of values can be passed in which case the output is a matrix
%          evaluated at each point.
%        N The total number of cells that contribute to the test region.
%          These are usually cells that are around the point being tested,
%          some distance outside of a guard region.
%        k The integer order to use k. That is, the kth largest sample in
%          the test region is used as the test statistic.
%
%OUTPUTS: PFA The false alarm probability associated with the given
%             parameters.
%
%This function implements Equation 37 of [1] with S=0. The function
%OSCFARThreshold4PFA is the inverse of this function.
%
%REFERENCES
%[1] P. P. Gandhi and S. A. Kassam, "Analysis of CFAR processors in
%    nonhomogeneous background," IEEE Transactions on Aerospace and
%    Electronic Systems, vol. 24, no. 4, pp. 427-445, Jul. 1988.
%
%February 2017  David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%Equation 37 with S=0.
PFA=ones(size(T));
for i=0:(k-1)
    curTerm=(N-i)./(N-i+T);
    PFA=PFA.*curTerm;
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
