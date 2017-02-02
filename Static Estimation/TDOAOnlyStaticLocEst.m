function sourceLoc=TDOAOnlyStaticLocEst(timeDelays,refRxLocs,nonRefRxLocs,c)
%%TDOAONLYSTATICLOCEST Given that there is an emitter simultaneously
%             detected by multiple receivers, and that one has the time
%             difference of arrival (TDOA) of a signal from the emitter to
%             various receivers with respect to one or more reference
%             receivers, determine the location of the emitter using a
%             closed-form, least-squares solution. A minimum of one
%             reference receiver and 4 TDOA measurements is needed for the
%             problem to be observable in three dimensions. In constrast,
%             the function TDOA2Cart can be used with a minimal number of
%             measurements.
%
%INPUTS: timeDelays If only one reference receiver is present, this can be
%               a numRxX1 vector of the time delays between the reach of
%               the receivers (in nonRefRxLocs) and the reference receiver.
%               In the general case, this is a numRefsX1 cell array where
%               the ith cell holds the numRxiX1 vector of delays for the
%               receivers associated with the ith reference. The type
%               (vector or cell array of vectors) of this parameter must
%               match the type of nonRefRxLocs.
%     refRxLocs If there is only one reference receiver (the receiver
%               that served as the reference for the time delays for all
%               other receivers), then this is its 3X1 Cartesian location.
%               If there are multiple reference receivers, then this is a
%               3XnumRefs matrix containing the 3X1 Cartesian
%               locations of all of the receivers.
%  nonRefRxLocs If there is only one reference receiver, this can be a
%               3XnumRx matrix containing the locations of the receivers
%               for which TDOA measurements are available. Otherwise, for
%               the general case of one or more reference receivers, this
%               is a numRefsX1 cell array where the ith cell contains a
%               3XnumRxi matrix of the locations of the receivers that took
%               TDOA measurements with the ith reference receiver. Each
%               reference receiver must have at least two receivers taking
%               TDOA measurements to contribute to solving for the position
%               of the target. A minimum of four total is required.
%             c The propagation speed in the medium in question. If this
%               parameter is omitted or an empty matrix is passed, the
%               default value of Constants.speedOfLight is used.
%
%OUTPUTS: sourceLoc The location of the emitter. In an error-free setting,
%               this is exact. In the presence of noise, this is a least
%               squares solution with respect to a non-standard cost
%               function.
%
%The algorithm implemented is taken from [1].
%
%REFERENCES:
%[1] M. D. Gillette and H. F. Silverman, "A linear closed-form algorithm
%    for source localization from time-differences of arrival," IEEE Signal
%    Processing Letters, vol. 15, pp. 1-4, 2008.
%
%April 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4||isempty(c))
   c=Constants.speedOfLight; 
end

if(isa(nonRefRxLocs,'double'))
%If only a single reference microphone is given as a receiver (cell arrays
%are not used), then directly build the matrices and solve the problem. 
    [A,w]=createAAndW(refRxLocs,nonRefRxLocs,timeDelays,c);
else
    %If multiple references are present (cell arrays are passed), then
    %building the matrix is a little bit more complciated.
    numRefs=size(refRxLocs,2);

    %count the total number of receivers for all of the references.
    totalNumRx=0;
    for curRef=1:numRefs
        totalNumRx=totalNumRx+size(nonRefRxLocs{curRef},2);
    end

    %Allocate space.
    A=zeros(totalNumRx,3+numRefs);
    w=zeros(totalNumRx,1);
    minIdx=1;
    for curRef=1:numRefs
        [ACur,wCur]=createAAndW(refRxLocs(:,curRef),nonRefRxLocs{curRef},timeDelays{curRef},c);
        
        numVals=size(ACur,1);
        idxSpan=minIdx:(minIdx+numVals-1);
        %Add the appropriate values to A and w
        A(idxSpan,1:3)=ACur(:,1:3);
        A(idxSpan,3+curRef)=ACur(:,4);
        
        w(idxSpan)=wCur;
        
        minIdx=minIdx+numVals;
    end
end
    if(size(A,1)<4)
       error('Not enough received signals to solve the problem. Consider TDOA2Cart for minimal systems.')
    end

    xs=pinv(A)*w;
    sourceLoc=xs(1:3);
end

function [A,w]=createAAndW(refRxLoc,nonRefRxLocs,timeDelays,c)
%%CREATEAANDW Create the two halves of the linear system of equations to
%             solve for the TDOA localization problem for a single
%             reference sensor. This is a separate subroutine so that it
%             can be called multiple times to construct the solution when
%             multiple reference sensors are present.

    numRx=size(nonRefRxLocs,2);

    dm0=timeDelays*c;

    %The final term in this sum is a scalar; the other two are numRxX1
    %vectors.
    w=0.5*(dm0.^2-sum(nonRefRxLocs.^2,1)'+sum(refRxLoc.^2));

    A=zeros(numRx,4);

    A(:,1:3)=bsxfun(@minus,refRxLoc',nonRefRxLocs');
    A(:,4)=dm0;
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
