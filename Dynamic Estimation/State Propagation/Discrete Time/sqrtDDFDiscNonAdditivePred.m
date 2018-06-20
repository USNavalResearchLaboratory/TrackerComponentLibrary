function [xPred,SPred]=sqrtDDFDiscNonAdditivePred(xPrev,SPrev,f,SQ,algorithm)
%%SQRTDDF2DISCNONADDITIVEPRED Perform the measurement update step in the
%               square-root version of the central difference filter (CDF),
%               or the first or second order divided difference filter
%               (DDF) with non-additive process noise.
%
%INPUTS: xPrev The xDim X 1 state estimate at the previous time-step.
%        SPrev The xDim X xDim lower-triangular square root of the state
%              covariance matrix at the previous time-step.
%            f A function handle for the state transition function This can
%              either be of the form f(xStacked), where xStacked is [x;w],
%              where x is the state and w is the process noise value, or
%              this can be of the form f(x,v), where the inputs are
%              separately specified. The function returns the propagated
%              state given the realization of the process noise.
%           SQ The xDimX xDim lower-triangular square root of the process
%              noise covariance matrix.
%    algorithm An optional parameter specifying which algorithm should be
%              used for the prediction. Possible values are:
%              0 Use the CDF of [3]. This is actually the same as a first-
%                order DDF of [1] and [2], but with a different step size
%                for the finite differences.
%              1 (The default if omitted or an empty matrix is passed) Use
%                the first-order divided difference filter of [1] and [2].
%              2 Use the second-order divided difference filter of [1] and
%                [2].
%
%OUTPUTS: xPred The xDim X 1 predicted state estimate.
%         SPred The xDim X xDim lower-triangular square root of the
%               predicted state covariance estimate.
%
%The divided difference filters are implemented as described in [1] and
%[2], where they are only derived for non-additive noise and the square
%root form is the most natural form to use. The central difference filter
%of [3] is the same as a first-order divided difference filter, but using a
%different distance for the finite differences.
%
%REFERENCES:
%[1] M. Nørgaard, N. K. Poulsen, and O. Ravn, "New developments in state
%    estimation for nonlinear systems," Automatica, vol. 36, no. 11, pp.
%    1627-1638, Nov. 2000.
%[2] M. Nørgaard, N. K. Poulsen, and O. Ravn, "Advances in derivative-free
%    state estimation for nonlinear systems (revised edition)," Department
%    of Informatics and Mathematical Modelling, Technical University of
%    Denmark, Tech. Rep. IMM-REP-1998-15, 29 Oct. 2004. [Online].
%    Available: http://www2.imm.dtu.dk/pubdb/views/publication details.php?id=2706
%[3] T. S. Schei, "A finite-difference method for linearization in
%    nonlinear estimation algorithms," Automatica, vol. 33, no. 11, pp.
%    2053-2058, Nov. 1997.
%
%October 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<5||isempty(algorithm))
    algorithm=1;
end

xDim=size(xPrev,1);

%This lets the user provide two different parameterization types for the
%measurement function.
if(nargin(f)<2)
    fSplit=@(x,v)f([x;v]);
else
    fSplit=f;
end

f0=fSplit(xPrev,zeros(xDim,1));

if(algorithm==0)
    %The central difference filter is the same as the first-order DDF,
    %except deltaH is just 1 not sqrt(3).
    h=1;
else
    %The term derived from the kurtosis of the assumed normalized 
    %distribution, as derived in Section 3.3 of [2].
    h=sqrt(3);
end

%These are to store values so that they do not have to be repeatedly
%recomputed.
fValsPlusX=zeros(xDim,xDim);
fValsMinusX=zeros(xDim,xDim);
fValsPlusV=zeros(xDim,xDim);
fValsMinusV=zeros(xDim,xDim);
for curDim=1:xDim
    fValsPlusX(:,curDim)=fSplit(xPrev+h*SPrev(:,curDim),zeros(xDim,1));
    fValsMinusX(:,curDim)=fSplit(xPrev-h*SPrev(:,curDim),zeros(xDim,1));
    
    fValsPlusV(:,curDim)=fSplit(xPrev,h*SQ(:,curDim));
    fValsMinusV(:,curDim)=fSplit(xPrev,-h*SQ(:,curDim));
end

%The loop fills in Sxx and Sxv according to Equation 54 in [1] (Equation 93
%in [2]).
Sxx=zeros(xDim,xDim);
Sxv=zeros(xDim,xDim);
for curDim=1:xDim
    %The zeros assumes that the mean process noise value is zero. If it is
    %not, then the f function could just be modified to reflect that.
    Sxx(:,curDim)=fValsPlusX(:,curDim)-fValsMinusX(:,curDim);
    Sxv(:,curDim)=fValsPlusV(:,curDim)-fValsMinusV(:,curDim);
end
Sxx=Sxx/(2*h);
Sxv=Sxv/(2*h);

if(algorithm==2)%Second order DDF
    %The loop fills in Sxx2 and Sxv2 according to the Equation in Section
    %4.3 of [1] (Section 4.3 of [2]). The loop also fills in xPred in
    %Equation 68 in [1] (Equation 107 in [2])
    Sxx2=zeros(xDim,xDim);
    Sxv2=zeros(xDim,xDim);
    xPred=zeros(xDim,1);
    for curDim=1:xDim
        Sxx2(:,curDim)=fValsPlusX(:,curDim)+fValsMinusX(:,curDim)-2*f0;
        Sxv2(:,curDim)=fValsPlusV(:,curDim)+fValsMinusV(:,curDim)-2*f0;
        xPred=xPred+fValsPlusX(:,curDim)+fValsMinusX(:,curDim)+fValsPlusV(:,curDim)+fValsMinusV(:,curDim);
    end
    Sxx2=(sqrt(h^2-1)/(2*h^2))*Sxx2;
    Sxv2=(sqrt(h^2-1)/(2*h^2))*Sxv2;
    xPred=xPred/(2*h^2)+((h^2-2*xDim)/(h^2))*f0;
    
    %Equation 108 in [2]
    SPred=tria([Sxx,Sxv,Sxx2,Sxv2]);
else
    %Equation 57 in [1] (Equation 96 in [2]).
    xPred=f0;
    %Equation 98 in [2].
    SPred=tria([Sxx,Sxv]);
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
