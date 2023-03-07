function [Rsp,U,V]=standardUVBeamPattern(T,xyPoints,normVals,Sigma,lineParams,bounds,numPoints)
%%STANDARDUVBEAMPATTERN Compute the beam pattern of a tapered array as a
%                       function of direction given in terms of direction
%                       cosines for a narrowband linear or planar array.
%                       The beam pattern is just the sum of all of the
%                       tapered elements taken for signals in different
%                       directions. This function can output the response
%                       over all (-1,+1) u-v values or in a rectangular
%                       subset. It can also output a linear cut across the
%                       values.
%
%INPUTS: T The numSubarraysXnumElements tapering matrix of the array. If
%          there are no subarrays, then this can be a numElementsX1
%          vector of weights for every element. The weights can be complex,
%          which means that difference beams can be formed and steering can
%          be taken into account. If an empty matrix is passed, then it is
%          assumed that no tapering is used so T will be a 1XnumElements
%          vector of all ones.
% xyPoints A 1XNumDim or 2XnumDim matrix of the [x;y] locations of the
%          points in the linear or planar array. The units of the distances
%          are in terms of wavelengths for the narrowband model.
% normVals This indicates how the sum beam value should be normalized.
%          Possible values are:
%          'ArrayGain' Return the ratio of the output power to the noise
%                      power. This is the default if this parameter is
%                      omitted or an empty matrix is passed.
%          'NormPowGain' Return the power of the output normalized such
%                        that the highest value is 1.
%          'AbsPowGain' Return the power of the output. This assumes that
%                       the tapering matrix contains the true gain or loss
%                       values for the tapering and is thus not scaled by
%                       any constant value.
%          'NormRealVal' Display the normalized real component of the
%                       output value. It is normalized such that the
%                       largest absolute value is one. Note that this is
%                       not squared, so it corresponds to an output
%                       voltage, not a power.
%           'RawOutput' Provide the raw sum beam output.
%   Sigma If the array gain is desired (normVals='ArrayGain'), then this
%         parameter is used. This is the numElsXnumEls covariance matrix of
%         the noise at the individual elements in the array. If this
%         parameter is omitted or an empty matrix is passed, then the
%         identity matrix will be used.
% lineParams If xyPoints is a 2XnumDims set of points (for a planar array),
%         then if this parameter is provided and is not an empty matrix,
%         Rsp will be a 1D cut across u and v values rather than all u and
%         v values. The equation for the line along with the beam pattern
%         will be evaluated is v=lineParams.intercept+lineParams.slope*u if
%         the value lineParams.vIndep=false or the vIndep component of
%         lineParams is omitted. Otherwise, the line is
%         u=lineParams.intercept+lineParams.slope*v .
%  bounds A 2X1 (or 1X2) or a 4X1 (or 1X4) vector with the bounds in u and
%         v of the plot. For 1-dimensional plots, which are the case if
%         xyPoints is 1XnumDim or lineParams is provided and xyPoints is
%         2XnumDim, then bounds=[minVal;maxVal] for the independent
%         variable. For two-dimensional plots, then
%         bounds=[minU;maxU;minV;maxV]. If this parameter is omitted or an
%         empty matrix is passed, then [-1;1;-1;1] is used to go over all
%         possible values in u and v (or just u for a 1D plot).
% numPoints A parameter determining the size of the output matrix. For 1D
%         plots, Rsp is numPointsX1. For 2D plots, Rsp is
%         numPointsXnumPoints. If this parameter is omitted or an empty
%         matrix is passed, then a default of numPoints=125 points is used.
%
%OUTPUTS: Rsp The array beam pattern over the selected region, normalized
%             as specified. The value in entry Rsp(i,j) corresponds to the
%             U and V values U(i,j) and V(i,j). For a 1D response for a
%             linear array, Rsp(i) corresponds to U(i) and if xyPoints is
%             2XnumDims, V(i) is the corresponding v value. The value of
%             Rsp for points outside of the visible region (u^2+v^2>1)
%             is set to zero.
%           U A matrix of points corresponding to the u values of the
%             responses in Rsp.
%           V A matrix of points corresponding to the v values of the
%             responses in Rsp. For 1D responses,
%
%The idea of a beam pattern for a narrowband array is discussed in Chapter
%2.2 of [1]. The array gain is discussed in Chapter 2.6.2.
%
%EXAMPLE 1:
%Here, we find the array gain beam pattern of a 20 element 1D linear array
%with Taylor tapering and half-wavelength element spacing. When consdiering
%the array gain, it does not matter that the elements are not provided
%centered about the origin. However, the points must be centered to get the
%proper tapering values from the TaylorLinearTapering function.
% xPoints=0:0.5:9.5;
% xPoints=xPoints-mean(xPoints);
% nBar=3;
% sidelobedB=-25;
% T=TaylorLinearTapering(nBar,sidelobedB,xPoints).';
% [Rsp,U]=standardUVBeamPattern(T,xPoints);
% figure(1)
% clf
% plot(U,10*log10(Rsp),'linewidth',2)
% axis([-1 1 -30 15])
% axis square
% h1=xlabel('u');
% h2=ylabel('Array Gain');
% title('Array Power Gain in Decibels')
% set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')
%
%EXAMPLE 2:
%In this example, we plot the normalized power gain of a 2D hexagonal array
%without any tapering.
% xyPoints=getShaped2DLattice([7;14],'hexagonal');
% [Rsp,U,V]=standardUVBeamPattern([],xyPoints,'NormPowGain');
% figure(2)
% clf
% surface(U,V,10*log10(Rsp),'EdgeColor','None')
% axis([-1 1 -1 1 -40 0])
% caxis([-40 0])
% colormap(jet(256))
% colorbar()
% view(45,45)
% light()
% h1=xlabel('u');
% h2=ylabel('v');
% h3=zlabel('Array Gain');
% title('Array Power Gain in Decibels')
% set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h3,'FontSize',14,'FontWeight','bold','FontName','Times')
%
%EXAMPLES 3:
%In this example, we plot the difference pattern of a circular array.
%First, we plot the normalized real value of the pattern (the pattern is
%real anyway so this is not an issue) in 2D and then we take 1D cuts of it
%in two different directions. Bayliss tapering are used to obtain the
%difference pattern.
% xyPoints=getShaped2DLattice([30;30],'circular');
% sidelobedB=-30;
% N=17;
% T=BaylissTapering(sidelobedB,N,xyPoints).';
% [Rsp,U,V]=standardUVBeamPattern(T,xyPoints,'NormRealVal');
% figure(3)
% clf
% surface(U,V,Rsp,'EdgeColor','None')
% axis([-1 1 -1 1 -1 1])
% caxis([-1 1])
% colormap(jet(256))
% colorbar()
% view(45,10)
% light()
% h1=xlabel('u');
% h2=ylabel('v');
% h3=zlabel('Real Response');
% title('Array Difference  Beam Pattern')
% set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h3,'FontSize',14,'FontWeight','bold','FontName','Times')
% 
% %Now, we take a cut along to v=0 line.
% lineParams=[];
% lineParams.intercept=0;
% lineParams.slope=0;
% lineParams.vIndep=false;
% 
% [Rsp,U]=standardUVBeamPattern(T,xyPoints,'NormRealVal',[],lineParams);
% figure(4)
% clf
% plot(U,Rsp,'linewidth',2)
% h1=xlabel('u');
% h2=ylabel('Real Response');
% title('Array Difference Beam Pattern')
% set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')
% 
% %And here we take a cut along the u=0 line.
% lineParams.vIndep=true;
% 
% [Rsp,~,V]=standardUVBeamPattern(T,xyPoints,'NormRealVal',[],lineParams);
% figure(5)
% clf
% plot(V,Rsp,'linewidth',2)
% h1=xlabel('v');
% h2=ylabel('Real Response');
% title('Array Difference Beam Pattern')
% set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')
%
%EXAMPLE 4:
%This is an example of using elements of a circular array with a fixed
%Taylor tapering that have been broken into subarrays. We then try to form
%the best approximation to a Bayliss difference pattern modifying only the
%subarray outputs. We also demonstrate how steering can be used to move the
%difference beam away from the center of the array.
% xyPoints=getShaped2DLattice([30;30],'circular');
% numEls=size(xyPoints,2);
% %Fixed element-level tapering.
% sidelobedB=-30;
% nBar=4;
% g=TaylorTapering(nBar,sidelobedB,xyPoints);
% %It is assumed the disjoint subarrays can be formed (no adjacency matrix
% %used).
% numLevels=3;
% N=17;
% T=double(findBaylissSubarrays(xyPoints,numLevels,sidelobedB,N));
% %We now apply the tapering to the matrix
% T=bsxfun(@times,T,g.');
% 
% %The desired element-level tapering
% g=BaylissTapering(sidelobedB,N,xyPoints);
% 
% %Next, we try to find the best subarray weights to approximate the desired
% %tapering. We specifically add a null at the boresight.
% g=findSubarrayWeights(T,g,ones(numEls,1));
% 
% %Apply the tapering to the elements
% T=bsxfun(@times,g,T);
% 
% %We also steer the array off boresight to u0. This is element-level
% %steering.
% u0=[-0.5;0.5];
% D=diag(exp(1j*2*pi*sum(bsxfun(@times,xyPoints,u0),1)));
% T=T*D;
% 
% [Rsp,U,V]=standardUVBeamPattern(T,xyPoints,'NormRealVal');
% figure(3)
% clf
% surface(U,V,Rsp,'EdgeColor','None')
% axis([-1 1 -1 1 -1 1])
% caxis([-1 1])
% colormap(jet(256))
% colorbar()
% view(45,10)
% light()
% h1=xlabel('u');
% h2=ylabel('v');
% h3=zlabel('Real Response');
% title('Array Difference  Beam Pattern')
% set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h3,'FontSize',14,'FontWeight','bold','FontName','Times')
%
%REFERENCES:
%[1] H. L. Van Trees, Optimum Array Processing. New York: Wiley-
%    Interscience, 2002.
%
%August 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numDim=size(xyPoints,1);
numEls=size(xyPoints,2);

if(isempty(T))
   T=ones(1,numEls);
end

if(nargin<7||isempty(numPoints))
   numPoints=125; 
end

if(nargin<6||isempty(bounds))
    if(numDim==1)
        bounds=[-1;1];
    else
        bounds=[-1;1;-1;1];
    end
end

if(nargin<5||isempty(lineParams))
    if(numDim>1)    
        %If a 2D plot is desired.
        uVals=linspace(bounds(1),bounds(2),numPoints);
        vVals=linspace(bounds(3),bounds(4),numPoints);
        [U,V]=meshgrid(uVals,vVals);
        Rsp=zeros(numPoints,numPoints);
        
        if(isvector(T))
            for curVal=1:(numPoints*numPoints)
                u=[U(curVal);V(curVal)];
                if(u'*u>1)
                    continue;
                end

                Rsp(curVal)=T*exp(-1j*2*pi*sum(bsxfun(@times,xyPoints,u),1)).';
            end
        else
            for curVal=1:(numPoints*numPoints)
                u=[U(curVal);V(curVal)];
                if(u'*u>1)
                    continue;
                end

                Rsp(curVal)=sum(sum(T*exp(-1j*2*pi*sum(bsxfun(@times,xyPoints,u),1)).'));
            end
        end
    else
        %If a 1D plot is required.
        U=linspace(bounds(1),bounds(2),numPoints)';
        V=[];
        
        if(isvector(T))
            Rsp=zeros(numPoints,1);
            for curVal=1:numPoints
                if(abs(U(curVal))>1)
                    continue;
                end

                Rsp(curVal)=T*exp(-1j*2*pi*xyPoints*U(curVal)).';
            end
        else
            Rsp=zeros(numPoints,1);
            for curVal=1:numPoints
                if(abs(U(curVal))>1)
                    continue;
                end

                Rsp(curVal)=sum(T*exp(-1j*2*pi*xyPoints*U(curVal)).');
            end
        end
    end
else
    %If a 1D cut across the 2D surface is desired. The
    slope=lineParams.slope;
    intercept=lineParams.intercept;
    if(isfield(lineParams,'vIndep'))
        vIndep=lineParams.vIndep;
    else
        vIndep=false;
    end
    
    if(vIndep==false)
        %The independent variable is u.
        U=linspace(bounds(1),bounds(2),numPoints)';
        V=zeros(numPoints,1);
        Rsp=zeros(numPoints,1);
        
        if(isvector(T))
            for curVal=1:numPoints
                V(curVal)=intercept+slope*U(curVal);
                uVec=[U(curVal);V(curVal)];

                if(uVec'*uVec>1)
                    continue;
                end

                Rsp(curVal)=T*exp(-1j*2*pi*sum(bsxfun(@times,xyPoints,uVec),1)).';
            end
        else
            for curVal=1:numPoints
                V(curVal)=intercept+slope*U(curVal);
                uVec=[U(curVal);V(curVal)];

                if(uVec'*uVec>1)
                    continue;
                end

                Rsp(curVal)=sum(sum(T*exp(-1j*2*pi*sum(bsxfun(@times,xyPoints,uVec),1)).'));
            end
        end
    else
        %The independent variable is v.
        V=linspace(bounds(1),bounds(2),numPoints)';
        U=zeros(numPoints,1);
        Rsp=zeros(numPoints,1);
        
        if(isvector(T))
             for curVal=1:numPoints
                U(curVal)=intercept+slope*V(curVal);
                uVec=[U(curVal);V(curVal)];

                if(uVec'*uVec>1)
                    continue;
                end

                Rsp(curVal)=T*exp(-1j*2*pi*sum(bsxfun(@times,xyPoints,uVec),1)).';
            end
        else
            for curVal=1:numPoints
                U(curVal)=intercept+slope*V(curVal);
                uVec=[U(curVal);V(curVal)];

                if(uVec'*uVec>1)
                    continue;
                end

                Rsp(curVal)=sum(sum(T*exp(-1j*2*pi*sum(bsxfun(@times,xyPoints,uVec),1)).'));
            end
        end
    end
end

if(nargin<3||isempty(normVals))
    normVals='ArrayGain';
end

switch(normVals)
    case 'ArrayGain'%Array power gain versus spatially white noise.
        if(nargin<4||isempty(Sigma))
           Sigma=eye(numEls,numEls);
        end
        
        if(isvector(T))
            Rsp=abs(Rsp).^2/sum(T*Sigma*T');
        else
        
            numSubarrays=size(T,1);
            e=ones(numSubarrays,1);
            %The sum sums up all of the subarray outputs.
            Rsp=abs(Rsp).^2/sum(e'*T*Sigma*T'*e);
        end
    case 'NormPowGain'%Normalized absolute array power gain
        Rsp=abs(Rsp);
        Rsp=(Rsp/max(Rsp(:))).^2;%Make it relative to the peak.
    case 'AbsPowGain'
        Rsp=abs(Rsp)^2;
    case 'NormRealVal'%Display the real component of the response.
        Rsp=real(Rsp);
        Rsp=Rsp/max(Rsp(:));%Make it relative to the peak.
    case 'RawOutput'
    otherwise
        error('Unknown plot type requested')
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
