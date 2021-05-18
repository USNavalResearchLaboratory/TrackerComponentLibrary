function demoCopulae()
%%DEMOCOPULAE Plots various examples of copulae for marginals drawn from a
%             Gaussian and Gumbel distribution.
%
%October 2020 Codie T. Lewis, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%Compute examples of copulae
numPoints=100;
x = linspace(0,1,numPoints);
var1 = GaussianD;
var2 = GumbelD;
[X,Y] = meshgrid(x);
Z1 = zeros(numPoints,numPoints);
Z2 = zeros(numPoints,numPoints);
Z3R1 = zeros(numPoints,numPoints);
Z3R2 = zeros(numPoints,numPoints);
Z3R3 = zeros(numPoints,numPoints);
Z3R4 = zeros(numPoints,numPoints);
Z4 = zeros(numPoints,numPoints);
Z5 = zeros(numPoints,numPoints);
Z6 = zeros(numPoints,numPoints);
Z7 = zeros(numPoints,numPoints);
for idx = 1:numPoints
    for iidx = 1:numPoints
        %Generate copulae values
        Z1(idx,iidx) = Gaussian2C.PDF([X(idx,iidx);Y(idx,iidx)],0.65);
        Z2(idx,iidx) = StudentT2C.PDF([X(idx,iidx);Y(idx,iidx)],0.65,3);
        Z3R1(idx,iidx) = Gumbel2C.PDF([X(idx,iidx);Y(idx,iidx)],2);
        Z3R2(idx,iidx) = Gumbel2C.PDF([X(idx,iidx);Y(idx,iidx)],2,90);
        Z3R3(idx,iidx) = Gumbel2C.PDF([X(idx,iidx);Y(idx,iidx)],2,180);
        Z3R4(idx,iidx) = Gumbel2C.PDF([X(idx,iidx);Y(idx,iidx)],2,270);
        Z4(idx,iidx) = Frank2C.PDF([X(idx,iidx);Y(idx,iidx)],0.65);
        Z5(idx,iidx) = Clayton2C.PDF([X(idx,iidx);Y(idx,iidx)],0.65);
        Z6(idx,iidx) = NakagamiM2C.PDF([X(idx,iidx);Y(idx,iidx)],2,0.65);
        Z7(idx,iidx) = Exponential2C.PDF([X(idx,iidx);Y(idx,iidx)],0.65);
    end
end

figure(1)
tiledlayout('flow')
sgtitle('Copulae Densities on the Unit Square')

nexttile
mesh(X,Y,Z1)
title('Gaussian')

nexttile
mesh(X,Y,Z2)
title('StudentT')

nexttile
mesh(X,Y,Z3R1)
title('Gumbel (Standard)')

nexttile
mesh(X,Y,Z3R2)
title('Gumbel (90 Deg.)')

nexttile
mesh(X,Y,Z3R3)
title('Gumbel (180 Deg.)')

nexttile
mesh(X,Y,Z3R4)
title('Gumbel (270 Deg.)')

nexttile
mesh(X,Y,Z4)
title('Frank')

nexttile
mesh(X,Y,Z5)
title('Clayton')

nexttile
mesh(X,Y,Z6)
title('Nakagami-M (M=2)')

nexttile
mesh(X,Y,Z7)
title('Exponential')

Z0 = zeros(numPoints,numPoints);
for idx = 1:numPoints
    for iidx = 1:numPoints
        X(idx,iidx) = var2.invCDF(X(idx,iidx),1,2);
        Y(idx,iidx) = var1.invCDF(Y(idx,iidx));
        Z0(idx,iidx) = var2.PDF(X(idx,iidx),1,2).*var1.PDF(Y(idx,iidx));
        Z1(idx,iidx) = Z1(idx,iidx).*var2.PDF(X(idx,iidx),1,2).*var1.PDF(Y(idx,iidx));
        Z2(idx,iidx) = Z2(idx,iidx).*var2.PDF(X(idx,iidx),1,2).*var1.PDF(Y(idx,iidx));
        Z3R1(idx,iidx) = Z3R1(idx,iidx).*var2.PDF(X(idx,iidx),1,2).*var1.PDF(Y(idx,iidx));
        Z3R2(idx,iidx) = Z3R2(idx,iidx).*var2.PDF(X(idx,iidx),1,2).*var1.PDF(Y(idx,iidx));
        Z3R3(idx,iidx) = Z3R3(idx,iidx).*var2.PDF(X(idx,iidx),1,2).*var1.PDF(Y(idx,iidx));
        Z3R4(idx,iidx) = Z3R4(idx,iidx).*var2.PDF(X(idx,iidx),1,2).*var1.PDF(Y(idx,iidx));
        Z4(idx,iidx) = Z4(idx,iidx).*var2.PDF(X(idx,iidx),1,2).*var1.PDF(Y(idx,iidx));
        Z5(idx,iidx) = Z5(idx,iidx).*var2.PDF(X(idx,iidx),1,2).*var1.PDF(Y(idx,iidx));
        Z6(idx,iidx) = Z6(idx,iidx).*var2.PDF(X(idx,iidx),1,2).*var1.PDF(Y(idx,iidx));
        Z7(idx,iidx) = Z7(idx,iidx).*var2.PDF(X(idx,iidx),1,2).*var1.PDF(Y(idx,iidx));
    end
end

% Delete X=0 and Y=0 columns for plotting marginals.
Xp = X(:,2:end);
Yp = Y(2:end,:);
Zp = Z0(2:end,2:end);

figure(2)
tiledlayout('flow')
sgtitle('Bivariate R.V. with Gaussian and Gumbel Components')

nexttile
plot(Yp(:,1),sum(Zp,2),'b','LineWidth',3)
title('R.V. 1 (Vertical)')

nexttile
plot(Xp(1,:),sum(Zp,1),'r','LineWidth',3)
title('R.V. 2 (Horizontal)')

nexttile
contourf(X,Y,Z0)
title('Independent')

nexttile
contourf(X,Y,Z1)
title('Gaussian')

nexttile
contourf(X,Y,Z2)
title('StudentT')

nexttile
contourf(X,Y,Z3R1)
title('Gumbel (Standard)')

nexttile
contourf(X,Y,Z3R2)
title('Gumbel (90 Deg.)')

nexttile
contourf(X,Y,Z3R3)
title('Gumbel (180 Deg.)')

nexttile
contourf(X,Y,Z3R4)
title('Gumbel (270 Deg.)')

nexttile
contourf(X,Y,Z4)
title('Frank')

nexttile
contourf(X,Y,Z5)
title('Clayton')

nexttile
contourf(X,Y,Z6)
title('Nakagami-M (M=2)')

nexttile
contourf(X,Y,Z7)
title('Exponential')
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
