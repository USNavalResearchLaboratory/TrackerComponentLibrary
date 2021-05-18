%%COMPARE2DASSIGN Compile the C and C++ implementations of the assignment
%                 algorithms and run the various 2D assignment algorithms
%                 on random cost  matrices of differing sizes to compare
%                 the execution times. This function demonstrates how fast
%                 2D assignment can be, but also how the speed is
%                 influenced by some programmatic options such as whether
%                 one uses a direct Matlab implementation or a compiled C
%                 or C++ implementation. Additionally, whether one scans
%                 by row or by column makes a difference, because Matlab
%                 stores data by column, so scanning across columns
%                 increases the likelihood of a cache miss, slowing things
%                 down.
%
%This function requires that the CompileCLibraries function has been run so
%that the necessary functions have been compiled. This function comes with
%the additional supporting files
%assign2DAlt.cpp
%assign2DByCol.c
%assign2DByCol.m
%which all mirror the functionality of assign2D.m/assign2D.c, but which
%differ in how they are implemented so as to demonstrate the effects on
%speed of the differences. The type of comparison done is here essentially
%that done in [1] and the second half of [2].
%
%When recording execution times, the functions tic and toc are used.
%However, to speed up the Monte Carlo runs, if installed, the parallel
%processing toolbox will be used (by using parfor instead of for) so that
%the Monte Carlo runs can occur in parallel across processors/ cores. This
%means that the loop has to keep a local variable, named ticLoc here, so
%that the timer across processors can properly match starting times to
%stopping times.
%
%REFERENCES:
%[1] D. F. Crouse, "On implementing 2D rectangular assignment algorithms,"
%    IEEE Transactions on Aerospace and Electronic Systems, accepted 2016.
%[2] D. F. Crouse, "Advances in displaying uncertain estimates of multiple
%    targets," in Proceedings of SPIE: Signal Processing, Sensor Fusion, 
%    and Target Recognition XXII, Baltimore, MD, Apr. 2013.
%
%November 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%Change this value to 1000 to redo the simulations used in the paper.
numRuns=2;
fprintf('%i Monte Carlo runs shall be used used\n',numRuns)

%The search path and current working directories will be modified. Thus,
%save the old search path and working directory so that they can be
%restored when the script exits.
oldPath=path();
curDir=pwd;
ScriptPath=mfilename('fullpath');
ScriptFolder=fileparts(ScriptPath);
%Set the current working directory to the folder in which this script is
%located.
cd(ScriptFolder)

%%Test the Matlab implementations of the 2D assignment code
%Remove the compiled code directory from the search path, so that the
%Matlab implementations are used and the C and C++ implementationa are
%ignored.
compiledCodeFolder=[fileparts(fileparts(ScriptFolder)),'/0_Compiled_Code'];
rmpath(compiledCodeFolder)

disp('Running the Monte Carlo runs for the 2D assignment code in Matlab')
MatlabByRow100=zeros(numRuns,1);
MatlabByRow200=zeros(numRuns,1);
MatlabByRow500=zeros(numRuns,1);
MatlabByRow500b1000=zeros(numRuns,1);
MatlabByRow3000=zeros(numRuns,1);

MatlabByCol100=zeros(numRuns,1);
MatlabByCol200=zeros(numRuns,1);
MatlabByCol500=zeros(numRuns,1);
MatlabByCol500b1000=zeros(numRuns,1);
MatlabByCol3000=zeros(numRuns,1);
for curRun=1:numRuns
    curRun
    %100X100 Matrix
    C=rand(100,100);
    %Matlab Implementation, where the inner loop scans by row.
    ticLoc=tic;
    [col4row,row4col,gain]=assign2D(C,false);
    MatlabByRow100(curRun)=toc(ticLoc);
    %Matlab implementation where the inner loop scans by column.
    ticLoc=tic;
    [col4row,row4col,gain]=assign2DByCol(C,false);
    MatlabByCol100(curRun)=toc(ticLoc);
    
    %200X200 matrix
    C=rand(200,200);
    ticLoc=tic;
    [col4row,row4col,gain]=assign2D(C,false);
    MatlabByRow200(curRun)=toc(ticLoc);
    ticLoc=tic;
    [col4row,row4col,gain]=assign2DByCol(C,false);
    MatlabByCol200(curRun)=toc(ticLoc);
    
    %500X500 Matrix
    C=rand(500,500);
    ticLoc=tic;
    [col4row,row4col,gain]=assign2D(C,false);
    MatlabByRow500(curRun)=toc(ticLoc);
    ticLoc=tic;
    [col4row,row4col,gain]=assign2DByCol(C,false);
    MatlabByCol500(curRun)=toc(ticLoc);
    
    %500X1000 Matrix
    C=rand(500,1000);
    ticLoc=tic;
    [col4row,row4col,gain]=assign2D(C,false);
    MatlabByRow500b1000(curRun)=toc(ticLoc);
    ticLoc=tic;
    [col4row,row4col,gain]=assign2DByCol(C,false);
    MatlabByCol500b1000(curRun)=toc(ticLoc);
    
    %3000X3000 Matrix
    C=rand(3000,3000);
    ticLoc=tic;
    [col4row,row4col,gain]=assign2D(C,false);
    MatlabByRow3000(curRun)=toc(ticLoc);
    ticLoc=tic;
    [col4row,row4col,gain]=assign2DByCol(C,false);
    MatlabByCol3000(curRun)=toc(ticLoc);
end

%%Test the Matlab implementation of the k-best 2D assignment code
disp('Running the Monte Carlo runs for the k-best 2D assignment code in Matlab')
MatlabkB10b20h2=zeros(numRuns,1);
MatlabkB10b20h25=zeros(numRuns,1);
MatlabkB10b20h50=zeros(numRuns,1);

MatlabkB50b100h2=zeros(numRuns,1);
MatlabkB50b100h25=zeros(numRuns,1);
MatlabkB50b100h50=zeros(numRuns,1);

MatlabkB100b100h2=zeros(numRuns,1);
MatlabkB100b100h25=zeros(numRuns,1);
MatlabkB100b100h50=zeros(numRuns,1);
for curRun=1:numRuns
    %10X20 Matrix
    C=rand(10,20);
    k=2;
    ticLoc=tic;
    [col4row,row4col,gain]=kBest2DAssign(C,k,false);
    MatlabkB10b20h2(curRun)=toc(ticLoc);
    k=25;
    ticLoc=tic;
    [col4row,row4col,gain]=kBest2DAssign(C,k,false);
    MatlabkB10b20h25(curRun)=toc(ticLoc);
    k=50;
    ticLoc=tic;
    [col4row,row4col,gain]=kBest2DAssign(C,k,false);
    MatlabkB10b20h50(curRun)=toc(ticLoc);
    
    %50X100 Matrix
    C=rand(50,100);
    k=2;
    ticLoc=tic;
    [col4row,row4col,gain]=kBest2DAssign(C,k,false);
    MatlabkB50b100h2(curRun)=toc(ticLoc);
    k=25;
    ticLoc=tic;
    [col4row,row4col,gain]=kBest2DAssign(C,k,false);
    MatlabkB50b100h25(curRun)=toc(ticLoc);
    k=50;
    ticLoc=tic;
    [col4row,row4col,gain]=kBest2DAssign(C,k,false);
    MatlabkB50b100h50(curRun)=toc(ticLoc);
    
    %100X100 Matrix
    C=rand(100,100);
    k=2;
    ticLoc=tic;
    [col4row,row4col,gain]=kBest2DAssign(C,k,false);
    MatlabkB100b100h2(curRun)=toc(ticLoc);
    k=25;
    ticLoc=tic;
    [col4row,row4col,gain]=kBest2DAssign(C,k,false);
    MatlabkB100b100h25(curRun)=toc(ticLoc);
    k=50;
    ticLoc=tic;
    [col4row,row4col,gain]=kBest2DAssign(C,k,false);
    MatlabkB100b100h50(curRun)=toc(ticLoc);
end

%%Test the C and C++ implementations of the 2D assignment code

%Add the compiled code directory to the search path so that the Matlab
%versions of the algorithms will no longer be used
addpath(compiledCodeFolder)

disp('Running the Monte Carlo runs for the 2D assignment code in C and C++')
CPPByRow100=zeros(numRuns,1);
CPPByRow200=zeros(numRuns,1);
CPPByRow500=zeros(numRuns,1);
CPPByRow500b1000=zeros(numRuns,1);
CPPByRow3000=zeros(numRuns,1);

CByRow100=zeros(numRuns,1);
CByRow200=zeros(numRuns,1);
CByRow500=zeros(numRuns,1);
CByRow500b1000=zeros(numRuns,1);
CByRow3000=zeros(numRuns,1);

CByCol100=zeros(numRuns,1);
CByCol200=zeros(numRuns,1);
CByCol500=zeros(numRuns,1);
CByCol500b1000=zeros(numRuns,1);
CByCol3000=zeros(numRuns,1);

for curRun=1:numRuns
    %100X100 Matrix
    C=rand(100,100);
    %C++ Implementation, where the inner loop scans by row.
    ticLoc=tic;
    [col4row,row4col,gain]=assign2DAlt(C,false);
    CPPByRow100(curRun)=toc(ticLoc);
    %C implementation, where the inner loop scans by row.
    ticLoc=tic;
    [col4row,row4col,gain]=assign2D(C,false);
    CByRow100(curRun)=toc(ticLoc);
    %C implementation where the inner loop scans by column.
    ticLoc=tic;
    [col4row,row4col,gain]=assign2DByCol(C,false);
    CByCol100(curRun)=toc(ticLoc);

    %200X200 Matrix
    C=rand(200,200);
    ticLoc=tic;
    [col4row,row4col,gain]=assign2DAlt(C,false);
    CPPByRow200(curRun)=toc(ticLoc);
    ticLoc=tic;
    [col4row,row4col,gain]=assign2D(C,false);
    CByRow200(curRun)=toc(ticLoc);
    ticLoc=tic;
    [col4row,row4col,gain]=assign2DByCol(C,false);
    CByCol200(curRun)=toc(ticLoc);
    
    %500X500 Matrix
    C=rand(500,500);
    ticLoc=tic;
    [col4row,row4col,gain]=assign2DAlt(C,false);
    CPPByRow500(curRun)=toc(ticLoc);
    ticLoc=tic;
    [col4row,row4col,gain]=assign2D(C,false);
    CByRow500(curRun)=toc(ticLoc);
    ticLoc=tic;
    [col4row,row4col,gain]=assign2DByCol(C,false);
    CByCol500(curRun)=toc(ticLoc);
    
    %500X1000 Matrix
    C=rand(500,1000);
    ticLoc=tic;
    [col4row,row4col,gain]=assign2DAlt(C,false);
    CPPByRow500b1000(curRun)=toc(ticLoc);
    ticLoc=tic;
    [col4row,row4col,gain]=assign2D(C,false);
    CByRow500b1000(curRun)=toc(ticLoc);
    ticLoc=tic;
    [col4row,row4col,gain]=assign2DByCol(C,false);
    CByCol500b1000(curRun)=toc(ticLoc);
    
    %3000X3000 Matrix
    C=rand(3000,3000);
    ticLoc=tic;
    [col4row,row4col,gain]=assign2DAlt(C,false);
    CPPByRow3000(curRun)=toc(ticLoc);
    ticLoc=tic;
    [col4row,row4col,gain]=assign2D(C,false);
    CByRow3000(curRun)=toc(ticLoc);
    ticLoc=tic;
    [col4row,row4col,gain]=assign2DByCol(C,false);
    CByCol3000(curRun)=toc(ticLoc);
end

%%Test the C++ implementation of the k-best 2D assignment code
disp('Running the Monte Carlo runs for the k-best 2D assignment code in C++')
CPPkB10b20h2=zeros(numRuns,1);
CPPkB10b20h25=zeros(numRuns,1);
CPPkB10b20h50=zeros(numRuns,1);

CPPkB50b100h2=zeros(numRuns,1);
CPPkB50b100h25=zeros(numRuns,1);
CPPkB50b100h50=zeros(numRuns,1);

CPPkB100b100h2=zeros(numRuns,1);
CPPkB100b100h25=zeros(numRuns,1);
CPPkB100b100h50=zeros(numRuns,1);
for curRun=1:numRuns
    curRun
    %10X20 Matrix
    C=rand(10,20);
    k=2;
    ticLoc=tic;
    [col4row,row4col,gain]=kBest2DAssign(C,k,false);
    CPPkB10b20h2(curRun)=toc(ticLoc);
    k=25;
    ticLoc=tic;
    [col4row,row4col,gain]=kBest2DAssign(C,k,false);
    CPPkB10b20h25(curRun)=toc(ticLoc);
    k=50;
    ticLoc=tic;
    [col4row,row4col,gain]=kBest2DAssign(C,k,false);
    CPPkB10b20h50(curRun)=toc(ticLoc);
    
    %50X100 Matrix
    C=rand(50,100);
    k=2;
    ticLoc=tic;
    [col4row,row4col,gain]=kBest2DAssign(C,k,false);
    CPPkB50b100h2(curRun)=toc(ticLoc);
    k=25;
    ticLoc=tic;
    [col4row,row4col,gain]=kBest2DAssign(C,k,false);
    CPPkB50b100h25(curRun)=toc(ticLoc);
    k=50;
    ticLoc=tic;
    [col4row,row4col,gain]=kBest2DAssign(C,k,false);
    CPPkB50b100h50(curRun)=toc(ticLoc);
    
    %100X100 Matrix
    C=rand(100,100);
    k=2;
    ticLoc=tic;
    [col4row,row4col,gain]=kBest2DAssign(C,k,false);
    CPPkB100b100h2(curRun)=toc(ticLoc);
    k=25;
    ticLoc=tic;
    [col4row,row4col,gain]=kBest2DAssign(C,k,false);
    CPPkB100b100h25(curRun)=toc(ticLoc);
    k=50;
    ticLoc=tic;
    [col4row,row4col,gain]=kBest2DAssign(C,k,false);
    CPPkB100b100h50(curRun)=toc(ticLoc);
end

%Now, compute and display the median and maximum run times.
disp('2D Assignment')
disp('Problem Size, Matlab by Row Median/Worst Case, Matlab by Column Median/Worst Case')
fprintf('100 X 100   %f/%f %f/%f\n', median(MatlabByRow100),      max(MatlabByRow100),      median(MatlabByCol100),      max(MatlabByCol100))
fprintf('200 X 200   %f/%f %f/%f\n', median(MatlabByRow200),      max(MatlabByRow200),      median(MatlabByCol200),      max(MatlabByCol200))
fprintf('500 X 500   %f/%f %f/%f\n', median(MatlabByRow500),      max(MatlabByRow500),      median(MatlabByCol500),      max(MatlabByCol500))
fprintf('500 X 1000  %f/%f %f/%f\n', median(MatlabByRow500b1000), max(MatlabByRow500b1000), median(MatlabByCol500b1000), max(MatlabByCol500b1000))
fprintf('3000 X 3000 %f/%f %f/%f\n', median(MatlabByRow3000),     max(MatlabByRow3000),     median(MatlabByCol3000),     max(MatlabByCol3000))

disp('2D Assignment')
disp('Problem Size, C++ by Row Median/Worst Case, C by Row Median/Worst Case, C By Column Median/Worst Case')
fprintf('100 X 100   %f/%f %f/%f %f/%f\n', median(CPPByRow100),      max(CPPByRow100),      median(CByRow100),      max(CByRow100),      median(CByCol100),      max(CByCol100))
fprintf('200 X 200   %f/%f %f/%f %f/%f\n', median(CPPByRow200),      max(CPPByRow200),      median(CByRow200),      max(CByRow200),      median(CByCol200),      max(CByCol200))
fprintf('500 X 500   %f/%f %f/%f %f/%f\n', median(CPPByRow500),      max(CPPByRow500),      median(CByRow500),      max(CByRow500),      median(CByCol500),      max(CByCol500))
fprintf('500 X 1000  %f/%f %f/%f %f/%f\n', median(CPPByRow500b1000), max(CPPByRow500b1000), median(CByRow500b1000), max(CByRow500b1000), median(CByCol500b1000), max(CByCol500b1000))
fprintf('3000 X 3000 %f/%f %f/%f %f/%f\n', median(CPPByRow3000),     max(CPPByRow3000),     median(CByRow3000),     max(CByRow3000),     median(CByCol3000),     max(CByCol3000))

disp('k-Best 2D Assignment')
disp('Problem Size, Matlab Median/Worst Case, C++ Median/Worst Case')
fprintf('10 X 20,   2  Hyp %f/%f %f/%f\n', median(MatlabkB10b20h2),    max(MatlabkB10b20h2),    median(CPPkB10b20h2),    max(CPPkB10b20h2))
fprintf('10 X 20,   25 Hyp %f/%f %f/%f\n', median(MatlabkB10b20h25),   max(MatlabkB10b20h25),   median(CPPkB10b20h25),   max(CPPkB10b20h25))
fprintf('10 X 20,   50 Hyp %f/%f %f/%f\n', median(MatlabkB10b20h50),   max(MatlabkB10b20h50),   median(CPPkB10b20h50),   max(CPPkB10b20h50))
fprintf('50 X 100,  2  Hyp %f/%f %f/%f\n', median(MatlabkB50b100h2),   max(MatlabkB50b100h2),   median(CPPkB50b100h2),   max(CPPkB50b100h2))
fprintf('50 X 100,  25 Hyp %f/%f %f/%f\n', median(MatlabkB50b100h25),  max(MatlabkB50b100h25),  median(CPPkB50b100h25),  max(CPPkB50b100h25))
fprintf('50 X 100,  50 Hyp %f/%f %f/%f\n', median(MatlabkB50b100h50),  max(MatlabkB50b100h50),  median(CPPkB50b100h50),  max(CPPkB50b100h50))
fprintf('100 X 100, 2  Hyp %f/%f %f/%f\n', median(MatlabkB100b100h2),  max(MatlabkB100b100h2),  median(CPPkB100b100h2),  max(CPPkB100b100h2))
fprintf('100 X 100, 25 Hyp %f/%f %f/%f\n', median(MatlabkB100b100h25), max(MatlabkB100b100h25), median(CPPkB100b100h25), max(CPPkB100b100h25))
fprintf('100 X 100, 50 Hyp %f/%f %f/%f\n', median(MatlabkB100b100h50), max(MatlabkB100b100h50), median(CPPkB100b100h50), max(CPPkB100b100h50))

%Restore the old working directory.
cd(curDir);
%Restore the old path
path(oldPath);

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