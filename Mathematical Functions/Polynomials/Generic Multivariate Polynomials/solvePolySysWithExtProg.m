function [xSol,otherInfo]=solvePolySysWithExtProg(polyStrings,varNames,algorithm,opts,scratchFolderPath,execPath)
%%SOLVEPOLYSYSWITHEXTPROG Solve a system of simultaneous multivariate
%            polynomials using an external program (which must be installed
%            on this computer). Programs supported are Bertini, PHCpack and
%            the certified solver in NAG4M2 in Macaulay2. All of the
%            programs read input files and then write their results to
%            disk. This function formats the input, and then reads the
%            output.
%
%INPUTS: polyStrings A numPolyX1 or 1XnumPoly cell array containing the
%           polynomials given as strings. The variable names should begin
%           with letters but may containg numbers. Imaginary coefficients
%           should be indicated by i in Bertini and PHCpack and by ii when
%           using NAG4M2 in Macaulay2. An example of a valid equation for
%           Bertini and PHCpack is 'x^2+(1+2*i)*y^2'. Variable names should
%           not contain numbers when using Bertini.
%  varNames A numVarX1 or 1XnumVar cell array containing the names
%           used for the different variables in polyStrings. For example,
%           'x'. Generally, numVar=numPoly. No variable is allowed to be
%           named 'i', 'ii' or be named 'eqn' followed by a number.
%           Variable names should not contain numbers when using Bertini.
% algorithm An optional parameter specifying which algorithm to use. The
%          solvers for simultaneous multivariate polynomials are external
%          programs. Possible values are:
%          0 (The default if omitted or an empty matrix is passed) Use
%            Bertini.
%          1 Use PHCpack.
%          2 Use the certified homotopy algorithm that is built into
%            Macaulay2 (in NAG4M2). This only uses the normalized solver
%            with the default options.
%      opts An optional input specifying options for the solver. This is
%           not used with Macaulay2; with PHCpack, it is a string of extra
%           options to pass to the solver. For example, '-t8' makes the
%           solver use eight threads to solve the problem, speeding up the
%           solution on multicore systems. For Bertini, this is a structure
%           configuration containing parameters to set. The name of the
%           entry in the structure elects the option and the value is the
%           value associated with the option. For example, if one wishes to
%           change the final tolerance, one can just use
%           opts.EndpointFiniteThreshold=1e-14. If Bertini has been
%           compiled with support for the message passing interface and the
%           mpirun command exists on this computer, then one can include
%           the parameter MPIRunProcs set to an integer above 0 and that
%           many processes will be used to parallelize Bertini. Setting the
%           number too high will result in slower performance. This
%           parameter can be omitted or an empty matrix passed if default
%           values should not be changed.
% scratchFolderPath An optional parameter specifying the folder in which
%           temporary files will be written. This is needed, because
%           all solvers only works through files. If this parameter is
%           omitted or an empty matrix is passed, then a folder named temp
%           in the folder enclosing this function is used. Note that if an
%           error occurs while executing this function, then temporary
%           files might be left in the temp folder after this function
%           ends.
%  execPath The command line command to use to execute the solver. The 
%           default if this parameter is omitted or an empty matrix is 
%           passed is just the standard name on the command line: bertini,
%           phc and M2 for each of the algorithms. The default assumes that
%           the program is already in the default search path. For example,
%           on *NIX systems, one can usually put the executable in
%           /usr/local/bin for it to be in the default search path.
%           However, for some reason, that directory is often not
%           included in Matlab's search path. Matlab 2016a's help under
%           "Run External Commands, Scripts, and Programs" specifically
%           says how to add that folder to the default search path.
%
%OUTPUTS: xSol The numVarXnumSol isolated (complex) solutions found to the
%              polynomial. The variables are in the same order as they
%              appear in varNames.
%    otherInfo Other info returned by the solver. For algorithm 0, this is
%              the text returned by the solver on the command line. For
%              algorithm 2, this is not used and an empty matrix is
%              returned. For algorithm 1, this is a structure with elements
%              multiplicity A numSolX1 vector of the multiplicity of each
%                  solution.
%              err A numSolX1 vector of magnitudes of last correction term 
%                  used in Newton's method when solving for each root.
%              rco A numSolX1 vector of the estimated inverse condition 
%                  numbers of the roots. This relates to numerical
%                  stability.
%              res A numSolX1 vector of magnitudes of the polynomial vector
%                  evaluated at each root.
%
%Note that if the external solver fails, then this function will usually
%have an error when trying to read the solutions and the temporary files
%will not always be deleted. Note that one can place calls to this function
%in try-catch statements to keep it from having an error if the solver
%fails.
%
%Binaries and source code for Bertini can be downloaded from
%https://bertini.nd.edu
%The primary source of documentation for Bertini is the book [1]. If one
%wishes to parallelize Bertini, then one can needs to install MPICH (or a
%similar message passing interface program). obtain source code and
%compiler instructions for MPICH from
%http://www.mpich.org
%
%Binaries and source code for PHCPack can be downloaded from
%http://homepages.math.uic.edu/~jan/download.html
%A description of the implementation of the original algorithm is given in
%[2].
%
%Binaries and source code for Macaulay2 may be downloaded from
%http://www.math.uiuc.edu/Macaulay2/
%The algorithm implemented in Macaulay2 that is called here is described in
%[3].
%
%EXAMPLE:
% polyStrings=cell(2,1);
% polyStrings{1}='-x1^2+2*x1*x2+x2^2+5*x1-3*x2-4';
% polyStrings{2}='x1^2+2*x1*x2+x2^2-1';
% varNames={'x1','x2'};
% [xSol,otherInfo]=solvePolySysWithExtProg(polyStrings,varNames)
% %One should get the solutions (0,-1), (1,0), (3,-2), and (4,-5).
%
%REFERENCES:
%[1] D. J. Bates, A. J. Sommese, J. D. Hauenstein, and C. W. Wampler,
%    Numerically Solving Polynomial Systems with Bertini. Philadelphia:
%    Society for Industrial and Applied Mathematics, 2013.
%[2] J. Verschelde, "Algorithm 795: PHCpack: A general-purpose solver for
%    polynomial systems by homotopy continuation," ACM Transactions on
%    Mathematical Software, vol. 25, no. 2, pp. 251-276, Jun. 1999.
%[3] C. Beltrán and A. Leykin, "Certified numerical homotopy tracking,"
%    Experimental Mathematics, vol. 21, no. 1, 2012.
%
%May 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(algorithm))
    algorithm=0;%Bertini is the default.
end

if(nargin<4)
   opts=[]; 
end

if(nargin<5)
    scratchFolderPath=[];
end

if(nargin<6)
    execPath=[];
end

switch(algorithm)
    case 0%Bertini
        if(~isempty(opts))
            if(isfield(opts,'MPIRunProcs'))
                MPIRunProcs=opts.MPIRunProcs;
                opts=rmfield(opts,'MPIRunProcs');
            else
                MPIRunProcs=[];
            end
        else
            MPIRunProcs=[];
        end

        [xSol,solverText]=solvePolySysUsingBertini(polyStrings,varNames,opts,MPIRunProcs,scratchFolderPath,execPath);
        otherInfo=solverText;
    case 1%PHCpack
        [xSol,multiplicity,err,rco,res]=solvePolySysUsingPHCpack(polyStrings,varNames,opts,scratchFolderPath,execPath);
        otherInfo.multiplicity=multiplicity;
        otherInfo.err=err;
        otherInfo.rco=rco;
        otherInfo.res=res;
    case 2%NAG4M@ in Macaulay2
        xSol=solvePolySysUsingMacaulay2(polyStrings,varNames,scratchFolderPath,execPath);
        otherInfo=[];
    otherwise
        error('Invalid algorithm specified')
end
end

%%%%%%%%%%%%%%%%%%%%%%%
%%FUNCTIONS FOR BERTINI
%%%%%%%%%%%%%%%%%%%%%%%
function [xSol,solverText]=solvePolySysUsingBertini(polyStrings,varNames,opts,MPIRunProcs,scratchFolderPath,BertiniPath)
%%SOLVEPOLYSYSUSINGBERTINI Solve a system of simultaneous multivariate
%           polynomials using the external program Bertini. This function
%           just formats the input, properly passes it to Bertini and reads
%           the output. Bertini only takes inputs as files on disk and
%           provides outputs as files. It is assumed that the polynomial
%           system has a finite number of soluitons. Note that this
%           function cannot be used in a parfor loop unless care is taken
%           to allow different function instances to use different scartch
%           folders.
%
%INPUTS: polyStrings A numPolyX1 or 1XnumPoly cell array containing the
%           polynomials given as strings. The variable names should begin
%           with letters but may containg numbers. Imaginary coefficients
%           should be indicated by i. An example of a valid equation is
%           'x^2+(1+2*i)*y^2'.
%  varNames A numVarX1 or 1XnumVar cell array containing the names
%           used for the different variables in polyStrings. For example,
%           'x'. Generally, numVar=numPoly. No variable is allowed to be
%           named'i' or be named 'eqn' followed by a number.
%      opts An optional  structure containg options to pass to Bertini.
%           These are the same options that one would use when making an
%           input file on the command line. The name of the entry in the
%           structure elects the option and the value is the value
%           associated with the option. For example, if one wishes to
%           change the final tolerance, one can just use
%           opts.EndpointFiniteThreshold=1e-14. THis parameter can be
%           omitted or an empty matrix passed if default values should not
%           be changed.
%  MPIRunProcs Bertini can be compiled to run multiple processes using the
%           message passing interface standard (MPI). If that is the case,
%           then one can provide this parameter to specify how many threads
%           should be run. If this parameter is omitted, an empty matrix is
%           passed, or if MPIRunProcs=0, then Bertini is called without
%           using the mpi program (the mpirun command on the command line).
% scratchFolderPath An optional parameter specifying the folder in which
%           temporary files will be written. This is needed, because
%           Bertini only works through files. If this parameter is omitted
%           or an empty matrix is passed, then a folder named temp in the
%           folder enclosing this function is used. Note that if an error
%           occurs while executing this function, then temporary files
%           might  be left in the temp folder.
% BertiniPath The command line command to use to execute Bertini. The 
%           default if this parameter is omitted or an empty matrix is 
%           passed is just bertini, which assumes that the program is
%           already in the default search path. For example, on *NIX
%           systems, one can usually put the executable in /usr/local/bin
%           for it to be in the default search path.
%
%OUTPUTS: xSol The numVarXnumSol isolated (complex) solutions found to the
%              polynomial. The variables are in the same order as they
%              appear in varNames.
%   solverText The text of the status message that is printed by bertini on
%              the command line when it exits.
%
%This function work by writing the commands for Bertini to a file and then
%calling Bertini to solve the problem. Bertini does not allow the output
%files to be given non-default names. Thus, one cannot run multiple
%instances of this file in a parfor loop in Matlab without making each
%instance use a different scratch folder.
%
%Binaries and ource code for Bertini can be downloaded from
%https://bertini.nd.edu
%The primary source of documentation for Bertini is the book [1].
%
%EXAMPLE:
% polyStrings=cell(2,1);
% polyStrings{1}='-x1^2+2*x1*x2+x2^2+5*x1-3*x2-4';
% polyStrings{2}='x1^2+2*x1*x2+x2^2-1';
% varNames={'x1','x2'};
% [xSol]=solvePolySysUsingBertini(polyStrings,varNames)
%
%REFERENCES:
%[1] D. J. Bates, A. J. Sommese, J. D. Hauenstein, and C. W. Wampler,
%    Numerically Solving Polynomial Systems with Bertini. Philadelphia:
%    Society forIndustrial and Applied Mathematics, 2013.
%
%May 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3)
   opts=[]; 
end

if(nargin<4)
   MPIRunProcs=[]; 
end

if(nargin<5||isempty(scratchFolderPath))
    %The default is to put the temporary files into a folder called temp
    %that is located in the same folder as this script.
    ScriptPath=mfilename('fullpath');
    scratchFolderPath=[fileparts(ScriptPath),'/temp/'];
end

if(nargin<6||isempty(BertiniPath))
    %The default assumes that the program is already in the user's path and
    %can be called directly.
    BertiniPath='bertini';
end

if(~isempty(MPIRunProcs)&&MPIRunProcs~=0)
    commandStart=['mpirun -np ',num2str(MPIRunProcs),' '];
else
    commandStart=[];
end

oldFolder=pwd();
cd(scratchFolderPath);

writePolynomialForBertini('input',polyStrings,opts,varNames);
%Call Bertini; this created a bunch of output files with names that cannot
%be programmatically changed.
str=[commandStart,BertiniPath,' input'];
[~,solverText]=system(str);

xSol=solutionsFromFileBertini('finite_solutions',length(varNames));

%Delete files created. bertini makes lots of file.
delete('input');
delete('finite_solutions');
delete('failed_paths');
delete('main_data');
delete('midpath_data');
delete('nonsingular_solutions');
delete('output');
delete('raw_data');
delete('raw_solutions');
delete('real_finite_solutions');
delete('singular_solutions');
delete('start');

cd(oldFolder);

end

function writePolynomialForBertini(inFilename,polyStrings,opts,varNames)
numVar=length(varNames);
numEq=length(polyStrings);

id = fopen(inFilename,'wt');
if(~isempty(opts))
    fprintf(id,'CONFIG\n');
    names=fieldnames(opts);
    numOpts=length(names);
    
    for curOpt=1:numOpts
        fprintf(id,[names{curOpt},': ',num2str(opts.(names{curOpt}),16),';\n']);
    end
    fprintf(id,'END;\n\n');
end

fprintf(id,'INPUT\nvariable_group');
fprintf(id,[' ',varNames{1}]);
for curVar=2:numVar
    fprintf(id,[', ',varNames{curVar}]);
end
fprintf(id,';\n');

fprintf(id,'function');
fprintf(id,[' eqn',num2str(1)]);
for curEq=2:numEq
    fprintf(id,[', eqn',num2str(curEq)]);
end
fprintf(id,';\n');

for curPoly=1:numEq
    fprintf(id,['eqn',num2str(curPoly),'=',polyStrings{curPoly},';\n']);
end
fprintf(id,'END;\n');
fclose(id);
end

function xSol=solutionsFromFileBertini(outfileName,numVar)
id = fopen(outfileName, 'rt');
solText = fread(id,'*char')';
fclose(id);

%Read all of the numbers. The first number is the number of solutions. All
%others are actual values of the solutions.
A=sscanf(solText,'%f');

numSol=A(1);

%Combine real and imaginary parts.
xSol=reshape(A(2:2:end)+1i*A(3:2:(end)),numVar,numSol);
end

%%%%%%%%%%%%%%%%%%%%%%%
%%FUNCTIONS FOR PHCPACK
%%%%%%%%%%%%%%%%%%%%%%%
function [xSol,multiplicity,err,rco,res]=solvePolySysUsingPHCpack(polyStrings,varNames,opts,scratchFolderPath,PHCPath)
%%SOLVEPOLYSYSUSINGPHCPACK Solve a system of simultaneous multivariate
%           polynomials using the external program PHCpack. This function
%           just formats the input, properly passes it to PHCpack and reads
%           the output. PHCpack only takes inputs as files on disk and
%           provides outputs as files. It is assumed that the polynomial
%           system has a finite number of soluitons.
%
%INPUTS: polyStrings A numPolyX1 or 1XnumPoly cell array containing the
%           polynomials given as strings. The variable names should begin
%           with letters but may containg numbers. Imaginary coefficients
%           should be indicated by i. An example of a valid equation is
%           'x^2+(1+2*i)*y^2'.
%  varNames A numVarX1 or 1XnumVar cell array containing the names
%           used for the different variables in polyStrings. For example,
%           'x'. Generally, numVar=numPoly. No variable is allowed to be
%           named 'I' or 'i'.
%      opts An optional string of additional inputs that should be passed
%           to PHCpack on the command line. For example, '-t8' makes the
%           solver use eight threads to solve the problem, speeding up the
%           solution on multicore systems. If omitted, an empty matrix is
%           used.
% scratchFolderPath An optional parameter specifying the folder in which
%           temporary files will be written. This is needed, because
%           PHCPack only works through files. If this parameter is omitted
%           or an empty matrix is passed, then a folder named temp in the
%           folder enclosing this function is used. Note that if an error
%           occurs while executing this function, then temporary files
%           might  be left in the temp folder.
%   PHCPath The command line command to use to execute PHCpack. The default
%           if this parameter is omitted or an empty matrix is passed is
%           just phc, which assumes that the program is already in the
%           default search path. For example, on *NIX systems, one can
%           usually put the executable in /usr/local/bin for it to be in
%           the default search path.
%
%OUTPUTS: xSol The numVarXnumSol isolated (complex) solutions found to the
%              polynomial. The variables are in the same order as they
%              appear in varNames.
% multiplicity A numSolX1 vector of the multiplicity of each solution.
%          err A numSolX1 vector of magnitudes of last correction term used
%              in Newton's method when solving for each root.
%          rco A numSolX1 vector of the estimated inverse condition numbers
%              of the roots. This relates to numerical stability.
%          res A numSolX1 vector of magnitudes of the polynomial vector
%              evaluated at each root.
%
%This function work by writing the commands for PHCPack to a file and then
%calling PHCpack to solve the problem. The files are given filenames with
%UUIDs in the names so that multiple instances of PHCpack can be run at the
%same time if one uses a parfor loop in Matlab (assuming that the parallel
%processing toolbox is installed).
%
%Binaries and source code for PHCPack can be downloaded from
%http://homepages.math.uic.edu/~jan/download.html
%A description of the implementation of the original algorithm is given in
%[1].
%
%EXAMPLE:
% polyStrings=cell(2,1);
% polyStrings{1}='-x1^2+2*x1*x2+x2^2+5*x1-3*x2-4';
% polyStrings{2}='x1^2+2*x1*x2+x2^2-1';
% varNames={'x1','x2'};
% [xSol,multiplicity,err,rco,res]=solvePolySysUsingPHCpack(polyStrings,varNames)
% %One should get the solutions (0,-1), (1,0), (3,-2), and (4,-5).
%
%REFERENCES:
%[1] J. Verschelde, "Algorithm 795: PHCpack: A general-purpose solver for
%    polynomial systems by homotopy continuation," ACM Transactions on
%    Mathematical Software, vol. 25, no. 2, pp. 251-276, Jun. 1999.
%
%May 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3)
   opts=[]; 
end

if(nargin<4||isempty(scratchFolderPath))
    %The default is to put the temporary files into a folder called temp
    %that is located in the same folder as this script.
    ScriptPath=mfilename('fullpath');
    scratchFolderPath=[fileparts(ScriptPath),'/temp/'];
end

if(nargin<5||isempty(PHCPath))
    %The default assumes that the program is already in the user's path and
    %can be called directly.
    PHCPath='phc';
end

%A random number for the filenames. This should help avoid name conflicts
%when doing runs in parallel. We use the call to genUUID to try to make
%sure that the same number of not used if this function is called multiple
%times in parallel in a parfor loop.
randNum = genUUID();

inFileName=[scratchFolderPath,'/in',randNum,'.txt'];
outFileName=[scratchFolderPath,'/out',randNum,'.txt'];
stripFileName=[scratchFolderPath,'/strip',randNum,'.txt'];

%Create the input file.
writePolynomialForPHCPack(inFileName,polyStrings,varNames);

%Call PHCpack; this creates the output files and modified the input file.
str=[PHCPath,' -b ',opts,' "',inFileName,'" "',outFileName,'"'];
system(str);

%Strip the output into a format that is easy to read. Output is also put
%into the input file.
str = [PHCPath,' -z "',inFileName,'" "', stripFileName,'"'];
system(str);

[xSol,multiplicity,err,rco,res]=solutionsFromFilePHCPack(stripFileName,varNames);

%Get rid of the temporary files.
delete(inFileName);
delete(outFileName);
delete(stripFileName);
end

function writePolynomialForPHCPack(inFilename,polyStrings,varNames)
numVar=length(varNames);
numEq=length(polyStrings);

id = fopen(inFilename,'wt');
fprintf(id,[num2str(numEq),' ',num2str(numVar),'\n']);
for curPoly=1:numEq
    fprintf(id,[polyStrings{curPoly},';\n']);
end
fclose(id);
end

function [xSol,multiplicity,err,rco,res]=solutionsFromFilePHCPack(stripFileName,varNames)
numVar=length(varNames);
varNameLengths=zeros(numVar,1);

for curVar=1:varNameLengths
    varNameLengths(curVar)=length(varNames{curVar});
end

id = fopen(stripFileName, 'rt');
solText = fread(id,'*char')';
fclose(id);

%Replace the I for a complex number in the output with i for a complex
%number in Matlab.
sel=strfind(solText,'I');
solText(sel)='i';

startBraces=strfind(solText,'[');
endBraces=strfind(solText,']');
numSol=length(startBraces)-1;

xSol=zeros(numVar,numSol);
multiplicity=zeros(numSol,1);
err=zeros(numSol,1);
rco=zeros(numSol,1);
res=zeros(numSol,1);
for curSol=1:numSol
    idxStart=startBraces(curSol+1)+1;
    idxEnd=endBraces(curSol)-1;
    curSolText=solText(idxStart:idxEnd);
    
    startIdx=strfind(curSolText,'=');
    endIdx=strfind(curSolText,',');
    
    %The first entry is always "time". We will skip that.
    sel=(startIdx(2)+1):(endIdx(2)-1);
    %The second entry is multiplicity. Get that.
    multiplicity(curSol)=eval(curSolText(sel));
    
    %The next numVar entries are the variables. The variables are given in
    %the same order for all of the solutions. Thus, we will extract them in
    %that order and then figure out the associated names.
    for curVar=1:numVar
        sel=(startIdx(2+curVar)+1):(endIdx(2+curVar)-1);
        xSol(curVar,curSol)=eval(curSolText(sel));
    end
    %Get the final three entries.
    
    %Magnitude of last correction term used in Newton's method.
    sel=(startIdx(2+numVar+1)+1):(endIdx(2+numVar+1)-1);
    err(curSol)=eval(curSolText(sel));
    %Estimated inverse condition number of the root
    sel=(startIdx(2+numVar+2)+1):(endIdx(2+numVar+2)-1);
    rco(curSol)=eval(curSolText(sel));
    %Magnitude of the polynomial vector evaluated at the root.
    res(curSol)=eval(curSolText((startIdx(2+numVar+3)+1):end));
end

%Now, we must find the names of the variables in order. They are in the
%same order for all solutions, so we will just get them for the first
%solution.
varIdx=zeros(numVar,1);

idxStart=startBraces(1+1)+1;
idxEnd=endBraces(1)-1;
curSolText=solText(idxStart:idxEnd);%The first solution
%The variable names begin after two spaces.
spaceList=strfind(curSolText,'  ');
%The first set of two spaces is the multiplicity and the last is err. We
%will just extract the names
for curVar=1:numVar
    startIdx=spaceList(curVar+1)+2;
    endIdx=startIdx;
    while(curSolText(endIdx+1)~=' ')
        endIdx=endIdx+1;
    end
    curVarName=curSolText(startIdx:endIdx);

    %Now, find which input variable has this name.
    for k=1:numVar
        if(strcmp(varNames{k},curVarName))
            varIdx(curVar)=k;
            break;
        end
    end
end
%Reorder the estimates.
xSol=xSol(varIdx,:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%FUNCTIONS FOR MACAULAY2
%%%%%%%%%%%%%%%%%%%%%%%%%
function xSol=solvePolySysUsingMacaulay2(polyStrings,varNames,scratchFolderPath,M2Path)
%%SOLVEPOLYSYSUSINGMACAULAY2 This function uses the certified multivariate
%           polynomial rooting algorithm in the NumericalAlgebraicGeometry
%           package of Macaulay2 to solve a system of simultaneous
%           multivariate polynomials. This function just formats the input
%           properly and passes it to the M2 function.
%
%INPUTS: polyStrings A numPolyX1 or 1XnumPoly cell array containing the
%           polynomials given as strings. The variable names should begin
%           with letters but may containg numbers. Imaginary coefficients
%           should be indicated by ii. An example of a valid equation is
%           'x^2+(1+2*ii)*y^2'. Note that 'ii' is used for a complex
%           variable instead of just i.
%  varNames A numVarX1 or 1XnumVar cell array containing the names
%           used for the different variables in polyStrings. For example,
%           'x'. Generally, numVar=numPoly.
% scratchFolderPath An optional parameter specifying the folder in which
%           temporary files will be written. This is needed, because
%           Macaulay2 only works through files. If this parameter is 
%           omitted or an empty matrix is passed, then a folder named temp
%           in the folder enclosing this function is used. Note that if an
%           error occurs while executing this function, then temporary
%           files might  be left in the temp folder.
%    M2Path The command line command to use to execute Macaulay2 The
%           default if this parameter is omitted or an empty matrix is
%           passed is just M2, which assumes that the program is already
%           in the default search path. For example, on *NIX systems, one
%           can usually put the executable in /usr/local/bin for it to be
%           in the default search path. On Mac OS X, the executable is
%           often placed in a folder like
%           '/Applications/Macaulay2-1.8.2/bin/M2', which is not in the
%           default search path for the system unless one adds it.
%
%OUTPUTS: xSol The numVarXnumSol isolated (complex) solutions found to the
%              polynomial. Note that solutions marked as failing in the
%              refinement step of the algorithm are not marked as such
%              here.
%
%This function works by writing the commands to call Macaulay2 to a
%temporary file, calling the function, can then reading the output file.
%The files are given filenames with UUIDs in the names so that multiple
%instances of M2 can be run at the same time if one uses a parfor loop in
%Matlab (assuming that the parallel processing toolbox is installed).
%
%Binaries and source code for Macaulay2 may be downloaded from
%http://www.math.uiuc.edu/Macaulay2/
%The algorithm implemented in Macaulay2 that is called here is described in
%[1].
%
%EXAMPLE:
% polyStrings=cell(2,1);
% polyStrings{1}='-x1^2+2*x1*x2+x2^2+5*x1-3*x2-4';
% polyStrings{2}='x1^2+2*x1*x2+x2^2-1';
% varNames={'x1','x2'};
% xSol=solvePolySysUsingMacaulay2(polyStrings,varNames)
% %One should get the solutions (0,-1), (1,0), (3,-2), and (4,-5).
%
%REFERENCES:
%[1] C. Beltrán and A. Leykin, "Certified numerical homotopy tracking,"
%    Experimental Mathematics, vol. 21, no. 1, 2012.
%
%March 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(scratchFolderPath))
    %The default is to put the temporary files into a folder called temp
    %that is located in the same folder as this script.
    ScriptPath=mfilename('fullpath');
    scratchFolderPath=[fileparts(ScriptPath),'/temp/'];
end

if(nargin<4||isempty(M2Path))
    %The default assumes that the program is already in the user's path and
    %can be called directly.
    M2Path='M2';
end

%A random number for the filenames. This should help avoid name conflicts
%when doing runs in parallel. We use the call to genUUID to try to make
%sure that the same number of not used if this function is called multiple
%times in parallel in a parfor loop.
randNum = genUUID();

inFileName=[scratchFolderPath,'/in',randNum,'.m2'];
outFileName=[scratchFolderPath,'/out',randNum,'.m2'];

writePolynomialForMacaulay2(inFileName,outFileName,polyStrings,varNames);

str=[M2Path,' --script "',inFileName,'"'];
system(str);%Call Macaulay2

numVar=length(varNames);
xSol=solutionsFromFileMacaulay2(outFileName,numVar);

%Get rid of the temporary files.
delete(inFileName);
delete(outFileName);
end

function writePolynomialForMacaulay2(inFilename,outfileName,polyStrings,varNames)
id = fopen(inFilename,'wt');
fprintf(id,'loadPackage \"NumericalAlgebraicGeometry"\n');
fprintf(id,'R=CC[%s];\n',strjoin(varNames,', '));
fprintf(id,'F={%s};\n',strjoin(polyStrings,', '));
fprintf(id,'printingPrecision = 53;\nsols=solveSystem(F,Predictor=>Certified,Normalize=>true);\n');
fprintf(id,['"',outfileName,'"','<< sols / coordinates << close\n']);
fclose(id);
end

function xSol=solutionsFromFileMacaulay2(outfileName,numVar)
id = fopen(outfileName, 'rt');
solText = fread(id,'*char');
fclose(id);

ii=1i;%Necessary so that complex numbers are properly identified in the
%eval statement below as the output uses ii for complex numbers.
solCells=eval(solText);
numSol=length(solCells);

xSol=cell2mat(reshape([solCells{:}],numVar,numSol));
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
