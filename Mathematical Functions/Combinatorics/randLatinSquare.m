function LatinSquare=randLatinSquare(n,method,numProperMoves)
%%RANDLATINSQUARE Generate a random Latin square. A Latin square is a
%          square where each row and each column contains the number 1:n
%          and no number is repeated in any row or column.
%
%INPUTS: n The size of the Latin square.
%   method An optional parameter specifying the approach for generating
%          random Latin squares. Possible values are:
%          0 (The default if omitted or an empty matrix is passed) Use the
%            algorithm of [1], implemented in a manner similar to the
%            "clear but not inefficient" way described in [2]. This
%            starts with a deterministic Latin square and then makes a
%            series of modifications going down a Markov chain of Latin
%            squares. The stationary distribution of the chain is uniform.
%          1 This is a simple heuristic approach that produces random Latin
%            squares within a fixed isotopy class. The heuristic is
%            Permute the symbols 1:n, place them into a circulant matrix
%            via circulantMatrix (circulant matrices are Latin squares),
%            permute the rows, then permute the columns of the matrix.
%            For n=2,3 this will produce uniformly distributed random Latin
%            squares. However, for n=4 and higher, there exists more than
%            one isotopy class, so this will exclude an increasingly large
%            portion of the possible Latin squares. This type of random
%            generation is very fast and can be sufficient for certain
%            applications where true randomness is not needed.
% numProperMoves This value is only used if method 0 is chosen. This is the
%            number of proper moves to go down the Markov chain before
%            declaring the result sufficiently random. The default if
%            omitted or an empty matrix is passed is n^3.
%
%OUTPUTS: LatinSquare An nXn random Latin square.
%
%REFERENCES:
%[1] M. T. Jacobson and P. Matthews, "Generating uniformly distributed
%    random latin squares," Journal of Combinatorial Designs, vol. 4, no.
%    6, pp. 405-437, 1996.
%[2] G. Sagastume, Ignacio, "Generation of random latin squares step by
%    step graphically," in XX Congreso Argentino de Ciencias de la
%    Computación, Buenos Aires, 20-24 Oct. 2014.
%
%October 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(method))
    method=0;
end

if(n==1)
    LatinSquare=1;
    return;
end

switch(method)
    case 0
        if(nargin<3||isempty(numProperMoves))
            numIter=n^3;
        else
            numIter=numProperMoves;
        end

        IncMat=zeros(n,n,n);
        %Initialize the incidence matrix of the cube as a circulant matrix.
        for i=0:(n-1)
            for j=0:(n-1)
                k=mod(i+j,n);
                IncMat(i+1,j+1,k+1)=1;
            end
        end

        improperCell=[];

        curProper=0;
        while(1)
            if(isempty(improperCell))
                %If the latin square is proper.

                %If we have visited enough proper latin squares on the random walk.
                curProper=curProper+1;
                if(curProper==numIter)
                    break;
                end

                %First, select a zero cell in IncMat at random.
                iR=randi(n);
                jR=randi(n);
                kR=randi(n);
                while(IncMat(iR,jR,kR)~=0)
                    iR=randi(n);
                    jR=randi(n);
                    kR=randi(n);
                end

                %Find the coordinate in i that has the 1 for a fixed jR and kR.
                i1=1;
                while(IncMat(i1,jR,kR)~=1)
                    i1=i1+1; 
                end

                %Find the coordinate in j that has the 1 for a fixed iR and kR.
                j1=1;
                while(IncMat(iR,j1,kR)~=1)
                    j1=j1+1; 
                end

                %Find the coordinate in k that has the 1 for a fixed iR and jR.
                k1=1;
                while(IncMat(iR,jR,k1)~=1)
                    k1=k1+1; 
                end

                %Do the move.
                IncMat(iR,jR,kR)=IncMat(iR,jR,kR)+1;
                IncMat(iR,j1,k1)=IncMat(iR,j1,k1)+1;
                IncMat(i1,j1,kR)=IncMat(i1,j1,kR)+1;
                IncMat(i1,jR,k1)=IncMat(i1,jR,k1)+1;
                IncMat(iR,jR,k1)=IncMat(iR,jR,k1)-1;
                IncMat(iR,j1,kR)=IncMat(iR,j1,kR)-1;
                IncMat(i1,jR,kR)=IncMat(i1,jR,kR)-1;
                %This decriment might make the value 0 or -1.
                IncMat(i1,j1,k1)=IncMat(i1,j1,k1)-1;

                %Determine whether the incidence matrix has been made improper.
                %Mark the -1 cell if the latin square is not proper.
                if(IncMat(i1,j1,k1)==-1)
                    improperCell=[i1;j1;k1];
                end
            else%If the latin square is not proper.
                %The coordinate sof the improper cell.
                iI=improperCell(1);
                jI=improperCell(2);
                kI=improperCell(3);

                %Decide whether we take the first or the second occurrence of a 1
                %in the i coordinate when fixing jI and kI. Then, find the
                %coordinate.
                num2Take=randi(2);
                curSeen=0;
                i1=0;
                while(curSeen<num2Take)
                    i1=i1+1;
                    if(IncMat(i1,jI,kI)==1)
                        curSeen=curSeen+1;
                    end
                end

                %Decide whether we take the first or the second occurrence of a 1
                %in the j coordinate when fixing iI and kI. Then, find the
                %coordinate.
                num2Take=randi(2);
                curSeen=0;
                j1=0;
                while(curSeen<num2Take)
                    j1=j1+1;
                    if(IncMat(iI,j1,kI)==1)
                        curSeen=curSeen+1;
                    end
                end

                %Decide whether we take the first or the second occurrence of a 1
                %in the k coordinate when fixing iI and jI. Then, find the
                %coordinate.
                num2Take=randi(2);
                curSeen=0;
                k1=0;
                while(curSeen<num2Take)
                    k1=k1+1;
                    if(IncMat(iI,jI,k1)==1)
                        curSeen=curSeen+1;
                    end
                end

                %Do the move.
                IncMat(iI,jI,kI)=IncMat(iI,jI,kI)+1;
                IncMat(iI,j1,k1)=IncMat(iI,j1,k1)+1;
                IncMat(i1,j1,kI)=IncMat(i1,j1,kI)+1;
                IncMat(i1,jI,k1)=IncMat(i1,jI,k1)+1;
                IncMat(iI,jI,k1)=IncMat(iI,jI,k1)-1;
                IncMat(iI,j1,kI)=IncMat(iI,j1,kI)-1;
                IncMat(i1,jI,kI)=IncMat(i1,jI,kI)-1;
                %This decriment might make the value 0 or -1.
                IncMat(i1,j1,k1)=IncMat(i1,j1,k1)-1;

                if(IncMat(i1,j1,k1)==-1)
                    improperCell=[i1;j1;k1];
                else
                    improperCell=[];
                end
            end
        end

        %Convert the incidence matrix into a Latin square
        LatinSquare=zeros(n,n);
        for i=1:n
           for j=1:n
               for k=1:n
                   if(IncMat(i,j,k)==1)
                       LatinSquare(i,j)=k;
                   end
               end
           end
        end
    case 1%The heuristic that limits things to one isotopy class.
        letters=1:n;
        letters=letters(randperm(n));
        
        %Create a cyclic latin square.
        LatinSquare=circulantMatrix(letters);
        %Apple a random permutation to the rows.
        LatinSquare=LatinSquare(randperm(n),:);
        %Apply a random permutation to the columns.
        LatinSquare=LatinSquare(:,randperm(n));
    otherwise
        error('Unknown method specified.')
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
