function [editDist,editCount,T]=findEditDistance(string1,string2,costs)
%%FINDEDITDISTANCE Determine the distance between two strings in terms of
%               the minimum cost sequence of "edit operations". Edit
%               operations consists of changing a symbol, inserting a
%               symbol and deleting a symbol.
%
%INPUTS: string1 An N1X1 or 1XN1 array. It could consists of characters or
%                numbers. N1>=1.
%        string2 An N2X1 or 1XN2 array. It could consists of characters or
%                numbers. N2>=1.
%          costs An optional array containing the (positive) costs of the
%                edit operations. costs(1) is the cost of changing a
%                symbol. costs(2) is the cost of deleting a symbol, and 
%                costs(3) is the cost of inserting a symbol. If this
%                parameter is omitted, then a value of 1 is assigned to all
%                costs.
%
%OUTPUTS: editDist The edit distance between string1 and string2.
%        editCount A 3X1 vector where editCount(1) is the number of change
%                  operations, editCount(2) is the number of deletion
%                  options and editCount(3) is the number of insertion
%                  operations.
%                T The trace of the edit. This is a 2XnumTrace matrix
%                  containing entries [i;j]=T(:,curEntry) which if
%                  string1(i)~=string2(j) indicates a change operations.
%                  Values of i that are missing from T indicate deletion
%                  operations. Values of j that are missing from T indicate
%                  insertion operations. The insertions occur between gaps
%                  in j for (i,j) pairs. The (i,j) pairs in T occur in an
%                  order of decreasing i and j.
%
%The Wagner-Fisher algorithm of [1] for determining the edit distance
%between two striungs is used. The algorithm has quadratic complexity,
%which according to [2], is optimal unless the Strong Exponential Time
%Hypothesis (SETH) is false. The memory used by the Wagner-Fisher algorithm
%scales quadratically. Memory usage in Hirschberg's Algorithm C of [3],
%scales linearly and might be good for a future upgrade of this algorithm.
%Compare this function to the alignSequences function, which implements a
%form of Hirschberg's linear-memory algorithm and produces a similar type
%of score, but which does not provide a trace matrix.
%
%EXAMPLE:
% string1='Mary had a little lamb its fleece was white as snow';
% string2='Mary hat hey lid tells ham ids fleas were white has known';
% [editDist,editCount,T]=findEditDistance(string1,string2)
%One obtains an edit count of 21, with 11 changes, 8 deletions and 2
%insertions.
%
%The trace can be difficult to understand. To help with that, the following
%code uses strings 1 and 2 to recreate string 2 using only operations in
%the trace. These include the insertions and the change operations.
%Performing the operations by directly modifying string 1 is tedious as the
%indices in T must be modified after every insertion and deletion.
%
% numInTrace=size(T,2);
% editedString=char(zeros(size(string2)));%Allocate space
% curJ=1;
% curEdit=0;
% for curTrace=1:numInTrace
%     i=T(1,curTrace);
%     j=T(2,curTrace);
% 
%     %Perform any insertion operations
%     for insertIdx=curJ:(j-1)
%         curEdit=curEdit+1;
%         editedString(curEdit)=string2(insertIdx);
%     end
%
%     %Deletion operations occur implicitly due to i indices being omitted.
%
%     %Perform any change operations.
%     curEdit=curEdit+1;
%     if(string1(i)~=string2(j))%The character is changed.
%         editedString(curEdit)=string2(j);
%     else%The character remains unchanged.
%         editedString(curEdit)=string1(i);
%     end
%     curJ=j+1;
% end
% 
% %Insertions at the end.
% for insertIdx=curJ:length(string2)
%     curEdit=curEdit+1;
%     editedString(curEdit)=string2(insertIdx);
% end
% editedString
%The editedString is now string2.
%
%REFERENCES:
%[1] R. A. Wagner and M. J. Fischer, "The string-to-string correction
%    problem," Journal of the Association for Computing Machinery, vol. 21,
%    no. 1, pp. 168-173, Jan. 1974.
%[2] A. Backurs and P. Indyk, "Edit distance cannot be computed in
%    strongly subquadratic time (unless SETH is false)," arXiv, 13 Apr.
%    2015. [Online]. Available: http://arxiv.org/abs/1412.0348
%[3] D. S. Hirschberg, "A linear space algorithm for computing maximal
%    common subsequences," Communications of the ACM, vol. 18, no. 6, pp.
%    341-343, Jun. 1975.
%
%July 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

N1=length(string1);
N2=length(string2);

d=zeros(N1+1,N2+1);

if(nargin<3||isempty(costs))
    costs(1,1)=1;%Cost of changing a symbol
    costs(2,1)=1;%Cost of deleting a symbol
    costs(3,1)=1;%Cost of inserting a symbol
end

%Initialize with costs of inserting symbols.
d(0+1,(0:N2)+1)=(0:N2)*costs(3);
%initialize with costs of deleting symbols.
d((0:N1)+1,0+1)=(0:N1)*costs(2);
    
for i=1:N1
    for j=1:N2
        %Cost of changing a symbol.
        m1=d(i-1+1,j-1+1)+(string1(i)~=string2(j))*costs(1);
        %Cost of inserting a symbol.
        m2=d(i-1+1,j+1)+costs(3);
        %Cost of deleting a symbol.
        m3=d(i+1,j-1+1)+costs(2);

        d(i+1,j+1) = min([m1,m2,m3]);
    end
end
editDist=d(end,end);

%If the number of each operation and/ or the trace of the edits is desired.
if(nargout>1)
    T=zeros(2,max(N1,N2));
    editCount=zeros(3,1);
    
    i=N1;
    j=N2;
    curEntry=0;
    while(i~=0&&j~=0)
        if(d(i+1,j+1)==d(i-1+1,j+1)+costs(3))%Symbol insertion
            editCount(3)=editCount(3)+1;
            
            i=i-1;
        elseif(d(i+1,j+1)==d(i+1,j-1+1)+costs(2))%Symbol deletion
            editCount(2)=editCount(2)+1;
            j=j-1;
        else
           curEntry=curEntry+1;           
           if(string1(i)~=string2(j))
               editCount(1)=editCount(1)+1;
           end
           
           T(:,curEntry)=[i;j];
           
           i=i-1;
           j=j-1;
        end
    end
    
    %Shrink the T to the actual length of the trace.
    T=T(:,1:curEntry);
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
