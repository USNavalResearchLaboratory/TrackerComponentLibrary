classdef ClusterSet < handle
%%CLUSTERSET A class to access a set of clusters by indexation. Given an
%            object CSObs, one can index the elements by cluster and by the
%            index of the item in the cluster as CSObj(curClust,curItem) as
%            described in the methods subsref and subsasgn, which overload
%            basic operators.
%
%Given a two-dimensional matrix A in Matlab, one can access the element in
%row r and column c as A(r,c). The ClusterSet class allows the analogous
%access of a set of data where each row can have a different number of
%columns. The data for all of the rows is stored in memory as one
%contiguous chunk, but is accessed in the same manner A(r,c) for element c
%in cluster r if A is an instance of the ClusterSet class. The insertion
%and deletion of data are not supported. However, one can overwrite values.
%
%The ClusterSet class is used, for example, to store coefficients of
%varying degrees and order for spherical harmonic synthesis. As the degree
%increases, the maximum order also increases. That is, the number of
%elements for the degree increases. Normally, if one wanted to access the
%elements with the matrix notation A(degree, order), one would need to
%allocate a square matrix of size (maximum degree X maximum degree) and
%about half of the elements would not be used. By using a ClusterSet class,
%there are no wasted elements.
%
%The ClusterSet class inherits from the handle class, meaning that any
%function modifying the data in the ClusterSet class modifies the data
%everywhere. The duplicate method of the class can be used to make a new
%copy of an instance of the ClusterSet class to work on, if a function
%wants to change the values of the elements in the class without modifying
%the original class.
%
%The ClusterSet class has a parallel C++ class named ClusterSetCPP. The C++
%class overloads the [] operators to provide functionality similar to what
%the Matlab ClusterSet class provides. However, the C++ version is
%significantly faster. The Matlab implementation of the ClusterSet class
%provides a useful data type to pass to mex files so that pointers to the
%data can be extracted and inserted into a ClusterSetCPP class.
%
%November 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    properties
        clusterEls
        clusterSizes
        offsetArray
    end

    methods
        function newClust=ClusterSet(clusterEls,clusterSizes,offsetArray)
        %%CLUSTERSET The constructor for a ClusterSet class. This is called
        %            by a DisjointSet class when creating a ClusterSet
        %            instance from the data in the disjoint set.
        %
        %November 2013 David F. Crouse, Naval Research Laboratory,
        %Washington D.C.
        
            if(nargin==0)
            %The default constructor creates an empty instance of the
            %ClusterSet class.
                newClust.clusterEls=0;
                newClust.clusterSizes=[];
                newClust.offsetArray=0;
            elseif(nargin==1)
                newClust.clusterEls=clusterEls;
                newClust.clusterSizes=length(clusterEls);
                newClust.offsetArray=0;
            elseif(nargin==2)
                newClust.clusterEls=clusterEls;
                newClust.clusterSizes=clusterSizes;
                if(size(clusterSizes,1)==1)
                    newClust.offsetArray=cumsum([0;clusterSizes(1:(end-1))']);
                else
                    newClust.offsetArray=cumsum([0;clusterSizes(1:(end-1))]);
                end
            else
                newClust.clusterEls=clusterEls;
                newClust.clusterSizes=clusterSizes;
                newClust.offsetArray=offsetArray;
            end
        end
        
        function val=getEl(CSObj,clusterIdx,itemIdx)
        %%GETEL This is the same as calling the subsref function, which has
        %       the syntax theSet(clusterIdx,itemIdx), to get a signle
        %       element. This function has the syntax
        %       theSet.getEl(clusterIdx,itemIdx). The advantage to using
        %       this function is that it is faster. The subsref function
        %       checks the type of the input and looks for parentheses and
        %       stuff. This one just directly accesses the element and
        %       doesn't do type checking on the input. The indices can't be
        %       different data types.
        %
        %December 2020 David F. Crouse, Naval Research Laboratory,
        %Washington D.C.
            
            val=CSObj.clusterEls(CSObj.offsetArray(clusterIdx)+itemIdx);
        end
        
        function val=subsref(CSObj,S)
        %%SUBSREF Access the elements of a set of clusters of indices by
        %        cluster number and by the number of the item in the
        %        cluster. The syntax is CSObj(clusterIdx,itemIdx); One can
        %        use CSObj(clusterIdx,:) to return all elements of a
        %        cluster, or CSObj(clusterIdx,a:b) to return the elements
        %        between a and b, but the indexation CSObj(:,itemIdx) has
        %        no meaning. Using CSObj(:) will return all of the elements
        %        of all of the clusters as a single array. Using
        %        CSObj(clusterIdx) will return all of the elements in a
        %        particular cluster. CSObj(:,:) will generate an error.
        %
        % Indexation is as for arrays in Matlab. Note that this does not
        % check whether a zero index is given. In some instances, passing a
        % zero index will not result in an error, but will erroneously
        % return a value. For example, CSObj(2,0) might return the last
        % element in cluster 1 is cluster 1 if not empty. Code could be
        % added to check for that, but it would make this class less
        % efficient.
        %
        %November 2013 David F. Crouse, Naval Research Laboratory,
        %Washington D.C.
            
            %If it is not just indexation, then pass it along to call the
            %appropriate class method.
            if(length({S.type})~=1||strcmp(S.type,'.')==true)
                val = builtin('subsref',CSObj,S);
                return
            end
            if(strcmp(S.type,'()')==true)
                switch(length(S.subs))
                    case 1
                        if(islogical(S.subs{1}))
                            if(length(S.subs{1})~=length(CSObj.clusterEls))
                                error('Invalid Indexation Given')
                            end
                    
                            val=CSObj.clusterEls(S.subs{1});
                        else
                            if(strcmp(S.subs{1},':')==true)
                                val=CSObj.clusterEls;
                                return
                            elseif(isnumeric(S.subs{1}))
                                if(isscalar(S.subs{1}))%Return all of the elements in the selected cluster.

                                    idx=CSObj.offsetArray(S.subs{1})+(1:CSObj.clusterSizes(S.subs{1}));
                                    val=CSObj.clusterEls(idx);
                                    return
                                else%Return elements from multiple clusters.
                                    elementCount=sum(CSObj.clusterSizes(S.subs{1}));
                                    numClust=length(S.subs{1});

                                    val=zeros(elementCount,1);
                                    valStartIdx=1;
                                    for clustIdx=1:numClust
                                        curClust=S.subs{1}(clustIdx);
                                        clustStart=1+CSObj.offsetArray(curClust);
                                        clustEnd=clustStart+CSObj.clusterSizes(curClust)-1;

                                        valEndIdx=valStartIdx+CSObj.clusterSizes(curClust)-1;
                                        val(valStartIdx:valEndIdx)=CSObj.clusterEls(clustStart:clustEnd);
                                        valStartIdx=valEndIdx+1;
                                    end
                                    return
                                end
                            else
                                error('Invalid Indexation Given')
                            end
                        end
                    case 2%Two indices are given
                        if(~isnumeric(S.subs{1})||~isscalar(S.subs{1}))
                            error('Invalid Indexation Given')
                        end

                        if(isnumeric(S.subs{2}))
                            if(S.subs{2}>CSObj.clusterSizes(S.subs{1}))
                                error('Attempt to read past the end of a cluster.')
                            end
                            
                            %The double conversions are necessary to avoid
                            %errors if, for example, a subscript vector
                            %passed is doubles but the offsets are an
                            %integer type.
                            idx=double(CSObj.offsetArray(S.subs{1}))+double(S.subs{2});
                            val=CSObj.clusterEls(idx);
                        else
                            if(strcmp(S.subs{2},':')==true)%Return all of the elements in the selected cluster.
                                idx=CSObj.offsetArray(S.subs{1})+(1:CSObj.clusterSizes(S.subs{1}));
                                val=CSObj.clusterEls(idx);
                            else
                                error('Invalid Indexation Given')
                            end
                        end
                        return
                    otherwise
                        error('Invalid Indexation Given')
                end
            else
                error('Invalid Indexation Given')
            end
        end
        
        function CSObj=subsasgn(CSObj,S,x)
        %SUBSASGN Set the elements of a set of clusters of indices by
        %         cluster number and by the number of the item in the
        %         cluster. The syntax is CSObj(clusterIdx,itemIdx)=x;
        %
        %The same indexation issues exist as for the subsref method.
        %
        %November 2013 David F. Crouse, Naval Research Laboratory,
        %Washington D.C.

            %If it is not just indexation, then pass it along to call the
            %appropriate class method.
            if(length({S.type})~=1||strcmp(S.type,'.')==true)
                CSObj=builtin('subsasgn',CSObj,S,x);
                return
            end

            if(strcmp(S.type,'()')==true)
                switch(length(S.subs))
                    case 1
                        if(islogical(S.subs{1}))
                            if(length(S.subs{1})~=length(CSObj.clusterEls))
                                error('Invalid Indexation Given')
                            end
                    
                            CSObj.clusterEls(S.subs{1})=x;
                            return
                        elseif(strcmp(S.subs{1},':')==true)
                            CSObj.clusterEls(:)=x;
                            return
                        elseif(isnumeric(S.subs{1}))
                            if(isscalar(S.subs{1}))%Set all of the elements in the selected cluster.
                                idx=CSObj.offsetArray(S.subs{1})+(1:CSObj.clusterSizes(S.subs{1}));
                                CSObj.clusterEls(idx)=x;
                                return
                            else%Setting elements in multiple clusters is not supported.
                                error('Invalid Indexation Given')
                            end
                        end
                        return
                    case 2
                        if(isnumeric(S.subs{2}))
                            if(S.subs{2}>CSObj.clusterSizes(S.subs{1}))
                                error('Attempt to read past the end of a cluster.')
                            end
                            
                            %The double conversions are necessary to avoid
                            %errors if, for example, a subscript vector
                            %passed is doubles but the offsets are an
                            %integer type.
                            idx=double(CSObj.offsetArray(S.subs{1}))+double(S.subs{2});
                            CSObj.clusterEls(idx)=x;
                        else
                            if(strcmp(S.subs{2},':')==true)%Set all of the elements in the selected cluster.
                                idx=CSObj.offsetArray(S.subs{1})+(1:CSObj.clusterSizes(S.subs{1}));
                                CSObj.clusterEls(idx)=x;
                            else
                                error('Invalid Indexation Given')
                            end
                        end
                        return
                    otherwise
                        error('Invalid Indexation Given')
                end
            else
                error('Invalid Indexation Given')
            end
        end
        
        function val=numClusters(CSObj)
        %NUMCLUSTERS Determine the number of clusters present.
        %
        %November 2013 David F. Crouse, Naval Research Laboratory,
        %Washington D.C.
        
            val=length(CSObj.clusterSizes);
        end
        
        function newClust=duplicate(CSObj)
        %DUPLICATE Make a copy of the current cluster set. Since the
        %          ClusterSet class is a subclass of the handle class,
        %          passing a ClusterSet object to a function only passes a
        %          pointer. If the function wishes to modify the class
        %          without changing the original, then it must duplicate
        %          the original class.
        %
        %November 2013 David F. Crouse, Naval Research Laboratory,
        %Washington D.C.
        
            newClust=ClusterSet(CSObj.clusterEls,CSObj.clusterSizes,CSObj.offsetArray);
        end
        
        function val=clustSize(CSObj,clustIdx)
        %CLUSTSIZE Determine the size of the cluster of index clustIdx. If
        %          no index is specified, then all of the cluster sizes are
        %          returned.
        %
        %November 2013 David F. Crouse, Naval Research Laboratory,
        %Washington D.C.
        
            if(nargin==1)
                val=CSObj.clusterSizes;
            else
                val=CSObj.clusterSizes(clustIdx);
            end
        end
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
