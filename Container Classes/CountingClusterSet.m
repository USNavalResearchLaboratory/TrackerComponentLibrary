classdef CountingClusterSet < handle
%%COUNTINGCLUSTERSET A class to access a set of clusters of elements by
%           indexation. The number of items in the clusters increases
%           linearly. Cluster 1 has 1 item, 2 has 2 items, etc. Given an
%           object CSObs, one can index the elements by cluster and by the
%           index of the item in the cluster as CSObj(curClust,curItem) as
%           described in the methods subsref and subsasgn, which overload
%           basic operators. A more version of this class that can handle
%           different sizes of clusters is the ClusterSet class. This class
%           can handle multiple sets of data at once.
%
%Given a two-dimensional matrix A in Matlab, one can access the element in
%row r and column c as A(r,c). The CountingClusterSet class allows the
%analogous access of a set of data where each row has the same number of
%columns as the row number. The data for all of the rows is stored in
%memory as one contiguous chunk, but is accessed in the same manner A(r,c)
%for element c in cluster r if A is an instance of the ClusterSet class.
%The insertion and deletion of data are not supported. However, one can
%overwrite values.
%
%It is possible to have multiple sets of data at once. In this instance,
%A(r,c) returns a 1XnumSets vector of points. Similarly, A(r,:) returns a
%numElsXnumSets matrix of values.
%
%The CountingClusterSet class inherits from the handle class, meaning that
%any function modifying the data in the CountingClusterSet class modifies
%the data everywhere. The duplicate method of the class can be used to make
%a new copy of an instance of the CountingClusterSet class to work on, if a
%function wants to change the values of the elements in the class without
%modifying the original class.
%
%July 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.   

properties
    clusterEls
    numClust
end

methods
    function newClust=CountingClusterSet(clusterEls)
    %%COUNTINGCLUSTERSET The constructor for a CountingClusterSet
    %        class. This should generally be handed all of the data
    %        (or a zero vector to reserve space) when allocated. The
    %        data is ordered by rows, so all elements for the first
    %        row, followed by those for the second row, etc. To have
    %        have M rows, clusterEls must have M*(M+1)/2 rows.
    %
    %INPUTS: clusterEls A numElsXnumSets collection of the elements of
    %                   numSets sets of data.
    %
    %OUTPUTS: newClust The newly created CountingClusterSet object.
    %
    %July 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
    
        if(nargin>0)
            newClust.clusterEls=clusterEls; 
            totalEls=size(clusterEls,1);

            MTest=(1/2)*(sqrt(1+8*totalEls)-1);
            if(MTest~=fix(MTest))
                error('The length of clusterEls must be of the form M*(M+1)/2 with M being a positive integer.');
            end
            newClust.numClust=MTest;
        else
            newClust.clusterEls=[];
            newClust.numClust=0;
        end
    end

    function val=subsref(CSObj,S)
    %SUBSREF Access the elements of a set of clusters of indices by
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
    %Indexation is as for arrays in Matlab, except if multiple sets of
    %data are passed, then multiple columns of values will be returned.
    %
    %July 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
    
        %If it is not just indexation, then pass it along to call the
        %appropriate class method.
        if(length({S.type})~=1||strcmp(S.type,'.')==true)
            val=builtin('subsref',CSObj,S);
            return
        end

        if(strcmp(S.type,'()')==true)
            switch(length(S.subs))
                case 1
                    if(islogical(S.subs{1}))
                        if(length(S.subs{1})~=length(CSObj.clusterEls))
                            error('Invalid Indexation Given')
                        end

                        val=CSObj.clusterEls(S.subs{1},:);
                    else
                        if(strcmp(S.subs{1},':')==true)
                            val=CSObj.clusterEls;
                            return
                        elseif(isnumeric(S.subs{1}))
                            if(isscalar(S.subs{1}))
                            %Return all of the elements in the selected
                            %cluster.

                                clustStartIdx=(1/2)*S.subs{1}*(S.subs{1}-1)+1;
                                clustEndIdx=clustStartIdx+S.subs{1}-1;
                                idx=clustStartIdx:clustEndIdx;

                                val=CSObj.clusterEls(idx,:);
                                return
                            else%Return elements from multiple
                                %clusters.
                                elementCount=sum(S.subs{1});
                                numIdx=length(S.subs{1});

                                val=zeros(elementCount,1);
                                valStartIdx=1;
                                for clustIdx=1:numIdx
                                    curClust=S.subs{1}(clustIdx);
                                    clustStartIdx=(1/2)*curClust*(curClust-1)+1;
                                    clustEndIdx=clustStartIdx+curClust-1;
                                    idx=clustStartIdx:clustEndIdx;

                                    valEndIdx=valStartIdx+curClust-1;
                                    val(valStartIdx:valEndIdx)=CSObj.clusterEls(idx,:);
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
                        if(S.subs{2}>S.subs{1})
                            error('Attempt to read past the end of a cluster.')
                        end

                        idx=(1/2)*double(S.subs{1}*(S.subs{1}-1))+double(S.subs{2});
                        val=CSObj.clusterEls(idx,:);
                    else
                        if(strcmp(S.subs{2},':')==true)%Return all of the elements in the selected cluster.
                            clustStartIdx=(1/2)*S.subs{1}*(S.subs{1}-1)+1;
                            clustEndIdx=clustStartIdx+S.subs{1}-1;
                            idx=clustStartIdx:clustEndIdx;

                            val=CSObj.clusterEls(idx,:);
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
    %July 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
    
        %If it is not just indexation, then pass it along to call the
        %appropriate class method.
        if(length({S.type})~=1||strcmp(S.type,'.')==true)
            CSObj = builtin('subsasgn',CSObj,S,x);
            return
        end

        if(strcmp(S.type,'()')==true)
            switch(length(S.subs))
                case 1
                    if(islogical(S.subs{1}))
                        if(length(S.subs{1})~=size(CSObj.clusterEls,1))
                            error('Invalid Indexation Given')
                        end

                        CSObj.clusterEls(S.subs{1},:)=x;
                        return
                    elseif(strcmp(S.subs{1},':')==true)
                        CSObj.clusterEls(:)=x;
                        return
                    elseif(isnumeric(S.subs{1}))
                        if(isscalar(S.subs{1}))%Set all of the elements in the selected cluster.
                            clustStartIdx=(1/2)*S.subs{1}*(S.subs{1}-1)+1;
                            clustEndIdx=clustStartIdx+S.subs{1}-1;
                            idx=clustStartIdx:clustEndIdx;

                            CSObj.clusterEls(idx,:)=x;
                            return
                        else%Setting elements in multiple clusters is not supported.
                            error('Invalid Indexation Given')
                        end
                    end
                    return
                case 2
                    if(isnumeric(S.subs{2}))
                        if(S.subs{2}>S.subs{1})
                            error('Attempt to read past the end of a cluster.')
                        end

                        %The double conversions are necessary to avoid
                        %errors if, for example, a subscript vector
                        %passed is doubles but the offsets are an
                        %integer type.
                        idx=(1/2)*double(S.subs{1}*(S.subs{1}-1))+double(S.subs{2});
                        CSObj.clusterEls(idx,:)=x;
                    else
                        if(strcmp(S.subs{2},':')==true)%Set all of the elements in the selected cluster.
                            clustStartIdx=(1/2)*S.subs{1}*(S.subs{1}-1)+1;
                            clustEndIdx=clustStartIdx+S.subs{1}-1;
                            idx=clustStartIdx:clustEndIdx;
                            CSObj.clusterEls(idx,:)=x;
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
    %July 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
    
        val=CSObj.numClust;
    end

    function newClust=duplicate(CSObj)
    %DUPLICATE Make a copy of the current cluster set. Since the
    %          CountingClusterSet class is a subclass of the handle
    %          class, passing a CountingClusterSet object to a function
    %          only passes a pointer. If the function wishes to modify
    %          the class without changing the original, then it must
    %          duplicate the original class.
    %
    %July 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
    
        newClust=CountingClusterSet(CSObj.clusterEls);
    end

    function val=clustSize(CSObj,clustIdx)
    %CLUSTSIZE Determine the size of the cluster of index clustIdx. If
    %          no index is specified, then all of the cluster sizes are
    %          returned.
    %
    %July 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
    
        if(nargin==1)
            val=1:CSObj.numClust;
        else
            val=clustIdx;
        end
    end
    
    function val=numSets(CSObj)
    %NUMSETS Return the number of sets of data that are present.
    %
    %July 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.

        val=size(CSObj.clusterEls,2);
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
