%before running this script, you must have run the
%CM_socialnetwork_analysis.R script in R and then saved comm4 or commL4 as
%a csv
%write.csv(comm4,file='sorted_communities_4.csv')
%then you need to sort by N number and remove headers and the first column
%you will also have to replace all 'N's with ''s
cd('F:\0_parcellation_analysis\scripts-data-sharing\')
comm=dlmread('sorted_communities_L4.csv');
x=1:size(comm,1);
y=1:size(comm,1);
[xx,yy]=meshgrid(x,y);
module_similarity=[xx(:),yy(:)];
equality=nan(size(module_similarity,1),1);
count=0;
for ii=1:size(comm,1)
    for jj=1:size(comm,1)
        count=count+1;
        subj1=comm(ii,2);
        subj2=comm(jj,2);
        equality(count,1)=isequal(subj1,subj2);
    end
end
    output=cat(2,module_similarity,equality);
%     indices=find(output(:,2)==output(:,1));
%     output(indices,:)=[];
    
vec_comm_similarity=reshape(output(:,3),size(comm,1),[]);
vec_comm_similarity=vec_comm_similarity.*~eye(size(vec_comm_similarity));
[con1,h]=contourf(vec_comm_similarity,2);
h.LineStyle='none';
csvwrite('communities_L4.csv',vec_comm_similarity);


comm=dlmread('sorted_communities_4.csv');
x=1:size(comm,1);
y=1:size(comm,1);
[xx,yy]=meshgrid(x,y);
module_similarity=[xx(:),yy(:)];
equality=nan(size(module_similarity,1),1);
count=0;
for ii=1:size(comm,1)
    for jj=1:size(comm,1)
        count=count+1;
        subj1=comm(ii,2);
        subj2=comm(jj,2);
        equality(count,1)=isequal(subj1,subj2);
    end
end
    output=cat(2,module_similarity,equality);
%     indices=find(output(:,2)==output(:,1));
%     output(indices,:)=[];
    
vec_comm_similarity=reshape(output(:,3),size(comm,1),[]);
vec_comm_similarity=vec_comm_similarity.*~eye(size(vec_comm_similarity));
figure
[con2,h]=contourf(vec_comm_similarity,2);
h.LineStyle='none';
csvwrite('communities_4.csv',vec_comm_similarity);

comm=dlmread('sorted_communities_700.csv');
x=1:size(comm,1);
y=1:size(comm,1);
[xx,yy]=meshgrid(x,y);
module_similarity=[xx(:),yy(:)];
equality=nan(size(module_similarity,1),1);
count=0;
for ii=1:size(comm,1)
    for jj=1:size(comm,1)
        count=count+1;
        subj1=comm(ii,2);
        subj2=comm(jj,2);
        equality(count,1)=isequal(subj1,subj2);
    end
end
    output=cat(2,module_similarity,equality);
%     indices=find(output(:,2)==output(:,1));
%     output(indices,:)=[];
    
vec_comm_similarity=reshape(output(:,3),size(comm,1),[]);
vec_comm_similarity=vec_comm_similarity.*~eye(size(vec_comm_similarity));
figure
[con2,h]=contourf(vec_comm_similarity,2);
h.LineStyle='none';
csvwrite('communities_700s.csv',vec_comm_similarity);

