%Carolyn McNabb - 2019
%This script takes the timeseries from each of 272 parcels and correlates it with every
%other parcel. then compares correlation matrices between participants and
%writes to a csv file called corr_all
%It will also use the brain connectivity toolbox 
% https://sites.google.com/site/bctnet/
%(which you will need to install) to get connectomics measurements

%add BCT to path before running this script
addpath(genpath('F:\0_parcellation_analysis\scripts-data-sharing'))
addpath(genpath('F:\BCT\'))
cd F:\0_parcellation_analysis\scripts-data-sharing
Ns={201	204	210	211	213	216	217	218	219	220	222	223	224	226	227	232	233	240	244	245	246	247	256	300	302	305	307	308	309	312	315	317	318	321	322	325	332	333	334	335	340	344	346	347	350	352	355	361	366	367	369 702 704 706 708 709 711 713 715 720 722 737 739 742 743 745 748 751};
all_corr=zeros(length(Ns),272^2);
all_zs=zeros(length(Ns),272^2);
pos_strength=zeros(length(Ns),272);
neg_strength=zeros(length(Ns),272);
pos_total=zeros(length(Ns),1);
neg_total=zeros(length(Ns),1);
for i=1:length(Ns)
    ppt=(Ns{i});
    name=sprintf('timeseries_rest_%d.txt', ppt);
    ts=dlmread(name);
    corr_ts=corr(ts');
    corr_ts=corr_ts.*~eye(size(corr_ts));
    [Spos, Sneg, Wpos, Wneg]=strengths_und_sign(corr_ts);
    pos_strength(i,:)=Spos;
    neg_strength(i,:)=Sneg;
    pos_total(i,:)=Wpos;
    neg_total(i,:)=Wneg;
    eval(['corr_ts' num2str(ppt) '=corr_ts']);
    vec=reshape(corr_ts',1,[]);
    eval(['vec_' num2str(ppt) '=vec']);
    z=atanh(vec);%RtoZ transformation
    all_corr(i,:)=vec;
    all_zs(i,:)=z;
end
% correlations=corr(all_corr');
correlations=corr(all_zs');
Unique_corr=tril(correlations);
Unique_corr=Unique_corr.*~eye(size(Unique_corr));
csvwrite('F:\0_parcellation_analysis\scripts-data-sharing\corr_all.csv',Unique_corr);
    
correlations=corr(neg_strength');
Unique_corr=tril(correlations);
Unique_corr=Unique_corr.*~eye(size(Unique_corr));
csvwrite('F:\0_parcellation_analysis\scripts-data-sharing\neg_strength.csv',Unique_corr);

correlations=corr(pos_strength');
Unique_corr=tril(correlations);
Unique_corr=Unique_corr.*~eye(size(Unique_corr));
csvwrite('F:\0_parcellation_analysis\scripts-data-sharing\pos_strength.csv',Unique_corr);


%% part II - total positive strength similarities between participants

total_str_L4=nan(23,2);
total_str_4=nan(28,2);
total_str_700=nan(17,2);
for i=1:length(total_str_L4)
total_str_L4(i,1)=(Ns{i});
end
for i=length(total_str_L4)+1:length(pos_total)-length(total_str_700)
total_str_4(i-length(total_str_L4),1)=(Ns{i});
end
for i=length(total_str_L4)+length(total_str_4)+1:length(pos_total)
total_str_700(i-(length(total_str_L4)+length(total_str_4)),1)=(Ns{i});
end
total_str_L4(:,2)=pos_total(1:length(total_str_L4),1);
total_str_4(:,2)=pos_total(length(total_str_L4)+1:length(pos_total)-length(total_str_700),1);
total_str_700(:,2)=pos_total(length(total_str_L4)++length(total_str_4)+1:length(pos_total),1);

%difference between individuals in each dyad L4
x=1:size(total_str_L4,1);
y=1:size(total_str_L4,1);
[xx,yy]=meshgrid(x,y);
grid_L4=[xx(:),yy(:)];
difference_L4=nan(size(grid_L4,1),1);
count=0;
for ii=1:size(total_str_L4,1)
    for jj=1:size(total_str_L4,1)
        count=count+1;
        subj1=total_str_L4(ii,2);
        subj2=total_str_L4(jj,2);
        difference_L4(count,1)=abs(subj1-subj2);
    end
end
output_L4=cat(2,grid_L4,difference_L4);

%difference between individuals in each dyad year 4
x=1:size(total_str_4,1);
y=1:size(total_str_4,1);
[xx,yy]=meshgrid(x,y);
grid_4=[xx(:),yy(:)];
difference_4=nan(size(grid_4,1),1);
count=0;
for ii=1:size(total_str_4,1)
    for jj=1:size(total_str_4,1)
        count=count+1;
        subj1=total_str_4(ii,2);
        subj2=total_str_4(jj,2);
        difference_4(count,1)=abs(subj1-subj2);
    end
end
output_4=cat(2,grid_4,difference_4);

%difference between individuals in each dyad 700s
x=1:size(total_str_700,1);
y=1:size(total_str_700,1);
[xx,yy]=meshgrid(x,y);
grid_4=[xx(:),yy(:)];
difference_700=nan(size(grid_4,1),1);
count=0;
for ii=1:size(total_str_700,1)
    for jj=1:size(total_str_700,1)
        count=count+1;
        subj1=total_str_700(ii,2);
        subj2=total_str_700(jj,2);
        difference_700(count,1)=abs(subj1-subj2);
    end
end
output_700=cat(2,grid_4,difference_700);

%now get rid of self-self connections and plot
figure
pos_total_similarity_L4=reshape(output_L4(:,3),size(total_str_L4,1),[]);
pos_total_similarity_L4=pos_total_similarity_L4.*~eye(size(pos_total_similarity_L4));
[con1,h]=contourf(pos_total_similarity_L4,2);
h.LineStyle='none';
figure
pos_total_similarity_4=reshape(output_4(:,3),size(total_str_4,1),[]);
pos_total_similarity_4=pos_total_similarity_4.*~eye(size(pos_total_similarity_4));
[con2,h]=contourf(pos_total_similarity_4,2);
h.LineStyle='none';
figure
pos_total_similarity_700=reshape(output_700(:,3),size(total_str_700,1),[]);
pos_total_similarity_700=pos_total_similarity_700.*~eye(size(pos_total_similarity_700));
[con3,h]=contourf(pos_total_similarity_700,2);
h.LineStyle='none';
csvwrite('F:\0_parcellation_analysis\scripts-data-sharing\pos_total_L4.csv',pos_total_similarity_L4);
csvwrite('F:\0_parcellation_analysis\scripts-data-sharing\pos_total_4.csv',pos_total_similarity_4);
csvwrite('F:\0_parcellation_analysis\scripts-data-sharing\pos_total_700.csv',pos_total_similarity_700);
  
%% part III - total negative strength similarities between participants

total_str_L4=nan(23,2);
total_str_4=nan(28,2);
total_str_700=nan(17,2);
for i=1:length(total_str_L4)
total_str_L4(i,1)=(Ns{i});
end
for i=length(total_str_L4)+1:length(neg_total)-length(total_str_700)
total_str_4(i-length(total_str_L4),1)=(Ns{i});
end
for i=length(total_str_L4)+length(total_str_4)+1:length(neg_total)
total_str_700(i-(length(total_str_L4)+length(total_str_4)),1)=(Ns{i});
end
total_str_L4(:,2)=neg_total(1:length(total_str_L4),1);
total_str_4(:,2)=neg_total(length(total_str_L4)+1:length(neg_total)-length(total_str_700),1);
total_str_700(:,2)=neg_total(length(total_str_L4)++length(total_str_4)+1:length(neg_total),1);

%difference between individuals in each dyad L4
x=1:size(total_str_L4,1);
y=1:size(total_str_L4,1);
[xx,yy]=meshgrid(x,y);
grid_L4=[xx(:),yy(:)];
difference_L4=nan(size(grid_L4,1),1);
count=0;
for ii=1:size(total_str_L4,1)
    for jj=1:size(total_str_L4,1)
        count=count+1;
        subj1=total_str_L4(ii,2);
        subj2=total_str_L4(jj,2);
        difference_L4(count,1)=abs(subj1-subj2);
    end
end
output_L4=cat(2,grid_L4,difference_L4);

%difference between individuals in each dyad year 4
x=1:size(total_str_4,1);
y=1:size(total_str_4,1);
[xx,yy]=meshgrid(x,y);
grid_4=[xx(:),yy(:)];
difference_4=nan(size(grid_4,1),1);
count=0;
for ii=1:size(total_str_4,1)
    for jj=1:size(total_str_4,1)
        count=count+1;
        subj1=total_str_4(ii,2);
        subj2=total_str_4(jj,2);
        difference_4(count,1)=abs(subj1-subj2);
    end
end
output_4=cat(2,grid_4,difference_4);

%difference between individuals in each dyad 700s
x=1:size(total_str_700,1);
y=1:size(total_str_700,1);
[xx,yy]=meshgrid(x,y);
grid_4=[xx(:),yy(:)];
difference_700=nan(size(grid_4,1),1);
count=0;
for ii=1:size(total_str_700,1)
    for jj=1:size(total_str_700,1)
        count=count+1;
        subj1=total_str_700(ii,2);
        subj2=total_str_700(jj,2);
        difference_700(count,1)=abs(subj1-subj2);
    end
end
output_700=cat(2,grid_4,difference_700);

%now get rid of self-self connections and plot
figure
neg_total_similarity_L4=reshape(output_L4(:,3),size(total_str_L4,1),[]);
neg_total_similarity_L4=neg_total_similarity_L4.*~eye(size(neg_total_similarity_L4));
[con1,h]=contourf(neg_total_similarity_L4,2);
h.LineStyle='none';
figure
neg_total_similarity_4=reshape(output_4(:,3),size(total_str_4,1),[]);
neg_total_similarity_4=neg_total_similarity_4.*~eye(size(neg_total_similarity_4));
[con2,h]=contourf(neg_total_similarity_4,2);
h.LineStyle='none';
figure
neg_total_similarity_700=reshape(output_700(:,3),size(total_str_700,1),[]);
neg_total_similarity_700=neg_total_similarity_700.*~eye(size(neg_total_similarity_700));
[con3,h]=contourf(neg_total_similarity_700,2);
h.LineStyle='none';
csvwrite('F:\0_parcellation_analysis\scripts-data-sharing\neg_total_L4.csv',neg_total_similarity_L4);
csvwrite('F:\0_parcellation_analysis\scripts-data-sharing\neg_total_4.csv',neg_total_similarity_4);
csvwrite('F:\0_parcellation_analysis\scripts-data-sharing\neg_total_700.csv',neg_total_similarity_700);
  