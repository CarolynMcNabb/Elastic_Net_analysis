%takes the timeseries from each of 272 parcels and correlates it with every
%other parcel. then runs community Louvain to identify modules in the
%network
addpath(genpath('F:\BCT\'));
addpath(genpath('F:\0_parcellation_analysis\'));
%addpath(genpath('F:\BrainNet-Viewer-master\'));
Ns={201	204	210	211	213	216	217	218	219	220	222	223	224	226	227 ...	
    232	233	240	244	245	246	247	256	300	302	305	307	308	309	312 ...	
    315	317	318	321	322	325	332	333	334	335	340	344	346	347	350 ... 	
    352	355	361	366	367	369 702 704 706 708 709 711 713 ...
    715 720 722 737 739 742 743 745 748 751};

%% part I - community structure for each ppt - sufficient for visualisation
modularity=zeros(length(Ns),1);
n_modules=zeros(length(Ns),1);
communities=zeros(length(Ns),272);
strength_pos=zeros(length(Ns),272);
strength_neg=zeros(length(Ns),272);
diversity_pos=zeros(length(Ns),272);
diversity_neg=zeros(length(Ns),272);
h_star = zeros(length(Ns),272);
s_star = zeros(length(Ns),272);
for i=1:length(Ns)
    ppt=(Ns{i});
    name=sprintf('timeseries_rest_%d.txt', ppt);
    ts=dlmread(name);
    corr_ts=corr(ts');
    corr_ts=corr_ts.*~eye(size(corr_ts));
    eval(['corr_ts' num2str(ppt) '=corr_ts']);
    %community detection (Ci will give module structure, Q will give modularity)
    %gamma is default (1),
    %M0 is empty (optional initial community affiliation vector)
    [Ci,Q]=community_louvain(corr_ts,[],[],'negative_asym');
    eval(['communities_' num2str(ppt) '=Ci']);
    eval(['modularity_' num2str(ppt) '=Q']);
    communities(i,:)=Ci;
    n_modules(i,1)=max(communities(i,:));
    modularity(i,1)=Q;
    %Diversity coefficient - provides a by-node diversity value for each parcel
    %Nodes with high diversity coefficients tend to facilitate global
    %intermodular integration; whereas nodes with low diversity coefficients
    %are characterised by more intramodular connections. Both positive and
    %negative diversity can be determined, providing information about
    %correlated and anticorrelated functional connectivity, respectively.
    [Hpos, Hneg] = diversity_coef_sign(corr_ts,Ci);
    eval(['diversity_neg_' num2str(ppt) '=Hneg']);
    eval(['diversity_pos_' num2str(ppt) '=Hpos']);
    diversity_pos(i,:)=Hpos;
    diversity_neg(i,:)=Hneg;
    %connection strength for positive and negative weights
    [Spos, Sneg] = strengths_und_sign(corr_ts);
    eval(['strength_pos_' num2str(ppt) '=Spos']);
    eval(['strength_neg_' num2str(ppt) '=Sneg']);
    strength_pos(i,:)=Spos;
    strength_neg(i,:)=Sneg;
    %find s* and h* using the formulas from Rubinov (2011) paper
    s_star(i,:) = Spos - (Sneg/(Spos + Sneg))*Sneg;
    h_star(i,:)=Hpos-(Sneg/(Spos + Sneg))*Hneg;
    
end
csvwrite('F:\0_parcellation_analysis\scripts-data-sharing\brain_modularity.csv',modularity);
csvwrite('F:\0_parcellation_analysis\scripts-data-sharing\brain_communities.csv',communities);
% csvwrite('F:\0_parcellation_analysis\scripts-data-sharing\diversity_pos.csv',diversity_pos);
% csvwrite('F:\0_parcellation_analysis\scripts-data-sharing\diversity_neg.csv',diversity_neg);
coords=dlmread('parcel_coords.csv');
parcels=reshape(coords,272,[]);
%change i in communities(i,:) to whichever ppt you would like to visualise
parcels(:,4)=(communities(12,:)');parcels(:,5)=ones(length(parcels),1);parcels(:,6)=(communities(1,:)');
dlmwrite('F:\0_parcellation_analysis\scripts-data-sharing\nodes.node',parcels,' ')
%now use BrainNet to visualise the modules in 3D brain space.

%correlations between participants for diversity of individual nodes
corr_div=corr(diversity_pos');
Unique_div=tril(corr_div);
Unique_div=Unique_div.*~eye(size(Unique_div));
csvwrite('F:\0_parcellation_analysis\scripts-data-sharing\corr_div_pos.csv',Unique_div);

corr_div=corr(diversity_neg');
Unique_div=tril(corr_div);
Unique_div=Unique_div.*~eye(size(Unique_div));
csvwrite('F:\0_parcellation_analysis\scripts-data-sharing\corr_div_neg.csv',Unique_div);

%correlations between participants for s* and h* - See Rubinov 2011 for
%definitions
corr_div=corr(h_star');
Unique_div=tril(corr_div);
Unique_div=Unique_div.*~eye(size(Unique_div));
csvwrite('F:\0_parcellation_analysis\scripts-data-sharing\corr_h_star.csv',Unique_div);

corr_str=corr(s_star');
Unique_str=tril(corr_str);
Unique_str=Unique_str.*~eye(size(Unique_str));
csvwrite('F:\0_parcellation_analysis\scripts-data-sharing\corr_s_star.csv',Unique_str);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% part II - modularity similarities between participants
%do not use this output to compare the similarity in functional
%connectivity and the similarity in modularity - this is kind of like
%double dipping. At least, it's not going to tell you anything useful.
modulesim_L4=nan(23,2);
modulesim_4=nan(28,2);
modulesim_700=nan(17,2);
for i=1:length(modulesim_L4)
modulesim_L4(i,1)=(Ns{i});
end
for i=length(modulesim_L4)+1:length(modularity)-length(modulesim_700)
modulesim_4(i-length(modulesim_L4),1)=(Ns{i});
end
for i=length(modulesim_L4)+length(modulesim_4)+1:length(modularity)
modulesim_700(i-(length(modulesim_L4)+length(modulesim_4)),1)=(Ns{i});
end
modulesim_L4(:,2)=modularity(1:length(modulesim_L4),1);
modulesim_4(:,2)=modularity(length(modulesim_L4)+1:length(modularity)-length(modulesim_700),1);
modulesim_700(:,2)=modularity(length(modulesim_L4)++length(modulesim_4)+1:length(modularity),1);

%difference between individuals in each dyad L4
x=1:size(modulesim_L4,1);
y=1:size(modulesim_L4,1);
[xx,yy]=meshgrid(x,y);
grid_L4=[xx(:),yy(:)];
difference_L4=nan(size(grid_L4,1),1);
count=0;
for ii=1:size(modulesim_L4,1)
    for jj=1:size(modulesim_L4,1)
        count=count+1;
        subj1=modulesim_L4(ii,2);
        subj2=modulesim_L4(jj,2);
        difference_L4(count,1)=abs(subj1-subj2);
    end
end
output_L4=cat(2,grid_L4,difference_L4);

%difference between individuals in each dyad year 4
x=1:size(modulesim_4,1);
y=1:size(modulesim_4,1);
[xx,yy]=meshgrid(x,y);
grid_4=[xx(:),yy(:)];
difference_4=nan(size(grid_4,1),1);
count=0;
for ii=1:size(modulesim_4,1)
    for jj=1:size(modulesim_4,1)
        count=count+1;
        subj1=modulesim_4(ii,2);
        subj2=modulesim_4(jj,2);
        difference_4(count,1)=abs(subj1-subj2);
    end
end
output_4=cat(2,grid_4,difference_4);

%difference between individuals in each dyad 700s
x=1:size(modulesim_700,1);
y=1:size(modulesim_700,1);
[xx,yy]=meshgrid(x,y);
grid_4=[xx(:),yy(:)];
difference_700=nan(size(grid_4,1),1);
count=0;
for ii=1:size(modulesim_700,1)
    for jj=1:size(modulesim_700,1)
        count=count+1;
        subj1=modulesim_700(ii,2);
        subj2=modulesim_700(jj,2);
        difference_700(count,1)=abs(subj1-subj2);
    end
end
output_700=cat(2,grid_4,difference_700);

%now get rid of self-self connections and plot
figure
module_similarity_L4=reshape(output_L4(:,3),size(modulesim_L4,1),[]);
module_similarity_L4=module_similarity_L4.*~eye(size(module_similarity_L4));
[con1,h]=contourf(module_similarity_L4,2);
h.LineStyle='none';
figure
module_similarity_4=reshape(output_4(:,3),size(modulesim_4,1),[]);
module_similarity_4=module_similarity_4.*~eye(size(module_similarity_4));
[con2,h]=contourf(module_similarity_4,2);
h.LineStyle='none';
figure
module_similarity_700=reshape(output_700(:,3),size(modulesim_700,1),[]);
module_similarity_700=module_similarity_700.*~eye(size(module_similarity_700));
[con3,h]=contourf(module_similarity_700,2);
h.LineStyle='none';
csvwrite('F:\0_parcellation_analysis\scripts-data-sharing\modularity_sim_L4.csv',module_similarity_L4);
csvwrite('F:\0_parcellation_analysis\scripts-data-sharing\modularity_sim_4.csv',module_similarity_4);
csvwrite('F:\0_parcellation_analysis\scripts-data-sharing\modularity_sim_700.csv',module_similarity_700);
%% 
