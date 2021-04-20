%Prediction modelling for rs-fMRI data from MOTIVC study

%writes to a csv file called MODELLING

%It will also use the brain connectivity toolbox 
% https://sites.google.com/site/bctnet/
%(which you will need to install) to get connectomics measurements

%add BCT to path before running this script
addpath(genpath('F:\0_parcellation_analysis\scripts-data-sharing\'))
addpath(genpath('F:\BCT\'))
%cd E:\FSLNets_stuff\deltacon\deltacon_matlab\GRAPHS
cd F:\0_parcellation_analysis\scripts-data-sharing\
%% part II - Get difference in none-node correlation for each connection
%between members of a dyad
yL4=nan(23,1);
y4=nan(28,1);
y700=nan(17,1);
ppts=cell2mat(Ns)';
% pts_zs=cat(2,ppts,all_zs);
  
zs_L4=all_zs(1:size(yL4,1),:);
zs_4=all_zs((size(yL4,1)+1):((size(yL4,1))+size(y4,1)),:);
zs_700=all_zs(((size(yL4,1))+size(y4,1)+1):size(ppts,1),:);

% %difference between individuals in each dyad L4
x=1:size(yL4,1);
y=1:size(yL4,1);
[xx,yy]=meshgrid(x,y);
grid_L4=[xx(:),yy(:)];
difference_L4=nan(size(grid_L4,1),size(all_zs,2));
count=0;

for ii=1:size(yL4,1)
    for jj=1:size(yL4,1)
        count=count+1;
        parfor z_col=1:length(all_zs)
        subj1=zs_L4(ii,z_col);
        subj2=zs_L4(jj,z_col);
        difference_L4(count,z_col)=abs(subj1-subj2);
        end
    end
end
output_L4=cat(2,grid_L4,difference_L4);

%difference between individuals in each dyad year 4
x=1:size(y4,1);
y=1:size(y4,1);
[xx,yy]=meshgrid(x,y);
grid_4=[xx(:),yy(:)];
difference_4=nan(size(grid_4,1),size(all_zs,2));
count=0;
tic
for ii=1:size(y4,1)
    for jj=1:size(y4,1)
        count=count+1;
        parfor z_col=1:length(all_zs)
        subj1=zs_4(ii,z_col);
        subj2=zs_4(jj,z_col);
        difference_4(count,z_col)=abs(subj1-subj2);
        end
    end
end
toc
output_4=cat(2,grid_4,difference_4);

%difference between individuals in each dyad 700s
x=1:size(y700,1);
y=1:size(y700,1);
[xx,yy]=meshgrid(x,y);
grid_700=[xx(:),yy(:)];
difference_700=nan(size(grid_700,1),size(all_zs,2));
count=0;

for ii=1:size(y700,1)
    for jj=1:size(y700,1)
        count=count+1;
        parfor z_col=1:length(all_zs)
        subj1=zs_700(ii,z_col);
        subj2=zs_700(jj,z_col);
        difference_700(count,z_col)=abs(subj1-subj2);
        end
    end
end
output_700=cat(2,grid_700,difference_700);
%save variables so that you don't have to run the top two sections
%everytime
save pred_modelling_var

%% EXCLUDE ALL SELF-SELF AND REPEATED (person TO person) CONNECTIONS
% addpath(genpath('F:\0_parcellation_analysis\scripts-data-sharing\'))
% addpath(genpath('F:\BCT\'))
% 
% cd F:\0_parcellation_analysis\scripts-data-sharing\
load pred_modelling_var
%PREPARE MATRICES FOR DATA INPUT IN FOLLOWING TWO LOOPS
layer_L4 = nan(size(yL4,1),size(yL4,1),size(zs_L4,2));
layer_4 = nan(size(y4,1),size(y4,1),size(zs_4,2));
layer_700 = nan(size(y700,1),size(y700,1),size(zs_700,2));

DATA_L4=nan(((size(yL4,1)^2)-size(yL4,1))/2);
DATA_4=nan(((size(y4,1)^2)-size(y4,1))/2);
DATA_700=nan(((size(y700,1)^2)-size(y700,1))/2);

%CREATE A MATRIX FOR EVERY NODE-NODE CONNECTION DIFFERENCE((272*272-272)/2 CONNECTIONS)
parfor l=1:size(layer_L4,3)
   % layer_L4(:,:,l)=vec2mat(difference_L4(:,l),size(yL4,1));vec2mat no
   % longer recommended
    layer_L4(:,:,l)=reshape(difference_L4(:,l),size(yL4,1),[]);
end
parfor l=1:size(layer_4,3)
%     layer_4(:,:,l)=vec2mat(difference_4(:,l),size(y4,1));
    layer_4(:,:,l)=reshape(difference_4(:,l),size(y4,1),[]);
end
parfor l=1:size(layer_700,3)
%     layer_700(:,:,l)=vec2mat(difference_700(:,l),size(y700,1));
    layer_700(:,:,l)=reshape(difference_700(:,l),size(y700,1),[]);
end

% TAKE THE MATRICES CREATED ABOVE AND MASK THE LOWER TRIANGLE ONLY. THEN
% COPY THAT INFORMATION INTO A FILE CALLED DATA_yeargroup
for l=1:size(layer_L4,3)
    mat=layer_L4(:,:,l);
    mask=tril(true(size(mat)),-1);
    vec=mat(mask);
    DATA_L4(:,l)=vec;
end
for l=1:size(layer_4,3)
    mat=layer_4(:,:,l);
    mask=tril(true(size(mat)),-1);
    vec=mat(mask);
    DATA_4(:,l)=vec;
end
for l=1:size(layer_700,3)
    mat=layer_700(:,:,l);
    mask=tril(true(size(mat)),-1);
    vec=mat(mask);
    DATA_700(:,l)=vec;
end
% Separate the output file from R containing all the distance 
%and community data in it so that it can be used in the model (as DVs)
%THEN ADD ATTRIBUTE NAMES SO THAT PREDICTIVE MODELLING WORKS
attr_names=1:length(all_zs)+1;

dist=dlmread('yL4_dist.csv');
DATA_L4d=cat(2,dist,DATA_L4);
DATA_L4d=cat(1,attr_names,DATA_L4d);
comm=dlmread('yL4_comm.csv');
DATA_L4c=cat(2,comm,DATA_L4);
DATA_L4c=cat(1,attr_names,DATA_L4c);

dist=dlmread('y4_dist.csv');
DATA_4d=cat(2,dist,DATA_4);
DATA_4d=cat(1,attr_names,DATA_4d);
comm=dlmread('y4_comm.csv');
DATA_4c=cat(2,comm,DATA_4);
DATA_4c=cat(1,attr_names,DATA_4c);

dist=dlmread('y700_dist.csv');
DATA_700d=cat(2,dist,DATA_700);
DATA_700d=cat(1,attr_names,DATA_700d);
comm=dlmread('y700_comm.csv');
DATA_700c=cat(2,comm,DATA_700);
DATA_700c=cat(1,attr_names,DATA_700c);

%SAVE DATA FILES
csvwrite('F:\0_parcellation_analysis\scripts-data-sharing\MODELLING_L4.csv',DATA_L4);
csvwrite('F:\0_parcellation_analysis\scripts-data-sharing\MODELLING_4.csv',DATA_4);
csvwrite('F:\0_parcellation_analysis\scripts-data-sharing\MODELLING_700.csv',DATA_700);
csvwrite('F:\0_parcellation_analysis\scripts-data-sharing\MODELLING_dist_L4.csv',DATA_L4d);
csvwrite('F:\0_parcellation_analysis\scripts-data-sharing\MODELLING_dist_4.csv',DATA_4d);
csvwrite('F:\0_parcellation_analysis\scripts-data-sharing\MODELLING_dist_700.csv',DATA_700d);
csvwrite('F:\0_parcellation_analysis\scripts-data-sharing\MODELLING_comm_L4.csv',DATA_L4c);
csvwrite('F:\0_parcellation_analysis\scripts-data-sharing\MODELLING_comm_4.csv',DATA_4c);
csvwrite('F:\0_parcellation_analysis\scripts-data-sharing\MODELLING_comm_700.csv',DATA_700c);
%


