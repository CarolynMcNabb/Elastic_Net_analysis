%Carolyn McNabb Jan 2019
%This script takes the timeseries from each of 272 parcels and correlates it with every
%other parcel. Then selects only those parcels that are included in a RSN
%(either DMN, FPN or exec control network). Any parcels with less than 50
%voxels included in the network are excluded.
%then compares correlation matrices between participants and
%writes to a csv file for each network called "RSN"_all
Ns={201	204	210	211	213	216	217	218	219	220	222	223	224	226	227	232	233	240	244	245	246	247	256	300	302	305	307	308	309	312	315	317	318	321	322	325	332	333	334	335	340	344	346	347	350	352	355	361	366	367	369 702 704 706 708 709 711 713 715 720 722 737 739 742 743 745 748 751};
rFPN=[2	4	6	11	12	17	19	21	23	25	27	30	32	34	36	38	40	42	47	55	58	60	62	66	68	79	81	83	85	88	92	94	98	105	107	133	135	137	139	144	146	147	148	149	150	151	152	154	156	165	167	169	171	180	189	193	199	202	223	225	272];
lFPN=[1	3	5	7	11	12	15	17	20	22	24	26	28	30	32	34	37	39	41	45	53	56	58	64	66	77	81	85	89	92	96	102	104	113	128	132	134	136	140	142	144	146	147	149	151	153	155	159	161	168	176	185	186	189	221	240	268	269];
DMN=[14	15	45	46	51	53	54	95	97	149	150	151	152	155	158	159	162	163	166	168	169	170	193	194	200	201	205	206	207	211	212	215	217	218];
salience=[1	2	9	10	11	13	14	17	18	23	24	26	28	35	37	40	42	43	45	46	47	60	61	72	73	74	75	78	86	87	92	164	166	167	168	172	173	178	187	197	198	199	200	201	202	204	205	206	210	211	212	213	217	218	220];
exec=[1	2	3	4	5	6	8	12	13	14	15	16	17	21	22	23	24	25	26	28	29	31	35	42	45	46	57	82	85	166	167	168	183	185	186	187	191	193	195	196	197	198	199	201	203	204	206	207	208	211	218	219	229	230	243	244	246	247	248	249	250	251	253	255	256	257	259	263	264	265	266	268	270	271	272];
all_corr=zeros(length(Ns),(272^2-272)/2);
all_zs=zeros(length(Ns),(272^2-272)/2);
for i=1:length(Ns)
    ppt=(Ns{i});
    name=sprintf('timeseries_rest_%d.txt', ppt);
    ts=dlmread(name);
    corr_ts=corr(ts');
    corr_ts=corr_ts.*~eye(size(corr_ts));
    net_DMN=corr_ts(DMN,DMN);
    net_salience=corr_ts(salience,salience);
    net_lFPN=corr_ts(lFPN, lFPN);
    net_rFPN=corr_ts(rFPN, rFPN);
    
    
    %remove upper half of matrix
    mask=tril(true(size(corr_ts)),-1);
    vec=corr_ts(mask)'; 


    z=atanh(vec);%RtoZ transformation
    all_corr(i,:)=vec;
    all_zs(i,:)=z;
    %DMN correlations
    mask=tril(true(size(net_DMN)),-1);
    vec=corr_ts(mask)';

    z=atanh(vec);%RtoZ transformation
    DMN_corr(i,:)=vec;
    DMN_zs(i,:)=z;
    %salience network correlations
        mask=tril(true(size(net_salience)),-1);
    vec=corr_ts(mask)';

    z=atanh(vec);%RtoZ transformation
    salience_corr(i,:)=vec;
    salience_zs(i,:)=z;
    %left FPN correlations
    mask=tril(true(size(net_lFPN)),-1);
    vec=corr_ts(mask)';

    z=atanh(vec);%RtoZ transformation
    lFPN_corr(i,:)=vec;
    lFPN_zs(i,:)=z;
    %right FPN correlations
     mask=tril(true(size(net_rFPN)),-1);
    vec=corr_ts(mask)';

    z=atanh(vec);%RtoZ transformation
    rFPN_corr(i,:)=vec;
    rFPN_zs(i,:)=z;

    
end
% correlations=corr(all_corr');
correlations=corr(all_zs');
Unique_corr=tril(correlations);
Unique_corr=Unique_corr.*~eye(size(Unique_corr));
csvwrite('F:\0_parcellation_analysis\scripts-data-sharing\corr_all.csv',Unique_corr);

DMNcorrelations=corr(DMN_zs');
Unique_DMNcorr=tril(DMNcorrelations);
Unique_DMNcorr=Unique_DMNcorr.*~eye(size(Unique_DMNcorr));
csvwrite('F:\0_parcellation_analysis\scripts-data-sharing\corrDMN_all.csv',Unique_DMNcorr);

saliencecorrelations=corr(salience_zs');
Unique_saliencecorr=tril(saliencecorrelations);
Unique_saliencecorr=Unique_saliencecorr.*~eye(size(Unique_saliencecorr));
csvwrite('F:\0_parcellation_analysis\scripts-data-sharing\corrsalience_all.csv',Unique_saliencecorr);

lFPNcorrelations=corr(lFPN_zs');
Unique_lFPNcorr=tril(lFPNcorrelations);
Unique_lFPNcorr=Unique_lFPNcorr.*~eye(size(Unique_lFPNcorr));
csvwrite('F:\0_parcellation_analysis\scripts-data-sharing\corrlFPN_all.csv',Unique_lFPNcorr);

rFPNcorrelations=corr(rFPN_zs');
Unique_rFPNcorr=tril(rFPNcorrelations);
Unique_rFPNcorr=Unique_rFPNcorr.*~eye(size(Unique_rFPNcorr));
csvwrite('F:\0_parcellation_analysis\scripts-data-sharing\corrrFPN_all.csv',Unique_rFPNcorr);

%     eval(['ts' num2str(ppt) '=dlmread(name)']);
%     eval(['corr_ts' num2str(ppt) '=corr(eval([''ts'' num2str(ppt)])'')']);
