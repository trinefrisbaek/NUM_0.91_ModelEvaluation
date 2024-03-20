clear;
allparameters={'r^*_d','y','r^*_l','aF','c_{passive}','\alpha_{Max}','\alpha_R','mortHTL',...
    'rhoC:N','fracHTL:N','\epsilon_L','\alpha_N','\epsilon_F','cF','\beta','\sigma',...
    'remin2','reminF','rho','mort2','rem_{POM}','vel_1','vel_2','mHTL'};

maxdepth=200;
%Plot Taylor diagram for alle som densitetsplot?
makeDensityPlot=true;
%Plot taylor diagram of best ones:
plotBestOnes=true;
%plot histogram of all stats for all simulations
plotHistOfStats=true;
% Define acceptable boundaries:
lim_COR=[0.9 1];
lim_RMS=[0 1];
lim_STD=[0.5 1.5];

% Minimum allowed overlap
% minoverlap=0;
minoverlap_CCE=0;% 7i phd
minoverlap_HOT=6; %6;

%Define the Parameters to use statistics on
theversions={'TOTpico','TOTnano','ACpico','ACnano','ACmicro'};
theversions_name={'TOT_{pico}','TOT_{nano}','AC_{pico}','AC_{nano}','AC_{micro}'};

thestats={'cRMS' 'STD' 'COR' 'RMSd'};
% define reference data values:
CORref=1;STDref=1;RMSref=0;

%% Choose dataset
% ----------------------- full range wout. nudging: -------------------------%

HOTfile={'stattable_HOT_nudge_0_ACv3_dmax200_time_final_both1.mat'}; % de 100000
% HOTfile={'stattable_HOT_nudge_0_ACv3_dmax200_time_final_CCEtoHOT.mat'}; %genkørsel af de 9 gode

CCEfile={'stattable_CCE_nudge_0_ACv3_dmax200_time_final_both1.mat'}; % de 100000
% CCEfile={'stattable_CCE_nudge_0_ACv3_dmax200_time_final_CCEtoHOT.mat'}; %genkørsel af de 9 gode
% CCEfile={'stattable_CCE_nudge_0_ACv3_dmax200_time_final_reduced_first.mat'}; % de 10000 restricted





%% Load and combine model results
for i=1:length(CCEfile)
    if i==1
        CCE=load(CCEfile{i});
        CCEfilenames=string(repmat(CCEfile{i},[size(CCE.randParam,1) 1]));
    else
        CCE_extra=load(CCEfile{i});
        [CCE]=CombineDatasets(CCE,CCE_extra);
        CCEfilenames=[CCEfilenames;string(repmat(CCEfile{i},[size(CCE_extra.randParam,1) 1]))];
    end
    
end

% for i=1:length(CCE_Nfile)
%     if i==1
%         CCE_N=load(CCE_Nfile{i});
%     else
%         CCE_N_extra=load(CCE_Nfile{i});
%         [CCE_N]=CombineDatasetsN(CCE_N,CCE_N_extra);
%     end
% end

for i=1:length(HOTfile)
    if i==1
        HOT=load(HOTfile{i});
        HOTfilenames=string(repmat(HOTfile{i},[size(HOT.randParam,1) 1]));
    else
        HOT_extra=load(HOTfile{i});
        [HOT]=CombineDatasets(HOT,HOT_extra);
        HOTfilenames=[HOTfilenames;string(repmat(HOTfile{i},[size(HOT_extra.randParam,1) 1]))];
    end
    
end
% for i=1:length(HOT_Nfile)
%     if i==1
%         HOT_N=load(HOT_Nfile{i});
%     else
%         HOT_N_extra=load(HOT_Nfile{i});
%         [HOT_N]=CombineDatasetsN(HOT_N,HOT_N_extra);
%     end
% end
%% Load data on HOT and CCE 
HOT.TayLanFig1a=load(fullfile('..','..','External_data','Taylor and Landry 2018','Taylor and Landry 2018 data for HOT',['HOT_TayLanFig1a_',num2str(maxdepth),'m_wSTD.mat']));
HOT.GasolFig3=load(fullfile('..','..','External_data','Taylor and Landry 2018','Taylor and Landry 2018 data for HOT',['HOT_GasolFig3_',num2str(maxdepth),'m_wSTD.mat']));
HOT.TayLanFig7=load(fullfile('..','..','External_data','Taylor and Landry 2018','Taylor and Landry 2018 data for HOT',['HOT_TayLanFig7_',num2str(maxdepth),'m_wSTD.mat']));
HOT.TOTpn=load(fullfile('..','..','External_data','Taylor and Landry 2018','Taylor and Landry 2018 data for HOT',['HOT_totPN_',num2str(maxdepth),'m_wSTD.mat']));
CCE.TayLanFig1a=load(fullfile('..','..','External_data','Taylor and Landry 2018','Landry CEE 2004-2011',['CCE_TayLanFig1a_',num2str(maxdepth),'m_wSTD.mat']));
CCE.GasolFig3=load(fullfile('..','..','External_data','Taylor and Landry 2018','Landry CEE 2004-2011',['CCE_GasolFig3_',num2str(maxdepth),'m_wSTD.mat']));
CCE.TayLanFig7=load(fullfile('..','..','External_data','Taylor and Landry 2018','Landry CEE 2004-2011',['CCE_TayLanFig7_',num2str(maxdepth),'m_wSTD.mat']));
CCE.TOTpn=load(fullfile('..','..','External_data','Taylor and Landry 2018','Landry CEE 2004-2011',['CCE_totPN_',num2str(maxdepth),'m_wSTD.mat']));


%% extract info:
% ----------------------------- filenames: -------------------------------%
CCE_names=CCE.stattable.filename;
HOT_names=HOT.stattable.filename;
% -------------------------- random paramters: ---------------------------%
CCE_randParam=CCE.randParam;
HOT_randParam=HOT.randParam;
% -------------------------- acceptable ranges: --------------------------%
[max_cRMS_CCE,max_std_CCE,min_cor_CCE,max_RMSd_CCE]=defineAcceptableRange('CCE',theversions);
[max_cRMS_HOT,max_std_HOT,min_cor_HOT,max_RMSd_HOT]=defineAcceptableRange('HOT',theversions);
% ------------------------ overlap with data: ----------------------------%
for ii=1:length(theversions)
    overlapCCE(:,ii)=CCE.stattable.(['nOverlap',theversions{ii}]);
    overlapHOT(:,ii)=HOT.stattable.(['nOverlap',theversions{ii}]);
end

rowname={};
for i=1:11
    rowname{i}=['nr overlapping bins: ',num2str(i-1)];
    overlap_vector(i,:)=sum(overlapHOT==i-1)./size(overlapHOT,1).*100;
end
overlap_vectorHOT=overlap_vector;
disp('Overlaps in % for HOT')
T = array2table(overlap_vector,'VariableNames',theversions,'RowName',rowname);
disp(T)

rowname={};
for i=1:12
    rowname{i}=['nr overlapping bins: ',num2str(i-1)];
    overlap_vector(i,:)=sum(overlapCCE==i-1)./size(overlapCCE,1).*100;
end
overlap_vectorCCE=overlap_vector;
disp('Overlaps in % for CCE')
T = array2table(overlap_vector,'VariableNames',theversions,'RowName',rowname);
disp(T)

%% Load, normalize and sort data
[CCE_CORs,CCE_cRMSs,CCE_STDs,CCE_RMSds]=loadTheData(CCE,theversions);
[HOT_CORs,HOT_cRMSs,HOT_STDs,HOT_RMSds]=loadTheData(HOT,theversions);


%% Find data with minimal overlap to exclude
idxremove_CCE=overlapCCE<minoverlap_CCE;
idxremove_HOT=overlapHOT<minoverlap_HOT;

% Plot histogram of the "bad" ones
disp('This is the number that is sorted out for CCE:')
disp(num2str(sum(idxremove_CCE)))
% i=2;
% plotHistParamCompare(CCE_randParam,idxremove_CCE(:,i),20,'yes')
% sgtitle(['Parameters for frasorterede CCE ',theversions_name{i}])

disp('This is the number that is sorted out for HOT:')
disp(num2str(sum(idxremove_HOT)))
% i=2;
% plotHistParamCompare(HOT_randParam,~idxremove_HOT(:,i),20,'yes')
% sgtitle(['Parameters for frasorterede HOT ',theversions_name{i}])
%
largeCCE=(~isnan(CCE.AC_pico_mean_bin_save(:,13)));
largeHOT=(~isnan(HOT.AC_pico_mean_bin_save(:,13)));

%% Set results with too few overlapping points to NaN;
HOT_CORs(idxremove_HOT)=NaN;HOT_cRMSs(idxremove_HOT)=NaN;HOT_STDs(idxremove_HOT)=NaN;HOT_RMSds(idxremove_HOT)=NaN;
CCE_CORs(idxremove_CCE)=NaN;CCE_cRMSs(idxremove_CCE)=NaN;CCE_STDs(idxremove_CCE)=NaN;CCE_RMSds(idxremove_CCE)=NaN;

%% Plot figures on ALL data (above min overlap)
xx=0.05; %max ratio on CLim
% plotRandomDensity(CCE,'CCE',maxdepth,xx)
%%
xx=0.04;
% plotRandomDensity(HOT,'HOT',maxdepth,xx)
%%
% ** Plot taylor diagram with density **
% plotTaylorDensity(theversions,theversions_name,CCE_CORs,CCE_cRMSs,CCE_STDs);
%%
% plotTaylorDensity(theversions,theversions_name,HOT_CORs,HOT_cRMSs,HOT_STDs);
%%
% plotTaylorDensityWPOINTS(theversions,theversions_name,HOT_CORs,HOT_cRMSs,HOT_STDs);

%% Which data is best?
load('mDelta.mat');
mass=CCE.mass;
biospecter_CCE=CCE.sumB./mDelta(3:12);

biospecter_HOT=CCE.sumB./mDelta(3:12);
% plotHistParamCompare(CCE_randParam,idx,nbins,'yes')
%%
idx=WhichAreTheBest('HOT',HOTfilenames,theversions,theversions_name,HOT_cRMSs,HOT_STDs,HOT_CORs,HOT_RMSds,HOT_names,thestats,max_cRMS_HOT,max_std_HOT,min_cor_HOT,max_RMSd_HOT,HOT_randParam,allparameters,biospecter_HOT,mass,idxremove_HOT,HOT,largeHOT);
%%
idx=WhichAreTheBest('CCE',CCEfilenames,theversions,theversions_name,CCE_cRMSs,CCE_STDs,CCE_CORs,CCE_RMSds,CCE_names,thestats,max_cRMS_CCE,max_std_CCE,min_cor_CCE,max_RMSd_CCE,CCE_randParam,allparameters,biospecter_CCE,mass,idxremove_CCE,CCE,largeCCE);

theRandParam=HOT_randParam(idx,:);
V=[];
thepct=0.5;
for i=1:24
    withinRange{i}=find(HOT_randParam(:,i)>theRandParam(i)*(1-thepct) & HOT_randParam(:,i)<theRandParam(i)*(1+thepct));
    V=[V;withinRange{i}];
end
nrWithin=accumarray(V,1);
find(nrWithin==24); 

idx=WhichAreTheBest('HOT',theversions,theversions_name,HOT_cRMSs,HOT_STDs,HOT_CORs,HOT_RMSds,HOT_names,thestats,max_cRMS_HOT,max_std_HOT,min_cor_HOT,max_RMSd_HOT,HOT_randParam,allparameters,HOT)


%%

idx=WhichAreTheBest('CCE',theversions,theversions_name,CCE_cRMSs,CCE_STDs,CCE_CORs,CCE_RMSds,CCE_names,thestats,max_cRMS_CCE,max_std_CCE,min_cor_CCE,max_RMSd_CCE,CCE_randParam,allparameters)
%% Look at the "good" data
% plot taylor for the 10 best COR for HOT for stat nr {mypar}:
mypar=3;
best=find(HOT_cRMSs(:,mypar)<1);
[RMS_origposit,RMS_idx]=sortrows(HOT_cRMSs(:,mypar));
RMS_idx=RMS_idx(1);
HOT_names(RMS_idx,1)
disp(['RMS: ',num2str(HOT_cRMSs(RMS_idx,mypar))])
disp(['STD: ',num2str(HOT_STDs(RMS_idx,mypar))])
disp(['COR: ',num2str(HOT_CORs(RMS_idx,mypar))])


% best=find(HOT_RMSs_rating(:,mypar)<11);
figure
subplot(1,2,1)
taylordiag([STDref;HOT_STDs(best,mypar)],[RMSref;HOT_cRMSs(best,mypar)],[CORref;HOT_CORs(best,mypar)],'plotScatter',0,'limSTD',3,'titleCOR',0,'titleSTD',0,'titleRMS',0,'tickRMS',[0.5:0.5:3]);

%%



% Define acceptable boundaries;
lim_COR=[0.9 1];lim_RMS=[0 3];lim_STD=[-3 3];
% lim_COR=[0.9 1];lim_RMS=[0 1.5];lim_STD=[0.5 1.5];
lim_COR=repmat(lim_COR,[length(theversions) 1]);
lim_RMS=repmat(lim_RMS,[length(theversions) 1]);
lim_STD=repmat(lim_STD,[length(theversions) 1]);
% lim_STD(1,:)=[0 2];

% find the indices for iterations with results within this area
[isolate_HOT_idx,maxrep_HOT,idxMaxRep_HOT]=findFittingpoints(HOT_CORs,HOT_cRMSs,HOT_STDs,lim_COR,lim_RMS,lim_STD,theversions);
[isolate_CCE_idx,maxrep_CCE,idxMaxRep_CCE]=findFittingpoints(CCE_CORs,CCE_cRMSs,CCE_STDs,lim_COR,lim_RMS,lim_STD,theversions);
disp ('*** Simulations within boundaries for HOT: ***')
for ii=1:length(theversions)
    nr_in=sum(~isnan(isolate_HOT_idx(:,ii)));
    if nr_in>0
        disp(['Categori ',theversions{ii},':'])
        for i=1:nr_in
            disp(['STD: ',num2str(HOT_STDs(isolate_HOT_idx(i,ii),ii)),', COR: ',num2str(HOT_CORs(isolate_HOT_idx(i,ii),ii)),', RMS: ',num2str(HOT_cRMSs(isolate_HOT_idx(i,ii),ii)),', ',char(HOT_names(isolate_HOT_idx(i,ii)))])
        end
    else
        disp(['Categori ',theversions{ii},' has none that is within this range'])
    end
end

plotFromName(HOT_names(isolate_HOT_idx(1,1)),'wide',thissetup)
plotFromName("p_point_1_restore_a_nrRanditer_497_2023_01_18_16_27_rand0.025445.mat",'wide',thissetup)

theta=real(acos(CORs(idx,myversion)));
thex=STDs(idx,myversion).*cos(theta);
they=STDs(idx,myversion).*sin(theta);



disp(['For HOT max overlap is ',num2str(maxrep_HOT),' for files:'])
HOT_thebest=idxMaxRep_HOT{maxrep_HOT}
disp(['For CCE max overlap is ',num2str(maxrep_CCE),' for files:'])

for ii=1:6;plotFromName(CCE_names(isolate_CCE_idx(ii,1)),'wide',thissetup);pause;end




% find the ones that are good in the first catagory
[val,pos]=intersect(isolate_CCE_idx(:,1),cce_thebest)
CCE_Names_to_plot=CCE_names(idxMaxRep_CCE{maxrep_CCE});

plotFromName(CCE_Names_to_plot(1),'wide',thissetup)



% Define acceptable boundaries:
lim_COR=[0.9 1];lim_RMS=[0 1.5];lim_STD=[0.5 1.5];
% find the indices for iterations with results within this area

disp(['For CCE max overlap is ',num2str(maxrep_CCE),' for files:'])
cce_thebest=idxMaxRep_CCE{maxrep_CCE};
% find the ones that are good in the first catagory
[val,pos]=intersect(isolate_CCE_idx(:,1),cce_thebest)
CCE_Names_to_plot=CCE_names(idxMaxRep_CCE{maxrep_CCE});

plotFromName(CCE_Names_to_plot(1),'wide',thissetup)

[isolate_HOT_idx,maxrep_HOT,idxMaxRep_HOT]=findFittingpoints(HOT_CORs,HOT_cRMSs,HOT_STDs,lim_COR,lim_RMS,lim_STD,theversions);

%plot them on a taylor diagram
plotTaylor_orig(theversions,CCE_CORs,CCE_cRMSs,CCE_STDs,HOT_CORs,HOT_cRMSs,HOT_STDs,...
    isolate_CCE_idx,isolate_HOT_idx);

% %plot histogram of their parameters for all 5 categories
% nbins=20;
% normalization='yes';
% for ii=1:length(theversions)
%     plotHistParamCompare(HOT_randParam,isolate_HOT_idx{ii},nbins,normalization)
%     sgtitle(['HOT, the best in category ',theversions{ii}])
%
plotHistParamCompare(CCE_randParam,isolate_CCE_idx{ii},nbins,normalization)
%     sgtitle(['CCE, the best in category ',theversions{ii}])
% end


%% Find certain data
myversion=3;
idx_inside=findDatapoints(myversion,theversions,CCE_CORs,CCE_cRMSs,CCE_STDs);

plotHistParamCompare(CCE_randParam,idx_inside,20,'yes')

CCE_Names_to_plot=CCE_names(idx_inside);
CCE_Names_to_plot=CCE_Names_to_plot(1:2);





figure('Color','w','Name','Mean biomass');
[AC_bin,~,~,~]=createACbin;
AC_x=(AC_bin(1:14-1)+AC_bin(2:14))./2;

CCE.AC_pico_mean_bin_save(idxRemove_CCE,:)=[];
CCE.AC_nano_mean_bin_save(idxRemove_CCE,:)=[];
CCE.AC_micro_mean_bin_save(idxRemove_CCE,:)=[];
HOT.AC_pico_mean_bin_save(idxRemove_CCE,:)=[];
HOT.AC_nano_mean_bin_save(idxRemove_CCE,:)=[];
HOT.AC_micro_mean_bin_save(idxRemove_CCE,:)=[];

figure('Color','w','Name','Mean biomass');
subplot(2,3,1)
hold on
tmp=mean(CCE.AC_pico_mean_bin_save,1,'omitnan');
mmax= prctile(CCE.AC_pico_mean_bin_save,100); mmin= prctile(CCE.AC_pico_mean_bin_save,0); mmin(mmin<0.001)=0.001;
fill([AC_x fliplr(AC_x)],[mmax(2:end) fliplr(mmin(2:end))],'k','EdgeColor','none','FaceAlpha',0.2)
set(gca,'xscale','log','yscale','log')
loglog(AC_x,CCE.TayLanFig1a.AC_pico_mean_bin(2:14),'-ok'); hold on;
loglog(AC_x,tmp(2:14),'-or'); hold on;
tmp=median(CCE.AC_pico_mean_bin_save,1,'omitnan');
loglog(AC_x,tmp(2:14),'-ob'); hold on;
subtitle('CCE AC_{Pico}')
ylabel('A_{pico}biomass (\mug C/l)');legend('','Data','Mean data','Median data')

subplot(2,3,2)
hold on
tmp=mean(CCE.AC_nano_mean_bin_save,1,'omitnan');
mmax= prctile(CCE.AC_nano_mean_bin_save,100); mmin= prctile(CCE.AC_nano_mean_bin_save,0); mmin(mmin<0.001)=0.001;
fill([AC_x fliplr(AC_x)],[mmax(2:end) fliplr(mmin(2:end))],'k','EdgeColor','none','FaceAlpha',0.2)
set(gca,'xscale','log','yscale','log')
loglog(AC_x,CCE.TayLanFig1a.AC_nano_mean_bin(2:14),'-ok'); hold on;
loglog(AC_x,tmp(2:14),'-or'); hold on;
tmp=median(CCE.AC_nano_mean_bin_save,1,'omitnan');
loglog(AC_x,tmp(2:14),'-ob'); hold on;
subtitle('CCE AC_{Nano}')
ylabel('A_{nano}biomass (\mug C/l)');legend('','Data','Mean data','Median data')

subplot(2,3,3)
hold on
tmp=mean(CCE.AC_micro_mean_bin_save,1,'omitnan');
mmax= prctile(CCE.AC_micro_mean_bin_save,100); mmin= prctile(CCE.AC_micro_mean_bin_save,0); mmin(mmin<0.001)=0.001;
fill([AC_x fliplr(AC_x)],[mmax(2:end) fliplr(mmin(2:end))],'k','EdgeColor','none','FaceAlpha',0.2)
set(gca,'xscale','log','yscale','log')
loglog(AC_x,CCE.TayLanFig1a.AC_micro_mean_bin(2:14),'-ok'); hold on;
loglog(AC_x,tmp(2:14),'-or'); hold on;
tmp=median(CCE.AC_micro_mean_bin_save,1,'omitnan');
loglog(AC_x,tmp(2:14),'-ob'); hold on;
subtitle('CCE AC_{Micro}')
ylabel('A_{micro}biomass (\mug C/l)');legend('','Data','Mean data','Median data')

subplot(2,3,4)
hold on
tmp=mean(HOT.AC_pico_mean_bin_save,1,'omitnan');
mmax= prctile(HOT.AC_pico_mean_bin_save,100); mmin= prctile(HOT.AC_pico_mean_bin_save,0); mmin(mmin<0.001)=0.001;
fill([AC_x fliplr(AC_x)],[mmax(2:end) fliplr(mmin(2:end))],'k','EdgeColor','none','FaceAlpha',0.2)
set(gca,'xscale','log','yscale','log')
loglog(AC_x,HOT.TayLanFig1a.AC_pico_mean_bin(2:14),'-ok'); hold on;
loglog(AC_x,tmp(2:14),'-or'); hold on;
tmp=median(HOT.AC_pico_mean_bin_save,1,'omitnan');
loglog(AC_x,tmp(2:14),'-ob'); hold on;
subtitle('HOT AC_{Pico}')
ylabel('A_{pico}biomass (\mug C/l)');legend('','Data','Mean data','Median data')

subplot(2,3,5)
hold on
tmp=mean(HOT.AC_nano_mean_bin_save,1,'omitnan');
mmax= prctile(HOT.AC_nano_mean_bin_save,100); mmin= prctile(HOT.AC_nano_mean_bin_save,0); mmin(mmin<0.001)=0.001;
fill([AC_x fliplr(AC_x)],[mmax(2:end) fliplr(mmin(2:end))],'k','EdgeColor','none','FaceAlpha',0.2)
set(gca,'xscale','log','yscale','log')
loglog(AC_x,HOT.TayLanFig1a.AC_nano_mean_bin(2:14),'-ok'); hold on;
loglog(AC_x,tmp(2:14),'-or'); hold on;
tmp=median(HOT.AC_nano_mean_bin_save,1,'omitnan');
loglog(AC_x,tmp(2:14),'-ob'); hold on;
subtitle('HOT AC_{Nano}')
ylabel('A_{nano}biomass (\mug C/l)');legend('','Data','Mean data','Median data')

subplot(2,3,6)
hold on
tmp=mean(HOT.AC_micro_mean_bin_save,1,'omitnan');
mmax= prctile(HOT.AC_micro_mean_bin_save,100); mmin= prctile(HOT.AC_micro_mean_bin_save,0); mmin(mmin<0.001)=0.001;
fill([AC_x fliplr(AC_x)],[mmax(2:end) fliplr(mmin(2:end))],'k','EdgeColor','none','FaceAlpha',0.2)
set(gca,'xscale','log','yscale','log')
loglog(AC_x,HOT.TayLanFig1a.AC_micro_mean_bin(2:14),'-ok'); hold on;
loglog(AC_x,tmp(2:14),'-or'); hold on;
tmp=median(HOT.AC_micro_mean_bin_save,1,'omitnan');
loglog(AC_x,tmp(2:14),'-ob'); hold on;
subtitle('HOT AC_{Micro}')
ylabel('A_{micro}biomass (\mug C/l)');legend('','Data','Mean data','Median data')

axH = findall(gcf,'type','axes');
set(axH,'Xlim',[2 1000],'Ylim',[10^-1 10^3]);
xlabel(axH,'AC_{tot} (\mug C/l)')

%% Find datapoints within defined acceptable COR, RMS og STD range
% ... and plot them
figure('Color','w','Name','Taylor ones within limits');
idx_CCE=nan(size(CCE_CORs,1),length(theversions));
idx_HOT=nan(size(HOT_CORs,1),length(theversions));

% find de simulationer hvor COR, RMS og STD er gode nok i mindst én
% katagori
for ii=1:length(theversions)
    tmp_CCE=find((CCE_CORs(:,ii)>=lim_COR(1) & CCE_CORs(:,ii)<=lim_COR(2)) & (CCE_cRMSs(:,ii)>=lim_RMS(1) & CCE_cRMSs(:,ii)<=lim_RMS(2)) & (CCE_STDs(:,ii)>=lim_STD(1) & CCE_STDs(:,ii)<=lim_STD(2)))';
    idx_CCE(1:length(tmp_CCE),ii)=tmp_CCE;
    tmp_HOT=find((HOT_CORs(:,ii)>=lim_COR(1) & HOT_CORs(:,ii)<=lim_COR(2)) & (HOT_RMSs(:,ii)>=lim_RMS(1) & HOT_cRMSs(:,ii)<=lim_RMS(2)) & (HOT_STDs(:,ii)>=lim_STD(1) & HOT_STDs(:,ii)<=lim_STD(2)))';
    idx_HOT(1:length(tmp_HOT),ii)=tmp_HOT;

    subplot(2,5,ii)
    %     taylordiag([STDref;CCE_STDs(tmp_CCE,ii)],[RMSref;CCE_RMSs(tmp_CCE,ii)],[CORref;CCE_CORs(tmp_CCE,ii)],'theLabels',tmp_CCE);
    taylordiag_density([STDref;CCE_STDs(tmp_CCE,ii)],[RMSref;CCE_cRMSs(tmp_CCE,ii)],[CORref;CCE_CORs(tmp_CCE,ii)],...
        'plotScatter',1,'limSTD',1.5,'titleCOR',0,'titleSTD',0,'titleRMS',0,'showlabelsCOR',0);
    subtitle(theversions{ii})

    subplot(2,5,length(theversions)+ii)
    %     taylordiag([STDref;HOT_STDs(tmp_HOT,ii)],[RMSref;HOT_RMSs(tmp_HOT,ii)],[CORref;HOT_CORs(tmp_HOT,ii)],'theLabels',tmp_HOT);
    taylordiag_density([STDref;HOT_STDs(tmp_HOT,ii)],[RMSref;HOT_cRMSs(tmp_HOT,ii)],[CORref;HOT_CORs(tmp_HOT,ii)],...
        'plotScatter',1,'limSTD',1.5,'titleCOR',0,'titleSTD',0,'titleRMS',0,'showlabelsCOR',0);
    subtitle(theversions{ii})
end
subplot(2,length(theversions),1);hold on;ylabel('CCE')
subplot(2,length(theversions),length(theversions)+1);hold on;ylabel('HOT')
sgtitle('Simulations within defined limits')
%% Plot the Histograms for these ones for CCE
for ii=1:length(theversions)
    acceptablePoints_CCE=find((CCE_CORs(:,ii)>=lim_COR(1) & CCE_CORs(:,ii)<=lim_COR(2)) & (CCE_cRMSs(:,ii)>=lim_RMS(1) & CCE_cRMSs(:,ii)<=lim_RMS(2)) & (CCE_STDs(:,ii)>=lim_STD(1) & CCE_STDs(:,ii)<=lim_STD(2)))';
    figure;
    isnormal=false(1,size(CCE_randParam,2));
    allparameters={'r^*_d','y','r^*_l','aF','c_{passive}','\alpha_{Max}','\alpha_R','mortHTL',...
        'rhoC:N','fracHTL:N','\epsilon_L','\alpha_N','\epsilon_F','cF','\beta','\sigma',...
        'remin2','reminF','rho','mort2','rem_{POM}','vel_1','vel_2'};
    isnormal([9,10,11,12,13,14,15,16,17,19,20,23])=true;
    isnormal(:)=true;
    nbins=10;
    thecolor={'b','r'};
    for i=1:23
        subplot(6,4,i)
        if isnormal(i)==true
            bins=linspace(floor(min(CCE_randParam(:,i))), ceil(max(CCE_randParam(:,i))), nbins);
            yyaxis left
            h=histogram(CCE_randParam(:,i),bins,'Normalization','probability','EdgeColor','w');
            binlim=h.BinLimits;
            hold on
            yyaxis right
            h=histogram(CCE_randParam(acceptablePoints_CCE,i),bins,'Normalization','probability','Edgecolor','none','FaceColor','none');
            plot((h.BinEdges(2:end)+h.BinEdges(1:end-1))./2,h.Values, ['-',thecolor{2}]);
            xlim(binlim)

        else
            bins=logspace(log10(min(CCE_randParam(:,i))),log10(max(CCE_randParam(:,i))),nbins);
            yyaxis left
            h=histogram(CCE_randParam(:,i),bins,'Normalization','probability','EdgeColor','w');
            binlim=h.BinLimits;
            hold on
            h=histogram(CCE_randParam(acceptablePoints_CCE,i),bins,'Normalization','probability','Edgecolor','none','FaceColor','none');
            plot((h.BinEdges(2:end)+h.BinEdges(1:end-1))./2,h.Values, ['-',thecolor{2}]);
            set(gca,'xscale','log')
            xlim(binlim)
        end
        subtitle(allparameters{i})
    end
    % thetext={'RMDS sum','RMSD HvsA','RMSD H:A vs A','RMSD Apico','RMSD Anano','RMSD Amicro'};
    sgtitle(['Distribution of parameters for CCE the best in ',theversions{ii}]);
end




%% Plot the Histograms for these ones for HOT

for ii=1:length(theversions)
    acceptablePoints_HOT=find((HOT_CORs(:,ii)>=lim_COR(1) & HOT_CORs(:,ii)<=lim_COR(2)) & (HOT_cRMSs(:,ii)>=lim_RMS(1) & HOT_cRMSs(:,ii)<=lim_RMS(2)) & (HOT_STDs(:,ii)>=lim_STD(1) & HOT_STDs(:,ii)<=lim_STD(2)))';

    figure;
    isnormal=false(1,size(HOT_randParam,2));
    allparameters={'r^*_d','y','r^*_l','aF','c_{passive}','\alpha_{Max}','\alpha_R','mortHTL',...
        'rhoC:N','fracHTL:N','\epsilon_L','\alpha_N','\epsilon_F','cF','\beta','\sigma',...
        'remin2','reminF','rho','mort2','rem_{POM}','vel_1','vel_2'};
    isnormal([9,10,11,12,13,14,15,16,17,19,20,23])=true;
    isnormal(:)=true;
    nbins=10;
    thecolor={'b','r'};
    for i=1:23
        subplot(6,4,i)
        if isnormal(i)==true
            bins=linspace(floor(min(HOT_randParam(:,i))), ceil(max(HOT_randParam(:,i))), nbins);
            yyaxis left
            h=histogram(HOT_randParam(:,i),bins,'Normalization','probability','EdgeColor','w');
            binlim=h.BinLimits;
            hold on
            yyaxis right
            h=histogram(HOT_randParam(acceptablePoints_HOT,i),bins,'Normalization','probability','Edgecolor','none','FaceColor','none');
            plot((h.BinEdges(2:end)+h.BinEdges(1:end-1))./2,h.Values, ['-',thecolor{2}]);
            xlim(binlim)

        else
            bins=logspace(log10(min(HOT_randParam(:,i))),log10(max(HOT_randParam(:,i))),nbins);
            yyaxis left
            h=histogram(HOT_randParam(:,i),bins,'Normalization','probability','EdgeColor','w');
            binlim=h.BinLimits;
            hold on
            h=histogram(HOT_randParam(acceptablePoints_HOT,i),bins,'Normalization','probability','Edgecolor','none','FaceColor','none');
            plot((h.BinEdges(2:end)+h.BinEdges(1:end-1))./2,h.Values, ['-',thecolor{2}]);
            set(gca,'xscale','log')
            xlim(binlim)
        end
        subtitle(allparameters{i})
    end
    % thetext={'RMDS sum','RMSD HvsA','RMSD H:A vs A','RMSD Apico','RMSD Anano','RMSD Amicro'};
    sgtitle(['Distribution of parameters for HOT the best in ',theversions{ii}]);
end



















%% Find out if any of these are good in multiple catagories
% Find ud af hvilke af disse simulationer der er gode nok for alle 5 katagorier:
tmp=idx_CCE(:);
nrRep=accumarray(tmp(~isnan(tmp)),1); % Hvor mange er der af hver. rækkenr i nrRep er simulationsnummer
disp(['Max for CCE is: ',num2str(max(nrRep))]);
themin=3;
idxMaxRep_CCE=find(nrRep>=themin); %idxMaxRep er tal svarende til simulationsnummer
%check that they are all good before using them (formulation from taylerdiag.m)
notworking=fix([repmat(RMSref,[1 length(theversions)]);CCE_cRMSs(idxMaxRep_CCE,:)].*100)./100-fix(sqrt([repmat(STDref,[1 length(theversions)]);CCE_STDs(idxMaxRep_CCE,:)].^2 + repmat(STDref(1),[1 length(theversions)]).^2 - 2.*[repmat(STDref,[1 length(theversions)]);CCE_STDs(idxMaxRep_CCE,:)].*repmat(STDref,[1 length(theversions)]).*[repmat(CORref,[1 length(theversions)]);CCE_CORs(idxMaxRep_CCE,:)]).*100)./100;
[x,~]=find(notworking~=0);
if ~isempty(x)
    disp('removing CCE file number :')
    disp(idxMaxRep_CCE(x-1))
    idxMaxRep_CCE(x-1)=[];
end

tmp=idx_HOT(:);
nrRep=accumarray(tmp(~isnan(tmp)),1);
disp(['Max for HOT is: ',num2str(max(nrRep))]);
idxMaxRep_HOT=find(nrRep>=themin);
% %check that they are all good before using them (formulation from taylerdiag.m)
notworking=fix([repmat(RMSref,[1 length(theversions)]);HOT_cRMSs(idxMaxRep_HOT,:)].*100)./100-fix(sqrt([repmat(STDref,[1 length(theversions)]);HOT_STDs(idxMaxRep_HOT,:)].^2 + repmat(STDref(1),[1 length(theversions)]).^2 - 2.*[repmat(STDref,[1 length(theversions)]);HOT_STDs(idxMaxRep_HOT,:)].*repmat(STDref,[1 length(theversions)]).*[repmat(CORref,[1 length(theversions)]);HOT_CORs(idxMaxRep_HOT,:)]).*100)./100;
[x,~]=find(notworking~=0);
if ~isempty(x)
    disp('removing HOT file number :')
    disp(idxMaxRep_HOT(x-1))
    idxMaxRep_HOT(x-1)=[];
end



%% Plot taylor diagram of best ones
if plotBestOnes
    figure('Color','w','Name','Taylor best ones');
    for ii=1:length(theversions)
        subplot(2,5,ii)
        taylordiag([STDref;CCE_STDs(idxMaxRep_CCE,ii)],[RMSref;CCE_cRMSs(idxMaxRep_CCE,ii)],[CORref;CCE_CORs(idxMaxRep_CCE,ii)],'theLabels',idxMaxRep_CCE);
        subtitle(theversions{ii})
        if ii==1
            ylabel('CCE')
        end
        subplot(2,5,length(theversions)+ii)
        taylordiag([STDref;HOT_STDs(idxMaxRep_HOT,ii)],[RMSref;HOT_cRMSs(idxMaxRep_HOT,ii)],[CORref;HOT_CORs(idxMaxRep_HOT,ii)],'theLabels',idxMaxRep_HOT);
        subtitle(theversions{ii})
        if ii==1
            ylabel('HOT')
        end

    end
end

% disp plotting
disp('** Plotting for CCE: **')
for i=1:length(idxMaxRep_CCE)
    CCEplottetOnes(i)=extractBetween(CCE_names(idxMaxRep_CCE(i)),'nrRanditer_','.mat');
    disp(strjoin(['     T',num2str(idxMaxRep_CCE(i)),': ',CCE_names(idxMaxRep_CCE(i))]))
end
disp('** Plotting for HOT: **')
for i=1:length(idxMaxRep_HOT)
    HOTplottetOnes(i)=extractBetween(HOT_names(idxMaxRep_HOT(i)),'nrRanditer_','.mat');
    disp(strjoin(['     T',num2str(idxMaxRep_HOT(i)),': ',HOT_names(idxMaxRep_HOT(i))]))
end

%% Check if any of the good ones are repeated in both HOT and CCE
[D,~,X]=unique([CCEplottetOnes';HOTplottetOnes']);
Y = hist(X,unique(X));
idx=Y==2;

if sum(idx)==0
    disp('No files are both in CCE and HOT')
else
    disp(['file ', D(idx),'is in both'])
end
CCE_further=CCE_names(contains(CCE_names,D(idx)));
HOT_further=HOT_names(contains(HOT_names,D(idx)));


Plot_H_vs_AC=true;
Plot_ACpnm=true;
plotOneResult(CCE_further,CCE.thissetup1,CCE.thissetup,CCE.ACversion,CCE.AC_fraction,CCE.maxdepth,CCE.tstart,Plot_H_vs_AC,Plot_ACpnm)
plotOneResult(HOT_further,HOT.thissetup1,HOT.thissetup,HOT.ACversion,HOT.AC_fraction,HOT.maxdepth,HOT.tstart,Plot_H_vs_AC,Plot_ACpnm)

%% Plot histogram of all stats
nbest=1000;
if plotHistOfStats
    fig1=figure('Color','w','Name','Normalized stats for CCE');
    nbins=100;
    minpct=0;
    maxpct=99;
    for ii=1:length(theversions)
        subplot(3,length(theversions),ii)
        yyaxis left
        bins=linspace(-1,1,nbins);
        histogram(CCE_CORs(:,ii),bins,'Edgecolor','none','FaceColor','b')
        hold on
        histogram(CCE_CORs(CCE_CORs_rating(1:nbest,ii),ii),bins,'Edgecolor','none','FaceColor','r')
        yyaxis right
        histogram(CCE_CORs(idxMaxRep_CCE,ii),bins,'Normalization','Probability','Edgecolor','none','FaceColor','y')
        subtitle(theversions{ii})
        if ii==1
            ylabel('CORRELATION')
        end

        subplot(3,length(theversions),length(theversions)+ii)
        bins=linspace(prctile(CCE_cRMSs(:,ii),minpct),prctile(CCE_cRMSs(:,ii),maxpct),nbins);
        yyaxis left
        histogram(CCE_cRMSs(:,ii),bins,'Edgecolor','none','FaceColor','b')
        hold on
        histogram(CCE_cRMSs(CCE_RMSs_rating(1:nbest,ii),ii),bins,'Edgecolor','none','FaceColor','r')
        yyaxis right
        histogram(CCE_cRMSs(idxMaxRep_CCE,ii),bins,'Normalization','Probability','Edgecolor','none','FaceColor','y')
        if ii==1
            ylabel('cRMS')
        end
        subplot(3,length(theversions),2*length(theversions)+ii)
        yyaxis left
        bins=linspace(prctile(CCE_STDs(:,ii),minpct),prctile(CCE_STDs(:,ii),maxpct),nbins);
        histogram(CCE_STDs(:,ii),bins,'Edgecolor','none','FaceColor','b')
        hold on
        histogram(CCE_STDs(CCE_STDs_rating(1:nbest,ii),ii),bins,'Edgecolor','none','FaceColor','r')
        yyaxis right
        histogram(CCE_STDs(idxMaxRep_CCE,ii),bins,'Normalization','Probability','Edgecolor','none','FaceColor','y')
        if ii==1
            ylabel('STD')
        end
    end
    sgtitle('Normalized statistics for CCE')

    fig2=figure('Color','w','Name','Normalized stats for HOT');
    for ii=1:length(theversions)
        subplot(3,length(theversions),ii)
        yyaxis left
        bins=linspace(-1,1,nbins);
        histogram(HOT_CORs(:,ii),bins,'Edgecolor','none','FaceColor','b')
        hold on
        histogram(HOT_CORs(HOT_CORs_rating(1:nbest,ii),ii),bins,'Edgecolor','none','FaceColor','r')
        yyaxis right
        histogram(HOT_CORs(idxMaxRep_HOT,ii),bins,'Normalization','Probability','Edgecolor','none','FaceColor','y')
        subtitle(theversions{ii})
        if ii==1
            ylabel('CORRELATION')
        end
        subplot(3,length(theversions),length(theversions)+ii)
        yyaxis left
        bins=linspace(prctile(HOT_cRMSs(:,ii),minpct),prctile(HOT_cRMSs(:,ii),maxpct),nbins);
        histogram(HOT_cRMSs(:,ii),bins,'Edgecolor','none','FaceColor','b')
        hold on
        histogram(HOT_cRMSs(HOT_RMSs_rating(1:nbest,ii),ii),bins,'Edgecolor','none','FaceColor','r')
        yyaxis right
        histogram(HOT_cRMSs(idxMaxRep_HOT,ii),bins,'Normalization','Probability','Edgecolor','none','FaceColor','y')

        if ii==1
            ylabel('cRMS')
        end
        subplot(3,length(theversions),2*length(theversions)+ii)
        yyaxis left
        bins=linspace(prctile(HOT_STDs(:,ii),minpct),prctile(HOT_STDs(:,ii),maxpct),nbins);
        histogram(HOT_STDs(:,ii),bins,'Edgecolor','none','FaceColor','b')
        hold on
        histogram(HOT_STDs(HOT_STDs_rating(1:nbest,ii),ii),bins,'Edgecolor','none','FaceColor','r')
        yyaxis right
        histogram(HOT_STDs(idxMaxRep_HOT,ii),bins,'Normalization','Probability','Edgecolor','none','FaceColor','y')
        if ii==1
            ylabel('STD')
        end
    end
    sgtitle('Normalized statistics for HOT')
end

%% FUNCTION Load the data

function [CORs,cRMSs,STDs,RMSds]=loadTheData(The,theversions)

for ii=1:length(theversions)
    corr=(The.stattable.(['r_',theversions{ii}]));
    crms=(sqrt(The.stattable.(['cRMSd_',theversions{ii}])));
    std_t=(The.stattable.(['std_',theversions{ii}]));
    std_r=(The.stattable.(['std_',theversions{ii},'_ref']));
    rmsd=(The.stattable.(['RMSd_',theversions{ii}]));

    %     corr=real(The.stattable.(['r_',theversions{ii}]));
    %     crms=real(sqrt(The.stattable.(['cRMSd_',theversions{ii}])));
    %     std_t=real(The.stattable.(['std_',theversions{ii}]));
    %     std_r=real(The.stattable.(['std_',theversions{ii},'_ref']));

    % normalize data for plot
    CORs(:,ii)=corr;
    cRMSs(:,ii)=crms./std_r; %normalized
    STDs(:,ii)=std_t./std_r; %normalized
    RMSds(:,ii)=rmsd;

end
end
%% FUNCTION
function [max_cRMS,max_std,min_cor,max_RMSd]=defineAcceptableRange(whichone,theversions)
switch whichone
    case 'CCE'
        % Define acceptable range for CCE:
        max_cRMS_HCpico=3.6;    max_std_HCpico=3.9;  min_cor_HCpico=-0.4;   max_RMSd_HCpico=1.6;
        max_cRMS_HCnano=1;      max_std_HCnano=1.6;  min_cor_HCnano=0.5;   max_RMSd_HCnano=0.4;
        max_cRMS_ACpico=1.7;    max_std_ACpico=2.1;   min_cor_ACpico=0;   max_RMSd_ACpico=2.7;
        max_cRMS_ACnano=0.7;    max_std_ACnano=1.4;   min_cor_ACnano=0.7;   max_RMSd_ACnano=0.4;
        max_cRMS_ACmicro=0.5;   max_std_ACmicro=1.4;  min_cor_ACmicro=0.96;  max_RMSd_ACmicro=1.0;
        max_cRMS_TOTpico=2;     max_std_TOTpico=2.0;    min_cor_TOTpico=0.2;    max_RMSd_TOTpico=0.4;
        max_cRMS_TOTnano=0.7;   max_std_TOTnano=1.4;    min_cor_TOTnano=0.7;    max_RMSd_TOTnano=0.3;


    case 'HOT'
        % Define acceptable range for HOT:
        max_cRMS_HCpico=0.6;    max_std_HCpico=1.5;  min_cor_HCpico=0.88;   max_RMSd_HCpico=0.03;
        max_cRMS_HCnano=1.1;    max_std_HCnano=1.4;  min_cor_HCnano=0.2;   max_RMSd_HCnano=0.77;
        max_cRMS_ACpico=0.4;    max_std_ACpico=1.4;   min_cor_ACpico=0.96;   max_RMSd_ACpico=0.21;
        max_cRMS_ACnano=0.8;    max_std_ACnano=1.5;   min_cor_ACnano=0.9;   max_RMSd_ACnano=0.39;
        max_cRMS_ACmicro=1.7;   max_std_ACmicro=2.6;  min_cor_ACmicro=0.7;  max_RMSd_ACmicro=1.03;
        max_cRMS_TOTpico=0.4;   max_std_TOTpico=1.2;    min_cor_TOTpico=0.97;    max_RMSd_TOTpico=0.04;
        max_cRMS_TOTnano=0.7;   max_std_TOTnano=1.6;    min_cor_TOTnano=0.9;    max_RMSd_TOTnano=0.47;
end

if isequal({'TOTpico','ACnano','ACmicro'},theversions)
    % TOTpico, TOTnano, ACmicro
    max_cRMS=[max_cRMS_TOTpico,max_cRMS_ACnano,max_cRMS_ACmicro];
    max_std=[max_std_TOTpico,max_std_ACnano,max_std_ACmicro];
    min_cor=[min_cor_TOTpico,min_cor_ACnano,min_cor_ACmicro];
    max_RMSd=[max_RMSd_TOTpico,max_RMSd_ACnano,max_RMSd_ACmicro];
elseif isequal({'TOTpico','TOTnano','ACmicro'},theversions)
    % TOTpico, TOTnano, ACmicro
    max_cRMS=[max_cRMS_TOTpico,max_cRMS_TOTnano,max_cRMS_ACmicro];
    max_std=[max_std_TOTpico,max_std_TOTnano,max_std_ACmicro];
    min_cor=[min_cor_TOTpico,min_cor_TOTnano,min_cor_ACmicro];
    max_RMSd=[max_RMSd_TOTpico,max_RMSd_TOTnano,max_RMSd_ACmicro];
elseif isequal({'TOTpico','TOTnano','ACpico','ACnano','ACmicro'},theversions)
    % TOTpico, TOTnano, ACmicro
    max_cRMS=[max_cRMS_TOTpico,max_cRMS_TOTnano,max_cRMS_ACpico,max_cRMS_ACnano,max_cRMS_ACmicro];
    max_std=[max_std_TOTpico,max_std_TOTnano,max_std_ACpico,max_std_ACnano,max_std_ACmicro];
    min_cor=[min_cor_TOTpico,min_cor_TOTnano,min_cor_ACpico,min_cor_ACnano,min_cor_ACmicro];
    max_RMSd=[max_RMSd_TOTpico,max_RMSd_TOTnano,max_RMSd_ACpico,max_RMSd_ACnano,max_RMSd_ACmicro];
elseif isequal({'ACpico','ACnano','ACmicro'},theversions)
    % TOTpico, TOTnano, ACmicro
    max_cRMS=[max_cRMS_ACpico,max_cRMS_ACnano,max_cRMS_ACmicro];
    max_std=[max_std_ACpico,max_std_ACnano,max_std_ACmicro];
    min_cor=[min_cor_ACpico,min_cor_ACnano,min_cor_ACmicro];
    max_RMSd=[max_RMSd_ACpico,max_RMSd_ACnano,max_RMSd_ACmicro];
elseif isequal({'TOTpico','HC_Hnano','ACnano','ACmicro'},theversions)
    % TOTpico, HCnano, ACnano, ACmicro
    max_cRMS=[max_cRMS_TOTpico,max_cRMS_HCnano,max_cRMS_ACnano,max_cRMS_ACmicro];
    max_std=[max_std_TOTpico,max_std_HCnano,max_std_ACnano,max_std_ACmicro];
    min_cor=[min_cor_TOTpico,min_cor_HCnano,min_cor_ACnano,min_cor_ACmicro];
    max_RMSd=[max_RMSd_TOTpico,max_RMSd_HCnano,max_RMSd_ACnano,max_RMSd_ACmicro];
elseif isequal({'HC_Hbac','HC_Hnano','ACpico','ACnano','ACmicro'},theversions)
    % % HCpico, HCnano, ACpico, ACnano, ACmicro
    max_cRMS=[max_cRMS_HCpico,max_cRMS_HCnano,max_cRMS_ACpico,max_cRMS_ACnano,max_cRMS_ACmicro];
    max_std=[max_std_HCpico,max_std_HCnano,max_std_ACpico,max_std_ACnano,max_std_ACmicro];
    min_cor=[min_cor_HCpico,min_cor_HCnano,min_cor_ACpico,min_cor_ACnano,min_cor_ACmicro];
    max_RMSd=[max_RMSd_HCpico,max_RMSd_HCnano,max_RMSd_ACpico,max_RMSd_ACnano,max_RMSd_ACmicro];
else
    disp('boundaries not defined for this case')
    return
end
end

%% FUNCTION plot Histogram of params to compare

function plotHistParamCompare(randParam,thesePoints,nbins,normalization,plotall)
allparameters={'r^*_d','y','r^*_l','aF','c_{passive}','\alpha_{Max}','\alpha_R','mortHTL',...
    'rhoC:N','fracHTL:N','\epsilon_L','\alpha_N','\epsilon_F','cF','\beta','\sigma',...
    'remin2','reminF','rho','mort2','rem_{POM}','vel_1','vel_2','mHTL'};
% figure;
switch normalization
    case 'no'
        mynorm='count';
    case 'yes'
        mynorm='probability';
end

for i=1:24
    subplot(5,5,i)
    h=histogram(randParam(:,i),nbins,'Normalization',mynorm,'EdgeColor','w');
    binedges=h.BinEdges;
    myxlim=xlim;
    if plotall
    hold on
    end
    for ii=1:size(thesePoints,2)
        a=thesePoints(:,ii)>0;
        histogram(randParam(thesePoints(a,ii),i),binedges,'Normalization',mynorm,'EdgeColor','w');
        xlim(myxlim)
        hold on;
    end

    subtitle(allparameters{i})
end
end

%% FUNCTION sort the data for good or bad

function [CORs_rating,RMSs_rating,STDs_rating]=sortbest(CORs,RMSs,STDs)
% Find the rating for each of the 3 catagories. Remove the ones that are
% NaN.

% correlation coefficient: best is 1, worst is -1
% (here NaN is on top so they need to be sorted out)
[tmp_val,CORs_rating]=sort(CORs,'descend');
rm=isnan(tmp_val);
CORs_rating(rm)=NaN;

% centered RMS: best is 0, worst is high
[tmp_val,RMSs_rating]=sort(RMSs,'ascend');
rm=isnan(tmp_val);
RMSs_rating(rm)=NaN;

% standard deviation: best is 1, worst is high or low
[tmp_val,STDs_rating]=sort(abs(1-STDs),'ascend');
rm=isnan(tmp_val);
STDs_rating(rm)=NaN;
end
%% FUNCTION Plot taylor density map

function plotTaylorDensity(theversions,theversions_name,CCE_CORs,CCE_RMSs,CCE_STDs)
CORref=1;STDref=1;RMSref=0;
load('BlueYellow.mat')
tilelist=[1,2,4,5,6];
fig=figure('Name','Taylor density plots','Color','w','units','centimeters','position',[15,6,11.2,12]);
tiledlayout(3,2,'TileIndexing','columnmajor','TileSpacing','none','Padding','tight')
rmax=4;
for ii=1:length(theversions)
    nexttile(tilelist(ii))
    %find the ones that are not nan:
    idx=find(~isnan(CCE_CORs(:,ii)) & ~isnan(CCE_RMSs(:,ii)) & ~isnan(CCE_STDs(:,ii)));
    taylordiag_density([STDref;CCE_STDs(idx,ii)],[RMSref;CCE_RMSs(idx,ii)],[CORref;CCE_CORs(idx,ii)],...
        'plotScatter',1,'limSTD',rmax,'titleCOR',0,'titleSTD',0,'titleRMS',0,'tickRMS',1:1:rmax,'tickSTD',1:1:rmax);%,'showlabelsCOR',0);
    %     subtitle(theversions_name{ii})
    xlim([-4.5 4.5])
    ylim([0 5])
    colormap(BlueYellow)
    %     plotBoundaries(ii,'CCE',theversions,rmax)
end
exportgraphics(fig,'TaylorDensity_all.png','Resolution',1200)

% create smaller legend:
% fig2=figure('Name','Taylor density plots','Color','w','units','centimeters','position',[15,6,7,8]);
% tiledlayout(3,2,'TileIndexing','columnmajor','TileSpacing','none','Padding','tight')
% nexttile(3)
% idx=[];
% taylordiag_density([STDref;CCE_STDs(idx,ii)],[RMSref;CCE_RMSs(idx,ii)],[CORref;CCE_CORs(idx,ii)],...
%     'plotScatter',1,'limSTD',rmax,'titleCOR',0,'titleSTD',0,'titleRMS',0,'tickRMS',1:1:rmax,'tickSTD',1:1:rmax);%,'showlabelsCOR',0);
% xlim([-4.5 4.5])
% ylim([0 5])
% colormap(BlueYellow)
% exportgraphics(fig2,'TaylorDensityLegend.pdf')
end
%%

function plotTaylorDensityWPOINTS(theversions,theversions_name,CCE_CORs,CCE_RMSs,CCE_STDs)
idxCORs=CCE_CORs(end-4:end,:);
idxSTDs=CCE_STDs(end-4:end,:);
CORref=1;STDref=1;RMSref=0;
load('BlueYellow.mat')
fig=figure('Name','Taylor density plots','Color','w','units','centimeters','position',[15,5,11.2,17]);
tiledlayout(length(theversions),1,'TileSpacing','tight','TileIndexing','columnmajor')
rmax=3;
for ii=1:length(theversions)
    nexttile
    if ii==3
        nexttile
    end
    %find the ones that are not nan:
    idx=find(~isnan(CCE_CORs(:,ii)) & ~isnan(CCE_RMSs(:,ii)) & ~isnan(CCE_STDs(:,ii)));
    taylordiag_density([STDref;CCE_STDs(idx,ii)],[RMSref;CCE_RMSs(idx,ii)],[CORref;CCE_CORs(idx,ii)],...
        'plotScatter',1,'limSTD',rmax,'titleCOR',0,'titleSTD',0,'titleRMS',0,'tickRMS',0.5:0.5:rmax);%,'showlabelsCOR',0);
%     subtitle(theversions_name{ii})
    xlim([-6 6])
    ylim([0 5])
    colormap(BlueYellow)
    %     plotBoundaries(ii,'CCE',theversions,rmax)
    hold on
    theta = real(acos(idxCORs(:,ii)));
    x=idxSTDs(:,ii).*cos(theta);
    y=idxSTDs(:,ii).*sin(theta);
    scatter(x,y,'.','r')

end
% exportgraphics(fig,'mygif.png')
end

function plotTaylorDensity_singleVertical(theversions,theversions_name,HOT_CORs,HOT_RMSs,HOT_STDs)
CORref=1;STDref=1;RMSref=0;
figure('Color','w','Name','Taylor density plots','units','centimeters','position',[15,5,11.2,22]);
t=tiledlayout(length(theversions),1);
t.TileSpacing = 'compact';
t.Padding = 'compact';
for ii=1:length(theversions)
    nexttile
    idx=find(~isnan(HOT_CORs(:,ii)) & ~isnan(HOT_RMSs(:,ii)) & ~isnan(HOT_STDs(:,ii))); % er denne linje nødvendig?
    taylordiag_density([STDref;HOT_STDs(idx,ii)],[RMSref;HOT_RMSs(idx,ii)],[CORref;HOT_CORs(idx,ii)],...
        'plotScatter',1,'limSTD',3,'titleCOR',0,'titleSTD',0,'titleRMS',0,'tickRMS',0.5:0.5:3);%,'showlabelsCOR',0);
    subtitle(theversions_name{ii})
end
end

function plotTaylorDensity_single(theversions,theversions_name,HOT_CORs,HOT_RMSs,HOT_STDs)
CORref=1;STDref=1;RMSref=0;
figure('Color','w','Name','Taylor density plots');
for ii=1:length(theversions)
    subplot(2,length(theversions),ii)
    subtitle(theversions_name{ii})

    subplot(2,length(theversions),length(theversions)+ii)
    idx=find(~isnan(HOT_CORs(:,ii)) & ~isnan(HOT_RMSs(:,ii)) & ~isnan(HOT_STDs(:,ii))); % er denne linje nødvendig?
    taylordiag_density([STDref;HOT_STDs(idx,ii)],[RMSref;HOT_RMSs(idx,ii)],[CORref;HOT_CORs(idx,ii)],...
        'plotScatter',1,'limSTD',3,'titleCOR',0,'titleSTD',0,'titleRMS',0,'tickRMS',0.5:0.5:3);%,'showlabelsCOR',0);
end
subplot(2,length(theversions),length(theversions)+1);hold on;ylabel('HOT')
end

%% FUNCTION Plot taylor with RMSd colored
function plotTaylorRMS(theversions,theversions_name,CCE_CORs,CCE_RMSs,CCE_STDs,HOT_CORs,HOT_RMSs,HOT_STDs,HOT_RMSds,CCE_RMSds)
CORref=1;STDref=1;RMSref=0;
figure('Color','w','Name','Taylor density plots');
for ii=1:length(theversions)
    subplot(2,length(theversions),ii)
    %find the ones that are not nan:
    idx=find(~isnan(CCE_CORs(:,ii)) & ~isnan(CCE_RMSs(:,ii)) & ~isnan(CCE_STDs(:,ii)));
    taylordiag_rms([STDref;CCE_STDs(idx,ii)],[RMSref;CCE_RMSs(idx,ii)],[CORref;CCE_CORs(idx,ii)],...
        'plotScatter',1,'limSTD',3,'titleCOR',0,'titleSTD',0,'titleRMS',0,'tickRMS',0.5:0.5:3,'RMSds',CCE_RMSds(idx,ii));%,'showlabelsCOR',0);
    subtitle(theversions_name{ii})

    subplot(2,length(theversions),length(theversions)+ii)
    idx=find(~isnan(HOT_CORs(:,ii)) & ~isnan(HOT_RMSs(:,ii)) & ~isnan(HOT_STDs(:,ii))); % er denne linje nødvendig?
    taylordiag_rms([STDref;HOT_STDs(idx,ii)],[RMSref;HOT_RMSs(idx,ii)],[CORref;HOT_CORs(idx,ii)],...
        'plotScatter',1,'limSTD',3,'titleCOR',0,'titleSTD',0,'titleRMS',0,'tickRMS',0.5:0.5:3,'RMSds',HOT_RMSds(idx,ii));%,'showlabelsCOR',0);
end
subplot(2,length(theversions),1);hold on;ylabel('CCE')
subplot(2,length(theversions),length(theversions)+1);hold on;ylabel('HOT')
end


%% FUNCTION Plot taylor as it is original
function plotTaylor_orig(theversions,CCE_CORs,CCE_RMSs,CCE_STDs,HOT_CORs,HOT_RMSs,HOT_STDs,isolate_CCE_idx,isolate_HOT_idx)
CORref=1;STDref=1;RMSref=0;
figure('Color','w','Name','Taylor density plots');
for ii=1:length(theversions)
    subplot(2,length(theversions),ii)
    %find the ones that are not nan:
    idx=~isnan(isolate_CCE_idx(:,1));
    taylordiag([STDref;CCE_STDs(idx,ii)],[RMSref;CCE_RMSs(idx,ii)],[CORref;CCE_CORs(idx,ii)],...
        'plotScatter',0,'limSTD',3,'titleCOR',0,'titleSTD',0,'titleRMS',0,'tickRMS',[0.5:0.5:3]);
    subtitle(theversions{ii})

    subplot(2,length(theversions),length(theversions)+ii)
    idx=~isnan(isolate_HOT_idx(:,1));
    taylordiag([STDref;HOT_STDs(idx,ii)],[RMSref;HOT_RMSs(idx,ii)],[CORref;HOT_CORs(idx,ii)],...
        'plotScatter',0,'limSTD',3,'titleCOR',0,'titleSTD',0,'titleRMS',0,'tickRMS',[0.5:0.5:3]);
end
subplot(2,length(theversions),1);hold on;ylabel('CCE')
subplot(2,length(theversions),length(theversions)+1);hold on;ylabel('HOT')
end
%% FUNCTION
function plotTaylor_ONE(CCE_CORs,CCE_RMSs,CCE_STDs)
CORref=1;STDref=1;RMSref=0;
taylordiag([STDref;CCE_STDs(:)],[RMSref;CCE_RMSs(:)],[CORref;CCE_CORs(:)],...
    'plotScatter',0,'limSTD',3,'titleCOR',0,'titleSTD',0,'titleRMS',0,'tickRMS',[0.5:0.5:3]);
end
%% FUNCTION
function plotHistogramStat(theversions,thename,CORs,RMSs,STDs,lim_COR,lim_RMS,lim_STD,minpct,maxpct,nbins)
for ii=1:length(theversions)
    subplot(3,length(theversions),ii)
    bins=linspace(-1,1,nbins);
    h=histogram(CORs(:,ii),bins,'Edgecolor','none','FaceColor','b');
    hold on;
    plot([lim_COR(1) lim_COR(1)],[h.Parent.YLim],':k',[lim_COR(2) lim_COR(2)],[h.Parent.YLim],':k')
    subtitle(theversions{ii})

    subplot(3,length(theversions),length(theversions)+ii)
    bins=linspace(prctile(RMSs(:,ii),minpct),prctile(RMSs(:,ii),maxpct),nbins);
    h=histogram(RMSs(:,ii),bins,'Edgecolor','none','FaceColor','b');
    hold on;plot([lim_RMS(1) lim_RMS(1)],[h.Parent.YLim],':k',[lim_RMS(2) lim_RMS(2)],[h.Parent.YLim],':k')

    subplot(3,length(theversions),2*length(theversions)+ii)
    bins=linspace(prctile(STDs(:,ii),minpct),prctile(STDs(:,ii),maxpct),nbins);
    h=histogram(STDs(:,ii),bins,'Edgecolor','none','FaceColor','b');
    hold on;plot([lim_STD(1) lim_STD(1)],[h.Parent.YLim],':k',[lim_STD(2) lim_STD(2)],[h.Parent.YLim],':k')
end
subplot(3,length(theversions),1);hold on;ylabel('CORRELATION')
subplot(3,length(theversions),length(theversions)+1);hold on;ylabel('cRMS')
subplot(3,length(theversions),2*length(theversions)+1);hold on;ylabel('STD')
sgtitle(['Normalized statistics for ',thename])
end
%% FUNCTION
function [theones,maxrep,idxMaxRep]=findFittingpoints(CORs,RMSs,STDs,lim_COR,lim_RMS,lim_STD,theversions)
% find interations where COR, RMS and STD is within acceptable range
theones=nan(size(CORs));
for ii=1:length(theversions)
    i=find((CORs(:,ii)>=lim_COR(ii,1) & CORs(:,ii)<=lim_COR(ii,2)) & (RMSs(:,ii)>=lim_RMS(ii,1) & RMSs(:,ii)<=lim_RMS(ii,2)) & (STDs(:,ii)>=lim_STD(ii,1) & STDs(:,ii)<=lim_STD(ii,2)))';
    theones(1:length(i),ii)=i;
end
tmp=theones(:);
nrRep=accumarray(tmp(~isnan(tmp)),1); % Hvor mange er der af hver. rækkenr i nrRep er simulationsnummer
maxrep=max(nrRep);
idxMaxRep=cell(1,maxrep);
for i=1:maxrep
    idxMaxRep{i}=find(nrRep==i); %idxMaxRep er tal svarende til simulationsnummer
end

end
%% FUNCTION
function idx_inside=findDatapoints(myversion,theversions,CORs,RMSs,STDs)
CORref=1;STDref=1;RMSref=0;
figure('Color','w','Name','Taylor density plots');
subplot(1,2,1)
%find the ones that are not nan:
idx=find(~isnan(CORs(:,myversion)) & ~isnan(RMSs(:,myversion)) & ~isnan(STDs(:,myversion))); % er denne linje nødvendig?
taylordiag_density([STDref;STDs(idx,myversion)],[RMSref;RMSs(idx,myversion)],[CORref;CORs(idx,myversion)],...
    'plotScatter',1,'limSTD',3,'titleCOR',0,'titleSTD',0,'titleRMS',0,'tickRMS',[0.5:0.5:3]);%,'showlabelsCOR',0);
subtitle(theversions{myversion})

roi = drawcircle;
theta=real(acos(CORs(idx,myversion)));
thex=STDs(idx,myversion).*cos(theta);
they=STDs(idx,myversion).*sin(theta);
[in,on] = inpolygon(thex,they,roi.Vertices(:,1),roi.Vertices(:,2));
a=find(in==1);
subplot(1,2,2)
taylordiag_density([STDref;STDs(idx(a),myversion)],[RMSref;RMSs(idx(a),myversion)],[CORref;CORs(idx(a),myversion)],...
    'plotScatter',1,'limSTD',3,'titleCOR',0,'titleSTD',0,'titleRMS',0,'tickRMS',[0.5:0.5:3]);%,'showlabelsCOR',0);
subtitle(theversions{myversion})
idx_inside=idx(a);
end
%% FUNCTION
function HOT_CCE_corr(HOT_names,CCE_names,theversions,theversions_name,HOT_RMSs,CCE_RMSs,HOT_STDs,CCE_STDs,HOT_CORs,CCE_CORs)
tmp_RMSs_HOT=nan(max(length(HOT_names),length(CCE_names)),length(theversions));
tmp_RMSs_CCE=tmp_RMSs_HOT;
tmp_STDs_HOT=tmp_RMSs_HOT;
tmp_STDs_CCE=tmp_RMSs_HOT;
tmp_CORs_HOT=tmp_RMSs_HOT;
tmp_CORs_CCE=tmp_RMSs_HOT;

for ii=1:length(HOT_names)
    name_CCE_tmp=['p_point_',num2str(2),char(extractAfter(HOT_names(ii),'point_1'))];
    CCE_idx_tmp=find(strcmp(name_CCE_tmp,CCE_names)==1);
    if ~isempty(CCE_idx_tmp)
        tmp_RMSs_HOT(ii,:)=HOT_RMSs(ii,:);
        tmp_RMSs_CCE(ii,:)=CCE_RMSs(CCE_idx_tmp,:);
        tmp_STDs_HOT(ii,:)=HOT_STDs(ii,:);
        tmp_STDs_CCE(ii,:)=CCE_STDs(CCE_idx_tmp,:);
        tmp_CORs_HOT(ii,:)=HOT_CORs(ii,:);
        tmp_CORs_CCE(ii,:)=CCE_CORs(CCE_idx_tmp,:);
    end
end
%%
figure('Color','w','units','normalized','outerposition',[0 0 1 1]);
for ii=1:length(theversions)
    subplot(3,length(theversions),ii);hold on;
    myidx=find(all(~isnan([tmp_RMSs_HOT(:,ii),tmp_RMSs_CCE(:,ii)]),2));
    dscatter(tmp_RMSs_HOT(myidx,ii),tmp_RMSs_CCE(myidx,ii));
    set(gca,'DataAspectRatio',[1 1 1])
    %     plot(tmp_RMSs_HOT(:,ii),tmp_RMSs_CCE(:,ii),'.')
    subtitle(['cRMS: ',theversions_name{ii}])
    %     ax=gca;
    %     maxax=max(ax.XLim(2),ax.YLim(2));
    %     minax=min(ax.XLim(1),ax.YLim(1));
    %     ax.XLim=([minax maxax]);
    %     ax.YLim=([minax maxax]);
    %     axis('square')
    if ii==1
        ylabel('CCE')
    end

    subplot(3,length(theversions),length(theversions)+ii);hold on;
    myidx=find(all(~isnan([tmp_STDs_HOT(:,ii),tmp_STDs_CCE(:,ii)]),2));
    dscatter(tmp_STDs_HOT(myidx,ii),tmp_STDs_CCE(myidx,ii));
    set(gca,'DataAspectRatio',[1 1 1])
    %     plot(tmp_STDs_HOT(:,ii),tmp_STDs_CCE(:,ii),'.')
    subtitle(['STD: ',theversions_name{ii}])
    %     ax=gca;
    %     maxax=max(ax.XLim(2),ax.YLim(2));
    %     minax=min(ax.XLim(1),ax.YLim(1));
    %     ax.XLim=([minax maxax]);
    %     ax.YLim=([minax maxax]);
    %     axis('square')
    if ii==1
        ylabel('CCE')
    end

    subplot(3,length(theversions),2*length(theversions)+ii);hold on;
    myidx=find(all(~isnan([tmp_CORs_HOT(:,ii),tmp_CORs_CCE(:,ii)]),2));
    dscatter(tmp_CORs_HOT(myidx,ii),tmp_CORs_CCE(myidx,ii));
    set(gca,'DataAspectRatio',[1 1 1])
    %     plot(tmp_CORs_HOT(:,ii),tmp_CORs_CCE(:,ii),'.')
    subtitle(['COR: ',theversions_name{ii}])
    %     ax=gca;
    %     maxax=max(ax.XLim(2),ax.YLim(2));
    %     minax=min(ax.XLim(1),ax.YLim(1));
    %     ax.XLim=([minax maxax]);
    %     ax.YLim=([minax maxax]);
    %     axis('square')
    xlabel('HOT');
    if ii==1
        ylabel('CCE')
    end

end
sgtitle('How well does HOT and CCE statistics correlate?')
end
%% FUNCTION
function plotSubplotStat(CCE_RMSs,theversions,theversions_name)
xvector=repmat(1:length(theversions),[1 length(theversions)]);
yvector=repmat(1:length(theversions),[length(theversions) 1]);
yvector=yvector(:)';
nbins=100;minpct=0;maxpct=99;
for i=1:length(theversions)^2
    if yvector(i)>=xvector(i)
        subplot(length(theversions),length(theversions),i)
        myidx=find(all(~isnan([CCE_RMSs(:,xvector(i)),CCE_RMSs(:,yvector(i))]),2));
        dscatter(CCE_RMSs(myidx,xvector(i)),CCE_RMSs(myidx,yvector(i)));
        %     plot(CCE_RMSs(:,xvector(i)),CCE_RMSs(:,yvector(i)),'.')
        %     ax=gca;maxax=max(ax.XLim(2),ax.YLim(2));minax=min(ax.XLim(1),ax.YLim(1));
        %     ax.XLim=([minax maxax]);
        %     ax.YLim=([minax maxax]);
        set(gca,'DataAspectRatio',[1 1 1])
        %     axis('square')
        if xvector(i)==yvector(i)
            bins=linspace(prctile(CCE_RMSs(:,xvector(i)),minpct),prctile(CCE_RMSs(:,xvector(i)),maxpct),nbins);
            histogram(CCE_RMSs(:,xvector(i)),bins,'Edgecolor','none')
            axis('square')
        end
        if xvector(i)==1
            ylabel(theversions_name{yvector(i)})
        end
        if yvector(i)==5
            xlabel(theversions_name{xvector(i)})
            %         set(gca,'xaxisLocation','top')
        end
    end
end
end
%% FUNCTION
function dispStatforFile(sort_idx_Hpico_CCE_reduced,nrRepeat_Hpico_CCE,CCE_names,thestats)
idx=find(nrRepeat_Hpico_CCE==max(nrRepeat_Hpico_CCE));
for i=1:length(idx)
    [~,thecolumns]=find(sort_idx_Hpico_CCE_reduced==idx(i));
    disp(strjoin([CCE_names(idx(i)),'giver ',num2str(max(nrRepeat_Hpico_CCE)),': ',strjoin({thestats{thecolumns}})]))
end
end
%% FUNCTION
function statvsrand(HOT_STDs,titname,HOT_randParam,n,theversions,theversions_name,allparameters)
for ii=1:length(theversions)
    figure;
    [~,RMS_idx]=sortrows(HOT_STDs(:,ii));
    idxnr=1:n;
    for i=1:size(HOT_randParam,2)
        subplot(4,6,i)
        plot(HOT_randParam(RMS_idx(idxnr),i),HOT_STDs(RMS_idx(idxnr),ii),'.')
        xlabel(allparameters{i})
        ylabel(titname)
    end
    sgtitle([titname,' for ',theversions_name{ii}])
end
end
%% FUNCTION
function [idx]=WhichAreTheBest(thisone,CCEfilenames,theversions,theversions_name,cRMSs,STDs,CORs,RMSds,names,thestats,max_cRMS,max_std,min_cor,max_RMSd,randParam,variables,biospecter,mass,idxremove_CCE,CCE,largeCCE)
CORref=1;STDref=1;RMSref=0;
%% reduce to acceptable limits:
nrbads_size=nan(length(names),length(theversions)); % how many stats are bad for TOTpico or the others
nrbads_size(:,:)=0;
for i=1:length(names)
    thebad_table{i,1}=i;
    howbad{i,1}=i;
    howgood{i,1}=i;
    nrbads(i,1)=i;

    howbad_sum(i,1)=i;
    %howbad_sum(i,2)=sqrt(sum(((RMSds(i,:)./sum(RMSds(i,:)))-(max_RMSd./sum(max_RMSd))).^2)./length(theversions));
    howbad_sum(i,2)=sqrt(sum((RMSds(i,:)-max_RMSd).^2)./length(theversions));
    howbad_sum(i,3)=sqrt(sum((CORs(i,:)-min_cor).^2)./length(theversions));
    
    tmp=find(cRMSs(i,:)>max_cRMS | isnan(cRMSs(i,:)));
    if ~isempty(tmp)
        thebad_table{i,2}=strjoin([theversions_name(tmp)]);
        howbad{i,2}=max(cRMSs(i,tmp)-max_cRMS(tmp));
        howbad{i,2}=max(cRMSs(i,tmp)-max_cRMS(tmp));
        nrbads(i,2)=length(tmp);
        nrbads_size(i,tmp)=nrbads_size(i,tmp)+1;
    end
    howgood{i,2}=cRMSs(i,:)-max_cRMS(:)';
    tmp=find(STDs(i,:)>max_std | isnan(STDs(i,:)));
    if ~isempty(tmp)
        thebad_table{i,3}=strjoin([theversions_name(tmp)]);
        howbad{i,3}=max(STDs(i,tmp)-max_std(tmp));
        nrbads(i,3)=length(tmp);
        nrbads_size(i,tmp)=nrbads_size(i,tmp)+1;
    end
    howgood{i,3}=STDs(i,:)-max_std(:)';
    tmp=find(CORs(i,:)<min_cor | isnan(CORs(i,:)));
    if ~isempty(tmp)
        thebad_table{i,4}=strjoin([theversions_name(tmp)]);
        howbad{i,4}=min(CORs(i,tmp)-min_cor(tmp));
        nrbads(i,4)=length(tmp);
        nrbads_size(i,tmp)=nrbads_size(i,tmp)+1;
    end
    howgood{i,4}=CORs(i,:)-min_cor(:)';
    tmp=find(RMSds(i,:)>max_RMSd | isnan(RMSds(i,:)));
    if ~isempty(tmp)
        thebad_table{i,5}=strjoin([theversions_name(tmp)]);
        howbad{i,5}=max(RMSds(i,tmp)-max_RMSd(tmp));
        nrbads(i,5)=length(tmp);
        nrbads_size(i,tmp)=nrbads_size(i,tmp)+1;
    end
    howgood{i,5}=RMSds(i,:)-max_RMSd(:)';
    howbad{i,6}=names(i);
    howgood{i,6}=names(i);
%     nrbads{i,6}=names(i);

end
tf = cellfun('isempty',howbad); % true for empty cells
howbad(tf) = {0};
a=sum(double(cell2mat(howgood(:,5))<0),2);
[ix,~]=find(isnan(cell2mat(howbad(:,2:5))));
[ix1,~]=find(idxremove_CCE==1);
idelete=unique([ix;ix1]);
howbad(idelete,:)=[];
howgood=array2table(howgood,'VariableNames',[{'filenr'},strcat('The worst in ',thestats),{'filename'}]);
howbad=array2table(howbad,'VariableNames',[{'filenr'},strcat('The worst in ',thestats),{'filename'}]);
thebad_table=array2table(thebad_table,'VariableNames',[{'filenr_repeat'},thestats]);

nrbads_size(idelete,:)=NaN;

% --------- How many stats fail for each of the size catagories? ----------
% This is evaluated with nrbads_size that display excatly that number. How
% many of them have zero failing?
idx_all_good=find(sum(nrbads_size,2)==0);

% how many fails for pico?
idx_pico_good=find(nrbads_size(:,3)==0);
% TheGood_PlotResult(names(idx_pico_good),CCEfilenames(idx_pico_good))

% how many fails for nano?
idx_nano_good=find(nrbads_size(:,4)==0);
% tmp=RMSds(idx_nano_good,:)
% idx_nano_good=idx_nano_good(tmp(:,1)<1)
% TheGood_PlotResult(names(idx_nano_good),CCEfilenames(idx_nano_good))

idx_micro_good=find(nrbads_size(:,5)==0);
% TheGood_PlotResult(names(idx_nano_good),CCEfilenames(idx_nano_good))
theRMS=RMSds-max_RMSd;
theRMS(idelete,:)=NaN;
theCOR=CORs-min_cor;
theCOR(idelete,:)=NaN;

[val_cor1,idx_cor1]=sort(theCOR(:,1),'descend','MissingPlacement','last');
[val_cor2,idx_cor2]=sort(theCOR(:,2),'descend','MissingPlacement','last');
[val_cor3,idx_cor3]=sort(theCOR(:,3),'descend','MissingPlacement','last');
[val_cor4,idx_cor4]=sort(theCOR(:,4),'descend','MissingPlacement','last');
[val_cor5,idx_cor5]=sort(theCOR(:,5),'descend','MissingPlacement','last');
[val_rms1,idx_rms1]=sort(theRMS(:,1),'ascend');
[val_rms2,idx_rms2]=sort(theRMS(:,2),'ascend');
[val_rms3,idx_rms3]=sort(theRMS(:,3),'ascend');
[val_rms4,idx_rms4]=sort(theRMS(:,4),'ascend');
[val_rms5,idx_rms5]=sort(theRMS(:,5),'ascend');

% cor1=idx_cor1(1:32591);
% cor2=idx_cor2(1:19285);
% cor3=idx_cor3(1:45306);
% cor4=idx_cor4(1:39518);
% cor5=idx_cor5(1:56);
% 
% rms1=idx_rms1(1:2650);
% rms2=idx_rms2(1:2758);
% rms3=idx_rms3(1:4669);
% rms4=idx_rms4(1:424);
% rms5=idx_rms5(1:2466);

cor1=idx_cor1(1:find(val_cor1<0,1));
cor2=idx_cor2(1:find(val_cor2<0,1));
cor3=idx_cor3(1:find(val_cor3<0,1));
cor4=idx_cor4(1:find(val_cor4<0,1));
cor5=idx_cor5(1:find(val_cor5<0,1));

rms1=idx_rms1(1:find(val_rms1>0,1));
rms2=idx_rms2(1:find(val_rms2>0,1));
rms3=idx_rms3(1:find(val_rms3>0,1));
rms4=idx_rms4(1:find(val_rms4>0,1));
rms5=idx_rms5(1:find(val_rms5>0,1));

k=accumarray([cor1;cor2;cor3;cor4;cor5;rms1;rms2;rms3;rms4;rms5],1);max(k)
largeCCE=find(~isnan(CCE.AC_pico_mean_bin_save(:,10)) & ~isnan(CCE.AC_nano_mean_bin_save(:,10)));
idx=find(k>=6);
if strcmp(thisone,'CCE')
idx=intersect(idx,largeCCE);
end
%bedste corr3
[~,iA]=intersect(cor3,idx);
cor3_reduced=cor3(iA);
cor3sorted=sort(iA);
bestcor3=cor3(cor3sorted(1));
% TheGood_PlotResult(names(bestcor3),CCEfilenames(bestcor3))

%bedste rms3
[~,iA]=intersect(rms3,idx);
rms3_reduced=rms3(iA);
rms3sorted=sort(iA);
bestrms3=rms3(rms3sorted(1));
% TheGood_PlotResult(names(bestrms3),CCEfilenames(bestrms3))

%bedste rms1 
[~,iA]=intersect(rms1,idx);
rms1_reduced=rms1(iA);
rms1sorted=sort(iA);
bestrms1=rms1(rms1sorted(1));
% TheGood_PlotResult(names(bestrms1),CCEfilenames(bestrms1))

nrApp=find(accumarray([rms3_reduced;rms1_reduced],1)==2);
nrApp=find(accumarray([rms3_reduced;rms1_reduced;cor3_reduced],1)==3);
idx=nrApp


if strcmp(thisone,'CCE')
    idx(7:9)=[];
    %idx=find(k==8);
    %idx=idx(1:5);
    %idx=[375;8109;9603;14939;19817(;26071;35231;49271)]; % original
    idx=[1703;5801;39986;41402;49833;60815;85751;91224;99032]; % review

end
if strcmp(thisone,'HOT')
    %     idx=[50517;12879;14203;28022;13809;6431;11788]; % original
    idx=[29534;37705;69585;74294;75576;95236]; % review
end

TheGood_PlotResult(names(idx),CCEfilenames(idx))
colors=linspace(200,256,5);
for itile=[1 2 4 5 6]
    nexttile(itile)
    ax=gca;
    ax.Children(6).Color=BlueYellow(1,:);
    for i=1:5
        ax.Children(i).Color=BlueYellow(colors(i),:);
    end
end
fignew=gcf;
% exportgraphics(fignew,'GoodSolutions_CCE.pdf')



figParamDist=figure('Name','paramters for good','Color','w','units','centimeters','position',[15,5,9,21]);
tiledlayout(23,1,'TileSpacing','tight','Padding','tight')
for i=[1:22,24]
nexttile
plot([min(randParam(:,i)) max(randParam(:,i))],[0 0],'k','LineWidth',1.5)
hold on
plot(randParam(idx,i),0.1,'.r','MarkerSize',15)
plot([min(randParam(idx,i)) max(randParam(idx,i))],[0.1 0.1],'b','LineWidth',2)
xlim([min(randParam(:,i)) max(randParam(:,i))]);ylim([0 0.5])
ax=gca;
ax.XTick=linspace(min(randParam(:,i)),max(randParam(:,i)),3);
ax.XMinorTick='on';
ax.TickDir='out';
ylabel(variables{i})
ax.XAxis.Exponent=0;  % don't use exponent
ax.XAxis.TickLabelFormat='%.3f';
ax.YColor='k';
box off
end
exportgraphics(figParamDist,['figParamDist',thisone,'_redone.pdf'])

%% Nutrients:
if strcmp(thisone,'CCE')
    load('Nmat_CCE.mat')
else
    load('Nmat_HOT.mat')
end
WOA_z=[5 15 27.5 45 65 87.5 117.5 160 222.50 310 435 610 847.5 1160 1542.5 1975 2450 2950 3450 3950 4450 4950 5450]';

N=CCE.N;
N(:,idelete)=NaN;

new_depth=linspace(1,400,400);
for ii=1:size(N,2)
    tmp=interp1(WOA_z(1:size(N,1)),N(:,ii),new_depth)';
    N_new(:,ii)=tmp;
end

xlines=linspace(0,1000,1000);
h1=nan(size(N_new,1),length(xlines)-1);
for i=1:size(N_new,1)
h=histogram(N_new(i,:),xlines,'Normalization','probability');
h1(i,:)=h.Values;
end
[x,y]=meshgrid(xlines(1:end-1),new_depth');
load('BlueYellow.mat')
figNutrient=figure('Name','Nutrients','Color','w','units','centimeters','position',[15,5,6,15]);
tiledlayout(1,1,'TileSpacing','tight','Padding','tight')
contourf(x,-y,h1,200,'LineStyle','none');caxis([0 0.05])
ylim([-300 -1]);xlim([0 600])
colormap(BlueYellow)
hold on
plot(Nmat(1:size(N,1),1),-WOA_z(1:size(N,1)),'.:w')
box off
ax=gca;
ax.XAxisLocation='top';
ax.YTickLabel=[];
ax.XTickLabel=[];
ax.TickDir='out';
ax.TickLength=[0.005 0.005];
exportgraphics(figNutrient,['Nutrients',thisone,'.png'],'Resolution',1200)


figure;
subplot(2,2,1);
bins=-1:0.01:1;
histogram(CCE_N.r_N,bins,'Normalization','probability','LineStyle','none')
hold on
histogram(CCE_N.r_N(idx),bins,'Normalization','probability','LineStyle','none')
set(gca,'yscale','log')
subplot(2,2,3)
bins=0:0.1:400;
histogram(CCE_N.rmsd_N,bins,'Normalization','probability','LineStyle','none');








figure;
nexttile;plot(Relevant_DataMatrix(:,3),Relevant_DataMatrix(:,1),'ok');xlabel('r*l');ylabel('r*d');%r*d fra r*l
nexttile;plot(Relevant_DataMatrix(:,7),Relevant_DataMatrix(:,2),'ok');xlabel('alphaR');ylabel('y');%y fra mHTL
nexttile;plot(Relevant_DataMatrix(:,16),Relevant_DataMatrix(:,5),'ok');xlabel('sigma');ylabel('cpassive');%cpassive fra sigma
nexttile;plot(Relevant_DataMatrix(:,14),Relevant_DataMatrix(:,6),'ok');xlabel('cF');ylabel('alphaMax');%alphaMax fra cF
nexttile;plot(Relevant_DataMatrix(:,3),Relevant_DataMatrix(:,7),'ok');xlabel('r*l');ylabel('alphaR');%alphaR fra r*l
nexttile;plot(Relevant_DataMatrix(:,6),Relevant_DataMatrix(:,9),'ok');xlabel('alphaMax');ylabel('rhoC:N');%rhoC:N fra alphaMax
nexttile;plot(Relevant_DataMatrix(:,11),Relevant_DataMatrix(:,10),'ok');xlabel('epsilonL');ylabel('fracHTL:N');%fracHTL:N fra epsilonL
nexttile;plot(Relevant_DataMatrix(:,2),Relevant_DataMatrix(:,12),'ok');xlabel('y');ylabel('alphaN');%alphaN fra y
nexttile;plot(Relevant_DataMatrix(:,7),Relevant_DataMatrix(:,13),'ok');xlabel('alphaR');ylabel('epsilonF');%epsilonF fra alphaR
nexttile;plot(Relevant_DataMatrix(:,6),Relevant_DataMatrix(:,15),'ok');xlabel('alphaMax');ylabel('beta');%beta fra alphaMax
nexttile;plot(Relevant_DataMatrix(:,3),Relevant_DataMatrix(:,16),'ok');xlabel('r*l');ylabel('sigma');%sigma fra r*L
nexttile;plot(Relevant_DataMatrix(:,5),Relevant_DataMatrix(:,18),'ok');xlabel('cPassive');ylabel('reminF');%reminF fra cPassive
nexttile;plot(Relevant_DataMatrix(:,3),Relevant_DataMatrix(:,20),'ok');xlabel('r*l');ylabel('mort2');%mort2 fra r*L
nexttile;plot(Relevant_DataMatrix(:,18),Relevant_DataMatrix(:,22),'ok');xlabel('reminF');ylabel('vel1');%vel1 fra reminF
nexttile;plot(Relevant_DataMatrix(:,3),Relevant_DataMatrix(:,24),'ok');xlabel('r*l');ylabel('mHTL');%mHTL fra r*L
nexttile;plot(Relevant_DataMatrix(:,10),Relevant_DataMatrix(:,17),'ok');xlabel('fracHTL:N');ylabel('remin2');%remin2 fra fracHTL:N
nexttile;plot(Relevant_DataMatrix(:,16),Relevant_DataMatrix(:,4),'ok');xlabel('sigma');ylabel('aF');%remin2 fra fracHTL:N
nexttile;plot(Relevant_DataMatrix(:,4),Relevant_DataMatrix(:,19),'ok');xlabel('aF');ylabel('rho');%remin2 fra fracHTL:N
nexttile;plot(Relevant_DataMatrix(:,22),Relevant_DataMatrix(:,21),'ok');xlabel('vel');ylabel('remPOM');%remin2 fra fracHTL:N

%%
fig1=figure('Name','paramters for good','Color','w','units','centimeters','position',[15,5,11.2,21]);
tiledlayout(12,2,'TileSpacing','tight','Padding','tight')
for i=[1:22,24]
    nexttile
    plot([min(randParam(:,i)) max(randParam(:,i))],[0 0],'k','LineWidth',1.5)
    hold on
    plot([min(randParam(idx,i)) max(randParam(idx,i))],[0 0],'b','LineWidth',2)
    xlim([min(randParam(:,i)) max(randParam(:,i))])
    %ylim([0 0])
    ax=gca;
    ax.XTick=[min(randParam(:,i)) min(randParam(idx,i)) max(randParam(idx,i)) max(randParam(:,i))];
    title(variables{i})
    ax.XAxis.Exponent=0;  % don't use exponent
    ax.XAxis.TickLabelFormat='%.2f';
    ax.YColor='w';
    box off
end
%%
thisdir=pwd;
cd(fullfile('..','..','External_data','Taylor and Landry 2018','Landry CEE 2004-2011'))
thefig=openfig('biospecter.fig');
cd(thisdir)
figure(thefig);hold on
subplot(1,2,1);hold on
axlim=xlim;
for i=1:length(idx)
    plot(mass(3:10),log10(biospecter(idx(i),3:10)),'-*')
end
xlim(axlim);
subplot(1,2,2);hold on
ax=gca;
BinEdges=ax.Children.BinEdges;
P=nan(size(biospecter,1),2);
for i=1:size(biospecter,1)
    ii=find(~isinf(log10(biospecter(i,:))));
    P(i,:)=polyfit(mass(ii),log10(biospecter(i,ii)),1);
end
ix=find(imag(P(:,1))~=0);
P(ix,:)=NaN;
yyaxis right
histogram(P(:,1),BinEdges,'Normalization','probability')
histogram(P(idx,1),BinEdges,'Normalization','probability')
%%

for i=1:length(names)
    %     logspecter(i,:)=log10(biospecter(56170,:));
    ii=find(~isinf(log10(biospecter(i,:))));
    P(i,:)=polyfit(mass(ii),log10(biospecter(i,ii)),1);
end



% Christians: 
for this_descriptor=1:length(theversions)
[sorted_cRMS(:,this_descriptor),sorted_idx_cRMS(:,this_descriptor)]=sortrows(cRMSs(idx,this_descriptor),'ascend','MissingPlacement','last');
[sorted_STD(:,this_descriptor),sorted_idx_STD(:,this_descriptor)]=sortrows(abs(STDs(idx,this_descriptor)-STDref),'ascend','MissingPlacement','last');
[sorted_COR(:,this_descriptor),sorted_idx_COR(:,this_descriptor)]=sortrows(CORs(idx,this_descriptor),'descend','MissingPlacement','last');
[sorted_RMSd(:,this_descriptor),sorted_idx_RMSd(:,this_descriptor)]=sortrows(RMSds(idx,this_descriptor),'ascend','MissingPlacement','last');
end
% [~,sortidx_RMSds]=sortrows(cRMSs(idx,1),'ascend','MissingPlacement','last');
% Relevant_DataMatrix=randParam(idx(sortidx_RMSds),:);
Relevant_DataMatrix=randParam(idx,:);
figure
[S,AX,BigAx,H,HAx] = plotmatrix(Relevant_DataMatrix);
for i = 1:length(AX)
set(AX, 'YTick', [], 'XTick', []);
end
for i = 1:length(AX)
ylabel(AX(i,1),variables{i},'Rotation',45,'HorizontalAlignment','right', 'FontWeight', 'bold');
xlabel(AX(end,i),variables{i},'Rotation',45,'HorizontalAlignment','right', 'FontWeight', 'bold');
end
Relevant_DataMatrix=randParam(idx,:);
Relevant_DataMatrix(:,23)=[];
variables2=variables;
variables2(23:24)=[];
ndendro=2;
CGobj = clustergram(Relevant_DataMatrix,'ColorMap',redbluecmap);
CGobj.Standardize = 'column';                      % Standardizing data along columns
CGobj.RowPDist = 'spearman';       %passed to pdist; 'spearman' in geochem is often used; 'mahalanobis' in oceanography
CGobj.ColumnPDist = 'spearman';    %passed to pdist;  'spearman' in geochem is often used; 'mahalanobis' in oceanography
CGobj.Linkage = 'average';
CGobj.ColumnLabels = variables2;
CGobj.Dendrogram = ndendro;
% CGobj.Cluster = 'row';
set(CGobj,'Linkage','complete','Dendrogram',ndendro)
HFfig=plot(CGobj);
% cm = struct('GroupNumber',{1,2},'Annotation',{'Time1','Time2'},...
% 'Color',{[1 1 0],[0.1 0.6 1]});
% set(CGobj,'ColumnGroupMarker',cm)
% HFfig=plot(CGobj);



thisrmsd=[2.45197985665164 0.815553757742078 1.3751342751239];
thisstd=[2.32288535746887 1.8640219016168 1.15731629655333];
thiscor=[0.866523924688227 0.949216102544467 0.86766966132097];

thisrmsd=[2.5 0.82 1.4];
thisstd=[2.323 1.9 1.6];
thiscor=[0.9 0.95 0.87];



dx=find(RMSds(:,1)<=RMSds(56170,1) & RMSds(:,2)<=RMSds(56170,2) & RMSds(:,3)<=RMSds(56170,3) & sum(isnan(RMSds(:,:)),2)==0 &...
STDs(:,1)<=STDs(56170,1) & STDs(:,2)<=STDs(56170,2) & STDs(:,3)<=STDs(56170,3) & sum(isnan(STDs(:,:)),2)==0 &...
CORs(:,1)>=CORs(56170,1) & CORs(:,2)>=CORs(56170,2) & CORs(:,3)>=CORs(56170,3) & sum(isnan(CORs(:,:)),2)==0)

%how many bad ones do they have?
a=sum(nrbads(:,2:5),2);
%do any of them only have one bad?
idx=find(a==1);
% idx=find(k==max(k));
idx=[3038;5062;14922];
% 14922 er nogenlunde men skiller sig lidt ud på en højere rms i TOTpico
idx=[41426;3038;5062;26970;27317]; %GOOD for CCE !!!
% Taylor de gode:
figure('Color','w','Name','Taylor plot of acceptable files');
for ii=1:length(theversions)
    subplot(1,length(theversions),ii)
    %find the ones that are not nan:
    taylordiag([STDref;STDs(idx,ii)],[RMSref;cRMSs(idx,ii)],[CORref;CORs(idx,ii)],...
        'plotScatter',1,'limSTD',3,'titleCOR',0,'titleSTD',0,'titleRMS',0,'tickRMS',0.5:0.5:3,'theLabels',idx);%,'showlabelsCOR',0);
    subtitle(theversions_name{ii})
    plotBoundaries(ii,thisone,theversions,3)
end

for ii=1:22
    for iii=ii+1:24
        [a, sa, cov, r(ii,iii)] = linfit(randParam_reduced(:,ii),randParam_reduced(:,iii),0);
    end
end



% reduce the randParam to the good ones
randParam_reduced=randParam(idx,:);

figure;
plotall=true;
plotHistParamCompare(randParam,idx,30,'yes',plotall)


TheGood_PlotResult(names(idx),CCEfilenames(idx))

%%% check which is good for pico OR nano OR micro
mysize=2; %1==pico 2==nano 3==micro
[val_good,idx_goodrmsds]=sort(thegood_rmsds(:,mysize));
figure;plot(1:60,val_good,'-ok')
%brug figuren til at udvælge dem der er bedst og som skal plottes
mynum2plot=8;
idx2=idx;
idx2(1:mynum2plot,2)=idx(idx_goodrmsds(1:mynum2plot));
figure;
plotall=false;
plotHistParamCompare(randParam,idx2,30,'yes',plotall)


%%% check if there from parameters are sub areas:
myParam2use=2;
mymin=0.4;
a=find(randParam_reduced(:,myParam2use)<mymin);
idx2=idx;
idx2(1:length(a),2)=idx(a);

figure;
plotHistParamCompare(randParam./randParam(:,3),idx,30,'yes')


figure;
n=1;
for ii=1:22
    for iii=[1:22,24]
        nexttile
        plot(randParam_reduced(:,ii),randParam_reduced(:,iii),'.b');hold on;
        xlim([min(randParam(:,ii)) max(randParam(:,ii))])
        ylim([min(randParam(:,iii)) max(randParam(:,iii))])
        xlabel(variables{ii})
        ylabel(variables{iii})
        %                 lsline;
        n=n+1;
        if n==10
            figure;
            n=1;
        end
    end

end
%%
figure;
tiledlayout(6,6)
n=1;
thiIs={[2],[7 13],[9 14 16],[6 10],[15],[9 10],[11 12 13 15],[14 18],[14],[18],[12 13 15 20],[13 15 18],[15 20],[16 24],[17],[21],[22],[],[20 22],[21]}
for ii=[1:17,19,20]
for iii=thiIs{ii}
nexttile
plot(randParam_reduced(:,ii),randParam_reduced(:,iii),'ob');hold on;
hold on
plot(randParam_reduced2(:,ii),randParam_reduced2(:,iii),'or');hold on;
% xlim([min(randParam(:,ii)) max(randParam(:,ii))])
% ylim([min(randParam(:,iii)) max(randParam(:,iii))])
xlabel(variables{ii})
ylabel(variables{iii})
%                 lsline;
n=n+1;
% if n==10
% figure;
% n=1;
end
end

%%

ii=[1:17,19,20];
iii={[2],[7 13],[9 14 16],[6 10],[15],[9 10],[11 12 13 15],[14 18],[14],[18],[12 13 15 20],[13 15 18],[15 20],[16 24],[17],[21],[22],[20 22],[21]}




% Plot parameter krydsplots
subarea=false
if length(idx)>1
    % run through them and see if there is a correlation
    figure;
    n=1;
    
    for ii=1:22
        for iii=[ii+1:22 24]
                nexttile
                plot(randParam_reduced(:,ii),randParam_reduced(:,iii),'.b');hold on;
                if subarea
                plot(randParam_reduced(idx_goodrmsds(1:mynum2plot),ii),randParam_reduced(idx_goodrmsds(1:mynum2plot),iii),'.r');hold on;
                end
                xlim([min(randParam(:,ii)) max(randParam(:,ii))])
                ylim([min(randParam(:,iii)) max(randParam(:,iii))])
                xlabel(variables{ii})
                ylabel(variables{iii})
%                 lsline;
                n=n+1;
                if n==10
                    figure;
                    n=1;
                end
        end
        
    end
end

if length(idx)>1
    % run through them and see if there is a correlation
    figure;
    n=1;
    
    for ii=1:22
        for iii=[ii+1:22 24]
                nexttile
                plot(randParam_reduced(idx_goodrmsds(1:mynum2plot),ii),randParam_reduced(idx_goodrmsds(1:mynum2plot),iii),'.');hold on;
                xlim([min(randParam(:,ii)) max(randParam(:,ii))])
                ylim([min(randParam(:,iii)) max(randParam(:,iii))])
                xlabel(variables{ii})
                ylabel(variables{iii})
%                 lsline;
                n=n+1;
                if n==10
                    figure;
                    n=1;
                end
        end
        
    end
end






combinedTable=[howbad thebad_table];
t = sortrows(howbad,[5,4,3]);
% writetable(t,'howbadarethey.xlsx','Sheet',1);
figure;
nbins=30;
plotHistParamCompare(randParam,cell2mat(table2cell(t(1:100,'filenr'))),nbins,'yes')






idx=find(k~=length(thestats));
theothers=howbad(idx,:);
t = sortrows(theothers,[5,4,3])



ok_cRMS=cRMSs;
ok_STD=STDs;
ok_COR=CORs;
ok_RMSd=RMSds;

for i=1:length(theversions) %HC_bac, HC_nano, ACpico, ACnano, ACmicro
    acceptable_cRMSs{i}=find(cRMSs(:,i)<max_cRMS(i));
    ix_cRMSs=cRMSs(:,i)>max_cRMS(i);
    ok_cRMS(ix_cRMSs,i)=NaN;

    acceptable_STDs{i}=find(abs(STDs(:,i)-STDref)<max_std(i));
    ix_STDs=abs(STDs(:,i)-STDref)>max_std(i);
    ok_STD(ix_STDs,i)=NaN;

    acceptable_CORs{i}=find(CORs(:,i)>min_cor(i));
    ix_CORs=CORs(:,i)<min_cor(i);
    ok_COR(ix_CORs,i)=NaN;

    acceptable_RMSds{i}=find(RMSds(:,i)<max_RMSd(i));
    ix_RMSds=RMSds(:,i)>max_RMSd(i);
    ok_RMSd(ix_RMSds,i)=NaN;

end
if 1==0
    plottable(theversions,theversions_name,cRMSs,ok_cRMS,STDs,ok_STD,CORs,ok_COR,RMSds,ok_RMSd)
    plotHistwithinRange(theversions,theversions_name,cRMSs,ok_cRMS,STDs,ok_STD,CORs,ok_COR,RMSds,ok_RMSd)
end
% find out how many stats (cRMSs, STD, COR and RMS) a single file can be good at
disp('For a single catagory, do any simulation have all good stats?')
for i=1:length(theversions)
    tmp={acceptable_cRMSs{i},acceptable_STDs{i},acceptable_CORs{i},acceptable_RMSds{i}};
    tmp=cell2mat(tmp');
    nrRepeat{i}=accumarray(tmp,1);
    if max(nrRepeat{i}==4)
        disp([theversions_name{i},' : yes'])
    else
        disp([theversions_name{i},' : no'])
    end
end
disp('Are there any simulations that are good in several catagories?')
themins=4; %(minimum number of statistics that has to be within boundaries. Here, All)
thegoodones=[];nr_Hpico=[];

for i=1:length(theversions)
    thegoodones=[thegoodones; find(nrRepeat{i}>=themins)];
    nr_Hpico=[nr_Hpico;sum(nrRepeat{i}==max(nrRepeat{i}))];
end

k=accumarray(thegoodones,1);
for i=1:length(k)
howbad_crms(i,:)=cRMSs(i,:)-max_cRMS;
howbad_std(i,:)=STDs(i,:)-max_std;
howbad_cor(i,:)=CORs(i,:)-min_cor;
howbad_rmsd(i,:)=cRMSs(i,:)-max_RMSd;
end
idx=1:length(k);
%save info on which is bad
%     thebad=nan(length(idx),5);


disp (['Max catagory a file is good in is: ',num2str(max(k)),' out of ',num2str(length(theversions))])
a=max(k);
if a>1
    disp('the good files are: ')
    idx=find(k==max(k));
%     thebad = array2table(thebad,"VariableNames",["idx","bad cRMS idx","bad STD idx","bad COR idx","bad RMS idx"]);
    if length(idx)>10
        disp(['Too many to show all (',num2str(length(idx)),'). THe first one is:'])
        aa=find(thegoodones==idx(1));
        disp(strjoin([names(idx(1))])); %' for ',theversions_name(X(aa))]))
    else
        for i=1:length(idx)
            aa=find(thegoodones==idx(i));
            disp(strjoin([names(idx(i))])); %,' for ',theversions_name(X(aa))]))
        end
    end
end
%taylor alle:
figure('Color','w','Name','Taylor density plots');
for ii=1:length(theversions)
subplot(1,length(theversions),ii)
%find the ones that are not nan:
taylordiag([STDref;STDs(:,ii)],[RMSref;cRMSs(:,ii)],[CORref;CORs(:,ii)],...
'plotScatter',1,'limSTD',3,'titleCOR',0,'titleSTD',0,'titleRMS',0,'tickRMS',0.5:0.5:3,'theLabels',1:size(CORs,1));%,'showlabelsCOR',0);
subtitle(theversions_name{ii})
plotBoundaries(ii,thisone,theversions,3)
end

% Taylor de gode:
figure('Color','w','Name','Taylor density plots');
for ii=1:length(theversions)
    subplot(1,length(theversions),ii)
    %find the ones that are not nan:
    taylordiag([STDref;STDs(idx,ii)],[RMSref;cRMSs(idx,ii)],[CORref;CORs(idx,ii)],...
        'plotScatter',1,'limSTD',3,'titleCOR',0,'titleSTD',0,'titleRMS',0,'tickRMS',0.5:0.5:3,'theLabels',idx);%,'showlabelsCOR',0);
    subtitle(theversions_name{ii})
    plotBoundaries(ii,thisone,theversions,3)
end

nbins=20;
plotHistParamCompare(randParam,idx,nbins,'yes')

% reduce the randParam to the good ones
randParam_reduced=randParam(idx,:);
if length(idx)>1
    % run through them and see if there is a correlation
    figure
    for ii=1:23
        for iii=ii+1:24
            p_bootstrp = bootstrp(5000,'lsqfitma',randParam_reduced(:,ii),randParam_reduced(:,iii));
            [h,x1_pvalue_bootstrp]=ztest(mean(p_bootstrp(:,2)),0,std(p_bootstrp(:,2)),'alpha',0.01);
            if h==1
                nexttile
                plot(randParam_reduced(:,ii),randParam_reduced(:,iii),'.');hold on;
                [a,r,sa(2),sa(1),xbar,ybar]=lsqfitma(randParam_reduced(:,ii),randParam_reduced(:,iii));
                x=min(randParam_reduced(:,ii)):(max(randParam_reduced(:,ii))-min(randParam_reduced(:,ii)))/20:max(randParam_reduced(:,ii));
                plot(x,a(2).*x+a(1),'-*r');
                xlabel(variables{ii})
                ylabel(variables{iii})
                lsline;
            end
        end
    end
end
% disp(['HC_{bac} max repeat: ',num2str(max(nrRepeat_Hpico_CCE)),' out of 4: '])
% dispStatforFile(sort_idx_Hpico_CCE_reduced,nrRepeat_Hpico_CCE,names,thestats)
%
% disp(['HC_{nano} max repeat: ',num2str(max(nrRepeat_Hnano_CCE)),' out of 4: '])
% dispStatforFile(sort_idx_Hnano_CCE_reduced,nrRepeat_Hnano_CCE,names,thestats)

end
function plotBoundaries(ii,whichone,theversions,rmax)
[max_cRMS,max_std,min_cor,~]=defineAcceptableRange(whichone,theversions);

% Plot Max correlation line
th  = acos(min_cor(ii));
cst = cos(th); snt = sin(th);
x_lines=rmax*cst; y_lines=rmax*snt;
line([0 x_lines],[0 y_lines],'linestyle','-','color','r')
% plot max rms circle and crossing point
max_rms=max_cRMS(ii);c_rms_x=1;c_rms_y=0;
th = 0:pi/50:2*pi;
xunit = max_rms * cos(th) + c_rms_x;
yunit = max_rms * sin(th) + c_rms_y;
plot(xunit, yunit,'r');
% plot max std circle
max_std=max_std(ii);c_std_x=0;c_std_y=0;
xunit = max_std * cos(th) + c_std_x;
yunit = max_std * sin(th) + c_std_y;
plot(xunit, yunit,'b');
end

function plottable(theversions,theversions_name,cRMSs,ok_cRMS,STDs,ok_STD,CORs,ok_COR,RMSds,ok_RMSd)
% Plot max/min values in table
for i=1:length(theversions)
    max_cRMSs(i,1)=max(cRMSs(~isnan(ok_cRMS(:,i)),i));
    max_STDs(i,1)=max(STDs(~isnan(ok_STD(:,i)),i));
    min_CORs(i,1)=min(CORs(~isnan(ok_COR(:,i)),i));
    max_RMSds(i,1)=max(RMSds(~isnan(ok_RMSd(:,i)),i));
end

T = table(max_cRMSs, max_STDs,min_CORs,max_RMSds,...
    'VariableNames',{'max(cRMS)','max(abs(STD-1))','min(COR)','max(RMSd)'},...
    'RowName',theversions_name);

disp(T)
end

function plotHistwithinRange(theversions,theversions_name,cRMSs,ok_cRMS,STDs,ok_STD,CORs,ok_COR,RMSds,ok_RMSd)
figure;
nbins=50;
for i=1:length(theversions)
    subplot(4,length(theversions),i)
    h=histogram(cRMSs(:,i),nbins,'Normalization','count','EdgeColor','w');
    hold on
    ix=cRMSs(~isnan(ok_cRMS(:,i)),i);
    histogram(ix,h.BinEdges,'Normalization','count','EdgeColor','w')
    subtitle(theversions_name{i})
    if i==1;ylabel('cRMS');end

    subplot(4,length(theversions),length(theversions)+i)
    h=histogram(STDs(:,i),nbins,'Normalization','count','EdgeColor','w');
    hold on
    ix=STDs(~isnan(ok_STD(:,i)),i);
    histogram(ix,h.BinEdges,'Normalization','count','EdgeColor','w')
    if i==1;ylabel('abs(STD-1)');end

    subplot(4,length(theversions),2*length(theversions)+i)
    h=histogram(CORs(:,i),nbins,'Normalization','count','EdgeColor','w');
    hold on
    ix=CORs(~isnan(ok_COR(:,i)),i);
    histogram(ix,h.BinEdges,'Normalization','count','EdgeColor','w')
    if i==1;ylabel('COR');end

    subplot(4,length(theversions),3*length(theversions)+i)
    h=histogram(RMSds(:,i),nbins,'Normalization','count','EdgeColor','w');
    hold on
    ix=RMSds(~isnan(ok_RMSd(:,i)),i);
    histogram(ix,h.BinEdges,'Normalization','count','EdgeColor','w')
    if i==1;ylabel('RMSd');end
end
end

function [CCE]=CombineDatasets(CCE,CCEnext)
CCE.AC_micro_mean_bin_save=[CCE.AC_micro_mean_bin_save;CCEnext.AC_micro_mean_bin_save];
CCE.AC_nano_mean_bin_save=[CCE.AC_nano_mean_bin_save;CCEnext.AC_nano_mean_bin_save];
CCE.AC_pico_mean_bin_save=[CCE.AC_pico_mean_bin_save;CCEnext.AC_pico_mean_bin_save];
CCE.HC_Hbac_tot_mean_bin_save=[CCE.HC_Hbac_tot_mean_bin_save;CCEnext.HC_Hbac_tot_mean_bin_save];
CCE.HC_Hnano_tot_mean_bin_save=[CCE.HC_Hnano_tot_mean_bin_save;CCEnext.HC_Hnano_tot_mean_bin_save];
CCE.TOT_nano_mean_bin_save=[CCE.TOT_nano_mean_bin_save;CCEnext.TOT_nano_mean_bin_save];
CCE.TOT_pico_mean_bin_save=[CCE.TOT_pico_mean_bin_save;CCEnext.TOT_pico_mean_bin_save];
CCE.randParam=[CCE.randParam;CCEnext.randParam];
CCE.stattable=[CCE.stattable;CCEnext.stattable];
CCE.sumB=[CCE.sumB;CCEnext.sumB];
CCE.thisfile_name=[CCE.thisfile_name,CCEnext.thisfile_name];
CCEnext.thissetup=repmat(CCEnext.thissetup,[length(CCEnext.thisfile_name) 1]);
CCE.thissetup=[CCE.thissetup;CCEnext.thissetup];
end

function [CCE_N]=CombineDatasetsN(CCE_N,CCEnext_N)
CCE_N.N=[CCE_N.N,CCEnext_N.N];
CCE_N.cRMSd_N=[CCE_N.cRMSd_N,CCEnext_N.cRMSd_N];
CCE_N.r_N=[CCE_N.r_N,CCEnext_N.r_N];
CCE_N.rmsd_N=[CCE_N.rmsd_N,CCEnext_N.rmsd_N];
CCE_N.std_N=[CCE_N.std_N,CCEnext_N.std_N];
CCE_N.thisfile_name=[CCE_N.thisfile_name,CCEnext_N.thisfile_name];
end

function plotHistNY(CCE_randParam,randParam_reduced,allparameters)
figure;
n=1;
for ii=[1:22 24]
    for iii=[ii+1:22 25]
        nexttile
        plot(randParam_reduced(:,ii),randParam_reduced(:,iii),'.b');hold on;
        xlim([min(CCE_randParam(:,ii)) max(CCE_randParam(:,ii))])
        ylim([min(CCE_randParam(:,iii)) max(CCE_randParam(:,iii))])
        xlabel(allparameters{ii})
        ylabel(allparameters{iii})
        %                 lsline;
        n=n+1;
        if n==10
            figure;
            n=1;
        end
    end
end
end

function plotRandomDensity(CCE,thisone,maxdepth,xx)
load('BlueYellow.mat');
ylog=logspace(-2,3,100);
N_acpico=nan(99,14);
N_acnano=N_acpico;
N_acmicro=N_acpico;
N_totpico=N_acpico;
N_totnano=N_acpico;
for i=2:14
    [N_acpico(:,i),~] = histcounts(CCE.AC_pico_mean_bin_save(:,i),ylog,'Normalization','probability');
    [N_acnano(:,i),~] = histcounts(CCE.AC_nano_mean_bin_save(:,i),ylog,'Normalization','probability');
    [N_acmicro(:,i),~] = histcounts(CCE.AC_micro_mean_bin_save(:,i),ylog,'Normalization','probability');
    [N_totpico(:,i),~] = histcounts(CCE.TOT_pico_mean_bin_save(:,i),ylog,'Normalization','probability');
    [N_totnano(:,i),~] = histcounts(CCE.TOT_nano_mean_bin_save(:,i),ylog,'Normalization','probability');
end
AC_new=logspace(0,4,100);
AC_center=(CCE.TayLanFig1a.AC_bin(2:14)+CCE.TayLanFig1a.AC_bin(1:14-1))./2;
for i=1:99
    N_acpico_new(i,:) = exp(interp1(log(AC_center),log(N_acpico(i,2:14)), log(AC_new)));
    N_acnano_new(i,:) = exp(interp1(log(AC_center),log(N_acnano(i,2:14)), log(AC_new)));
    N_acmicro_new(i,:) = exp(interp1(log(AC_center),log(N_acmicro(i,2:14)), log(AC_new)));
    N_totpico_new(i,:) = exp(interp1(log(AC_center),log(N_totpico(i,2:14)), log(AC_new)));
    N_totnano_new(i,:) = exp(interp1(log(AC_center),log(N_totnano(i,2:14)), log(AC_new)));
end
ixuse=find(~isnan(N_acmicro_new(1,:)));

% x1 = exp(interp1(log(AC_center),log(N_totnano(1,2:14)), log(AC_new)))

[xold,~]=meshgrid((CCE.TayLanFig1a.AC_bin(2:14)+CCE.TayLanFig1a.AC_bin(1:14-1))./2,(ylog(2:end)+ylog(1:end-1)./2));
[x,y]=meshgrid(AC_new,(ylog(2:end)+ylog(1:end-1)./2));

fig_ACpnm=openfig(fullfile('..','..','External_data','Taylor and Landry 2018',['TOTpn_and_ACpnm_',thisone,'_',num2str(maxdepth),'m.fig']));
figure(fig_ACpnm)
nexttile(4)
ax=gca;hh=ax.Children(1);hh.Color='w';
hold on
contourf(x(:,ixuse),y(:,ixuse),N_acpico_new(:,ixuse),20,'LineStyle','none');
uistack(hh,'top');
xlim([min(min(xold)) max(max(xold))])
colormap(BlueYellow)
set(gca,'TickDir','out')
clim([0 0.02])

nexttile(5)
ax=gca;hh=ax.Children(1);hh.Color='w';
hold on
contourf(x(:,ixuse),y(:,ixuse),N_acnano_new(:,ixuse),20,'LineStyle','none');
uistack(hh,'top');
xlim([min(min(xold)) max(max(xold))])
colormap(BlueYellow)
set(gca,'TickDir','out')
clim([0 xx])

nexttile(6)
ax=gca;hh=ax.Children(1);hh.Color='w';
hold on
contourf(x(:,ixuse),y(:,ixuse),N_acmicro_new(:,ixuse),20,'LineStyle','none');
uistack(hh,'top');
xlim([min(min(xold)) max(max(xold))])
colormap(BlueYellow)
set(gca,'TickDir','out')
clim([0 xx])

nexttile(1)
ax=gca;hh=ax.Children(1);hh.Color='w';
hold on
contourf(x(:,ixuse),y(:,ixuse),N_totpico_new(:,ixuse),20,'LineStyle','none');
uistack(hh,'top');
xlim([min(min(xold)) max(max(xold))])
colormap(BlueYellow)
set(gca,'TickDir','out')
clim([0 0.02])


nexttile(2)
ax=gca;hh=ax.Children(1);hh.Color='w';
hold on
contourf(x(:,ixuse),y(:,ixuse),N_totnano_new(:,ixuse),20,'LineStyle','none');
uistack(hh,'top');
xlim([min(min(xold)) max(max(xold))])
colormap(BlueYellow)
set(gca,'TickDir','out')
clim([0 xx])

exportgraphics(fig_ACpnm,'BiomassDensityPlot_HOT.png','Resolution',1200)
end
