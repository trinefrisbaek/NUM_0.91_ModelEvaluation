
% thisParam=allrandParam{ii};
thisParam='mortHTL_rand';
% thisParam='test';
folder=fullfile('results',thisParam);
fileList = dir(fullfile(folder, '*.mat'));
f={fileList.name};

if exist([thisParam,'SUM.mat'],'file')
    disp('load summerized data')
    load([thisParam,'SUM.mat'])
else
    disp('load data for the first time')
    [B,jLreal,jDOC,jF,mass,randParam,dz]=total_carbon(f,folder);
    save([thisParam,'SUM.mat'],'B','jLreal','jDOC','jF','mass','randParam','dz')
end


AutoB1=B.*(jLreal./(jLreal+jDOC+jF));

AA=squeeze(B(:,1,:));
figure;hold on
for ii=1:length(f)
plot(mass,AA(ii,:))
end
set(gca,'xscale','log')

% Bsum=squeeze(sum(AutoB1,3,'omitnan'));
% dz_mat=repmat(dz',[length(f) 1]);
% Bsum=Bsum.*dz_mat; % mg/m3 --> mg/m2
% Bsum=squeeze(sum(Bsum,2,'omitnan')); % sum all size classes


% [c,ia,ic]=unique(randParam,'rows');
% 
% Relevant_DataMatrix=nan(12,length(ia));
% for ii=1:length(ia)
%     idx=find(ic==ii);
%     Relevant_DataMatrix(:,ii)=Bsum(idx);
% end
% 
% figure
% [S,AX,BigAx,H,HAx] = plotmatrix(Relevant_DataMatrix(:,1:10)); 


function [B,jLreal,jDOC,jF,mass,randParam,dz]=total_carbon(f,folder)

tmp=load(fullfile(folder,string(f(1))));
mass=tmp.simOutput.p.m(tmp.simOutput.p.ixStart(1):tmp.simOutput.p.ixEnd(1));
dz=tmp.simOutput.dznom;
B=nan(length(f),15,25);
jLreal=B;
jDOC=B;
jF=B;

randParam=nan(length(f),length(tmp.simOutput.rand_Param));
for filenr=1:length(f)
    thisfile=string(f(filenr));
    simOutput=load(fullfile(folder,thisfile));
    simOutput=simOutput.simOutput;

    B_tmp=simOutput.B(:,1:25,:);
    %load rates
    jLreal_tmp=simOutput.jLreal(:,1:25,:);
    jDOC_tmp=simOutput.jDOC(:,1:25,:);
    jF_tmp=simOutput.jF(:,1:25,:);
    % take only last year;
    idx=simOutput.t>2*365;
    B_tmp=B_tmp(:,:,idx);
    jLreal_tmp=jLreal_tmp(:,:,idx);
    jDOC_tmp=jDOC_tmp(:,:,idx);
    jF_tmp=jF_tmp(:,:,idx);
    %make annual mean:
    B_tmp=squeeze(mean(B_tmp,3,'omitnan'));
    jLreal_tmp=squeeze(mean(jLreal_tmp,3,'omitnan'));
    jDOC_tmp=squeeze(mean(jDOC_tmp,3,'omitnan'));
    jF_tmp=squeeze(mean(jF_tmp,3,'omitnan'));


    full_tmp=[B_tmp;nan(15-size(B_tmp,1),25)];
    B(filenr,:,:)=full_tmp;
    full_tmp=[jLreal_tmp;nan(15-size(jLreal_tmp,1),25)];
    jLreal(filenr,:,:)=full_tmp;
    full_tmp=[jDOC_tmp;nan(15-size(jDOC_tmp,1),25)];
    jDOC(filenr,:,:)=full_tmp;
    full_tmp=[jF_tmp;nan(15-size(jF_tmp,1),25)];
    jF(filenr,:,:)=full_tmp;

    randParam(filenr,:)=simOutput.rand_Param;
end
end