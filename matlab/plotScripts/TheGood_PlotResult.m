function TheGood_PlotResult(thegoodfiles,thegoodtable)
%% Calculate statistics and plot;
if length(thegoodtable)>1
    thisfile=load(thegoodtable{1});
else
    thisfile=load(thegoodtable);
end

Plot_ACpnm=true;
Plot_HbacHnano=true;
Plot_TOTpn=true;
Plot_gasol=false;


thissite=(extractBetween(thegoodfiles{1},'p_point_','_restore'));
if isequal(thissite{1},'1')
    mysite='HOT';
elseif isequal(thissite{1},'2')
    mysite='CCE';

end
maxdepth=thisfile.maxdepth;

switch mysite
    case 'HOT'
        sitenr = 1;
        TayLanFig1a=load(fullfile('..','..','External_data','Taylor and Landry 2018','Taylor and Landry 2018 data for HOT',['HOT_TayLanFig1a_',num2str(maxdepth),'m_wSTD.mat']));
        TayLanFig7=load(fullfile('..','..','External_data','Taylor and Landry 2018','Taylor and Landry 2018 data for HOT',['HOT_TayLanFig7_',num2str(maxdepth),'m.mat']));
        GasolFig3=load(fullfile('..','..','External_data','Taylor and Landry 2018','Taylor and Landry 2018 data for HOT',['HOT_GasolFig3_',num2str(maxdepth),'m.mat']));
        totPN=load(fullfile('..','..','External_data','Taylor and Landry 2018','Taylor and Landry 2018 data for HOT',['HOT_totPN_',num2str(maxdepth),'m_wSTD.mat']));
        load(fullfile('..','..','External_data','Taylor and Landry 2018','daily_insolation_HOT.mat'))

    case 'CCE'
        sitenr = 2;
        TayLanFig1a=load(fullfile('..','..','External_data','Taylor and Landry 2018','Landry CEE 2004-2011',['CCE_TayLanFig1a_',num2str(maxdepth),'m_wSTD.mat']));
        TayLanFig7=load(fullfile('..','..','External_data','Taylor and Landry 2018','Landry CEE 2004-2011',['CCE_TayLanFig7_',num2str(maxdepth),'m.mat']));
        GasolFig3=load(fullfile('..','..','External_data','Taylor and Landry 2018','Landry CEE 2004-2011',['CCE_GasolFig3_',num2str(maxdepth),'m.mat']));
        totPN=load(fullfile('..','..','External_data','Taylor and Landry 2018','Landry CEE 2004-2011',['CCE_totPN_',num2str(maxdepth),'m_wSTD.mat']));
        load(fullfile('..','..','External_data','Taylor and Landry 2018','daily_insolation_CCE.mat'))
end
% fig_ACpnm=openfig(fullfile('..','..','External_data','Taylor and Landry 2018',['ACpnm_',mysite,'_',num2str(maxdepth),'m.fig']));
% fig_TOTpn=openfig(fullfile('..','..','External_data','Taylor and Landry 2018',['TOTpn_',mysite,'_',num2str(maxdepth),'m.fig']));

fig_ACpnm=openfig(fullfile('..','..','External_data','Taylor and Landry 2018',['TOTpn_and_ACpnm_',mysite,'_',num2str(maxdepth),'m.fig']));
fig_HCpn=openfig(fullfile('..','..','External_data','Taylor and Landry 2018',['HCpn_',mysite,'_',num2str(maxdepth),'m.fig']));



ACversion=thisfile.ACversion;
AC_fraction=thisfile.AC_fraction;
AC_fraction_F=thisfile.AC_fraction_F;
AC_fraction_O=thisfile.AC_fraction_O;

EinConv=4.57;
PARfrac=0.4;

% depth of statistic zone
mindepth=0;

% start of statistic time span
tstart=thisfile.tstart;
%load('devon.mat')
load('BlueYellow.mat');
devon=BlueYellow;
allcolors_nr=floor(1:256/length(thegoodfiles):256);
allcolors=devon(allcolors_nr,:);

for i=1:length(thegoodfiles)
    runfile=string(thegoodfiles(i));
    % same bins as Taylor & Landry
    bin_nr=1:14;
    AC_bin(bin_nr)=0;
    AC_pico_mean_bin=AC_bin;
    AC_pico_mean_bin(:)=NaN;

    AC_nano_mean_bin=AC_pico_mean_bin; AC_micro_mean_bin=AC_pico_mean_bin; AC_total_mean_bin=AC_pico_mean_bin;
    HC_phago_pico_mean_bin=AC_pico_mean_bin; HC_phago_nano_mean_bin=AC_pico_mean_bin; HC_phago_micro_mean_bin=AC_pico_mean_bin;
    TOT_pico_mean_bin=AC_pico_mean_bin; TOT_nano_mean_bin=AC_pico_mean_bin;

    AC_bin(2)=5;
    for ii=3:14
        AC_bin(ii)=AC_bin(ii-1)+0.5*AC_bin(ii-1);
    end
    AC_bin(end)=AC_bin(end)*2;


    %% Data load
    thisfile=load(thegoodtable{i});
    AC_fraction_F=thisfile.AC_fraction_F;
    AC_fraction_O=thisfile.AC_fraction_O;
    load(fullfile(thisfile.thissetup1,thisfile.thissetup,runfile))
    lyears=(simOutput.p.tEnd-simOutput.p.tSaveFrom)/365;

    n=simOutput.p.ixEnd(1)-(simOutput.p.ixStart(1)-1);
    thedepths=find(simOutput.z>mindepth & simOutput.z<maxdepth);

    %% remove small biomass for comparison
    % everything below 1.5 um has to be removed
    m_min=mass_function(1.5/2);
    m_max=simOutput.p.m(2+n).*sqrt(simOutput.p.m(2+n)/simOutput.p.m(1+n));
    for tt=1:size(simOutput.B,3)
        for dd=1:size(simOutput.B,1)
            [~,f] = calcBiomassRange(squeeze(simOutput.B(dd,1:n,tt)),simOutput.p.m(3:2+n),m_min,m_max);
            simOutput.B(dd,1:n,tt)=simOutput.B(dd,1:n,tt).*f;
        end
    end
    %% Calculate Biomasses (AC and HC)
    jsum=(simOutput.jLreal(thedepths,1:n,tstart:end)+simOutput.jDOC(thedepths,1:n,tstart:end)+simOutput.jFreal(thedepths,1:n,tstart:end));
    switch ACversion
        case 'ACv1'
            % autotroph carbon is a fraction of the total that is the ratio
            % of the rates
            AC=simOutput.B(thedepths,1:n,tstart:end).*...
                simOutput.jLreal(thedepths,1:n,tstart:end)./jsum;
            HC=simOutput.B(thedepths,1:n,tstart:end).*...
                (simOutput.jDOC(thedepths,1:n,tstart:end)+simOutput.jFreal(thedepths,1:n,tstart:end))./jsum;
            HC_osmo=simOutput.B(thedepths,1:n,tstart:end).*simOutput.jDOC(thedepths,1:n,tstart:end)./jsum;
            HC_phago=simOutput.B(thedepths,1:n,tstart:end).*simOutput.jFreal(thedepths,1:n,tstart:end)./jsum;
        case 'ACv2'
            jAH_ratio=simOutput.jLreal(thedepths,1:n,tstart:end)./(simOutput.jLreal(thedepths,1:n,tstart:end)+simOutput.jFreal(thedepths,1:n,tstart:end));
            A_cells=jAH_ratio;
            A_cells(jAH_ratio>AC_fraction)=1; %phototroph
            A_cells(A_cells~=1)=0;
            F_cells = double(~A_cells);

            AC=simOutput.B(thedepths,1:n,tstart:end).*A_cells.*(simOutput.jLreal(thedepths,1:n,tstart:end)+simOutput.jFreal(thedepths,1:n,tstart:end))./jsum...
                +simOutput.B(thedepths,1:n,tstart:end).*F_cells.*(simOutput.jLreal(thedepths,1:n,tstart:end))./jsum;
            HC=simOutput.B(thedepths,1:n,tstart:end).*simOutput.jDOC(thedepths,1:n,tstart:end)./jsum...
                +simOutput.B(thedepths,1:n,tstart:end).*F_cells.*simOutput.jFreal(thedepths,1:n,tstart:end)./jsum;
            HC_osmo=simOutput.B(thedepths,1:n,tstart:end).*simOutput.jDOC(thedepths,1:n,tstart:end)./jsum;
            HC_phago=simOutput.B(thedepths,1:n,tstart:end).*F_cells.*simOutput.jFreal(thedepths,1:n,tstart:end)./jsum;
        case 'ACv3'
                        I0=ins(1,1:2:end);
                        % 1pr day for 2 years
                        I0=repmat(I0,[1 lyears]).*EinConv.*PARfrac;
                        I0_5pct=I0.*0.05;
                        deepcells=simOutput.L(thedepths,:)<=I0_5pct;
                        deepcells=repmat(deepcells,[1 1 10]);
                        deepcells=permute(deepcells,[1 3 2]);

                        %definer ratios af AC til Phago (F) og osmo (O)
                        jAF_ratio=simOutput.jLreal(thedepths,1:n,tstart:end)./(simOutput.jLreal(thedepths,1:n,tstart:end)+simOutput.jFreal(thedepths,1:n,tstart:end));
                        jAO_ratio=simOutput.jLreal(thedepths,1:n,tstart:end)./(simOutput.jLreal(thedepths,1:n,tstart:end)+simOutput.jDOC(thedepths,1:n,tstart:end));
                        %Find celler hvor A er større end ratioen:
                        AF_cells=zeros(size(jAF_ratio));
                        AF_cells(jAF_ratio>AC_fraction_F)=1; %phototroph i stedet for Phago
                        AO_cells=zeros(size(jAO_ratio));
                        AO_cells(jAO_ratio>AC_fraction_O & deepcells)=1; %phototroph i stedet for osmo i dybe celler
                        %Og marker dem som ikke er
                        HF_cells = double(~AF_cells);
                        HO_cells = double(~AO_cells);
                        % Find celler som både er HF og HO
                        H_cells = HF_cells==HO_cells & HF_cells==1;
                        % Find celler som både er AF og AO
                        AFO_cells = AO_cells==AF_cells & AO_cells==1;
                        % fjern dem fra de enkelte
                        AF_cells(AFO_cells==1)=0;
                        AO_cells(AFO_cells==1)=0;
                        % Beregn AC som summen af de celler hvor alt er auto, hvor
                        % noget er osmo, hvor noget er phago og hvor der er lidt af
                        % alle:
                        AC=simOutput.B(thedepths,1:n,tstart:end).*AFO_cells...
                            +simOutput.B(thedepths,1:n,tstart:end).*AF_cells.*(simOutput.jLreal(thedepths,1:n,tstart:end)+simOutput.jFreal(thedepths,1:n,tstart:end))./jsum...
                            +simOutput.B(thedepths,1:n,tstart:end).*AO_cells.*(simOutput.jLreal(thedepths,1:n,tstart:end)+simOutput.jDOC(thedepths,1:n,tstart:end))./jsum...
                            +simOutput.B(thedepths,1:n,tstart:end).*H_cells.*simOutput.jLreal(thedepths,1:n,tstart:end)./jsum;
                        HC_phago=simOutput.B(thedepths,1:n,tstart:end).*HF_cells.*simOutput.jFreal(thedepths,1:n,tstart:end)./jsum;
                        HC_osmo=simOutput.B(thedepths,1:n,tstart:end).*HO_cells.*simOutput.jDOC(thedepths,1:n,tstart:end)./jsum;

                        HC=HC_osmo+HC_phago;

    end

    % sum biomass from size classes
    AC_tot=squeeze(sum(AC,2,'omitnan'));
    HC_tot=squeeze(sum(HC,2,'omitnan'));

    %% Data for figure Taylor and Landry 2A
    if Plot_ACpnm==true
        % Divide AC and HC_phago into pico, nano and micro plankton
        AC_pnm=nan(size(AC,1),size(AC,3),3);
        for tt=1:size(AC,3)
            for dd=1:size(AC,1)
                AC_pnm(dd,tt,:) = calcPicoNanoMicro(squeeze(AC(dd,:,tt)),simOutput.p.m(3:2+n));

            end
        end
        AC_pico=squeeze(AC_pnm(:,:,1));
        AC_nano=squeeze(AC_pnm(:,:,2));
        AC_micro=squeeze(AC_pnm(:,:,3));

        for ii=2:14
            AC_pico_mean_bin(ii)=mean(AC_pico(AC_tot>AC_bin(ii-1)&AC_tot<=AC_bin(ii)),'omitnan');
            AC_nano_mean_bin(ii)=mean(AC_nano(AC_tot>AC_bin(ii-1)&AC_tot<=AC_bin(ii)),'omitnan');
            AC_micro_mean_bin(ii)=mean(AC_micro(AC_tot>AC_bin(ii-1)&AC_tot<=AC_bin(ii)),'omitnan');

        end
        %make negative data 0:
        AC_pico_mean_bin(AC_pico_mean_bin<0)=0;
        AC_nano_mean_bin(AC_nano_mean_bin<0)=0;
        AC_micro_mean_bin(AC_micro_mean_bin<0)=0;

        AC_x=(AC_bin(1:14-1)+AC_bin(2:14))./2;
        figure(fig_ACpnm)
        nexttile(4);hold on
%         subplot(1,3,1);hold on
        loglog(AC_x,AC_pico_mean_bin(2:14),'-','Color',allcolors(i,:),'MarkerSize',4,'LineWidth',1,'MarkerFaceColor',allcolors(i,:))

%         subplot(1,3,2);hold on
        nexttile(5);hold on;
        loglog(AC_x,AC_nano_mean_bin(2:14),'-','Color',allcolors(i,:),'MarkerSize',4,'LineWidth',1,'MarkerFaceColor',allcolors(i,:))

%         subplot(1,3,3);hold on
        nexttile(6);hold on;
        loglog(AC_x,AC_micro_mean_bin(2:14),'-','Color',allcolors(i,:),'MarkerSize',4,'LineWidth',1,'MarkerFaceColor',allcolors(i,:))

    end

    %% Data for figure 7
    if Plot_HbacHnano==true
        % divide phagotrophs into pico, nano, and micro
        HC_phago_pnm=nan(size(HC_phago,1),size(HC_phago,3),3);
        for tt=1:size(AC,3)
            for dd=1:size(AC,1)
                HC_phago_pnm(dd,tt,:) = calcPicoNanoMicro(squeeze(HC_phago(dd,:,tt)),simOutput.p.m(3:2+n));
            end
        end
        HC_phago_pico=squeeze(HC_phago_pnm(:,:,1));
        HC_phago_nano=squeeze(HC_phago_pnm(:,:,2));
        HC_phago_micro=squeeze(HC_phago_pnm(:,:,3));
        % Define Hbac as osmo + pico phago
        HC_Hbac_tot=squeeze(sum(HC_osmo,2,'omitnan'))+HC_phago_pico;
        % Define Hnano as phago nano
        HC_Hnano_tot=HC_phago_nano;

        % Find the mean biomass within each AC bin
        AC_tot_mean_bin=nan(size(AC_bin));
        HC_Hbac_tot_mean_bin=nan(size(AC_bin));
        HC_Hnano_tot_mean_bin=nan(size(AC_bin));

        for ii=2:14
            AC_tot_mean_bin(ii)=mean(AC_tot(AC_tot>AC_bin(ii-1)&AC_tot<=AC_bin(ii)),'omitnan');
            HC_Hbac_tot_mean_bin(ii)=mean(HC_Hbac_tot(AC_tot>AC_bin(ii-1)&AC_tot<=AC_bin(ii)),'omitnan');
            HC_Hnano_tot_mean_bin(ii)=mean(HC_Hnano_tot(AC_tot>AC_bin(ii-1)&AC_tot<=AC_bin(ii)),'omitnan');
        end
        %make negative data 0:
        HC_Hbac_tot_mean_bin(HC_Hbac_tot_mean_bin<0)=0;
        HC_Hnano_tot_mean_bin(HC_Hnano_tot_mean_bin<0)=0;


        %% Plot figure fig 7
        figure(fig_HCpn)
        subplot(1,3,1);hold on;
        loglog((AC_bin(1:14-1)+AC_bin(2:14))./2, HC_Hbac_tot_mean_bin(2:14),'-o','Color',allcolors(i,:),'MarkerSize',4,'LineWidth',1,'MarkerFaceColor',allcolors(i,:));
        subplot(1,3,2);hold on;
        loglog((AC_bin(1:14-1)+AC_bin(2:14))./2, HC_Hnano_tot_mean_bin(2:14),'-o','Color',allcolors(i,:),'MarkerSize',4,'LineWidth',1,'MarkerFaceColor',allcolors(i,:));

    end
    %% calculate sum pico and sum nano
    if Plot_TOTpn==true
        TotB=simOutput.B(thedepths,1:n,tstart:end);


        TOT_pnm=nan(size(TotB,1),size(TotB,3),3);
        for tt=1:size(TotB,3)
            for dd=1:size(TotB,1)
                TOT_pnm(dd,tt,:) = calcPicoNanoMicro(squeeze(TotB(dd,:,tt)),simOutput.p.m(3:2+n));

            end
        end

        TOT_pico=squeeze(TOT_pnm(:,:,1));
        TOT_nano=squeeze(TOT_pnm(:,:,2));

        for ii=2:14
            TOT_pico_mean_bin(ii)=mean(TOT_pico(AC_tot>AC_bin(ii-1)&AC_tot<=AC_bin(ii)),'omitnan');
            TOT_nano_mean_bin(ii)=mean(TOT_nano(AC_tot>AC_bin(ii-1)&AC_tot<=AC_bin(ii)),'omitnan');
        end

        %make negative data 0:
        TOT_pico_mean_bin(TOT_pico_mean_bin<0)=0;
        TOT_nano_mean_bin(TOT_nano_mean_bin<0)=0;


        %% Plot figure Taylor and Landry 2A
%         figure(fig_TOTpn)
        figure(fig_ACpnm)
        AC_x=(AC_bin(1:14-1)+AC_bin(2:14))./2;
%         subplot(1,3,1);hold on;
        nexttile(1);hold on;
        loglog(AC_x,TOT_pico_mean_bin(2:14),'-','Color',allcolors(i,:),'MarkerSize',4,'LineWidth',1,'MarkerFaceColor',allcolors(i,:))

%         subplot(1,3,2);hold on;
        nexttile(2);hold on
        loglog(AC_x,TOT_nano_mean_bin(2:14),'-','Color',allcolors(i,:),'MarkerSize',4,'LineWidth',1,'MarkerFaceColor',allcolors(i,:))

    end
    if Plot_gasol==true
        %% Data for figure Gasol fig 3
        Find the mean biomass within each AC bin
        AC_tot_mean_bin=nan(size(AC_bin));
        HC_tot_mean_bin=nan(size(AC_bin));

        for ii=2:14
            AC_tot_mean_bin(ii)=mean(AC_tot(AC_tot>AC_bin(ii-1)&AC_tot<=AC_bin(ii)),'omitnan');
            HC_tot_mean_bin(ii)=mean(HC_tot(AC_tot>AC_bin(ii-1)&AC_tot<=AC_bin(ii)),'omitnan');
        end


        %% Plot figure Gasol fig 3
        figure;
        subplot(2,1,1)
        h=loglog(AC_tot,HC_tot,'ok');
        hold on
        h1=loglog((AC_bin(1:14-1)+AC_bin(2:14))./2, HC_tot_mean_bin(2:14),'-or');
        h2=loglog(GasolFig3.AC_tot,GasolFig3.HC_tot,'oy');
        h3=loglog((GasolFig3.AC_bin(1:14-1)+GasolFig3.AC_bin(2:14))./2, GasolFig3.HC_tot_mean_bin(2:14),'-ob');
        set(gca,'xscale','log','yscale','log')
        xlabel('AC_{tot} (\mugC/l)')
        ylabel('H-biomass (\mugC/l)')
        h=h(1);h(2)=h1(1);h(3)=h2(1);h(4)=h3(1);
        legend(h([1 2 3 4]),{'model result','binned mean','data','HOT binned mean'})

        subplot(2,1,2)
        h=plot(AC_tot,(HC_tot./AC_tot),'ok');
        hold on
        h1=loglog( (AC_bin(1:14-1)+AC_bin(2:14))./2, HC_tot_mean_bin(2:14)./AC_tot_mean_bin(2:14),'-or');
        h2=plot(GasolFig3.AC_tot,(GasolFig3.HC_tot./GasolFig3.AC_tot),'oy');
        h3=loglog( (GasolFig3.AC_bin(1:14-1)+GasolFig3.AC_bin(2:14))./2, GasolFig3.HC_tot_mean_bin(2:14)./GasolFig3.AC_tot_mean_bin(2:14),'-ob');
        set(gca,'xscale','log','yscale','log')
        xlabel('AC_{tot} (\mugC/l)')
        ylabel('Biomass ratio(H:A)')
        h=h(1);h(2)=h1(1);h(3)=h2(1);h(4)=h3(1);
        legend(h([1 2 3 4]),{'model result','binned mean','data','HOT binned mean'})


    end
    figure(fig_ACpnm)
%     figure(fig_TOTpn)
%     figure(fig_HCpn)
end

end




