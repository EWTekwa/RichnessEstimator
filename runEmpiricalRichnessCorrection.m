set(0,'DefaultAxesFontSize',22)
scrsz = get(0,'ScreenSize');
rng(6); %set random number generator seed

% Data2012=importdata('community_2012.csv');
% Data2013=importdata('community_2013.csv');
% Data2014=importdata('community_2014.csv');
% Data2015=importdata('community_2015.csv');
% Data2016=importdata('community_2016.csv');
% Data2017=importdata('community_2017.csv');
% Data2018=importdata('community_2018.csv');
% Data2019=importdata('community_2019.csv');
 
load('Data_all.mat'); %load .mat file containing the cell Data_all with species % cover data from Martone seaweed survey


%collect data from individual years (note: columns are species, rows are
%quadrats sequenced according to transects (10 consecutive rows are
%quadrats belonging to one transect)
Data_all={Data2012,Data2013,Data2014,Data2015,Data2016,Data2017,Data2018,Data2019};

NumYears=length(Data_all);
YearLabels=[12:19];
numResample=20; %number of resamples for subsampled estimates (20)
numBoot=500; %number of bootstramps for CI (500)


ksub=9; %number of transects
ksub_sub=10; %number of quadrats in each transect
Richness_raw=zeros(1,NumYears);
Chao1_all=zeros(1,NumYears);
Richness_apx=zeros(1,NumYears);
min_cover=zeros(1,NumYears);
Samples=zeros(1,NumYears);
for yr=1:NumYears
    Richness_raw(yr)=sum(sum(Data_all{yr}.data)>0);
    min_cover(yr)=min(Data_all{yr}.data(Data_all{yr}.data>0));
    Samples(yr)=size(Data_all{yr}.data,1);
end
mdl_Richness_raw=fitlm([1:NumYears],Richness_raw);
RichnessSlope_raw=mdl_Richness_raw.Coefficients.Estimate(2);
RichnessSlopeSD_raw=mdl_Richness_raw.Coefficients.SE(2)*sqrt(NumYears);
RichnessSlopeP_raw=mdl_Richness_raw.Coefficients.pValue(2);

expectedRichness_raw=zeros(numBoot,NumYears);
expectedChao1=zeros(numBoot,NumYears);
expectedChao2=zeros(numBoot,NumYears);
expectedACE=zeros(numBoot,NumYears);
expectedS_aj2=zeros(numBoot,NumYears);
expectedS_ij2=zeros(numBoot,NumYears);
expectedRichness_apx=zeros(numBoot,NumYears);
for yr=1:NumYears
    Data_all{yr}.data(Data_all{yr}.data<0.5 & Data_all{yr}.data>0)=0.5; %set nonzerocovers to minimum value of 0.5
    SpeciesIDs=find(sum(Data_all{yr}.data)>0);
    TransectAbundance=zeros(Samples(yr)/ksub_sub,Richness_raw(yr));
    for species=1:Richness_raw(yr)
        for transect=1:ksub %counting transect as a sample of the community
            TransectAbundance(transect,species)=round(sum(Data_all{yr}.data((transect-1)*10+1:(transect-1)*10+10,SpeciesIDs(species)))/0.5); %convert cover to count
        end
    end
    %correct bootstrapping accounting for non-independence of individuals belonging to same species and transect:
    [~,Chao1_all(yr),Chao2_all(yr),ACE_all(yr),S_aj2_all(yr),S_ij2_all(yr),Richness_apx(yr),expectedRichness_raw(:,yr),expectedChao1(:,yr),expectedChao2(:,yr),expectedACE(:,yr),expectedS_aj2(:,yr),expectedS_ij2(:,yr),expectedRichness_apx(:,yr)] = bootRichnessEsts_all(TransectAbundance,numBoot);
    %bootstrapping assuming individual presence is independent of other individuals
    %[~,Chao1_all(yr),Chao2_all(yr),ACE_all(yr),S_aj2_all(yr),S_ij2_all(yr),Richness_apx(yr),expectedRichness_raw(:,yr),expectedChao1(:,yr),expectedChao2(:,yr),expectedACE(:,yr),expectedS_aj2(:,yr),expectedS_ij2(:,yr),expectedRichness_apx(:,yr)] = bootRichnessEsts_all_indepIndivAssumed(TransectAbundance,numBoot);
end

figure('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3)/1.18 scrsz(3)/3.05]);
%test 1: plot all richness estimates from full dataset
subplot(2,4,1)
hold on
plot(Richness_raw,'k','LineWidth',2);
plot(Chao1_all,'b','LineWidth',2);
plot(Chao2_all,'--b','LineWidth',2);
plot(ACE_all,'g','LineWidth',2);
plot(S_aj2_all,'m','LineWidth',2);
plot(S_ij2_all,'--m','LineWidth',2);
plot(Richness_apx,'r','LineWidth',2);
legend({'raw','Chao1','Chao2','ACE','JK_a','JK_i','\Omega_o'},'location','NorthWest')
ylim([65 165])
xlabel 'year'
ylabel 'richness'
xticks([1:NumYears])
xticklabels(YearLabels)
xlim([1,NumYears])

subplot(2,4,2)
hold on
boundedline([1:NumYears]', mean(expectedRichness_apx)',[max(mean(expectedRichness_apx)'-prctile(expectedRichness_apx,2.5)',0),prctile(expectedRichness_apx,97.5)'-mean(expectedRichness_apx)'],'-r','alpha','transparency', 0.1);
boundedline([1:NumYears]', mean(expectedRichness_raw)',[max(mean(expectedRichness_raw)'-prctile(expectedRichness_raw,2.5)',0),prctile(expectedRichness_raw,97.5)'-mean(expectedRichness_raw)'],'-k','alpha','transparency', 0.1);
boundedline([1:NumYears]', mean(expectedChao1)',[max(mean(expectedChao1)'-prctile(expectedChao1,2.5)',0),prctile(expectedChao1,97.5)'-mean(expectedChao1)'],'-b','alpha','transparency', 0.1);
boundedline([1:NumYears]', mean(expectedChao2)',[max(mean(expectedChao2)'-prctile(expectedChao2,2.5)',0),prctile(expectedChao2,97.5)'-mean(expectedChao2)'],'--b','alpha','transparency', 0.1);
boundedline([1:NumYears]', mean(expectedACE)',[max(mean(expectedACE)'-prctile(expectedACE,2.5)',0),prctile(expectedACE,97.5)'-mean(expectedACE)'],'-g','alpha','transparency', 0.1);
boundedline([1:NumYears]', mean(expectedS_aj2)',[max(mean(expectedS_aj2)'-prctile(expectedS_aj2,2.5)',0),prctile(expectedS_aj2,97.5)'-mean(expectedS_aj2)'],'-m','alpha','transparency', 0.1);
boundedline([1:NumYears]', mean(expectedS_ij2)',[max(mean(expectedS_ij2)'-prctile(expectedS_ij2,2.5)',0),prctile(expectedS_ij2,97.5)'-mean(expectedS_ij2)'],'--m','alpha','transparency', 0.1);
ylim([65 165])
xlabel 'year'
ylabel 'richness'
xticks([1:NumYears])
xticklabels(YearLabels)
xlim([1,NumYears])

figure('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3)/1.18 scrsz(3)/3.05]);
%test 2: plot selected richness estimates across years on full dataset:
subplot(2,4,1)
hold on
plot(Richness_raw,'k','LineWidth',2);
plot(Chao1_all,'b','LineWidth',2);
% plot(Chao2_all,'--b','LineWidth',2);
% plot(ACE_all,'g','LineWidth',2);
% plot(S_aj2_all,'m','LineWidth',2);
% plot(S_ij2_all,'--m','LineWidth',2);
plot(Richness_apx,'r','LineWidth',2);
boundedline([1:NumYears]', mean(expectedRichness_apx)',[max(mean(expectedRichness_apx)'-prctile(expectedRichness_apx,2.5)',0),prctile(expectedRichness_apx,97.5)'-mean(expectedRichness_apx)'],'-r','alpha','transparency', 0.1);
boundedline([1:NumYears]', mean(expectedRichness_raw)',[max(mean(expectedRichness_raw)'-prctile(expectedRichness_raw,2.5)',0),prctile(expectedRichness_raw,97.5)'-mean(expectedRichness_raw)'],'-k','alpha','transparency', 0.1);
boundedline([1:NumYears]', mean(expectedChao1)',[max(mean(expectedChao1)'-prctile(expectedChao1,2.5)',0),prctile(expectedChao1,97.5)'-mean(expectedChao1)'],'-b','alpha','transparency', 0.1);
% boundedline([1:NumYears]', mean(expectedChao2)',[max(mean(expectedChao2)'-prctile(expectedChao2,2.5)',0),prctile(expectedChao2,97.5)'-mean(expectedChao2)'],'--b','alpha','transparency', 0.1);
% boundedline([1:NumYears]', mean(expectedACE)',[max(mean(expectedACE)'-prctile(expectedACE,2.5)',0),prctile(expectedACE,97.5)'-mean(expectedACE)'],'-g','alpha','transparency', 0.1);
% boundedline([1:NumYears]', mean(expectedS_aj2)',[max(mean(expectedS_aj2)'-prctile(expectedS_aj2,2.5)',0),prctile(expectedS_aj2,97.5)'-mean(expectedS_aj2)'],'-m','alpha','transparency', 0.1);
% boundedline([1:NumYears]', mean(expectedS_ij2)',[max(mean(expectedS_ij2)'-prctile(expectedS_ij2,2.5)',0),prctile(expectedS_ij2,97.5)'-mean(expectedS_ij2)'],'--m','alpha','transparency', 0.1);
%legend({'raw','Chao1','Chao2','ACE','JK_a','JK_i','\Omega_o'})

ylim([30 120])
xlabel 'year'
ylabel 'richness'
xticks([1:NumYears])
xticklabels(YearLabels)
xlim([1,NumYears])
title([num2str(ksub) '/9 transects, ' num2str(ksub_sub) '/10 quadrats'])

%test 3: get selected richness estimates on subsampled data:
ksub_perms=[9 3 2]; %pick out of 9 transects
ksub_sub_perms=[2 6 9]; %pick out of 10 quadrats in each transect

%figure('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3)/3 scrsz(3)]);
for testperm=1:length(ksub_perms)
    subplot(2,4,testperm+1) %plot estimated richness by year
    hold on
    ksub=ksub_perms(testperm);
    ksub_sub=ksub_sub_perms(testperm);
    
    Richness_raw_sub=zeros(numResample,NumYears); %point estimates
    Chao1_sub=zeros(numResample,NumYears); %point estimates
    Richness_apx_sub=zeros(numResample,NumYears); %point estimates
    Richness_raw_sub_mean=zeros(numResample,NumYears); %bootstrap means across resamples
    Chao1_sub_mean=zeros(numResample,NumYears);
    Richness_apx_sub_mean=zeros(numResample,NumYears);
    Richness_raw_sub_lo=zeros(numResample,NumYears); %bootstrap lower bound across resamples (95%)
    Chao1_sub_lo=zeros(numResample,NumYears);
    Richness_apx_sub_lo=zeros(numResample,NumYears);
    Richness_raw_sub_up=zeros(numResample,NumYears); %bootstrap upper bounds across resamples (95%)
    Chao1_sub_up=zeros(numResample,NumYears);
    Richness_apx_sub_up=zeros(numResample,NumYears);
    for resample=1:numResample
        transectsubID=randperm(9,ksub); %pick random transects for subsampling, same across years
        %QuadratSubID=randperm(10,ksub_sub);
        for yr=1:NumYears
            %             transectsubID=randperm(9,ksub); %pick random transects for subsampling, different across years
            SpeciesIDs=find(sum(Data_all{yr}.data)>0);
            TransectAbundance_sub=zeros(ksub,Richness_raw(yr));
            for species=1:Richness_raw(yr)
                for transect_sub=1:ksub %counting transect as a sample of the community
                    QuadratSubID=randperm(10,ksub_sub); %pick different quadrats, different across year and transects
                    TransectAbundance_sub(transect_sub,species)=round(sum(Data_all{yr}.data((transectsubID(transect_sub)-1)*10+QuadratSubID,SpeciesIDs(species)))/0.5);
                end
            end
            [Richness_raw_sub(resample,yr),Chao1_sub(resample,yr),Richness_apx_sub(resample,yr),expectedRichness_raw_boot,expectedChao1_boot,expectedRichness_apx_boot] = bootRichnessEsts(TransectAbundance_sub,numBoot);
            Richness_raw_sub_mean(resample,yr)=mean(expectedRichness_raw_boot);
            Chao1_sub_mean(resample,yr)=mean(expectedChao1_boot);
            Richness_apx_sub_mean(resample,yr)=mean(expectedRichness_apx_boot);
            Richness_raw_sub_lo(resample,yr)=prctile(expectedRichness_raw_boot,2.5);
            Chao1_sub_lo(resample,yr)=prctile(expectedChao1_boot,2.5);
            Richness_apx_sub_lo(resample,yr)=prctile(expectedRichness_apx_boot,2.5);
            Richness_raw_sub_up(resample,yr)=prctile(expectedRichness_raw_boot,97.5);
            Chao1_sub_up(resample,yr)=prctile(expectedChao1_boot,97.5);
            Richness_apx_sub_up(resample,yr)=prctile(expectedRichness_apx_boot,97.5);
        end
    end
    %         plot(Richness_apx_sub,'r','LineWidth',2);
    %         plot(Chao1_sub,'b','LineWidth',2);
    %plot(Richness_raw_sub,'k','LineWidth',2);
    %boundedline([1:NumYears]', mean(expectedRichness_raw)',std(expectedRichness_raw)','-k','alpha','transparency', 0.1);
    %boundedline([1:NumYears]', mean(expectedChao1)',std(expectedChao1)','-b','alpha','transparency', 0.1);
    %boundedline([1:NumYears]', mean(expectedRichness_apx)',std(expectedRichness_apx)','-r','alpha','transparency', 0.1);
%     boundedline([1:NumYears]', 100*mean(Richness_raw_sub./Richness_raw)',100*[max(mean(Richness_raw_sub./Richness_raw)'-prctile(Richness_raw_sub./Richness_raw,2.5)',0),prctile(Richness_raw_sub./Richness_raw,97.5)'-mean(Richness_raw_sub./Richness_raw)'],'-k','alpha','transparency', 0.1);
%     boundedline([1:NumYears]', 100*mean(Chao1_sub_mean./mean(expectedChao1))',100*[mean(Chao1_sub_mean./mean(expectedChao1))'-mean(Chao1_sub_lo./mean(expectedChao1))',mean(Chao1_sub_up./mean(expectedChao1))'-mean(Chao1_sub_mean./mean(expectedChao1))'],'-b','alpha','transparency', 0.1);
%     boundedline([1:NumYears]', 100*mean(Richness_apx_sub_mean./mean(expectedRichness_apx))',100*[mean(Richness_apx_sub_mean./mean(expectedRichness_apx))'-mean(Richness_apx_sub_lo./mean(expectedRichness_apx))',mean(Richness_apx_sub_up./mean(expectedRichness_apx))'-mean(Richness_apx_sub_mean./mean(expectedRichness_apx))'],'-r','alpha','transparency', 0.1);
    %bl1=boundedline([1:NumYears]', mean(Richness_raw_sub)',[max(mean(Richness_raw_sub)'-prctile(Richness_raw_sub,2.5)',0),prctile(Richness_raw_sub,97.5)'-mean(Richness_raw_sub)'],'-k','alpha','transparency', 0.1);
    boundedline([1:NumYears]', mean(Richness_raw_sub_mean)',[mean(Richness_raw_sub_mean)'-mean(Richness_raw_sub_lo)',mean(Richness_raw_sub_up)'-mean(Richness_raw_sub_mean)'],'-k','alpha','transparency', 0.1);
    boundedline([1:NumYears]', mean(Chao1_sub_mean)',[mean(Chao1_sub_mean)'-mean(Chao1_sub_lo)',mean(Chao1_sub_up)'-mean(Chao1_sub_mean)'],'-b','alpha','transparency', 0.1);
    boundedline([1:NumYears]', mean(Richness_apx_sub_mean)',[mean(Richness_apx_sub_mean)'-mean(Richness_apx_sub_lo)',mean(Richness_apx_sub_up)'-mean(Richness_apx_sub_mean)'],'-r','alpha','transparency', 0.1);
    plot(mean(Richness_raw_sub),'k','LineWidth',2);
    plot(mean(Chao1_sub),'b','LineWidth',2);
    plot(mean(Richness_apx_sub),'r','LineWidth',2);
    xlabel 'year'
    %ylim([0 120])
    %ylabel '% richness recovered'
    ylim([30 120])
    ylabel 'richness'
    title([num2str(ksub) '/9 transects, ' num2str(ksub_sub) '/10 quadrats'])
    xticks([1:NumYears])
    xticklabels(YearLabels)
    xlim([1,NumYears])
    %     text(1.3,97,['Taylor2 subsampled \Deltarichness/yr=' num2str(mean(RichnessSlope_apx_sub),2) '\pm' num2str(std(RichnessSlope_apx_sub),2) ' (p=' num2str(mean(RichnessSlopeP_apx_sub),2) '\pm' num2str(std(RichnessSlopeP_apx_sub),2) ')'],'color','r')
    %     text(1.3,92,['Chao1 subsampled \Deltarichness/yr=' num2str(mean(Chao1Slope_sub),2) '\pm' num2str(std(Chao1Slope_sub),2) ' (p=' num2str(mean(Chao1SlopeP_sub),2) '\pm' num2str(std(Chao1SlopeP_sub),2) ')'],'color','b')
    %     text(1.3,87,['raw subsampled \Deltarichness/yr=' num2str(mean(RichnessSlope_raw_sub),2) '\pm' num2str(std(RichnessSlope_raw_sub),2) ' (p=' num2str(mean(RichnessSlopeP_raw_sub),2) '\pm' num2str(std(RichnessSlopeP_raw_sub),2) ')'])
    
    subplot(2,4,4+testperm+1)
    %bl1=boundedline([1:NumYears]', 100*mean(Richness_raw_sub./Richness_raw)',100*[max(mean(Richness_raw_sub./Richness_raw)'-prctile(Richness_raw_sub./Richness_raw,2.5)',0),prctile(Richness_raw_sub./Richness_raw,97.5)'-mean(Richness_raw_sub./Richness_raw)'],'-k','alpha','transparency', 0.1);
    bl1=boundedline([1:NumYears]', 100*mean(Richness_raw_sub./Richness_raw)',100*[max(mean(Richness_raw_sub./Richness_raw)'-mean(Richness_raw_sub_lo./Richness_raw)',0),mean(Richness_raw_sub_up./Richness_raw)'-mean(Richness_raw_sub./Richness_raw)'],'-k','alpha','transparency', 0.1);
    bl2=boundedline([1:NumYears]', 100*mean(Chao1_sub./Chao1_all)',100*[mean(Chao1_sub./Chao1_all)'-mean(Chao1_sub_lo./Chao1_all)',mean(Chao1_sub_up./Chao1_all)'-mean(Chao1_sub./Chao1_all)'],'-b','alpha','transparency', 0.1);
    bl3=boundedline([1:NumYears]', 100*mean(Richness_apx_sub./Richness_apx)',100*[mean(Richness_apx_sub./Richness_apx)'-mean(Richness_apx_sub_lo./Richness_apx)',mean(Richness_apx_sub_up./Richness_apx)'-mean(Richness_apx_sub./Richness_apx)'],'-r','alpha','transparency', 0.1);
    set(bl1,'linewidth',2);
    set(bl2,'linewidth',2);
    set(bl3,'linewidth',2);
    h=refline(0,100)
    h.Color='k';
    h.LineStyle='--';
    xlabel 'year'
    ylabel '% richness recovered'
    ylim([40 110])
    %title(['transects subsampled=' num2str(ksub) '/9, quadrats subsampled=' num2str(ksub_sub) '/10'])
    xticks([1:NumYears])
    xticklabels(YearLabels)
    xlim([1,NumYears])
end

