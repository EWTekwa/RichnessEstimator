scrsz = get(0,'ScreenSize');
rng(6); %set random number generator seed


load('bci.tree_abundance.mat'); %load .mat file containing the cell Data_all with species % cover data from Martone seaweed survey

% %%%%% run following to import raw data from BCI and convert to
% %%%%% abundnace counts by quadrats (rows) and species (columns) per census
% %%%%% (third dimension)
% DataBCI1=importdata('bci.tree1.csv');
% DataBCI2=importdata('bci.tree2.csv');
% DataBCI3=importdata('bci.tree3.csv');
% DataBCI4=importdata('bci.tree4.csv');
% DataBCI5=importdata('bci.tree5.csv');
% DataBCI6=importdata('bci.tree6.csv');
% DataBCI7=importdata('bci.tree7.csv');
% DataBCI8=importdata('bci.tree8.csv');
% Data_all={DataBCI1,DataBCI2,DataBCI3,DataBCI4,DataBCI5,DataBCI6,DataBCI7,DataBCI8};
% clear DataBCI1 DataBCI2 DataBCI3 DataBCI4 DataBCI5 DataBCI6 DataBCI7 DataBCI8
%
% NumYears=length(Data_all);
%
% QuadratAbundance_m=zeros(500,1250,NumYears);
% for set=1:NumYears
%     speciesNames=unique(Data_all{set}.textdata(2:end,3)); %get species names
%     quadratNames=unique(Data_all{set}.textdata(2:end,4)); %get quadrat names
%
%     numQuadrats=length(quadratNames)-1;
%     maxRichness=length(speciesNames);
%     %QuadratAbundance_m=zeros(numQuadrats,maxRichness,NumYears);
%     for i=1:maxRichness
%         for j=1:numQuadrats
%             treeRow=find(strcmp(Data_all{set}.textdata(2:end,3),speciesNames{i}))+1; %find trees belonging to species j
%             numTreesAtJ=sum(strcmp(Data_all{set}.textdata(treeRow,4),quadratNames{j}) .* (strcmp(Data_all{set}.textdata(treeRow,13),'A') + strcmp(Data_all{set}.textdata(treeRow,13),'P')));
%             QuadratAbundance_m(i,j,set)=QuadratAbundance_m(i,j,set)+numTreesAtJ;
%         end
%     end
% end
%
% clear Data_all treeRow numTreesAtJ
% %%%%%


ksub_perms=[12 3]; %spatial subsampling experiments out of 1250 quadrats
m_perms=[1 0.1]; %downsampling experiments within each quadrats (fraction of individuals observed)

numYears=size(QuadratAbundance,3);
YearLabels=[1:numYears]; %these are census numbers, not consecutive years
numResample=40; %number of resamples for subsampled estimates (40)
numBoot=2000; %total number of bootstramps for CI (2000): numBoot/numResample per subsample
numQuad=size(QuadratAbundance,1); %number of quadrats
maxRich=550; %for display only
minRich=0;
diffRich=maxRich-minRich;

allSub_Chao1=zeros(length(m_perms)*(length(ksub_perms)+1),NumYears);
allSub_Chao2=zeros(length(m_perms)*(length(ksub_perms)+1),NumYears);
allSub_Richness_apx=zeros(length(m_perms)*(length(ksub_perms)+1),NumYears);
allSub_Richness_raw=zeros(length(m_perms)*(length(ksub_perms)+1),NumYears);
allSub_Chao1SD=zeros(length(m_perms)*(length(ksub_perms)+1),NumYears);
allSub_Chao2SD=zeros(length(m_perms)*(length(ksub_perms)+1),NumYears);
allSub_Richness_apxSD=zeros(length(m_perms)*(length(ksub_perms)+1),NumYears);
allSub_Richness_rawSD=zeros(length(m_perms)*(length(ksub_perms)+1),NumYears);

allSubSlope_Chao1=zeros(length(m_perms)*(length(ksub_perms)+1),1);
allSubSlopeSD_Chao1=zeros(length(m_perms)*(length(ksub_perms)+1),1);
allSubSlopeP_Chao1=zeros(length(m_perms)*(length(ksub_perms)+1),1);
allSubSlopePPt_Chao1=zeros(length(m_perms)*(length(ksub_perms)+1),1);
allSubSlopeSDP_Chao1=zeros(length(m_perms)*(length(ksub_perms)+1),1);
allSubSlopeBeta_Chao1=zeros(length(m_perms)*(length(ksub_perms)+1),1);
allSubSlopeBetaPt_Chao1=zeros(length(m_perms)*(length(ksub_perms)+1),1);

allSubSlope_Chao2=zeros(length(m_perms)*(length(ksub_perms)+1),1);
allSubSlopeSD_Chao2=zeros(length(m_perms)*(length(ksub_perms)+1),1);
allSubSlopeP_Chao2=zeros(length(m_perms)*(length(ksub_perms)+1),1);
allSubSlopePPt_Chao2=zeros(length(m_perms)*(length(ksub_perms)+1),1);
allSubSlopeSDP_Chao2=zeros(length(m_perms)*(length(ksub_perms)+1),1);
allSubSlopeBeta_Chao2=zeros(length(m_perms)*(length(ksub_perms)+1),1);
allSubSlopeBetaPt_Chao2=zeros(length(m_perms)*(length(ksub_perms)+1),1);

allSubSlope_apx=zeros(length(m_perms)*(length(ksub_perms)+1),1);
allSubSlopeSD_apx=zeros(length(m_perms)*(length(ksub_perms)+1),1);
allSubSlopeP_apx=zeros(length(m_perms)*(length(ksub_perms)+1),1);
allSubSlopePPt_apx=zeros(length(m_perms)*(length(ksub_perms)+1),1);
allSubSlopeSDP_apx=zeros(length(m_perms)*(length(ksub_perms)+1),1);
allSubSlopeBeta_apx=zeros(length(m_perms)*(length(ksub_perms)+1),1);
allSubSlopeBetaPt_apx=zeros(length(m_perms)*(length(ksub_perms)+1),1);

allSubSlope_raw=zeros(length(m_perms)*(length(ksub_perms)+1),1);
allSubSlopeSD_raw=zeros(length(m_perms)*(length(ksub_perms)+1),1);
allSubSlopeP_raw=zeros(length(m_perms)*(length(ksub_perms)+1),1);
allSubSlopePPt_raw=zeros(length(m_perms)*(length(ksub_perms)+1),1);
allSubSlopeSDP_raw=zeros(length(m_perms)*(length(ksub_perms)+1),1);
allSubSlopeBeta_raw=zeros(length(m_perms)*(length(ksub_perms)+1),1);
allSubSlopeBetaPt_raw=zeros(length(m_perms)*(length(ksub_perms)+1),1);

plotLabels=cell(1,length(m_perms)*(length(ksub_perms)+1));

%figure('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3)/2 scrsz(4)/2.4],'DefaultAxesFontSize',16,'defaultAxesColorOrder',[[0 0 0]; [1 0 1]]);
figure('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3)/2 scrsz(4)/1.2],'DefaultAxesFontSize',16,'defaultAxesColorOrder',[[0 0 0]; [0 0 0]]);

i=1; %subplot label
%prepare subsampling before simulations (so same set of resampled quadrats are used for paired downsampling experiments)
quadratsubIDs=cell(length(ksub_perms),numResample);
for testperm=1:length(ksub_perms)
    for resample=1:numResample
        quadratsubIDs{testperm,resample}=randperm(numQuad,ksub_perms(testperm)); %pick random quadrats without replacement for subsampling, same across years
    end
end

%run simulations for all subsampling and paired downsampling scenarios
for mperm=1:length(m_perms)
    QuadratAbundance_m=zeros(size(QuadratAbundance));
    Richness_raw=zeros(1,NumYears);
    Chao1_all=zeros(1,NumYears);
    Chao2_all=zeros(1,NumYears);
    Richness_apx=zeros(1,NumYears);
    Richness_taylor=zeros(1,NumYears);
    min_cover=zeros(1,NumYears);
    meanStates=zeros(5,NumYears); %mean mn,P, var mn,P, cov(mn,P)
    %Samples=zeros(1,NumYears);
    for yr=1:NumYears
        if m_perms(mperm)==1
            QuadratAbundance_m(:,:,yr)=QuadratAbundance(:,:,yr); %use original dataset
        else
            QuadratAbundance_m(:,:,yr)=poissrnd(QuadratAbundance(:,:,yr)*m_perms(mperm)); %downsample within each quadrat to a fraction of true abundances per species
        end
        Richness_raw(yr)=sum(sum(QuadratAbundance_m(:,:,yr))>0);
        %min_cover(yr)=min(Data_all{yr}.data(Data_all{yr}.data>0));
        %Samples(yr)=size(QuadratAbundance_m(:,:,yr),1);
    end
    
    expectedRichness_raw=zeros(numBoot,NumYears);
    expectedChao1=zeros(numBoot,NumYears);
    expectedChao2=zeros(numBoot,NumYears);
    expectedRichness_apx=zeros(numBoot,NumYears);
    expectedRichness_taylor=zeros(numBoot,NumYears);
    ke_yr=zeros(1,NumYears);
    for yr=1:NumYears
        %correct bootstrapping accounting for non-independence of individuals belonging to same species and quadrat:
        [~,Chao1_all(yr),~,Chao2_all(yr),~,~,~,Richness_apx(yr),Richness_taylor(yr),~,expectedRichness_raw(:,yr),expectedChao1(:,yr),~,expectedChao2(:,yr),~,~,~,expectedRichness_apx(:,yr),expectedRichness_taylor(:,yr),~,meanStates(:,yr)] = bootRichnessEsts(QuadratAbundance_m(:,:,yr),numBoot);
        %bootstrapping assuming individual presence is independent of other individuals
        %[~,Chao1_all(yr),Chao2_all(yr),ACE_all(yr),S_aj2_all(yr),S_ij2_all(yr),Richness_apx(yr),expectedRichness_raw(:,yr),expectedChao1(:,yr),expectedChao2(:,yr),expectedACE(:,yr),expectedS_aj2(:,yr),expectedS_ij2(:,yr),expectedRichness_apx(:,yr)] = bootRichnessEsts_indepIndivAssumed(QuadratAbundance_m(:,:,yr),numBoot);
    end
    
    %record estimates for full dataset:
    allSub_Chao1((mperm-1)*(length(ksub_perms)+1)+1,:)=Chao1_all;
    allSub_Chao2((mperm-1)*(length(ksub_perms)+1)+1,:)=Chao2_all;
    allSub_Richness_apx((mperm-1)*(length(ksub_perms)+1)+1,:)=Richness_apx;
    allSub_Richness_raw((mperm-1)*(length(ksub_perms)+1)+1,:)=Richness_raw;
    
    %get temporal richness slope for centred Chao1 and apx boots
    slopeBoot=numBoot;
    RawDiff=(Richness_raw-mean(expectedRichness_raw));
    Chao1Diff=(Chao1_all-mean(expectedChao1));
    Chao2Diff=(Chao2_all-mean(expectedChao2));
    ApxDiff=(Richness_apx-mean(expectedRichness_apx));
    TaylorDiff=(Richness_taylor-mean(expectedRichness_taylor));
    RichnessSlope_raw_boot=zeros(slopeBoot,1);
    RichnessSlopeSD_raw_boot=zeros(slopeBoot,1);
    RichnessSlopeP_raw_boot=zeros(slopeBoot,1);
    RichnessSlope_Chao1_boot=zeros(slopeBoot,1);
    RichnessSlopeSD_Chao1_boot=zeros(slopeBoot,1);
    RichnessSlopeP_Chao1_boot=zeros(slopeBoot,1);
    RichnessSlope_Chao2_boot=zeros(slopeBoot,1);
    RichnessSlopeSD_Chao2_boot=zeros(slopeBoot,1);
    RichnessSlopeP_Chao2_boot=zeros(slopeBoot,1);
    RichnessSlope_Apx_boot=zeros(slopeBoot,1);
    RichnessSlopeSD_Apx_boot=zeros(slopeBoot,1);
    RichnessSlopeP_Apx_boot=zeros(slopeBoot,1);
    RichnessSlope_Taylor_boot=zeros(slopeBoot,1);
    RichnessSlopeSD_Taylor_boot=zeros(slopeBoot,1);
    RichnessSlopeP_Taylor_boot=zeros(slopeBoot,1);
    
    %get temporal slopes:
    %regress based on point richness estimates:
    mdl_Richness_raw=fitlm([1:NumYears],Richness_raw);
    mdl_Chao1=fitlm([1:NumYears],Chao1_all);
    mdl_Chao2=fitlm([1:NumYears],Chao2_all);
    mdl_Apx=fitlm([1:NumYears],Richness_apx);
    mdl_Taylor=fitlm([1:NumYears],Richness_taylor);
    RichnessSlope_raw=mdl_Richness_raw.Coefficients.Estimate(2);
    RichnessSlopeSD_raw=mdl_Richness_raw.Coefficients.SE(2)*sqrt(NumYears);
    RichnessSlopeP_raw=mdl_Richness_raw.Coefficients.pValue(2);
    RichnessSlope_Chao1=mdl_Chao1.Coefficients.Estimate(2);
    RichnessSlopeSD_Chao1=mdl_Chao1.Coefficients.SE(2)*sqrt(NumYears);
    RichnessSlopeP_Chao1=mdl_Chao1.Coefficients.pValue(2);
    RichnessSlope_Chao2=mdl_Chao2.Coefficients.Estimate(2);
    RichnessSlopeSD_Chao2=mdl_Chao2.Coefficients.SE(2)*sqrt(NumYears);
    RichnessSlopeP_Chao2=mdl_Chao2.Coefficients.pValue(2);
    RichnessSlope_Apx=mdl_Apx.Coefficients.Estimate(2);
    RichnessSlopeSD_Apx=mdl_Apx.Coefficients.SE(2)*sqrt(NumYears);
    RichnessSlopeP_Apx=mdl_Apx.Coefficients.pValue(2);
    RichnessSlope_Taylor=mdl_Taylor.Coefficients.Estimate(2);
    RichnessSlopeSD_Taylor=mdl_Taylor.Coefficients.SE(2)*sqrt(NumYears);
    RichnessSlopeP_Taylor=mdl_Taylor.Coefficients.pValue(2);
    %regress based on bootstrapped richness estimates:
    for boot=1:slopeBoot
        mdl_Richness_raw_boot=fitlm([1:NumYears],expectedRichness_raw(boot,:)+RawDiff);
        mdl_Chao1_boot=fitlm([1:NumYears],expectedChao1(boot,:)+Chao1Diff);
        mdl_Chao2_boot=fitlm([1:NumYears],expectedChao2(boot,:)+Chao2Diff);
        mdl_Apx_boot=fitlm([1:NumYears],expectedRichness_apx(boot,:)+ApxDiff);
        mdl_Taylor_boot=fitlm([1:NumYears],expectedRichness_taylor(boot,:)+TaylorDiff);
        %end
        RichnessSlope_raw_boot(boot)=mdl_Richness_raw_boot.Coefficients.Estimate(2);
        RichnessSlopeSD_raw_boot(boot)=mdl_Richness_raw_boot.Coefficients.SE(2)*sqrt(NumYears);
        RichnessSlopeP_raw_boot(boot)=mdl_Richness_raw_boot.Coefficients.pValue(2);
        RichnessSlope_Chao1_boot(boot)=mdl_Chao1_boot.Coefficients.Estimate(2);
        RichnessSlopeSD_Chao1_boot(boot)=mdl_Chao1_boot.Coefficients.SE(2)*sqrt(NumYears);
        RichnessSlopeP_Chao1_boot(boot)=mdl_Chao1_boot.Coefficients.pValue(2);
        RichnessSlope_Chao2_boot(boot)=mdl_Chao2_boot.Coefficients.Estimate(2);
        RichnessSlopeSD_Chao2_boot(boot)=mdl_Chao2_boot.Coefficients.SE(2)*sqrt(NumYears);
        RichnessSlopeP_Chao2_boot(boot)=mdl_Chao2_boot.Coefficients.pValue(2);
        RichnessSlope_Apx_boot(boot)=mdl_Apx_boot.Coefficients.Estimate(2);
        RichnessSlopeSD_Apx_boot(boot)=mdl_Apx_boot.Coefficients.SE(2)*sqrt(NumYears);
        RichnessSlopeP_Apx_boot(boot)=mdl_Apx_boot.Coefficients.pValue(2);
        RichnessSlope_Taylor_boot(boot)=mdl_Taylor_boot.Coefficients.Estimate(2);
        RichnessSlopeSD_Taylor_boot(boot)=mdl_Taylor_boot.Coefficients.SE(2)*sqrt(NumYears);
        RichnessSlopeP_Taylor_boot(boot)=mdl_Taylor_boot.Coefficients.pValue(2);
    end
    
    %record estimated richness SD from full dataset:
    allSub_Chao1SD((mperm-1)*(length(ksub_perms)+1)+1,:)=std(expectedChao1);
    allSub_Chao2SD((mperm-1)*(length(ksub_perms)+1)+1,:)=std(expectedChao2);
    allSub_Richness_apxSD((mperm-1)*(length(ksub_perms)+1)+1,:)=std(expectedRichness_apx);
    allSub_Richness_rawSD((mperm-1)*(length(ksub_perms)+1)+1,:)=std(expectedRichness_raw);
    
    %record trend estimates from full dataset:
    allSubSlope_Chao1((mperm-1)*(length(ksub_perms)+1)+1,:)=mean(RichnessSlope_Chao1_boot(:));
    allSubSlopeSD_Chao1((mperm-1)*(length(ksub_perms)+1)+1,:)=mean(RichnessSlopeSD_Chao1_boot(:));
    allSubSlopeP_Chao1((mperm-1)*(length(ksub_perms)+1)+1,:)=mean(RichnessSlopeP_Chao1_boot(:));
    allSubSlopePPt_Chao1((mperm-1)*(length(ksub_perms)+1)+1,:)=mean(RichnessSlopeP_Chao1);
    allSubSlopeSDP_Chao1((mperm-1)*(length(ksub_perms)+1)+1,:)=std(RichnessSlopeP_Chao1_boot(:));
    allSubSlopeBeta_Chao1((mperm-1)*(length(ksub_perms)+1)+1,:)=sum(RichnessSlopeP_Chao1_boot(:)>0.05)/numBoot;
    allSubSlopeBetaPt_Chao1((mperm-1)*(length(ksub_perms)+1)+1,:)=sum(RichnessSlopeP_Chao1(:)>0.05)/numResample;
    
    allSubSlope_Chao2((mperm-1)*(length(ksub_perms)+1)+1,:)=mean(RichnessSlope_Chao2_boot(:));
    allSubSlopeSD_Chao2((mperm-1)*(length(ksub_perms)+1)+1,:)=mean(RichnessSlopeSD_Chao2_boot(:));
    allSubSlopeP_Chao2((mperm-1)*(length(ksub_perms)+1)+1,:)=mean(RichnessSlopeP_Chao2_boot(:));
    allSubSlopePPt_Chao2((mperm-1)*(length(ksub_perms)+1)+1,:)=mean(RichnessSlopeP_Chao2);
    allSubSlopeSDP_Chao2((mperm-1)*(length(ksub_perms)+1)+1,:)=std(RichnessSlopeP_Chao2_boot(:));
    allSubSlopeBeta_Chao2((mperm-1)*(length(ksub_perms)+1)+1,:)=sum(RichnessSlopeP_Chao2_boot(:)>0.05)/numBoot;
    allSubSlopeBetaPt_Chao2((mperm-1)*(length(ksub_perms)+1)+1,:)=sum(RichnessSlopeP_Chao2(:)>0.05)/numBoot;
    
    allSubSlope_apx((mperm-1)*(length(ksub_perms)+1)+1,:)=mean(RichnessSlope_Apx_boot(:));
    allSubSlopeSD_apx((mperm-1)*(length(ksub_perms)+1)+1,:)=mean(RichnessSlopeSD_Apx_boot(:));
    allSubSlopeP_apx((mperm-1)*(length(ksub_perms)+1)+1,:)=mean(RichnessSlopeP_Apx_boot(:));
    allSubSlopePPt_apx((mperm-1)*(length(ksub_perms)+1)+1,:)=mean(RichnessSlopeP_Apx);
    allSubSlopeSDP_apx((mperm-1)*(length(ksub_perms)+1)+1,:)=std(RichnessSlopeP_Apx_boot(:));
    allSubSlopeBeta_apx((mperm-1)*(length(ksub_perms)+1)+1,:)=sum(RichnessSlopeP_Apx_boot(:)>0.05)/numBoot;
    allSubSlopeBetaPt_apx((mperm-1)*(length(ksub_perms)+1)+1,:)=sum(RichnessSlopeP_Apx(:)>0.05)/numBoot;
    
    allSubSlope_raw((mperm-1)*(length(ksub_perms)+1)+1,:)=mean(RichnessSlope_raw_boot(:));
    allSubSlopeSD_raw((mperm-1)*(length(ksub_perms)+1)+1,:)=mean(RichnessSlopeSD_raw_boot(:));
    allSubSlopeP_raw((mperm-1)*(length(ksub_perms)+1)+1,:)=mean(RichnessSlopeP_raw_boot(:));
    allSubSlopePPt_raw((mperm-1)*(length(ksub_perms)+1)+1,:)=mean(RichnessSlopeP_raw);
    allSubSlopeSDP_raw((mperm-1)*(length(ksub_perms)+1)+1,:)=std(RichnessSlopeP_raw_boot(:));
    allSubSlopeBeta_raw((mperm-1)*(length(ksub_perms)+1)+1,:)=sum(RichnessSlopeP_raw_boot(:)>0.05)/numBoot;
    allSubSlopeBetaPt_raw((mperm-1)*(length(ksub_perms)+1)+1,:)=sum(RichnessSlopeP_raw(:)>0.05)/numBoot;

    
    %     figure('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3)/1.18 scrsz(3)/3.05]);
    %     %test 1: plot all richness estimates from full dataset
    %     subplot(2,4,mperm1)
    %     hold on
    %     plot(Chao1_all,'b','LineWidth',2);
    %     plot(Chao2_all,'--b','LineWidth',2);
    %     plot(ACE_all,'g','LineWidth',2);
    %     plot(S_aj2_all,'m','LineWidth',2);
    %     plot(S_ij2_all,'--m','LineWidth',2);
    %     plot(Richness_apx,'r','LineWidth',2);
    %     plot(Richness_raw,'--k','LineWidth',2);
    %     legend({'Chao1','Chao2','ACE','JK_a','JK_i','\Omega_o','raw'},'location','NorthWest')
    %     %ylim([65 165])
    %     xlabel 'year'
    %     ylabel 'richness'
    %     xticks([1:NumYears])
    %     xticklabels(YearLabels)
    %     xlim([1,NumYears])
    %
    %     subplot(2,4,2)
    %     hold on
    %     boundedline([1:NumYears]', mean(expectedRichness_apx)',[max(mean(expectedRichness_apx)'-prctile(expectedRichness_apx,2.5)',0),prctile(expectedRichness_apx,97.5)'-mean(expectedRichness_apx)'],'-r','alpha','transparency', 0.1);
    %     boundedline([1:NumYears]', mean(expectedChao1)',[max(mean(expectedChao1)'-prctile(expectedChao1,2.5)',0),prctile(expectedChao1,97.5)'-mean(expectedChao1)'],'-b','alpha','transparency', 0.1);
    %     boundedline([1:NumYears]', mean(expectedChao2)',[max(mean(expectedChao2)'-prctile(expectedChao2,2.5)',0),prctile(expectedChao2,97.5)'-mean(expectedChao2)'],'--b','alpha','transparency', 0.1);
    %     boundedline([1:NumYears]', mean(expectedACE)',[max(mean(expectedACE)'-prctile(expectedACE,2.5)',0),prctile(expectedACE,97.5)'-mean(expectedACE)'],'-g','alpha','transparency', 0.1);
    %     boundedline([1:NumYears]', mean(expectedS_aj2)',[max(mean(expectedS_aj2)'-prctile(expectedS_aj2,2.5)',0),prctile(expectedS_aj2,97.5)'-mean(expectedS_aj2)'],'-m','alpha','transparency', 0.1);
    %     boundedline([1:NumYears]', mean(expectedS_ij2)',[max(mean(expectedS_ij2)'-prctile(expectedS_ij2,2.5)',0),prctile(expectedS_ij2,97.5)'-mean(expectedS_ij2)'],'--m','alpha','transparency', 0.1);
    %     boundedline([1:NumYears]', mean(expectedRichness_raw)',[max(mean(expectedRichness_raw)'-prctile(expectedRichness_raw,2.5)',0),prctile(expectedRichness_raw,97.5)'-mean(expectedRichness_raw)'],'--k','alpha','transparency', 0.1);
    %
    %     %ylim([65 165])
    %     xlabel 'year'
    %     ylabel 'richness'
    %     xticks([1:NumYears])
    %     xticklabels(YearLabels)
    %     xlim([1,NumYears])
    %
    %    figure('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3)/1.18 scrsz(3)/3.05]);
    %test 2: plot selected richness estimates across years on full dataset:
    subplot(4,length(ksub_perms)+1,(mperm-1)*(length(ksub_perms)+1)+1)
    %yyaxis left
    hold on
%     boundedline([1:NumYears]', Richness_taylor',[max(mean(Richness_taylor_sub_mean-Richness_taylor)',0),prctile(expectedRichness_apx,97.5)'-mean(expectedRichness_apx)'],'-r','alpha','transparency', 0.1);
%     boundedline([1:NumYears]', mean(expectedRichness_raw)',[max(mean(expectedRichness_raw)'-prctile(expectedRichness_raw,2.5)',0),prctile(expectedRichness_raw,97.5)'-mean(expectedRichness_raw)'],'--k','alpha','transparency', 0.1);
%     boundedline([1:NumYears]', mean(expectedChao1)',[max(mean(expectedChao1)'-prctile(expectedChao1,2.5)',0),prctile(expectedChao1,97.5)'-mean(expectedChao1)'],'-b','alpha','transparency', 0.1);
    %centre confidence bounds to match bootstrapped mean to spot estimates
    %on the original dataset:
    %boundedline([1:NumYears]', Richness_apx',[max(mean(expectedRichness_apx)'-prctile(expectedRichness_apx,2.5)',0),prctile(expectedRichness_apx,97.5)'-mean(expectedRichness_apx)'],'-r','alpha','transparency', 0.1);
    boundedline([1:NumYears]', Chao2_all',[max(mean(expectedChao2)'-prctile(expectedChao2,2.5)',0),prctile(expectedChao1,97.5)'-mean(expectedChao1)'],'--c','alpha','transparency', 0.1);
    boundedline([1:NumYears]', Chao1_all',[max(mean(expectedChao1)'-prctile(expectedChao1,2.5)',0),prctile(expectedChao1,97.5)'-mean(expectedChao1)'],'-b','alpha','transparency', 0.1);
    boundedline([1:NumYears]', Richness_apx',[max(mean(expectedRichness_apx)'-prctile(expectedRichness_apx,2.5)',0),prctile(expectedRichness_apx,97.5)'-mean(expectedRichness_apx)'],'-r','alpha','transparency', 0.1);
    boundedline([1:NumYears]', Richness_raw',[max(mean(expectedRichness_raw)'-prctile(expectedRichness_raw,2.5)',0),prctile(expectedRichness_raw,97.5)'-mean(expectedRichness_raw)'],'--k','alpha','transparency', 0.1);
    
    plot(Chao1_all,'-b','LineWidth',2);
    plot(Chao2_all,'--','LineWidth',2,'Color',[0 0.7 0.7]);
    % plot(ACE_all,'g','LineWidth',2);
    % plot(S_aj2_all,'m','LineWidth',2);
    % plot(S_ij2_all,'--m','LineWidth',2);
    %plot(Richness_apx,'-r','LineWidth',2);
    plot(Richness_apx,'-r','LineWidth',2);
    plot(Richness_raw,'--k','LineWidth',2);
    
    %     plot(mean(expectedChao1),'b','LineWidth',2);
    %     plot(mean(expectedRichness_apx),'r','LineWidth',2);
    %     plot(mean(expectedRichness_raw),'--k','LineWidth',2);
    % boundedline([1:NumYears]', mean(expectedChao2)',[max(mean(expectedChao2)'-prctile(expectedChao2,2.5)',0),prctile(expectedChao2,97.5)'-mean(expectedChao2)'],'--b','alpha','transparency', 0.1);
    % boundedline([1:NumYears]', mean(expectedACE)',[max(mean(expectedACE)'-prctile(expectedACE,2.5)',0),prctile(expectedACE,97.5)'-mean(expectedACE)'],'-g','alpha','transparency', 0.1);
    % boundedline([1:NumYears]', mean(expectedS_aj2)',[max(mean(expectedS_aj2)'-prctile(expectedS_aj2,2.5)',0),prctile(expectedS_aj2,97.5)'-mean(expectedS_aj2)'],'-m','alpha','transparency', 0.1);
    % boundedline([1:NumYears]', mean(expectedS_ij2)',[max(mean(expectedS_ij2)'-prctile(expectedS_ij2,2.5)',0),prctile(expectedS_ij2,97.5)'-mean(expectedS_ij2)'],'--m','alpha','transparency', 0.1);
    %legend({'raw','Chao1','Chao2','ACE','JK_a','JK_i','\Omega_o'})
    
    %ylim([30 120])
    xlabel 'census'
    ylabel 'richness'
    xticks([1:NumYears])
    xticklabels(YearLabels)
    xlim([1,NumYears])
    title({['1250 quadrats'];[num2str(m_perms(mperm)*100) '% individuals observed']},'fontweight','normal')
    ylim([minRich maxRich])
%     text(1.2,minRich+diffRich*.95,['\DeltaRaw=' num2str(mean(RichnessSlope_raw_boot(:)),1) '\pm' num2str(mean(RichnessSlopeSD_raw_boot(:)),1) ', p=' num2str(mean(RichnessSlopeP_raw_boot(:)),1) '(' num2str(mean(RichnessSlopeP_raw),1) ')\pm' num2str(std(RichnessSlopeP_raw_boot(:)),1) ', \beta=' num2str(sum(RichnessSlopeP_raw_boot(:)>0.05)/numBoot,1)],'fontsize',14);
%     text(1.2,minRich+diffRich*.85,['\DeltaChao1=' num2str(mean(RichnessSlope_Chao1_boot(:)),1) '\pm' num2str(mean(RichnessSlopeSD_Chao1_boot(:)),1) ', p=' num2str(mean(RichnessSlopeP_Chao1_boot(:)),1) '(' num2str(mean(RichnessSlopeP_Chao1),1) ')\pm' num2str(std(RichnessSlopeP_Chao1_boot(:)),1) ', \beta=' num2str(sum(RichnessSlopeP_Chao1_boot(:)>0.05)/numBoot,1)],'fontsize',14,'color','b');
%     text(1.2,minRich+diffRich*.75,['\DeltaChao2=' num2str(mean(RichnessSlope_Chao2_boot(:)),1) '\pm' num2str(mean(RichnessSlopeSD_Chao2_boot(:)),1) ', p=' num2str(mean(RichnessSlopeP_Chao2_boot(:)),1) '(' num2str(mean(RichnessSlopeP_Chao2),1) ')\pm' num2str(std(RichnessSlopeP_Chao2_boot(:)),1) ', \beta=' num2str(sum(RichnessSlopeP_Chao2_boot(:)>0.05)/numBoot,1)],'fontsize',14,'Color',[0 0.7 0.7]);
%     text(1.2,minRich+diffRich*.63,['\Delta\Omega_T=' num2str(mean(RichnessSlope_Taylor_boot(:)),1) '\pm' num2str(mean(RichnessSlopeSD_Taylor_boot(:)),1) ', p=' num2str(mean(RichnessSlopeP_Taylor_boot(:)),1) '(' num2str(mean(RichnessSlopeP_Taylor),1) ')\pm' num2str(std(RichnessSlopeP_Taylor_boot(:)),1) ', \beta=' num2str(sum(RichnessSlopeP_Taylor_boot(:)>0.05)/numBoot,1)],'fontsize',14,'Color','r');
    %     text(1,2,maxRich*.95,['bootstrapped estimates:'],'fontsize',14)
    %     text(1.2,maxRich*.85,['raw slope=' num2str(mean(RichnessSlope_raw_boot),1) '\pm' num2str(mean(RichnessSlopeSD_raw_boot),1) ', p=' num2str(mean(RichnessSlopeP_raw_boot),1) '\pm' num2str(std(RichnessSlopeP_raw_boot),1) ', \beta=' num2str(sum(RichnessSlopeP_raw_boot(:)>0.05)/numBoot,1)],'fontsize',14);
    %     text(1.2,maxRich*.75,['Chao1 slope=' num2str(mean(RichnessSlope_Chao1_boot),1) '\pm' num2str(mean(RichnessSlopeSD_Chao1_boot),1) ', p=' num2str(mean(RichnessSlopeP_Chao1_boot),1) '\pm' num2str(std(RichnessSlopeP_Chao1_boot),1) ', \beta=' num2str(sum(RichnessSlopeP_Chao1_boot(:)>0.05)/numBoot,1)],'fontsize',14,'color','b');
    %     text(1.2,maxRich*.63,['\Omega_o slope=' num2str(mean(RichnessSlope_Apx_boot),1) '\pm' num2str(mean(RichnessSlopeSD_Apx_boot),1) ', p=' num2str(mean(RichnessSlopeP_Apx_boot),1) '\pm' num2str(std(RichnessSlopeP_Apx_boot),1) ', \beta=' num2str(sum(RichnessSlopeP_Apx_boot(:)>0.05)/numBoot,1)],'fontsize',14,'Color','r');
    %     text(1,2,maxRich*.45,['point estimates:'],'fontsize',14)
    %     text(1.2,maxRich*.35,['raw slope=' num2str(mean(RichnessSlope_raw),1) '\pm' num2str(mean(RichnessSlopeSD_raw),1) ', p=' num2str(mean(RichnessSlopeP_raw),1)],'fontsize',14);
    %     text(1.2,maxRich*.25,['Chao1 slope=' num2str(mean(RichnessSlope_Chao1),1) '\pm' num2str(mean(RichnessSlopeSD_Chao1),1) ', p=' num2str(mean(RichnessSlopeP_Chao1),1)],'fontsize',14,'color','b');
    %     text(1.2,maxRich*.13,['\Omega_o slope=' num2str(mean(RichnessSlope_Apx),1) '\pm' num2str(mean(RichnessSlopeSD_Apx),1) ', p=' num2str(mean(RichnessSlopeP_Apx),1)],'fontsize',14,'Color','r');
    text(1-(NumYears-1)*0.15,maxRich*1.15,char(64+i),'Fontsize',16)
    
    %plot state variable means and variances over years
    subplot(4,length(ksub_perms)+1,(mperm-1)*(length(ksub_perms)+1)+1+6)
    hold on
    plot(meanStates(1,:),'-','LineWidth',2,'color','b'); %E[mn]
    plot(meanStates(2,:),'-','LineWidth',2,'color',[255 139 37]/255); %E[P]
    plot(meanStates(3,:),'--','LineWidth',2,'color','b'); %var[mn]
    plot(meanStates(4,:),'--','LineWidth',2,'color',[255 139 37]/255); %var[P]
    plot(meanStates(5,:),'--','LineWidth',2,'color',[69 181 80]/255); %cov[mn,P]
    
%     plot(meanStates(1,:),'-','LineWidth',2,'color','b'); %E[mn]
%     plot(meanStates(2,:),'-','LineWidth',2,'color',[69 181 80]/255); %E[C]
%     plot(meanStates(3,:),'-','LineWidth',2,'color',[255 139 37]/255); %E[P]
%     plot(meanStates(4,:),'--','LineWidth',2,'color','b'); %var[mn]
%     plot(meanStates(5,:),'--','LineWidth',2,'color',[69 181 80]/255); %var[C]
%     plot(meanStates(6,:),'--','LineWidth',2,'color',[255 139 37]/255); %var[P]
    
    refline(0,1)
    ax=gca;
    ax.YAxis.Scale= 'log';
    ylim([0.01 2*10^5])
    yticks([0.01 1 10^2 10^4])
    ylabel('state')
    ylims=ylim;
    xlabel 'census'
    xticks([1:NumYears])
    xticklabels(YearLabels)
    xlim([1,NumYears])
    title({['1250 quadrats'];[num2str(m_perms(mperm)*100) '% individuals observed']},'fontweight','normal');
    text(1-(NumYears-1)*0.15,ylims(2)*3,char(64+i),'Fontsize',16)
    i=i+1;
    
    %     yyaxis right
    %     plot([1:NumYears],ke_yr,'m');
    %     ylabel 'E[k^e]'
    
    %test 3: get selected richness estimates on subsampled data:
    %ksub_perms=[500 20 5]; %pick out of 1250 quadrats
    %ksub_sub_perms=[2 6 9]; %pick out of 10 quadrats in each quadrat
    
    %figure('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3)/3 scrsz(3)]);
    for testperm=1:length(ksub_perms)
        subplot(4,length(ksub_perms)+1,(mperm-1)*(length(ksub_perms)+1)+testperm+1) %plot estimated richness by year
        %yyaxis left
        hold on
        ksub=ksub_perms(testperm);
        %ksub_sub=ksub_sub_perms(testperm);
        
        Richness_raw_sub=zeros(numResample,NumYears); %point estimates
        Chao1_sub=zeros(numResample,NumYears); %point estimates
        Chao2_sub=zeros(numResample,NumYears); %point estimates
        Richness_apx_sub=zeros(numResample,NumYears); %point estimates
        Richness_taylor_sub=zeros(numResample,NumYears); %point estimates
        Richness_raw_sub_mean=zeros(numResample,NumYears); %bootstrap means across resamples
        Chao1_sub_mean=zeros(numResample,NumYears);
        Chao2_sub_mean=zeros(numResample,NumYears);
        Richness_apx_sub_mean=zeros(numResample,NumYears);
        Richness_taylor_sub_mean=zeros(numResample,NumYears);
        Richness_raw_sub_lo=zeros(numResample,NumYears); %bootstrap lower bound across resamples (95%)
        Chao1_sub_lo=zeros(numResample,NumYears);
        Chao2_sub_lo=zeros(numResample,NumYears);
        Richness_apx_sub_lo=zeros(numResample,NumYears);
        Richness_taylor_sub_lo=zeros(numResample,NumYears);
        Richness_raw_sub_up=zeros(numResample,NumYears); %bootstrap upper bounds across resamples (95%)
        Chao1_sub_up=zeros(numResample,NumYears);
        Chao2_sub_up=zeros(numResample,NumYears);
        Richness_apx_sub_up=zeros(numResample,NumYears);
        Richness_taylor_sub_up=zeros(numResample,NumYears);
        meanStates=zeros(5,numResample,NumYears);
        ke_yr=zeros(1,NumYears);
        
        %get temporal richness slope for raw and centred Chao1 and apx
        %boots in subsampling experiments:
        subBoots=round(numBoot/numResample);
        RichnessSlope_raw=zeros(numResample,1);
        RichnessSlopeSD_raw=zeros(numResample,1);
        RichnessSlopeP_raw=zeros(numResample,1);
        RichnessSlope_Chao1=zeros(numResample,1);
        RichnessSlopeSD_Chao1=zeros(numResample,1);
        RichnessSlopeP_Chao1=zeros(numResample,1);
        RichnessSlope_Chao2=zeros(numResample,1);
        RichnessSlopeSD_Chao2=zeros(numResample,1);
        RichnessSlopeP_Chao2=zeros(numResample,1);
        RichnessSlope_Apx=zeros(numResample,1);
        RichnessSlopeSD_Apx=zeros(numResample,1);
        RichnessSlopeP_Apx=zeros(numResample,1);
        RichnessSlope_Taylor=zeros(numResample,1);
        RichnessSlopeSD_Taylor=zeros(numResample,1);
        RichnessSlopeP_Taylor=zeros(numResample,1);
        
        RichnessSlope_raw_boot=zeros(numResample,subBoots);
        RichnessSlopeSD_raw_boot=zeros(numResample,subBoots);
        RichnessSlopeP_raw_boot=zeros(numResample,subBoots);
        RichnessSlope_Chao1_boot=zeros(numResample,subBoots);
        RichnessSlopeSD_Chao1_boot=zeros(numResample,subBoots);
        RichnessSlopeP_Chao1_boot=zeros(numResample,subBoots);
        RichnessSlope_Chao2_boot=zeros(numResample,subBoots);
        RichnessSlopeSD_Chao2_boot=zeros(numResample,subBoots);
        RichnessSlopeP_Chao2_boot=zeros(numResample,subBoots);
        RichnessSlope_Apx_boot=zeros(numResample,subBoots);
        RichnessSlopeSD_Apx_boot=zeros(numResample,subBoots);
        RichnessSlopeP_Apx_boot=zeros(numResample,subBoots);
        RichnessSlope_Taylor_boot=zeros(numResample,subBoots);
        RichnessSlopeSD_Taylor_boot=zeros(numResample,subBoots);
        RichnessSlopeP_Taylor_boot=zeros(numResample,subBoots);
        
        expectedRichness_raw_yrs=zeros(numBoot,NumYears);
        expectedChao1_yrs=zeros(numBoot,NumYears);
        expectedChao2_yrs=zeros(numBoot,NumYears);
        expectedRichness_apx_yrs=zeros(numBoot,NumYears);
        expectedRichness_taylor_yrs=zeros(numBoot,NumYears);
        
        for resample=1:numResample
            expectedRichness_raw=zeros(subBoots,NumYears);
            expectedChao1=zeros(subBoots,NumYears);
            expectedChao2=zeros(subBoots,NumYears);
            expectedRichness_apx=zeros(subBoots,NumYears);
            expectedRichness_taylor=zeros(subBoots,NumYears);
            %quadratsubID=randperm(numQuad,ksub); %pick random quadrats without replacement for subsampling, same across years
            quadratsubID=quadratsubIDs{testperm,resample}; %load quadrat set
            for yr=1:NumYears
                quadratAbundance_sub=QuadratAbundance_m(quadratsubID,:,yr);
                [Richness_raw_sub(resample,yr),Chao1_sub(resample,yr),~,Chao2_sub(resample,yr),~,~,~,Richness_apx_sub(resample,yr),Richness_taylor_sub(resample,yr),~,expectedRichness_raw_boot,expectedChao1_boot,~,expectedChao2_boot,~,~,~,expectedRichness_apx_boot,expectedRichness_taylor_boot,~,meanStates(:,resample,yr)] = bootRichnessEsts(quadratAbundance_sub,subBoots);
                expectedRichness_raw(:,yr)=expectedRichness_raw_boot;
                expectedChao1(:,yr)=expectedChao1_boot;
                expectedChao2(:,yr)=expectedChao2_boot;
                expectedRichness_apx(:,yr)=expectedRichness_apx_boot;
                expectedRichness_taylor(:,yr)=expectedRichness_taylor_boot;
                expectedRichness_raw_yrs((resample-1)*subBoots+1:resample*subBoots,yr)=expectedRichness_raw_boot;
                expectedChao1_yrs((resample-1)*subBoots+1:resample*subBoots,yr)=expectedChao1_boot;
                expectedChao2_yrs((resample-1)*subBoots+1:resample*subBoots,yr)=expectedChao2_boot;
                expectedRichness_apx_yrs((resample-1)*subBoots+1:resample*subBoots,yr)=expectedRichness_apx_boot;
                expectedRichness_taylor_yrs((resample-1)*subBoots+1:resample*subBoots,yr)=expectedRichness_taylor_boot;
                Richness_raw_sub_mean(resample,yr)=mean(expectedRichness_raw_boot);
                Chao1_sub_mean(resample,yr)=mean(expectedChao1_boot);
                Chao2_sub_mean(resample,yr)=mean(expectedChao2_boot);
                Richness_apx_sub_mean(resample,yr)=mean(expectedRichness_apx_boot);
                Richness_taylor_sub_mean(resample,yr)=mean(expectedRichness_taylor_boot);
                Richness_raw_sub_lo(resample,yr)=prctile(expectedRichness_raw_boot,2.5);
                Chao1_sub_lo(resample,yr)=prctile(expectedChao1_boot,2.5);
                Chao2_sub_lo(resample,yr)=prctile(expectedChao2_boot,2.5);
                Richness_apx_sub_lo(resample,yr)=prctile(expectedRichness_apx_boot,2.5);
                Richness_taylor_sub_lo(resample,yr)=prctile(expectedRichness_taylor_boot,2.5);
                Richness_raw_sub_up(resample,yr)=prctile(expectedRichness_raw_boot,97.5);
                Chao1_sub_up(resample,yr)=prctile(expectedChao1_boot,97.5);
                Chao2_sub_up(resample,yr)=prctile(expectedChao2_boot,97.5);
                Richness_apx_sub_up(resample,yr)=prctile(expectedRichness_apx_boot,97.5);
                Richness_taylor_sub_up(resample,yr)=prctile(expectedRichness_taylor_boot,97.5);
            end
            RawDiff=(Richness_raw_sub(resample,:)-mean(expectedRichness_raw));
            Chao1Diff=(Chao1_sub(resample,:)-mean(expectedChao1));
            Chao2Diff=(Chao2_sub(resample,:)-mean(expectedChao2));
            ApxDiff=(Richness_apx_sub(resample,:)-mean(expectedRichness_apx));
            TaylorDiff=(Richness_taylor_sub(resample,:)-mean(expectedRichness_taylor));
            mdl_Richness_raw=fitlm([1:NumYears],Richness_raw_sub(resample,:));
            mdl_Chao1=fitlm([1:NumYears],Chao1_sub(resample,:));
            mdl_Chao2=fitlm([1:NumYears],Chao2_sub(resample,:));
            mdl_Apx=fitlm([1:NumYears],Richness_apx_sub(resample,:));
            mdl_Taylor=fitlm([1:NumYears],Richness_taylor_sub(resample,:));
            RichnessSlope_raw(resample)=mdl_Richness_raw.Coefficients.Estimate(2);
            RichnessSlopeSD_raw(resample)=mdl_Richness_raw.Coefficients.SE(2)*sqrt(NumYears);
            RichnessSlopeP_raw(resample)=mdl_Richness_raw.Coefficients.pValue(2);
            RichnessSlope_Chao1(resample)=mdl_Chao1.Coefficients.Estimate(2);
            RichnessSlopeSD_Chao1(resample)=mdl_Chao1.Coefficients.SE(2)*sqrt(NumYears);
            RichnessSlopeP_Chao1(resample)=mdl_Chao1.Coefficients.pValue(2);
            RichnessSlope_Chao2(resample)=mdl_Chao2.Coefficients.Estimate(2);
            RichnessSlopeSD_Chao2(resample)=mdl_Chao2.Coefficients.SE(2)*sqrt(NumYears);
            RichnessSlopeP_Chao2(resample)=mdl_Chao2.Coefficients.pValue(2);
            RichnessSlope_Apx(resample)=mdl_Apx.Coefficients.Estimate(2);
            RichnessSlopeSD_Apx(resample)=mdl_Apx.Coefficients.SE(2)*sqrt(NumYears);
            RichnessSlopeP_Apx(resample)=mdl_Apx.Coefficients.pValue(2);
            RichnessSlope_Taylor(resample)=mdl_Taylor.Coefficients.Estimate(2);
            RichnessSlopeSD_Taylor(resample)=mdl_Taylor.Coefficients.SE(2)*sqrt(NumYears);
            RichnessSlopeP_Taylor(resample)=mdl_Taylor.Coefficients.pValue(2);
            for boot=1:subBoots
                %if addBootUncert==0
                %else
                mdl_Richness_raw_boot=fitlm([1:NumYears],expectedRichness_raw(boot,:)+RawDiff);
                mdl_Chao1_boot=fitlm([1:NumYears],expectedChao1(boot,:)+Chao1Diff);
                mdl_Chao2_boot=fitlm([1:NumYears],expectedChao2(boot,:)+Chao2Diff);
                mdl_Apx_boot=fitlm([1:NumYears],expectedRichness_apx(boot,:)+ApxDiff);
                mdl_Taylor_boot=fitlm([1:NumYears],expectedRichness_taylor(boot,:)+TaylorDiff);
                %end
                RichnessSlope_raw_boot(resample,boot)=mdl_Richness_raw_boot.Coefficients.Estimate(2);
                RichnessSlopeSD_raw_boot(resample,boot)=mdl_Richness_raw_boot.Coefficients.SE(2)*sqrt(NumYears);
                RichnessSlopeP_raw_boot(resample,boot)=mdl_Richness_raw_boot.Coefficients.pValue(2);
                RichnessSlope_Chao1_boot(resample,boot)=mdl_Chao1_boot.Coefficients.Estimate(2);
                RichnessSlopeSD_Chao1_boot(resample,boot)=mdl_Chao1_boot.Coefficients.SE(2)*sqrt(NumYears);
                RichnessSlopeP_Chao1_boot(resample,boot)=mdl_Chao1_boot.Coefficients.pValue(2);
                RichnessSlope_Chao2_boot(resample,boot)=mdl_Chao2_boot.Coefficients.Estimate(2);
                RichnessSlopeSD_Chao2_boot(resample,boot)=mdl_Chao2_boot.Coefficients.SE(2)*sqrt(NumYears);
                RichnessSlopeP_Chao2_boot(resample,boot)=mdl_Chao2_boot.Coefficients.pValue(2);
                RichnessSlope_Apx_boot(resample,boot)=mdl_Apx_boot.Coefficients.Estimate(2);
                RichnessSlopeSD_Apx_boot(resample,boot)=mdl_Apx_boot.Coefficients.SE(2)*sqrt(NumYears);
                RichnessSlopeP_Apx_boot(resample,boot)=mdl_Apx_boot.Coefficients.pValue(2);
                RichnessSlope_Taylor_boot(resample,boot)=mdl_Taylor_boot.Coefficients.Estimate(2);
                RichnessSlopeSD_Taylor_boot(resample,boot)=mdl_Taylor_boot.Coefficients.SE(2)*sqrt(NumYears);
                RichnessSlopeP_Taylor_boot(resample,boot)=mdl_Taylor_boot.Coefficients.pValue(2);
            end
        end
        %record estimates for subsampled+downsampled dataset:
        allSub_Chao1((mperm-1)*(length(ksub_perms)+1)+testperm+1,:)=mean(Chao1_sub);
        allSub_Chao2((mperm-1)*(length(ksub_perms)+1)+testperm+1,:)=mean(Chao2_sub);
        allSub_Richness_apx((mperm-1)*(length(ksub_perms)+1)+testperm+1,:)=mean(Richness_apx_sub);
        allSub_Richness_raw((mperm-1)*(length(ksub_perms)+1)+testperm+1,:)=mean(Richness_raw_sub);
        allSub_Chao1SD((mperm-1)*(length(ksub_perms)+1)+testperm+1,:)=std(expectedChao1_yrs);
        allSub_Chao2SD((mperm-1)*(length(ksub_perms)+1)+testperm+1,:)=std(expectedChao2_yrs);
        allSub_Richness_apxSD((mperm-1)*(length(ksub_perms)+1)+testperm+1,:)=std(expectedRichness_apx_yrs);
        allSub_Richness_rawSD((mperm-1)*(length(ksub_perms)+1)+testperm+1,:)=std(expectedRichness_raw_yrs);
        
        
        
        %record trend estimates for subsampled+downsampled dataset:
        allSubSlope_Chao1((mperm-1)*(length(ksub_perms)+1)+testperm+1,:)=mean(RichnessSlope_Chao1_boot(:));
        allSubSlopeSD_Chao1((mperm-1)*(length(ksub_perms)+1)+testperm+1,:)=mean(RichnessSlopeSD_Chao1_boot(:));
        allSubSlopeP_Chao1((mperm-1)*(length(ksub_perms)+1)+testperm+1,:)=mean(RichnessSlopeP_Chao1_boot(:));
        allSubSlopePPt_Chao1((mperm-1)*(length(ksub_perms)+1)+testperm+1,:)=mean(RichnessSlopeP_Chao1);
        allSubSlopeSDP_Chao1((mperm-1)*(length(ksub_perms)+1)+testperm+1,:)=std(RichnessSlopeP_Chao1_boot(:));
        allSubSlopeBeta_Chao1((mperm-1)*(length(ksub_perms)+1)+testperm+1,:)=sum(RichnessSlopeP_Chao1_boot(:)>0.05)/numBoot;
        allSubSlopeBetaPt_Chao1((mperm-1)*(length(ksub_perms)+1)+testperm+1,:)=sum(RichnessSlopeP_Chao1(:)>0.05)/numResample;
        
        allSubSlope_Chao2((mperm-1)*(length(ksub_perms)+1)+testperm+1,:)=mean(RichnessSlope_Chao2_boot(:));
        allSubSlopeSD_Chao2((mperm-1)*(length(ksub_perms)+1)+testperm+1,:)=mean(RichnessSlopeSD_Chao2_boot(:));
        allSubSlopeP_Chao2((mperm-1)*(length(ksub_perms)+1)+testperm+1,:)=mean(RichnessSlopeP_Chao2_boot(:));
        allSubSlopePPt_Chao2((mperm-1)*(length(ksub_perms)+1)+testperm+1,:)=mean(RichnessSlopeP_Chao2);
        allSubSlopeSDP_Chao2((mperm-1)*(length(ksub_perms)+1)+testperm+1,:)=std(RichnessSlopeP_Chao2_boot(:));
        allSubSlopeBeta_Chao2((mperm-1)*(length(ksub_perms)+1)+testperm+1,:)=sum(RichnessSlopeP_Chao2_boot(:)>0.05)/numBoot;
        allSubSlopeBetaPt_Chao2((mperm-1)*(length(ksub_perms)+1)+testperm+1,:)=sum(RichnessSlopeP_Chao2(:)>0.05)/numResample;
        
        allSubSlope_apx((mperm-1)*(length(ksub_perms)+1)+testperm+1,:)=mean(RichnessSlope_Apx_boot(:));
        allSubSlopeSD_apx((mperm-1)*(length(ksub_perms)+1)+testperm+1,:)=mean(RichnessSlopeSD_Apx_boot(:));
        allSubSlopeP_apx((mperm-1)*(length(ksub_perms)+1)+testperm+1,:)=mean(RichnessSlopeP_Apx_boot(:));
        allSubSlopePPt_apx((mperm-1)*(length(ksub_perms)+1)+testperm+1,:)=mean(RichnessSlopeP_Apx);
        allSubSlopeSDP_apx((mperm-1)*(length(ksub_perms)+1)+testperm+1,:)=std(RichnessSlopeP_Apx_boot(:));
        allSubSlopeBeta_apx((mperm-1)*(length(ksub_perms)+1)+testperm+1,:)=sum(RichnessSlopeP_Apx_boot(:)>0.05)/numBoot;
        allSubSlopeBetaPt_apx((mperm-1)*(length(ksub_perms)+1)+testperm+1,:)=sum(RichnessSlopeP_Apx(:)>0.05)/numResample;
        
        allSubSlope_raw((mperm-1)*(length(ksub_perms)+1)+testperm+1,:)=mean(RichnessSlope_raw_boot(:));
        allSubSlopeSD_raw((mperm-1)*(length(ksub_perms)+1)+testperm+1,:)=mean(RichnessSlopeSD_raw_boot(:));
        allSubSlopeP_raw((mperm-1)*(length(ksub_perms)+1)+testperm+1,:)=mean(RichnessSlopeP_raw_boot(:));
        allSubSlopePPt_raw((mperm-1)*(length(ksub_perms)+1)+testperm+1,:)=mean(RichnessSlopeP_raw);
        allSubSlopeSDP_raw((mperm-1)*(length(ksub_perms)+1)+testperm+1,:)=std(RichnessSlopeP_raw_boot(:));
        allSubSlopeBeta_raw((mperm-1)*(length(ksub_perms)+1)+testperm+1,:)=sum(RichnessSlopeP_raw_boot(:)>0.05)/numBoot;
        allSubSlopeBetaPt_raw((mperm-1)*(length(ksub_perms)+1)+testperm+1,:)=sum(RichnessSlopeP_raw(:)>0.05)/numResample;
        
        %         plot(Richness_apx_sub,'r','LineWidth',2);
        %         plot(Chao1_sub,'b','LineWidth',2);
        %         plot(Richness_raw_sub,'--k','LineWidth',2);
        %boundedline([1:NumYears]', mean(expectedRichness_raw)',std(expectedRichness_raw)','-k','alpha','transparency', 0.1);
        %boundedline([1:NumYears]', mean(expectedChao1)',std(expectedChao1)','-b','alpha','transparency', 0.1);
        %boundedline([1:NumYears]', mean(expectedRichness_apx)',std(expectedRichness_apx)','-r','alpha','transparency', 0.1);
        %     boundedline([1:NumYears]', 100*mean(Richness_raw_sub./Richness_raw)',100*[max(mean(Richness_raw_sub./Richness_raw)'-prctile(Richness_raw_sub./Richness_raw,2.5)',0),prctile(Richness_raw_sub./Richness_raw,97.5)'-mean(Richness_raw_sub./Richness_raw)'],'-k','alpha','transparency', 0.1);
        %     boundedline([1:NumYears]', 100*mean(Chao1_sub_mean./mean(expectedChao1))',100*[mean(Chao1_sub_mean./mean(expectedChao1))'-mean(Chao1_sub_lo./mean(expectedChao1))',mean(Chao1_sub_up./mean(expectedChao1))'-mean(Chao1_sub_mean./mean(expectedChao1))'],'-b','alpha','transparency', 0.1);
        %     boundedline([1:NumYears]', 100*mean(Richness_apx_sub_mean./mean(expectedRichness_apx))',100*[mean(Richness_apx_sub_mean./mean(expectedRichness_apx))'-mean(Richness_apx_sub_lo./mean(expectedRichness_apx))',mean(Richness_apx_sub_up./mean(expectedRichness_apx))'-mean(Richness_apx_sub_mean./mean(expectedRichness_apx))'],'-r','alpha','transparency', 0.1);
        %bl1=boundedline([1:NumYears]', mean(Richness_raw_sub)',[max(mean(Richness_raw_sub)'-prctile(Richness_raw_sub,2.5)',0),prctile(Richness_raw_sub,97.5)'-mean(Richness_raw_sub)'],'-k','alpha','transparency', 0.1);
        %         boundedline([1:NumYears]', mean(Chao1_sub_mean)',[mean(Chao1_sub_mean)'-mean(Chao1_sub_lo)',mean(Chao1_sub_up)'-mean(Chao1_sub_mean)'],'-b','alpha','transparency', 0.1);
        %         boundedline([1:NumYears]', mean(Richness_apx_sub_mean)',[mean(Richness_apx_sub_mean)'-mean(Richness_apx_sub_lo)',mean(Richness_apx_sub_up)'-mean(Richness_apx_sub_mean)'],'-r','alpha','transparency', 0.1);
        %boundedline([1:NumYears]', mean(Richness_raw_sub_mean)',[mean(Richness_raw_sub_mean)'-mean(Richness_raw_sub_lo)',mean(Richness_raw_sub_up)'-mean(Richness_raw_sub_mean)'],'--k','alpha','transparency', 0.1);
        %centre confidence bounds to match bootstrapped mean to spot estimates
        %on the original dataset:
        boundedline([1:NumYears]', mean(Chao2_sub)',[max(mean(Chao2_sub_mean-Chao2_sub_lo)',0),mean(Chao2_sub_up-Chao2_sub_mean)'],'--c','alpha','transparency', 0.1);
        boundedline([1:NumYears]', mean(Chao1_sub)',[max(mean(Chao1_sub_mean-Chao1_sub_lo)',0),mean(Chao1_sub_up-Chao1_sub_mean)'],'-b','alpha','transparency', 0.1);
        %boundedline([1:NumYears]', mean(Richness_apx_sub)',[max(mean(Richness_apx_sub_mean)'-mean(Richness_apx_sub_lo)',0),mean(Richness_apx_sub_up)'-mean(Richness_apx_sub_mean)'],'-r','alpha','transparency', 0.1);
        boundedline([1:NumYears]', mean(Richness_apx_sub)',[max(mean(Richness_apx_sub_mean-Richness_apx_sub_lo)',0),mean(Richness_apx_sub_up-Richness_apx_sub_mean)'],'-r','alpha','transparency', 0.1);
        boundedline([1:NumYears]', mean(Richness_raw_sub)',[max(mean(Richness_raw_sub_mean-Richness_raw_sub_lo)',0),mean(Richness_raw_sub_up-Richness_raw_sub_mean)'],'--k','alpha','transparency', 0.1);
%         boundedline([1:NumYears]', mean(Richness_apx_sub)',[max(mean(expectedRichness_taylor_yrs)'-prctile(expectedRichness_taylor_yrs,2.5)',0),prctile(expectedRichness_taylor_yrs,97.5)'-mean(expectedRichness_taylor_yrs)'],'-r','alpha','transparency', 0.1);
%         boundedline([1:NumYears]', mean(Chao1_sub)',[max(mean(expectedChao1_yrs)'-prctile(expectedChao1_yrs,2.5)',0),prctile(expectedChao1_yrs,97.5)'-mean(expectedChao1_yrs)'],'-b','alpha','transparency', 0.1);
%         boundedline([1:NumYears]', mean(Chao2_sub)',[max(mean(expectedChao2_yrs)'-prctile(expectedChao2_yrs,2.5)',0),prctile(expectedChao1_yrs,97.5)'-mean(expectedChao1_yrs)'],'--c','alpha','transparency', 0.1);
%         boundedline([1:NumYears]', mean(Richness_raw_sub)',[max(mean(expectedRichness_taylor_yrs)'-prctile(expectedRichness_taylor_yrs,2.5)',0),prctile(expectedRichness_taylor_yrs,97.5)'-mean(expectedRichness_taylor_yrs)'],'--k','alpha','transparency', 0.1);
        
        plot(mean(Chao1_sub),'-b','LineWidth',2);
        plot(mean(Chao2_sub),'--','LineWidth',2,'Color',[0 0.7 0.7]);
        plot(mean(Richness_apx_sub),'-r','LineWidth',2);
        plot(mean(Richness_raw_sub),'--k','LineWidth',2);
        %         plot(mean(Chao1_sub_mean),'b','LineWidth',2);
        %         plot(mean(Richness_apx_sub_mean),'r','LineWidth',2);
        %         plot(mean(Richness_raw_sub_mean),'--k','LineWidth',2);
        xlabel 'census'
        %ylim([0 120])
        %ylabel '% richness recovered'
        ylim([minRich maxRich])
        ylabel 'richness'
        title({[num2str(ksub) ' quadrats'];[num2str(m_perms(mperm)*100) '% individuals observed']},'fontweight','normal')
        xticks([1:NumYears])
        xticklabels(YearLabels)
        xlim([1,NumYears])
%         text(1.2,minRich+diffRich*.95,['\DeltaChao1=' num2str(mean(RichnessSlope_Chao1_boot(:)),1) '\pm' num2str(mean(RichnessSlopeSD_Chao1_boot(:)),1) ', p=' num2str(mean(RichnessSlopeP_Chao1_boot(:)),1) '(' num2str(mean(RichnessSlopeP_Chao1),1) ')\pm' num2str(std(RichnessSlopeP_Chao1_boot(:)),1) ', \beta=' num2str(sum(RichnessSlopeP_Chao1_boot(:)>0.05)/numBoot,1) '(' num2str(sum(RichnessSlopeP_Chao1(:)>0.05)/numResample,1) ')'],'fontsize',14,'color','b');
%         text(1.2,minRich+diffRich*.85,['\DeltaChao2=' num2str(mean(RichnessSlope_Chao2_boot(:)),1) '\pm' num2str(mean(RichnessSlopeSD_Chao2_boot(:)),1) ', p=' num2str(mean(RichnessSlopeP_Chao2_boot(:)),1) '(' num2str(mean(RichnessSlopeP_Chao2),1) ')\pm' num2str(std(RichnessSlopeP_Chao2_boot(:)),1) ', \beta=' num2str(sum(RichnessSlopeP_Chao2_boot(:)>0.05)/numBoot,1) '(' num2str(sum(RichnessSlopeP_Chao2(:)>0.05)/numResample,1) ')'],'fontsize',14,'Color',[0 0.7 0.7]);
%         text(1.2,minRich+diffRich*.73,['\Delta\Omega_T=' num2str(mean(RichnessSlope_Taylor_boot(:)),1) '\pm' num2str(mean(RichnessSlopeSD_Taylor_boot(:)),1) ', p=' num2str(mean(RichnessSlopeP_Taylor_boot(:)),1) '(' num2str(mean(RichnessSlopeP_Taylor),1) ')\pm' num2str(std(RichnessSlopeP_Taylor_boot(:)),1) ', \beta=' num2str(sum(RichnessSlopeP_Taylor_boot(:)>0.05)/numBoot,1) '(' num2str(sum(RichnessSlopeP_Taylor(:)>0.05)/numResample,1) ')'],'fontsize',14,'Color','r');
%         text(1.2,minRich+diffRich*.63,['\DeltaRaw=' num2str(mean(RichnessSlope_raw_boot(:)),1) '\pm' num2str(mean(RichnessSlopeSD_raw_boot(:)),1) ', p=' num2str(mean(RichnessSlopeP_raw_boot(:)),1) '(' num2str(mean(RichnessSlopeP_raw),1) ')\pm' num2str(std(RichnessSlopeP_raw_boot(:)),1) ', \beta=' num2str(sum(RichnessSlopeP_raw_boot(:)>0.05)/numBoot,1) '(' num2str(sum(RichnessSlopeP_raw(:)>0.05)/numResample,1) ')' ],'fontsize',14);
    
        %         text(1,2,maxRich*.95,['bootstrapped estimates:'],'fontsize',14)
        %         text(1.2,maxRich*.85,['raw slope=' num2str(mean(RichnessSlope_raw_boot(:)),1) '\pm' num2str(mean(RichnessSlopeSD_raw_boot(:)),1) ', p=' num2str(mean(RichnessSlopeP_raw_boot(:)),1) '\pm' num2str(std(RichnessSlopeP_raw_boot(:)),1) ', \beta=' num2str(sum(RichnessSlopeP_raw_boot(:)>0.05)/numBoot,1)],'fontsize',14);
        %         text(1.2,maxRich*.75,['Chao1 slope=' num2str(mean(RichnessSlope_Chao1_boot(:)),1) '\pm' num2str(mean(RichnessSlopeSD_Chao1_boot(:)),1) ', p=' num2str(mean(RichnessSlopeP_Chao1_boot(:)),1) '\pm' num2str(std(RichnessSlopeP_Chao1_boot(:)),1) ', \beta=' num2str(sum(RichnessSlopeP_Chao1_boot(:)>0.05)/numBoot,1)],'fontsize',14,'color','b');
        %         text(1.2,maxRich*.63,['\Omega_o slope=' num2str(mean(RichnessSlope_Apx_boot(:)),1) '\pm' num2str(mean(RichnessSlopeSD_Apx_boot(:)),1) ', p=' num2str(mean(RichnessSlopeP_Apx_boot(:)),1) '\pm' num2str(std(RichnessSlopeP_Apx_boot(:)),1) ', \beta=' num2str(sum(RichnessSlopeP_Apx_boot(:)>0.05)/numBoot,1)],'fontsize',14,'Color','r');
        %         text(1,2,maxRich*.45,['point estimates:'],'fontsize',14)
        %         text(1.2,maxRich*.35,['raw slope=' num2str(mean(RichnessSlope_raw),1) '\pm' num2str(mean(RichnessSlopeSD_raw),1) ', p=' num2str(mean(RichnessSlopeP_raw),1) '\pm' num2str(std(RichnessSlopeP_raw),1) ', \beta=' num2str(sum(RichnessSlopeP_raw(:)>0.05)/numBoot,1)],'fontsize',14);
        %         text(1.2,maxRich*.25,['Chao1 slope=' num2str(mean(RichnessSlope_Chao1),1) '\pm' num2str(mean(RichnessSlopeSD_Chao1),1) ', p=' num2str(mean(RichnessSlopeP_Chao1),1) '\pm' num2str(std(RichnessSlopeP_Chao1),1) ', \beta=' num2str(sum(RichnessSlopeP_Chao1(:)>0.05)/numBoot,1)],'fontsize',14,'color','b');
        %         text(1.2,maxRich*.13,['\Omega_o slope=' num2str(mean(RichnessSlope_Apx),1) '\pm' num2str(mean(RichnessSlopeSD_Apx),1) ', p=' num2str(mean(RichnessSlopeP_Apx),1) '\pm' num2str(std(RichnessSlopeP_Apx),1) ', \beta=' num2str(sum(RichnessSlopeP_Apx(:)>0.05)/numBoot,1)],'fontsize',14,'Color','r');
        text(1-(NumYears-1)*0.15,maxRich*1.15,char(64+i),'Fontsize',16)
        
        %plot state variable means and variances over years (averaged over
        %replicates)
        subplot(4,length(ksub_perms)+1,(mperm-1)*(length(ksub_perms)+1)+testperm+1+6) %plot estimated richness by year
        hold on
        plot(reshape(mean(meanStates(1,:,:),2),[],1),'-','LineWidth',2,'color','b'); %E[mn]
        plot(reshape(mean(meanStates(2,:,:),2),[],1),'-','LineWidth',2,'color',[255 139 37]/255); %E[P]
        plot(reshape(mean(meanStates(3,:,:),2),[],1),'--','LineWidth',2,'color','b'); %var[mn]
        plot(reshape(mean(meanStates(4,:,:),2),[],1),'--','LineWidth',2,'color',[255 139 37]/255); %var[P]
        plot(reshape(mean(meanStates(5,:,:),2),[],1),'-','LineWidth',2,'color',[69 181 80]/255); %cov[mn,P]
        %         plot(reshape(mean(meanStates(1,:,:),2),[],1),'-','LineWidth',2,'color','b'); %E[mn]
        %         plot(reshape(mean(meanStates(2,:,:),2),[],1),'-','LineWidth',2,'color',[69 181 80]/255); %E[C]
        %         plot(reshape(mean(meanStates(3,:,:),2),[],1),'-','LineWidth',2,'color',[255 139 37]/255); %E[P]
        %         plot(reshape(mean(meanStates(4,:,:),2),[],1),'--','LineWidth',2,'color','b'); %var[mn]
        %         plot(reshape(mean(meanStates(5,:,:),2),[],1),'--','LineWidth',2,'color',[69 181 80]/255); %var[C]
        %         plot(reshape(mean(meanStates(6,:,:),2),[],1),'--','LineWidth',2,'color',[255 139 37]/255); %var[P]
        refline(0,1)
        ax=gca;
        ax.YAxis.Scale= 'log';
        ylim([0.01 2*10^5])
        yticks([0.01 1 10^2 10^4])
        ylabel('state')
        ylims=ylim;
        xlabel 'census'
        xticks([1:NumYears])
        xticklabels(YearLabels)
        xlim([1,NumYears])
        title({[num2str(ksub) ' quadrats'];[num2str(m_perms(mperm)*100) '% individuals observed']},'fontweight','normal')
        text(1-(NumYears-1)*0.15,ylims(2)*3,char(64+i),'Fontsize',16)
        if i==6
            lgd=legend('E[$$\hat{mn}$$]','E[$$\hat{P}$$]','var[$$\hat{mn}$$]','var[$$\hat{P}$$]','cov[$$\hat{mn}$$,$$\hat{P}$$]','Interpreter','Latex');
            %lgd=legend('E[$$\hat{mn}$$]','E[$$\hat{C}$$]','E[$$\hat{P}$$]','var[$$\hat{mn}$$]','var[$$\hat{C}$$]','var[$$\hat{P}$$]','Interpreter','Latex');
            lgd.AutoUpdate='off';
        end
        i=i+1;
        
        %         yyaxis right
        %         plot([1:NumYears],ke_yr,'m');
        %         ylabel 'E[k^e]'
        %     text(1.3,97,['Taylor2 subsampled \Deltarichness/yr=' num2str(mean(RichnessSlope_apx_sub),2) '\pm' num2str(std(RichnessSlope_apx_sub),2) ' (p=' num2str(mean(RichnessSlopeP_apx_sub),2) '\pm' num2str(std(RichnessSlopeP_apx_sub),2) ')'],'color','r')
        %     text(1.3,92,['Chao1 subsampled \Deltarichness/yr=' num2str(mean(Chao1Slope_sub),2) '\pm' num2str(std(Chao1Slope_sub),2) ' (p=' num2str(mean(Chao1SlopeP_sub),2) '\pm' num2str(std(Chao1SlopeP_sub),2) ')'],'color','b')
        %     text(1.3,87,['raw subsampled \Deltarichness/yr=' num2str(mean(RichnessSlope_raw_sub),2) '\pm' num2str(std(RichnessSlope_raw_sub),2) ' (p=' num2str(mean(RichnessSlopeP_raw_sub),2) '\pm' num2str(std(RichnessSlopeP_raw_sub),2) ')'])
        
        %     subplot(2,4,4+testperm+1)
        %     %bl1=boundedline([1:NumYears]', 100*mean(Richness_raw_sub./Richness_raw)',100*[max(mean(Richness_raw_sub./Richness_raw)'-prctile(Richness_raw_sub./Richness_raw,2.5)',0),prctile(Richness_raw_sub./Richness_raw,97.5)'-mean(Richness_raw_sub./Richness_raw)'],'-k','alpha','transparency', 0.1);
        %     bl2=boundedline([1:NumYears]', 100*mean(Chao1_sub./Chao1_all)',100*[mean(Chao1_sub./Chao1_all)'-mean(Chao1_sub_lo./Chao1_all)',mean(Chao1_sub_up./Chao1_all)'-mean(Chao1_sub./Chao1_all)'],'-b','alpha','transparency', 0.1);
        %     bl3=boundedline([1:NumYears]', 100*mean(Richness_apx_sub./Richness_apx)',100*[mean(Richness_apx_sub./Richness_apx)'-mean(Richness_apx_sub_lo./Richness_apx)',mean(Richness_apx_sub_up./Richness_apx)'-mean(Richness_apx_sub./Richness_apx)'],'-r','alpha','transparency', 0.1);
        %     bl1=boundedline([1:NumYears]', 100*mean(Richness_raw_sub./Richness_raw)',100*[max(mean(Richness_raw_sub./Richness_raw)'-mean(Richness_raw_sub_lo./Richness_raw)',0),mean(Richness_raw_sub_up./Richness_raw)'-mean(Richness_raw_sub./Richness_raw)'],'--k','alpha','transparency', 0.1);
        %
        % %     set(bl1,'linewidth',2);
        % %     set(bl2,'linewidth',2);
        % %     set(bl3,'linewidth',2);
        %     h=refline(0,100)
        %     h.Color='k';
        %     h.LineStyle='--';
        %     xlabel 'year'
        %     ylabel '% richness recovered'
        %     ylim([40 110])
        %     %title(['quadrats subsampled=' num2str(ksub) '/9, quadrats subsampled=' num2str(ksub_sub) '/10'])
        %     xticks([1:NumYears])
        %     xticklabels(YearLabels)
        %     xlim([1,NumYears])
        
    end
end
% 
% %figure('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3)/7 scrsz(4)/6],'DefaultAxesFontSize',16,'defaultAxesColorOrder',[[0 0 0]; [1 0 1]]);
% figure('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3)/7 scrsz(4)],'DefaultAxesFontSize',16,'defaultAxesColorOrder',[[0 0 0]; [1 0 1]]);
% subplot(3,1,1)
% hold on
% %plotOrder=[1 4 2 5 3 6];
% plotOrder=[1 2 3 4 5 6];
% reordered_Chao1=mean(allSub_Chao1(plotOrder,:),2);
% reordered_Chao2=mean(allSub_Chao2(plotOrder,:),2);
% reordered_Richness_apx=mean(allSub_Richness_apx(plotOrder,:),2);
% reordered_Richness_raw=mean(allSub_Richness_raw(plotOrder,:),2);
% plot(100*reordered_Chao1(1:3)./reordered_Chao1(1),'-b','LineWidth',5);
% plot(100*reordered_Chao2(1:3)./reordered_Chao2(1),'--','LineWidth',5,'Color',[0 0.7 0.7]);
% plot(100*reordered_Richness_apx(1:3)./reordered_Richness_apx(1),'-r','LineWidth',5);
% plot(100*reordered_Richness_raw(1:3)./reordered_Richness_raw(1),'--k','LineWidth',5);
% 
% plot(100*reordered_Chao1(4:6)./reordered_Chao1(1),'-b','LineWidth',2);
% plot(100*reordered_Chao2(4:6)./reordered_Chao2(1),'--','LineWidth',2,'Color',[0 0.7 0.7]);
% plot(100*reordered_Richness_apx(4:6)./reordered_Richness_apx(1),'-r','LineWidth',2);
% plot(100*reordered_Richness_raw(4:6)./reordered_Richness_raw(1),'--k','LineWidth',2);
% ylabel('% of full estimates')
% xlim([1 3]);
% xticks([1:3]);
% plotLabels=["1250" "12" "3"];
% xticklabels(plotLabels)
% xlabel('quadrats sampled')
% 
% subplot(3,1,2)
% hold on
% reordered_Chao1_slope=mean(allSubSlope_Chao1(plotOrder,:),2);
% reordered_Chao2_slope=mean(allSubSlope_Chao2(plotOrder,:),2);
% reordered_taylor_slope=mean(allSubSlope_taylor(plotOrder,:),2);
% reordered_raw_slope=mean(allSubSlope_raw(plotOrder,:),2);
% reordered_Chao1_slopeSD=mean(allSubSlopeSD_Chao1(plotOrder,:),2);
% reordered_Chao2_slopeSD=mean(allSubSlopeSD_Chao2(plotOrder,:),2);
% reordered_taylor_slopeSD=mean(allSubSlopeSD_taylor(plotOrder,:),2);
% reordered_raw_slopeSD=mean(allSubSlopeSD_raw(plotOrder,:),2);
% 
% plot(reordered_Chao1_slope(1:3),'-b','LineWidth',5);
% plot(reordered_Chao2_slope(1:3),'--','LineWidth',5,'Color',[0 0.7 0.7]);
% plot(reordered_taylor_slope(1:3),'-r','LineWidth',5);
% plot(reordered_raw_slope(1:3),'--k','LineWidth',5);
% 
% plot(reordered_Chao1_slope(4:6),'-b','LineWidth',2);
% plot(reordered_Chao2_slope(4:6),'--','LineWidth',2,'Color',[0 0.7 0.7]);
% plot(reordered_taylor_slope(4:6),'-r','LineWidth',2);
% plot(reordered_raw_slope(4:6),'--k','LineWidth',2);
% 
% ylabel('\Deltarichnes\\\DeltaT')
% xticks([1:3]);
% plotLabels=["1250" "12" "3"];
% xticklabels(plotLabels)
% xlabel('quadrats sampled')
% 
% subplot(3,1,3)
% hold on
% reordered_Chao1_slopeP=mean(allSubSlopePPt_Chao1(plotOrder,:),2);
% reordered_Chao2_slopeP=mean(allSubSlopePPt_Chao2(plotOrder,:),2);
% reordered_taylor_slopeP=mean(allSubSlopePPt_taylor(plotOrder,:),2);
% reordered_raw_slopeP=mean(allSubSlopePPt_raw(plotOrder,:),2);
% 
% plot(reordered_Chao1_slopeP(1:3),'-b','LineWidth',5);
% plot(reordered_Chao2_slopeP(1:3),'--','LineWidth',5,'Color',[0 0.7 0.7]);
% plot(reordered_taylor_slopeP(1:3),'-r','LineWidth',5);
% plot(reordered_raw_slopeP(1:3),'--k','LineWidth',5);
% 
% plot(reordered_Chao1_slopeP(4:6),'-b','LineWidth',2);
% plot(reordered_Chao2_slopeP(4:6),'--','LineWidth',2,'Color',[0 0.7 0.7]);
% plot(reordered_taylor_slopeP(4:6),'-r','LineWidth',2);
% plot(reordered_raw_slopeP(4:6),'--k','LineWidth',2);
% 
% ylabel('p_{\Deltarichness}')
% xticks([1:3]);
% plotLabels=["1250" "12" "3"];
% xticklabels(plotLabels)
% xlabel('quadrats sampled')


%figure('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3)/7 scrsz(4)/6],'DefaultAxesFontSize',16,'defaultAxesColorOrder',[[0 0 0]; [1 0 1]]);
%figure('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3) scrsz(4)/3],'DefaultAxesFontSize',16,'defaultAxesColorOrder',[[0 0 0]; [1 0 1]]);
figure('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3)/2 scrsz(4)/6],'DefaultAxesFontSize',16,'defaultAxesColorOrder',[[0 0 0]; [0 0 0]]);


%plotOrder=[1 4 2 5 3 6];
plotOrder=[1 2 3 4 5 6];
  boundedline([1:NumYears]', Richness_raw',[max(mean(expectedRichness_raw)'-prctile(expectedRichness_raw,2.5)',0),prctile(expectedRichness_raw,97.5)'-mean(expectedRichness_raw)'],'--k','alpha','transparency', 0.1);
  boundedline([1:NumYears]', mean(Richness_raw_sub)',[max(mean(expectedRichness_taylor_yrs)'-prctile(expectedRichness_taylor_yrs,2.5)',0),prctile(expectedRichness_taylor_yrs,97.5)'-mean(expectedRichness_taylor_yrs)'],'--k','alpha','transparency', 0.1);
  
reordered_Chao1=mean(allSub_Chao1(plotOrder,:),2);
reordered_Chao2=mean(allSub_Chao2(plotOrder,:),2);
reordered_apx=mean(allSub_Richness_apx(plotOrder,:),2);
reordered_raw=mean(allSub_Richness_raw(plotOrder,:),2);
reordered_Chao1SD=mean(allSub_Chao1SD(plotOrder,:),2);
reordered_Chao2SD=mean(allSub_Chao2SD(plotOrder,:),2);
reordered_apxSD=mean(allSub_Richness_apxSD(plotOrder,:),2);
reordered_rawSD=mean(allSub_Richness_rawSD(plotOrder,:),2);
reordered_Chao1CV=mean(allSub_Chao1SD(plotOrder,:)./allSub_Chao1(plotOrder,:),2);
reordered_Chao2CV=mean(allSub_Chao2SD(plotOrder,:)./allSub_Chao2(plotOrder,:),2);
reordered_apxCV=mean(allSub_Richness_apxSD(plotOrder,:)./allSub_Richness_apx(plotOrder,:),2);
reordered_rawCV=mean(allSub_Richness_rawSD(plotOrder,:)./allSub_Richness_raw(plotOrder,:),2);

subplot(1,5,1)
hold on
plot(reordered_Chao1(1:3),'-b','LineWidth',5);
plot(reordered_Chao2(1:3),'--','LineWidth',5,'Color',[0 0.7 0.7]);
plot(reordered_apx(1:3),'-r','LineWidth',5);
plot(reordered_raw(1:3),'--k','LineWidth',5);
% bl1=boundedline([1:3]', reordered_Chao1(1:3)',reordered_Chao1SD(1:3)','-b','alpha','transparency', 0.1);
% set(bl1,'LineWidth',5)
% bl2=boundedline([1:3]', reordered_Chao2(1:3)',reordered_Chao2SD(1:3)','--c','alpha','transparency', 0.1);
% set(bl2,'LineWidth',5)
% bl3=boundedline([1:3]', reordered_taylor(1:3)',reordered_taylorSD(1:3)','-r','alpha','transparency', 0.1);
% set(bl3,'LineWidth',5)
% bl4=boundedline([1:3]', reordered_raw(1:3)',reordered_rawSD(1:3)','--k','alpha','transparency', 0.1);
% set(bl4,'LineWidth',5)

plot(reordered_Chao1(4:6),'-b','LineWidth',2);
plot(reordered_Chao2(4:6),'--','LineWidth',2,'Color',[0 0.7 0.7]);
plot(reordered_apx(4:6),'-r','LineWidth',2);
plot(reordered_raw(4:6),'--k','LineWidth',2);
% bl5=boundedline([1:3]', reordered_Chao1(4:6)',reordered_Chao1SD(4:6)','-b','alpha','transparency', 0.1);
% set(bl5,'LineWidth',2)
% bl6=boundedline([1:3]', reordered_Chao2(4:6)',reordered_Chao2SD(4:6)','--c','alpha','transparency', 0.1);
% set(bl6,'LineWidth',2)
% bl7=boundedline([1:3]', reordered_taylor(4:6)',reordered_taylorSD(4:6)','-r','alpha','transparency', 0.1);
% set(bl7,'LineWidth',2)
% bl8=boundedline([1:3]', reordered_raw(4:6)',reordered_rawSD(4:6)','--k','alpha','transparency', 0.1);
% set(bl8,'LineWidth',2)

ylabel('mean richness')
xlim([1 3]);
xticks([1:3]);
plotLabels=["1250" "12" "3"];
xticklabels(plotLabels)
%xlabel('quadrats sampled')

subplot(1,5,2)
hold on
% plot(reordered_Chao1SD(1:3),'-b','LineWidth',5);
% plot(reordered_Chao2SD(1:3),'--','LineWidth',5,'Color',[0 0.7 0.7]);
% plot(reordered_taylorSD(1:3),'-r','LineWidth',5);
% plot(reordered_rawSD(1:3),'--k','LineWidth',5);
% plot(reordered_Chao1SD(4:6),'-b','LineWidth',2);
% plot(reordered_Chao2SD(4:6),'--','LineWidth',2,'Color',[0 0.7 0.7]);
% plot(reordered_taylorSD(4:6),'-r','LineWidth',2);
% plot(reordered_rawSD(4:6),'--k','LineWidth',2);
%ylabel('mean richness S.D.')

plot(reordered_Chao1CV(1:3),'-b','LineWidth',5);
plot(reordered_Chao2CV(1:3),'--','LineWidth',5,'Color',[0 0.7 0.7]);
plot(reordered_apxCV(1:3),'-r','LineWidth',5);
plot(reordered_rawCV(1:3),'--k','LineWidth',5);
plot(reordered_Chao1CV(4:6),'-b','LineWidth',2);
plot(reordered_Chao2CV(4:6),'--','LineWidth',2,'Color',[0 0.7 0.7]);
plot(reordered_apxCV(4:6),'-r','LineWidth',2);
plot(reordered_rawCV(4:6),'--k','LineWidth',2);
ylabel('mean richness C.V.')

xlim([1 3]);
xticks([1:3]);
plotLabels=["9x10" "9x2" "2x9"];
plotLabels=["1250" "12" "3"];
xticklabels(plotLabels)
%xlabel('quadrats sampled')

subplot(1,5,3)
hold on
%plotOrder=[1 4 2 5 3 6];
plotOrder=[1 2 3 4 5 6];
reordered_Chao1=mean(allSub_Chao1(plotOrder,:),2);
reordered_Chao2=mean(allSub_Chao2(plotOrder,:),2);
reordered_Richness_apx=mean(allSub_Richness_apx(plotOrder,:),2);
reordered_Richness_raw=mean(allSub_Richness_raw(plotOrder,:),2);
plot(100*reordered_Chao1(1:3)./reordered_Chao1(1),'-b','LineWidth',5);
plot(100*reordered_Chao2(1:3)./reordered_Chao2(1),'--','LineWidth',5,'Color',[0 0.7 0.7]);
plot(100*reordered_Richness_apx(1:3)./reordered_Richness_apx(1),'-r','LineWidth',5);
plot(100*reordered_Richness_raw(1:3)./reordered_Richness_raw(1),'--k','LineWidth',5);

plot(100*reordered_Chao1(4:6)./reordered_Chao1(1),'-b','LineWidth',2);
plot(100*reordered_Chao2(4:6)./reordered_Chao2(1),'--','LineWidth',2,'Color',[0 0.7 0.7]);
plot(100*reordered_Richness_apx(4:6)./reordered_Richness_apx(1),'-r','LineWidth',2);
plot(100*reordered_Richness_raw(4:6)./reordered_Richness_raw(1),'--k','LineWidth',2);
ylabel('% of full estimates')
ylim([0 100])
xlim([1 3]);
xticks([1:3]);
plotLabels=["1250" "12" "3"];
xticklabels(plotLabels)
xlabel('quadrats sampled')

subplot(1,5,4)
hold on
reordered_Chao1_slope=mean(allSubSlope_Chao1(plotOrder,:),2);
reordered_Chao2_slope=mean(allSubSlope_Chao2(plotOrder,:),2);
reordered_apx_slope=mean(allSubSlope_apx(plotOrder,:),2);
reordered_raw_slope=mean(allSubSlope_raw(plotOrder,:),2);
reordered_Chao1_slopeSD=mean(allSubSlopeSD_Chao1(plotOrder,:),2);
reordered_Chao2_slopeSD=mean(allSubSlopeSD_Chao2(plotOrder,:),2);
reordered_apx_slopeSD=mean(allSubSlopeSD_apx(plotOrder,:),2);
reordered_raw_slopeSD=mean(allSubSlopeSD_raw(plotOrder,:),2);

plot(reordered_Chao1_slope(1:3),'-b','LineWidth',5);
plot(reordered_Chao2_slope(1:3),'--','LineWidth',5,'Color',[0 0.7 0.7]);
plot(reordered_apx_slope(1:3),'-r','LineWidth',5);
plot(reordered_raw_slope(1:3),'--k','LineWidth',5);

plot(reordered_Chao1_slope(4:6),'-b','LineWidth',2);
plot(reordered_Chao2_slope(4:6),'--','LineWidth',2,'Color',[0 0.7 0.7]);
plot(reordered_apx_slope(4:6),'-r','LineWidth',2);
plot(reordered_raw_slope(4:6),'--k','LineWidth',2);

ylabel('\Deltarichness/\DeltaT')
plotLabels=["1250" "12" "3"];
xticklabels(plotLabels)
%xlabel('quadrats sampled')

subplot(1,5,5)
hold on
reordered_Chao1_slopeP=mean(allSubSlopePPt_Chao1(plotOrder,:),2);
reordered_Chao2_slopeP=mean(allSubSlopePPt_Chao2(plotOrder,:),2);
reordered_apx_slopeP=mean(allSubSlopePPt_apx(plotOrder,:),2);
reordered_raw_slopeP=mean(allSubSlopePPt_raw(plotOrder,:),2);

plot(reordered_Chao1_slopeP(1:3),'-b','LineWidth',5);
plot(reordered_Chao2_slopeP(1:3),'--','LineWidth',5,'Color',[0 0.7 0.7]);
plot(reordered_apx_slopeP(1:3),'-r','LineWidth',5);
plot(reordered_raw_slopeP(1:3),'--k','LineWidth',5);

plot(reordered_Chao1_slopeP(4:6),'-b','LineWidth',2);
plot(reordered_Chao2_slopeP(4:6),'--','LineWidth',2,'Color',[0 0.7 0.7]);
plot(reordered_apx_slopeP(4:6),'-r','LineWidth',2);
plot(reordered_raw_slopeP(4:6),'--k','LineWidth',2);

ylabel('p_{\Deltarichness}')
xticks([1:3]);
plotLabels=["1250" "12" "3"];
xticklabels(plotLabels)
%xlabel('quadrats sampled')