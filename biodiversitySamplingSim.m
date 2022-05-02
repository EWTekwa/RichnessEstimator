clear
set(0,'DefaultAxesFontSize',22)
scrsz = get(0,'ScreenSize');
%figure('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3)/2.5 scrsz(4)/3]);
figure('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3)/2.2 scrsz(4)/1.5]);

%define dataset
rng(1); %set random number generator seed (88)
simReps=2000; %number of communities to simulation (2000)
numTests=2000; %number of random pairwise community comparisons (2000)
minSpecies=1; %minimum number of species
maxSpecies=100; %maximum number of species
k_min=2; %minimum number of sampled sites within region
k_max=10; %maximum number of sampled sites within region
k_mean=(k_max-k_min)/2+k_min;
num_sites=100; %true number of sites within region
maxRichDiff_all=[2 10]; %maximum richness difference between community pairs to be compared

meanAbundance_var=[2 2]; %mean abundances

%case 1: main text 2 scenarios, a. low portion of inviduals observable, b. high spatial heterogeneitis
% m_var=[20 1]; %1/m is the portion of individuals observable, or equivalently m is the abundance needed for a species to be detected 63% of the time within a site
% c_var=[1 30]; %max clustering
% mvar_var=[10 0]; %variance in m across species
% sdC_var=[0 1];

% %case 2: 2 scenarios, a. lower portion of inviduals observable, b. higher spatial heterogeneitis
m_var=[60 1]; %1/m is the portion of individuals observable, or equivalently m is the abundance needed for a species to be detected 63% of the time within a site
c_var=[1 90]; %max clustering
mvar_var=[10 0]; %variance in m across species
sdC_var=[0 1];

% %case 3: 2 scenarios, a. bad mixed conditions, b. worse mixed conditions
% m_var=[20 60]; %1/m is the portion of individuals observable, or equivalently m is the abundance needed for a species to be detected 63% of the time within a site
% c_var=[30 90]; %max clustering
% mvar_var=[10 10]; %variance in m across species
% sdC_var=[1 1];

sdLogAbundance=0.5; %0.5
min_meanLogClustering=0.9; %0.5

n_m_ratio=meanAbundance_var./m_var
%write observational process model
syms nm C k F P
Dix=(1-exp(-C*nm))*P;
Dis=1-(1-Dix)^k;
d2Dis_dnm2=diff(Dis,'nm',2);
d2Dis_dC2=diff(Dis,'C',2);
d2Dis_dP2=diff(Dis,'P',2);
d2Dis_dnmC=diff(diff(Dis,'nm'),'C');
d2Dis_dnmP=diff(diff(Dis,'nm'),'P');
d2Dis_dCP=diff(diff(Dis,'C'),'P');

%run simulations
for scenario=1:length(meanAbundance_var)
    meanAbundance=meanAbundance_var(scenario);
    max_meanLogClustering=c_var(scenario);
    m_mean=m_var(scenario);
    sdLogClustering=sdC_var(scenario);
    all_k=zeros(1,simReps);
    all_Species=zeros(1,simReps);
    all_DsTrue=zeros(1,simReps);
    all_Ds_unbiased_apx=zeros(1,simReps);
    all_Ds_unbiased_apx_T=zeros(1,simReps);
    all_Chao1=zeros(1,simReps);
    all_Chao2=zeros(1,simReps);
    all_ACE=zeros(1,simReps);
    all_S_aj2=zeros(1,simReps);
    all_S_ij2=zeros(1,simReps);
    all_mean_n=zeros(1,simReps);
    all_var_n=zeros(1,simReps);
    all_mean_n_m=zeros(1,simReps);
    all_var_n_m=zeros(1,simReps);
    all_mean_C=zeros(1,simReps);
    all_var_C=zeros(1,simReps);
    all_mean_P=zeros(1,simReps);
    all_var_P=zeros(1,simReps);
    all_cov_nm_C=zeros(1,simReps);
    all_cov_nm_P=zeros(1,simReps);
    all_cov_C_P=zeros(1,simReps);
    all_mean_n_m_ratio=zeros(1,simReps);
    all_var_n_m_ratio=zeros(1,simReps);
    all_mean_C_ratio=zeros(1,simReps);
    all_var_C_ratio=zeros(1,simReps);
    all_mean_P_ratio=zeros(1,simReps);
    all_var_P_ratio=zeros(1,simReps);
    all_true_n_estRelDiff=zeros(1,simReps);
    all_detected_n_estRelDiff=zeros(1,simReps);
    all_Apx_detectP_terms=zeros(7,simReps);
    all_Apx_T_detectP_terms=zeros(7,simReps);
    
    for i=1:simReps
        m_true=1+rand*(m_mean-1);
        k=randi([k_min k_max]); %a different k per rep (community)
        k(k==0)=1;
        meanLogClustering=min_meanLogClustering+(max_meanLogClustering-min_meanLogClustering)*rand; %a different mean log C per rep (region)
        
        ni=random('lognormal',log(rand*meanAbundance),sdLogAbundance,1,randi([minSpecies,maxSpecies])); %generate random array of species densities
        nc=max(random('lognormal',log(rand*meanLogClustering),sdLogClustering,1,length(ni)),1); %generate random array of species clustering
        
        nk_true=zeros(k,length(ni)); %store true abundance at sites k
        nk=zeros(k,length(ni)); %store observed abundance at sites k
        pi=1-(1./(ni.*(nc-1))).^(1./(ni.*(nc-1)-1)); %calculate occupancy based on abundance and clustering
        for site=1:num_sites %create true species local abundances at sites
            nk_true(site,:)=(random('poisson',ni.*nc-1,1,length(ni))+1).*(rand(1,length(ni))<pi);
        end
        
        m_species=max(m_true-mvar_var(scenario)/2+rand(1,length(ni))*mvar_var(scenario),1); %generate independent random detection thresholds across species
        for site=1:k %sample from true sites
            rand_site=randi(num_sites); %pick random site
            nk(site,:)=random('poisson',nk_true(rand_site,:)./m_species); %generate independent random observed abundance per sample site k
        end
        
        %record true abundance data
        n_s=mean(nk_true(:,sum(nk_true,1)>0));
        nm_true=(nk_true(:,sum(nk_true,1)>0))./m_species(:,sum(nk_true,1)>0); %true observable abundance
        nm_s=mean(nm_true);
        num_species_true=sum(sum(nk_true)>0);
        sites_occupied_true=sum(nk_true(:,sum(nk_true,1)>0)>0);
        mean_n_true=mean(n_s);
        mean_n_m_true=mean(nm_s);
        var_n_true=var(n_s);
        var_n_m_true=var(nm_s);
        C_true=1+var(nk_true(:,sum(nk_true,1)>0),1)./(n_s.^2); %spatial variance normalized by N
        mean_C_true=mean(C_true);
        var_C_true=var(C_true);
        P_true=sites_occupied_true/num_sites;
        mean_P_true=mean(P_true);
        mean_1P_true=mean(1./P_true);
        var_P_true=var(P_true);
        var_1P_true=var(1./P_true);
        cov_nm_C_true=nancov(nm_s,C_true);
        cov_nm_P_true=nancov(nm_s,P_true);
        cov_C_P_true=nancov(C_true,P_true);
        if size(cov_nm_P_true,2)>1
            cov_nm_C_true=cov_nm_C_true(1,2);
            cov_nm_P_true=cov_nm_P_true(1,2);
            cov_C_P_true=cov_C_P_true(1,2);
        else
            cov_nm_C_true=0;
            cov_nm_P_true=0;
            cov_C_P_true=0;
        end
        
        %record observation data
        nm_s_detected=mean(nk(:,sum(nk,1)>0),1);
        sites_occupied_detected=sum(nk(:,sum(nk,1)>0)>0);
        num_detected=sum(sum(nk,1)>0);
        mean_detected=mean(nm_s_detected);
        var_detected=var(nm_s_detected);
        P_detected=sites_occupied_detected/k;
        mean_P_detected=mean(P_detected);
        mean_1P_detected=mean(1./P_detected);
        var_P_detected=var(P_detected);
        var_1P_detected=var(1./P_detected);
        if k==1
            C_detected=1+1./(mean(nk(:,sum(nk,1)>0))); %if only one sampled site, estimate spatial variance as poisson mean
        else
            C_detected=1+var(nk(:,sum(nk,1)>0),1)./(mean(nk(:,sum(nk,1)>0)).^2); %spatial variance normalized by N
        end
        mean_C_detected=mean(C_detected);
        var_C_detected=var(C_detected);
        cov_nm_C_detected=nancov(nm_s_detected,C_detected);
        cov_nm_P_detected=nancov(nm_s_detected,P_detected);
        cov_C_P_detected=nancov(C_detected,P_detected);
        if size(cov_nm_P_detected,2)>1
            cov_nm_C_detected=cov_nm_C_detected(1,2);
            cov_nm_P_detected=cov_nm_P_detected(1,2);
            cov_C_P_detected=cov_C_P_detected(1,2);
        else
            cov_nm_C_detected=0;
            cov_nm_P_detected=0;
            cov_C_P_detected=0;
        end

        
        [Richness_raw,Chao1,Chao2,ACE,S_aj2,S_ij2,Richness_apx,Apx_detectP_terms] = RichnessEstsCov(nk); %compute estimates
        [Richness_apx_T,Apx_T_detectP_terms] = RichnessEstsCov_ideal(nm_true,k,Richness_raw); %compute idealized approximation
        %compute proposed approximation using true C, P, and nm for
        %theoretical best performance projection:
%         nm=mean_n_m_true;
%         C=mean_C_true;
%         P=mean_P_true;
%  
%         find(all_Ds_unbiased_apx_T<-3000)
%         sp=1777
%         num_detected=all_DsTrue(sp)
%         num_species_true=all_Species(i)
%         k=all_k(sp)
%         nm=all_mean_n_m(sp)
%         C=all_mean_C(sp)
%         P=all_mean_P(sp)
%         var_n_m_true=all_var_n_m(sp)
%         var_C_true=all_var_C(sp)
%         var_P_true=all_var_P(sp)
%         cov_nm_C_true=all_cov_nm_C(sp)
%         cov_nm_P_true=all_cov_nm_P(sp)
%         cov_C_P_true=all_cov_C_P(sp)
%         
%         Ds_apx_T=eval(Dis)+eval(d2Dis_dnm2)*var_n_m_true/2+eval(d2Dis_dC2)*var_C_true/2+eval(d2Dis_dP2)*var_P_true/2+eval(d2Dis_dnmC)*cov_nm_C_true+eval(d2Dis_dnmP)*cov_nm_P_true+eval(d2Dis_dCP)*cov_C_P_true; %Approximated detection probability in community using true C, P, and nm values
%        % Ds_apx_T=eval(Dis)+eval(d2Dis_dnm2)*var_n_m_true/2+eval(d2Dis_dC2)*var_C_true/2+eval(d2Dis_dP2)*var_P_true/2; %Approximated detection probability in community using true C, P, and nm values
%         Richness_apx_T=num_detected/Ds_apx_T;
        
        all_k(i)=k;
        all_mean_n(i)=mean_n_true;
        all_var_n(i)=var_n_true;
        all_mean_n_m(i)=mean_n_m_true;
        all_var_n_m(i)=var_n_m_true;
        all_mean_C(i)=mean_C_true;
        all_var_C(i)=var_C_true;
        all_mean_P(i)=mean_P_true;
        all_var_P(i)=var_P_true;
        all_cov_nm_C(i)=cov_nm_C_true;
        all_cov_nm_P(i)=cov_nm_P_true;
        all_cov_C_P(i)=cov_C_P_true;
        
        all_Apx_detectP_terms(:,i)=Apx_detectP_terms;
        all_Apx_T_detectP_terms(:,i)=Apx_T_detectP_terms;
        
        all_mean_n_m_ratio(i)=mean_detected/mean_n_m_true;
        all_var_n_m_ratio(i)=var_detected/var_n_m_true;
        all_mean_C_ratio(i)=mean_C_detected/mean_C_true;
        all_var_C_ratio(i)=var_C_detected/var_C_true;
        all_mean_P_ratio(i)=mean_P_detected/mean_P_true;
        all_var_P_ratio(i)=var_P_detected/var_P_true;
        
        all_Species(i)=num_species_true;
        all_DsTrue(i)=num_detected;
        all_Ds_unbiased_apx(i)=Richness_apx;
        all_Ds_unbiased_apx_T(i)=Richness_apx_T;
        all_Chao1(i)=Chao1;
        all_Chao2(i)=Chao2;
        all_ACE(i)=ACE;
        all_S_aj2(i)=S_aj2;
        all_S_ij2(i)=S_ij2;
    end
    
    %pairwise richness comparison (between simulation repicates)
    numCorrect_detected=[0 0];
    numCorrect_Chao1=[0 0];
    numCorrect_Chao2=[0 0];
    numCorrect_ACE=[0 0];
    numCorrect_S_aj2=[0 0];
    numCorrect_S_ij2=[0 0];
    numCorrect_unbiased_apx=[0 0];
    numCorrect_unbiased_apx_T=[0 0];
    for test=1:numTests
        [~,temppos]=sort(all_Species);
        NonZeroPos=temppos(all_DsTrue(temppos)>0); %find indices of communities with detected species, ordered from low to high abundance
        CommunityOrder=randi(sum(all_DsTrue>0),1,2); %draw only from communities with detected species
        CommunityPair=NonZeroPos(CommunityOrder); %match picked community orders to indices
        rich1=all_Species(CommunityPair(1)); %richness of first random community
        for testvar=1:2
            pos2=find(all_Species(NonZeroPos)<rich1+maxRichDiff_all(testvar) & all_Species(NonZeroPos)>rich1-maxRichDiff_all(testvar) & all_Species(NonZeroPos)~=rich1); %pick a second community with richness within the defined range when compared to the first community, but does not have identical richness
            if testvar==1 && isempty(pos2)
                test=test-1;
                testvar=2;
            else
                CommunityPair(2)=NonZeroPos(pos2(randperm(length(pos2),1)));
                CommDiff_true=diff(all_Species(CommunityPair)); %difference (richness2-richness1)
                CommDiff_detected=diff(all_DsTrue(CommunityPair));
                CommDiff_Chao1=diff(all_Chao1(CommunityPair));
                CommDiff_Chao2=diff(all_Chao2(CommunityPair));
                CommDiff_ACE=diff(all_ACE(CommunityPair));
                CommDiff_S_aj2=diff(all_S_aj2(CommunityPair));
                CommDiff_S_ij2=diff(all_S_ij2(CommunityPair));
                CommDiff_unbiased_apx=diff(all_Ds_unbiased_apx(CommunityPair));
                CommDiff_unbiased_apx_T=diff(all_Ds_unbiased_apx_T(CommunityPair));
                numCorrect_detected(testvar)=numCorrect_detected(testvar)+((0^CommDiff_true)==(0^CommDiff_detected)); %add 1 if correct
                numCorrect_Chao1(testvar)=numCorrect_Chao1(testvar)+((0^CommDiff_true)==(0^CommDiff_Chao1)); %add 1 if correct
                numCorrect_Chao2(testvar)=numCorrect_Chao2(testvar)+((0^CommDiff_true)==(0^CommDiff_Chao2)); %add 1 if correct
                numCorrect_ACE(testvar)=numCorrect_ACE(testvar)+((0^CommDiff_true)==(0^CommDiff_ACE)); %add 1 if correct
                numCorrect_S_aj2(testvar)=numCorrect_S_aj2(testvar)+((0^CommDiff_true)==(0^CommDiff_S_aj2)); %add 1 if correct
                numCorrect_S_ij2(testvar)=numCorrect_S_ij2(testvar)+((0^CommDiff_true)==(0^CommDiff_S_ij2)); %add 1 if correct
                numCorrect_unbiased_apx(testvar)=numCorrect_unbiased_apx(testvar)+((0^CommDiff_true)==(0^CommDiff_unbiased_apx)); %add 1 if correct
                numCorrect_unbiased_apx_T(testvar)=numCorrect_unbiased_apx_T(testvar)+((0^CommDiff_true)==(0^CommDiff_unbiased_apx_T)); %add 1 if correct
            end
        end
    end
    
    scenario
    mean_n=nanmean(all_mean_n)
    mean_m=nanmean(all_mean_n_m./all_mean_n)
    mean_C=nanmean(all_mean_C)
    mean_P=nanmean(all_mean_P)
    
    %subplot(length(meanAbundance_var)/2,2,scenario)
    subplot(2,2,scenario)
    title(['E[m]=' num2str(mean_m,1) ', E[n]=' num2str(mean_n,1) ', E[C]=' num2str(mean_C,1) ', E[P]=' num2str(mean_P,1) ', E[K]=' num2str(k_mean,1)])
    hold on
    scatter(all_Species,all_DsTrue,18,'k','filled')
    scatter(all_Species,all_Chao1,18,'b')
    scatter(all_Species,all_Ds_unbiased_apx,18,'r')
    hline1=refline(1,0);
    hline1.Color = 'k';
    xlabel 'true species richness'
    ylabel 'estimated species richness'
    
    all_DsTrue(all_DsTrue==0)=NaN;
    mdl_detected=fitlm(all_Species,all_DsTrue);
    mdl_Chao1=fitlm(all_Species,all_Chao1);
    mdl_Chao2=fitlm(all_Species,all_Chao2);
    mdl_ACE=fitlm(all_Species,all_ACE);
    mdl_S_aj2=fitlm(all_Species,all_S_aj2);
    mdl_S_ij2=fitlm(all_Species,all_S_ij2);
    mdl_correctionApx=fitlm(all_Species,all_Ds_unbiased_apx);
    mdl_correctionApx_T=fitlm(all_Species,all_Ds_unbiased_apx_T);
    
    RSS_detected=nansum((all_Species-all_DsTrue).^2);
    RSS_Chao1=nansum((all_Species-all_Chao1).^2);
    RSS_Chao2=nansum((all_Species-all_Chao2).^2);
    RSS_ACE=nansum((all_Species-all_ACE).^2);
    RSS_S_aj2=nansum((all_Species-all_S_aj2).^2);
    RSS_S_ij2=nansum((all_Species-all_S_ij2).^2);
    RSS_correctionApx=nansum((all_Species-all_Ds_unbiased_apx).^2);
    RSS_correctionApx_T=nansum((all_Species-all_Ds_unbiased_apx_T).^2);
    TSS=nansum((all_Species-nanmean(all_Species)).^2);
    
    slope_detected=mdl_detected.Coefficients.Estimate(2);
    slope_Chao1=mdl_Chao1.Coefficients.Estimate(2);
    slope_Chao2=mdl_Chao2.Coefficients.Estimate(2);
    slope_ACE=mdl_ACE.Coefficients.Estimate(2);
    slope_S_aj2=mdl_S_aj2.Coefficients.Estimate(2);
    slope_S_ij2=mdl_S_ij2.Coefficients.Estimate(2);
    slope_correctionApx=mdl_correctionApx.Coefficients.Estimate(2);
    slope_correctionApx_T=mdl_correctionApx_T.Coefficients.Estimate(2);
    slopeSD_detected=mdl_detected.Coefficients.SE(2)*sqrt(simReps);
    slopeSD_Chao1=mdl_Chao1.Coefficients.SE(2)*sqrt(simReps);
    slopeSD_Chao2=mdl_Chao2.Coefficients.SE(2)*sqrt(simReps);
    slopeSD_ACE=mdl_ACE.Coefficients.SE(2)*sqrt(simReps);
    slopeSD_S_aj2=mdl_S_aj2.Coefficients.SE(2)*sqrt(simReps);
    slopeSD_S_ij2=mdl_S_ij2.Coefficients.SE(2)*sqrt(simReps);
    slopeSD_correctionApx=mdl_correctionApx.Coefficients.SE(2)*sqrt(simReps);
    slopeSD_correctionApx_T=mdl_correctionApx_T.Coefficients.SE(2)*sqrt(simReps);
    R2_detected=1-RSS_detected/TSS;
    R2_Chao1=1-RSS_Chao1/TSS;
    R2_Chao2=1-RSS_Chao2/TSS;
    R2_ACE=1-RSS_ACE/TSS;
    R2_S_aj2=1-RSS_S_aj2/TSS;
    R2_S_ij2=1-RSS_S_ij2/TSS;
    R2_correctionApx=1-RSS_correctionApx/TSS;
    R2_correctionApx_T=1-RSS_correctionApx_T/TSS;
    
    ylim([0 180])
    ylimits=ylim;
    text(minSpecies+5*(maxSpecies-minSpecies)/100,ylimits(2)*0.95,['raw slope=' num2str(slope_detected,2) '\pm' num2str(slopeSD_detected,2) ', R^2*=' num2str(R2_detected,2) ', pairwise_{' num2str(maxRichDiff_all(1)) ',' num2str(maxRichDiff_all(2)) '}=' num2str(numCorrect_detected(1)/numTests,2) ',' num2str(numCorrect_detected(2)/numTests,2)])
    text(minSpecies+5*(maxSpecies-minSpecies)/100,ylimits(2)*0.9,['Chao1 slope=' num2str(slope_Chao1,2) '\pm' num2str(slopeSD_Chao1,2) ', R^2*=' num2str(R2_Chao1,2) ', pairwise_{' num2str(maxRichDiff_all(1)) ',' num2str(maxRichDiff_all(2)) '}=' num2str(numCorrect_Chao1(1)/numTests,2) ',' num2str(numCorrect_Chao1(2)/numTests,2)],'Color','b')
    text(minSpecies+5*(maxSpecies-minSpecies)/100,ylimits(2)*0.85,['Chao2 slope=' num2str(slope_Chao2,2) '\pm' num2str(slopeSD_Chao2,2) ', R^2*=' num2str(R2_Chao2,2) ', pairwise_{' num2str(maxRichDiff_all(1)) ',' num2str(maxRichDiff_all(2)) '}=' num2str(numCorrect_Chao2(1)/numTests,2) ',' num2str(numCorrect_Chao2(2)/numTests,2)],'Color',[0.5 0.5 0.5])
    text(minSpecies+5*(maxSpecies-minSpecies)/100,ylimits(2)*0.8,['ACE slope=' num2str(slope_ACE,2) '\pm' num2str(slopeSD_ACE,2) ', R^2*=' num2str(R2_ACE,2) ', pairwise_{' num2str(maxRichDiff_all(1)) ',' num2str(maxRichDiff_all(2)) '}=' num2str(numCorrect_ACE(1)/numTests,2) ',' num2str(numCorrect_ACE(2)/numTests,2)],'Color',[0.5 0.5 0.5])
    text(minSpecies+5*(maxSpecies-minSpecies)/100,ylimits(2)*0.75,['JK_a slope=' num2str(slope_S_aj2,2) '\pm' num2str(slopeSD_S_aj2,2) ', R^2*=' num2str(R2_S_aj2,2) ', pairwise_{' num2str(maxRichDiff_all(1)) ',' num2str(maxRichDiff_all(2)) '}=' num2str(numCorrect_S_aj2(1)/numTests,2) ',' num2str(numCorrect_S_aj2(2)/numTests,2)],'Color',[0.5 0.5 0.5])
    text(minSpecies+5*(maxSpecies-minSpecies)/100,ylimits(2)*0.7,['JK_i slope=' num2str(slope_S_ij2,2) '\pm' num2str(slopeSD_S_ij2,2) ', R^2*=' num2str(R2_S_ij2,2) ', pairwise_{' num2str(maxRichDiff_all(1)) ',' num2str(maxRichDiff_all(2)) '}=' num2str(numCorrect_S_ij2(1)/numTests,2) ',' num2str(numCorrect_S_ij2(2)/numTests,2)],'Color',[0.5 0.5 0.5])
    text(minSpecies+5*(maxSpecies-minSpecies)/100,ylimits(2)*0.65,['Apx slope=' num2str(slope_correctionApx,2) '\pm' num2str(slopeSD_correctionApx,2) ', R^2*=' num2str(R2_correctionApx,2) ', pairwise_{' num2str(maxRichDiff_all(1)) ',' num2str(maxRichDiff_all(2)) '}=' num2str(numCorrect_unbiased_apx(1)/numTests,2) ',' num2str(numCorrect_unbiased_apx(2)/numTests,2)],'Color','r')
    text(minSpecies+5*(maxSpecies-minSpecies)/100,ylimits(2)*0.60,['Apx_T slope=' num2str(slope_correctionApx_T,2) '\pm' num2str(slopeSD_correctionApx_T,2) ', R^2*=' num2str(R2_correctionApx_T,2) ', pairwise_{' num2str(maxRichDiff_all(1)) ',' num2str(maxRichDiff_all(2)) '}=' num2str(numCorrect_unbiased_apx_T(1)/numTests,2) ',' num2str(numCorrect_unbiased_apx_T(2)/numTests,2)],'Color',[0.5 0.5 0.5])
    
    subplot(2,2,scenario+2)
    hold on
    xVal=repmat(1:size(all_Apx_detectP_terms,1)+1,1,simReps); %mean x-axis positions for correction terms
    xPos1=xVal+0.4*(rand(size(xVal))-1); %add jitter to x-axis positions (left)
    xPos2=xVal+0.4*(rand(size(xVal))); %add jitter to x-axis positions (right)
    scatter(xPos1,reshape([sum(all_Apx_detectP_terms); all_Apx_detectP_terms],1,[]),5,'o','MarkerEdgeColor','none','MarkerFaceColor','r','MarkerFaceAlpha',0.3);
    scatter(xPos2,reshape([sum(all_Apx_T_detectP_terms); all_Apx_T_detectP_terms],1,[]),5,'o','MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha',0.3);
    hline2=refline(0,0);
    hline2.Color = 'k';
    xline(1.5)
    xlabel 'correction factor'
    ylabel 'estimated observation probability'
    xticks(1:size(all_Apx_detectP_terms,1)+1)
    xticklabels({'sum','0^{th}','V_{nm}','V_C','V_P','C_{nm,C}','C_{nm,P}','C_{C,P}'})
    xlim([0.5 size(all_Apx_detectP_terms,1)+1.5])
end

% %%%%%%%%
% %diagnostic plots for last scenario or last community simulated:
% 
% diagfig=figure('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3)/3.5 scrsz(4)/6]);
% set(diagfig,'defaultAxesColorOrder',[[0 0 0]; [0 0 0]]);
% subplot(1,2,1)
% hold on
% [n_true_sorted,n_index]=sort(mean(nk_true(:,sum(nk_true,1)>0))./m_species(:,sum(nk_true,1)>0),'descend');
% scatter([1:num_species_true],n_true_sorted,'m','filled') %plot rank-true m/n curve
% scatter([1:num_detected],sort(mean(nk(:,sum(nk,1)>0),1),'descend'),'k','filled') %plot rank-observed m/n curve
% ax=gca;
% text(ax.XLim(2)*.1,ax.YLim(2)*.95,['mean true=' num2str(mean_n_m_true,2) ', observed=' num2str(mean_detected,2)])
% text(ax.XLim(2)*.1,ax.YLim(2)*0.9,['var true=' num2str(var_n_m_true,2) ', observed=' num2str(var_detected,2)])
% xlabel 'species abundance rank'
% ylabel 'observable abundance per site (nm)'
% legend('true','observed','location','east')
% 
% subplot(1,2,2)
% yyaxis left
% hold on
% [C_sorted,C_index]=sort(C_true,'descend');
% [C_detected_sorted,C_sorted_index]=sort(C_detected,'descend');
% NumSitesOccupied=sum(nk(:,sum(nk,1)>0)>0);
% scatter([1:num_species_true],C_sorted,'m','filled')
% scatter([num_species_true-num_detected+1:num_species_true],[C_detected_sorted],'k','filled')
% ax=gca;
% text(ax.XLim(2)*.1,ax.YLim(2)*.95,['mean(C) true=' num2str(mean_C_true,2) ', observed=' num2str(mean_C_detected,2)])
% text(ax.XLim(2)*.1,ax.YLim(2)*.9,['var(C) true=' num2str(var_C_true,2) ', observed=' num2str(var_C_detected,2)])
% text(ax.XLim(2)*.1,ax.YLim(2)*.85,['mean(P) true=' num2str(mean_P_true,2) ', observed=' num2str(mean_P_detected,2)])
% text(ax.XLim(2)*.1,ax.YLim(2)*.8,['var(P) true=' num2str(var_P_true,2) ', observed=' num2str(var_P_detected,2)])
% xlabel 'species clustering rank'
% ylabel 'clustering (C)'
% yyaxis right
% hold on
% scatter([1:num_species_true],sites_occupied_true(C_index)./num_sites,60,'xm')
% scatter([num_species_true-num_detected+1:num_species_true],sites_occupied_detected(C_sorted_index)./k,60,'xk')
% ylabel 'occupancy (P)'
% legend('C_T','C_{obs}','P_T','P_{obs}','location','best')
% 
% 
% figure('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3)/3.5 scrsz(4)/6]);
% subplot(1,2,1)
% hold on
% histogram(all_mean_n,80,'facecolor','w');
% xlim([0,8])
% xlabel 'mean species abundance per site per community'
% ylabel 'count'
% 
% subplot(1,2,2)
% hold on
% histogram(mean(nk(:,sum(nk,1)>0),1),'facecolor','w');
% xlabel 'observed species abundance per site within one community'
% ylabel 'count'