%biodiversitySamplingSim.m
%Eden Tekwa

clear
set(0,'DefaultAxesFontSize',22)
scrsz = get(0,'ScreenSize');
%figure('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3)/2.5 scrsz(4)/3]);
figure('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3)/2.2 scrsz(4)/1.5]);

%define dataset
rng(50); %set random number generator seed (10)
simReps=2000; %number of communities to simulation (2000)
numTests=2000; %number of random pairwise community comparisons (2000)
minSpecies=1; %minimum number of species
maxSpecies=100; %maximum number of species
k_min=2; %minimum number of sampled sites within region
k_max=10; %maximum number of sampled sites within region
k_mean=(k_max-k_min)/2+k_min;
num_sites=100; %true number of sites within region
maxRichDiff_all=[2 10 20]; %maximum richness difference between community pairs to be compared
numBoot=50; %number of bootstraps for each community estimate

%case 1: main text 2 scenarios, a. low portion of inviduals observable, b. high spatial heterogeneitis
m_scenarios=[.4 1]; %m is the portion of individuals observable, or equivalently 1/m is the abundance needed for a species to be detected 63% of the time within a site
P_scenarios=[1 0.4]; %mean occupancy (variance=mean)
meanAbundance_scenarios=[2 2]; %mean abundances

% %case 2: 2 scenarios, a. mixed conditions, b. worse mixed conditions
% m_scenarios=[0.4 0.2]; %m is the portion of individuals observable, or equivalently 1/m is the abundance needed for a species to be detected 63% of the time within a site
% P_scenarios=[0.4 0.2]; %mean occupancy (variance=mean)
% meanAbundance_scenarios=[2 2]; %mean abundances

sdLogAbundance=1; %0.5, 1
%min_meanLogClustering=0.9; %0.5
%min_Clustering=0.9;
%min_Abundance=0.1;

%n_m_ratio=meanAbundance_var./m_scenarios
%write observational process model
syms nm k P
%Dix=(1-exp(-C*nm))*P;
Dix=(1-exp(-nm/P))*P;
Dis=1-(1-Dix)^k;
% d2Dis_dnm2=diff(Dis,'nm',2);
% d2Dis_dC2=diff(Dis,'C',2);
% d2Dis_dP2=diff(Dis,'P',2);
% d2Dis_dnmC=diff(diff(Dis,'nm'),'C');
% d2Dis_dnmP=diff(diff(Dis,'nm'),'P');
% d2Dis_dCP=diff(diff(Dis,'C'),'P');
d2Dis_dnm2 = - (k*exp(-nm/P)*(P*(exp(-nm/P) - 1) + 1)^(k - 1))/P - k*exp(-(2*nm)/P)*(P*(exp(-nm/P) - 1) + 1)^(k - 2)*(k - 1);
d2Dis_dP2 = - k*(P*(exp(-nm/P) - 1) + 1)^(k - 2)*(k - 1)*(exp(-nm/P) + (nm*exp(-nm/P))/P - 1)^2 - (k*nm^2*exp(-nm/P)*(P*(exp(-nm/P) - 1) + 1)^(k - 1))/P^3;
d2Dis_dnmP = (k*nm*exp(-nm/P)*(P*(exp(-nm/P) - 1) + 1)^(k - 1))/P^2 + k*exp(-nm/P)*(P*(exp(-nm/P) - 1) + 1)^(k - 2)*(k - 1)*(exp(-nm/P) + (nm*exp(-nm/P))/P - 1);

%run simulations
ScenarioEvals=zeros(13,6,length(meanAbundance_scenarios)); %store 6 evaluation outcomes (columns) for all estimators (rows) and scenarios (3rd dimension)
ScenarioEvalsRel=zeros(13,6,length(meanAbundance_scenarios));
ScenarioEvalRankIndex=zeros(13,6,length(meanAbundance_scenarios));
ScenarioEvalRanks=zeros(13,7,length(meanAbundance_scenarios));
ScenarioEvalsRelDiff=zeros(13,7,length(meanAbundance_scenarios));
splot=1; %track subplot label
for scenario=1:length(meanAbundance_scenarios)
    meanAbundance=meanAbundance_scenarios(scenario);
    %    max_meanLogClustering=c_var(scenario);
    m_mean=m_scenarios(scenario);
    P_mean=P_scenarios(scenario);
    %sdLogClustering=sdC_var(scenario);
    all_k=zeros(1,simReps);
    all_Species=zeros(1,simReps);
    all_DsTrue=zeros(1,simReps);
    all_Ds_unbiased_apx=zeros(1,simReps);
    all_Ds_unbiased_T=zeros(1,simReps);
    all_Ds_unbiased_0=zeros(1,simReps);
    all_Ds_unbiased_apxC=zeros(1,simReps);
    all_Ds_unbiased_TC=zeros(1,simReps);
    all_Ds_unbiased_0C=zeros(1,simReps);
    all_Chao1=zeros(1,simReps);
    all_GP=zeros(1,simReps);
    all_Chao2=zeros(1,simReps);
    all_ACE=zeros(1,simReps);
    all_S_aj2=zeros(1,simReps);
    all_S_ij2=zeros(1,simReps);
    all_mean_n=zeros(1,simReps);
    all_var_n=zeros(1,simReps);
    all_mean_n_m=zeros(1,simReps);
    all_var_n_m=zeros(1,simReps);
    %all_mean_C=zeros(1,simReps);
    %all_var_C=zeros(1,simReps);
    all_mean_P=zeros(1,simReps);
    all_var_P=zeros(1,simReps);
    %all_cov_nm_C=zeros(1,simReps);
    all_cov_nm_P=zeros(1,simReps);
    %all_cov_C_P=zeros(1,simReps);
    all_mean_n_m_ratio=zeros(1,simReps);
    all_var_n_m_ratio=zeros(1,simReps);
    %all_mean_C_ratio=zeros(1,simReps);
    %all_var_C_ratio=zeros(1,simReps);
    all_mean_P_ratio=zeros(1,simReps);
    all_var_P_ratio=zeros(1,simReps);
    all_true_n_estRelDiff=zeros(1,simReps);
    all_detected_n_estRelDiff=zeros(1,simReps);
    all_Apx_T_detectP_terms=zeros(4,simReps);
    all_Apx_TC_detectP_terms=zeros(4,simReps);
    all_bootStd_raw=zeros(1,simReps);
    all_bootStd_Chao1=zeros(1,simReps);
    all_bootStd_Apx=zeros(1,simReps);
    
    for i=1:simReps %(community i)
        %set
        %m_community=1+rand*2*(m_mean-1); %mean m in community
        m_community=m_mean+(rand-0.5)*(m_mean*(1-m_mean));
        k=randi([k_min k_max]); %mean k in community
        k(k==0)=1;
        P_community=P_mean+(rand-0.5)*(P_mean*(1-P_mean));
        Abundance_community=rand*2*meanAbundance;
        %meanLogClustering=min_meanLogClustering+(max_meanLogClustering-min_meanLogClustering)*rand; %a different mean log C per rep (region)
        
        ni=random('lognormal',log(Abundance_community)-(sdLogAbundance^2)/2,sdLogAbundance,1,randi([minSpecies,maxSpecies]))/P_community; %generate random array of species densities
        %ni<1 will have inflated average abundance because at least one
        %individual has to be added to occupied patches (all species must
        %have average abundance of at least P). Thus the average abundance
        %is adjusted downward for all species:
        ni_totAdj=sum(ni(ni<1))-sum(ni<1);
        ni=ni*(sum(ni)+ni_totAdj)/sum(ni);
        
        %nc=max(random('lognormal',log(rand*meanLogClustering),sdLogClustering,1,length(ni)),1); %generate random array of species clustering
        
        nk_true=zeros(k,length(ni)); %store true abundance at sites k
        nk=zeros(k,length(ni)); %store observed abundance at sites k
        %pi=1-(1./(ni.*(nc-1))).^(1./(ni.*(nc-1)-1)); %calculate occupancy based on abundance and clustering
        %pi=min(max(random('normal',P_community,P_var,1,length(ni)),0),1); %limit occupancy to [0,1]
        pi=P_community+(rand(1,length(ni))-0.5)*(P_community*(1-P_community));
        for site=1:num_sites %create true species local abundances at sites
            %nk_true(site,:)=(random('poisson',ni.*nc-1,1,length(ni))+1).*(rand(1,length(ni))<pi);
            %nk_true(site,:)=(random('poisson',ni,1,length(ni))+1).*(rand(1,length(ni))<pi);
            nk_true(site,:)=(max(random('poisson',ni,1,length(ni)),1)).*(rand(1,length(ni))<pi);
        end
        m_species=m_community+(rand(1,length(ni))-0.5)*(m_community*(1-m_community));
        %m_species=1+rand(1,length(ni))*2*(m_community-1);
        %m_species=max(m_community-m_var/2+rand(1,length(ni))*m_var,1); %generate independent random detection thresholds across species
        rand_sites=randperm(num_sites,k); %pick k random site without replacement
        for site=1:k %sample from true sites
            %rand_site=randi(num_sites); %pick random site
            rand_site=rand_sites(site);
            nk(site,:)=random('poisson',nk_true(rand_site,:).*m_species); %generate independent random observed abundance per sample site k
        end
        
        %record true abundance data
        n_s=mean(nk_true(:,sum(nk_true,1)>0));
        nm_true=(nk_true(:,sum(nk_true,1)>0)).*m_species(:,sum(nk_true,1)>0); %true observable abundance
        nm_s=mean(nm_true);
        num_species_true=sum(sum(nk_true)>0);
        sites_occupied_true=sum(nk_true(:,sum(nk_true,1)>0)>0);
        mean_n_true=mean(n_s);
        mean_n_m_true=mean(nm_s);
        var_n_true=var(n_s);
        var_n_m_true=var(nm_s);
        %C_true=1+var(nk_true(:,sum(nk_true,1)>0),1)./(n_s.^2); %spatial variance normalized by N
        %mean_C_true=mean(C_true);
        %var_C_true=var(C_true);
        P_true=sites_occupied_true/num_sites;
        mean_P_true=mean(P_true);
        mean_1P_true=mean(1./P_true);
        var_P_true=var(P_true);
        var_1P_true=var(1./P_true);
        %cov_nm_C_true=nancov(nm_s,C_true);
        cov_nm_P_true=nancov(nm_s,P_true);
        %cov_C_P_true=nancov(C_true,P_true);
        if size(cov_nm_P_true,2)>1
            %cov_nm_C_true=cov_nm_C_true(1,2);
            cov_nm_P_true=cov_nm_P_true(1,2);
            %cov_C_P_true=cov_C_P_true(1,2);
        else
            %cov_nm_C_true=0;
            cov_nm_P_true=0;
            %cov_C_P_true=0;
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
        %         if k==1
        %             C_detected=1+1./(mean(nk(:,sum(nk,1)>0))); %if only one sampled site, estimate spatial variance as poisson mean
        %         else
        %             C_detected=1+var(nk(:,sum(nk,1)>0),1)./(mean(nk(:,sum(nk,1)>0)).^2); %spatial variance normalized by N
        %         end
        %         mean_C_detected=mean(C_detected);
        %         var_C_detected=var(C_detected);
        %         cov_nm_C_detected=nancov(nm_s_detected,C_detected);
        cov_nm_P_detected=nancov(nm_s_detected,P_detected);
        %cov_C_P_detected=nancov(C_detected,P_detected);
        if size(cov_nm_P_detected,2)>1
            %cov_nm_C_detected=cov_nm_C_detected(1,2);
            cov_nm_P_detected=cov_nm_P_detected(1,2);
            %cov_C_P_detected=cov_C_P_detected(1,2);
        else
            %cov_nm_C_detected=0;
            cov_nm_P_detected=0;
            %cov_C_P_detected=0;
        end
        
        
        %[Richness_raw,Chao1,Chao2,ACE,S_aj2,S_ij2,Richness_apx,Apx_detectP_terms] = RichnessEstsIncidence(nk); %compute estimates
        [Richness_raw,Chao1,GP,Chao2,ACE,S_aj2,S_ij2,Richness_apx,Richness_T,Richness_0,Apx_T_detectP_terms,~] = RichnessEsts(nk); %compute estimates
        if Richness_raw==0 || isnan(Richness_raw) %if no species observed, rerun simulation replicate
            i=i-1;
            continue
        end
        [~,~,~,~,~,~,~,~,~,~,expectedRichness_raw,expectedChao1,~,~,~,~,~,expectedRichness_apx,~,~,~] = bootRichnessEsts(nk,numBoot);
        [Richness_apxC,Richness_TC,Richness_0C,Apx_TC_detectP_terms] = RichnessEsts_ideal(nm_true,k,Richness_raw); %compute idealized approximation
        
        all_bootStd_raw(i)=nanstd(expectedRichness_raw);
        all_bootStd_Chao1(i)=nanstd(expectedChao1);
        all_bootStd_Apx(i)=nanstd(expectedRichness_apx);
        
        all_k(i)=k;
        all_mean_n(i)=mean_n_true;
        all_var_n(i)=var_n_true;
        all_mean_n_m(i)=mean_n_m_true;
        all_var_n_m(i)=var_n_m_true;
        %all_mean_C(i)=mean_C_true;
        %all_var_C(i)=var_C_true;
        all_mean_P(i)=mean_P_true;
        all_var_P(i)=var_P_true;
        %all_cov_nm_C(i)=cov_nm_C_true;
        all_cov_nm_P(i)=cov_nm_P_true;
        %all_cov_C_P(i)=cov_C_P_true;
        
        all_Apx_T_detectP_terms(:,i)=Apx_T_detectP_terms;
        all_Apx_TC_detectP_terms(:,i)=Apx_TC_detectP_terms;
        
        all_mean_n_m_ratio(i)=mean_detected/mean_n_m_true;
        all_var_n_m_ratio(i)=var_detected/var_n_m_true;
        %all_mean_C_ratio(i)=mean_C_detected/mean_C_true;
        %all_var_C_ratio(i)=var_C_detected/var_C_true;
        all_mean_P_ratio(i)=mean_P_detected/mean_P_true;
        all_var_P_ratio(i)=var_P_detected/var_P_true;
        
        all_Species(i)=num_species_true;
        all_DsTrue(i)=num_detected;
        all_Ds_unbiased_apx(i)=Richness_apx;
        all_Ds_unbiased_T(i)=Richness_T;
        all_Ds_unbiased_0(i)=Richness_0;
        all_Ds_unbiased_apxC(i)=Richness_apxC;
        all_Ds_unbiased_TC(i)=Richness_TC;
        all_Ds_unbiased_0C(i)=Richness_0C;
        all_Chao1(i)=Chao1;
        all_GP(i)=GP;
        all_Chao2(i)=Chao2;
        all_ACE(i)=ACE;
        all_S_aj2(i)=S_aj2;
        all_S_ij2(i)=S_ij2;
        
            %%%%%%%%
    %diagnostic plots for current community:
    %
%     diagfig=figure('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3)/3.5 scrsz(4)/4]);
%     %set(diagfig,'defaultAxesColorOrder',[[0 0 0]; [0 0 0]]);
%     subplot(1,2,1)
%     hold on
%     [nm_true_sorted,n_index]=sort(mean(nm_true(:,sum(nk_true,1)>0)),'descend');
%     scatter([1:num_species_true],nm_true_sorted,'m','filled') %plot rank-true m/n curve
%     scatter([1:num_detected],sort(mean(nk(:,sum(nk,1)>0),1),'descend'),'k','filled') %plot rank-observed m/n curve
%     ax=gca;
%     text(ax.XLim(2)*.1,ax.YLim(2)*.9,{['mean true=' num2str(mean_n_m_true,2) ', observed=' num2str(mean_detected,2)],['var true=' num2str(var_n_m_true,2) ', observed=' num2str(var_detected,2)]},'fontsize',16)
%     xlabel 'species abundance rank'
%     ylabel 'observable abundance (mn)'
%     legend('true','observed','location','east')
%     
%     subplot(1,2,2)
%     hold on
%     [P_sorted,P_index]=sort(P_true,'descend');
%     [P_detected_sorted,P_sorted_index]=sort(P_detected,'descend');
%     scatter([1:num_species_true],sites_occupied_true(P_index)./num_sites,'m','filled')
%     scatter([1:num_detected],sites_occupied_detected(P_sorted_index)./k,'k','filled')
%     ax=gca;
%     text(ax.XLim(2)*.1,ax.YLim(2)*.90,{['mean true=' num2str(mean_P_true,2) ', observed=' num2str(mean_P_detected,2)],['var true=' num2str(var_P_true,2) ', observed=' num2str(var_P_detected,2)]},'fontsize',16)
%     xlabel 'species occupancy rank'
%     ylabel 'occupancy (P)'
    
    end
    
    %pairwise richness comparison (between simulation repicates)
    numRichDiff=length(maxRichDiff_all);
    numCorrect_detected=zeros(1,numRichDiff);
    numCorrect_Chao1=zeros(1,numRichDiff);
    numCorrect_GP=zeros(1,numRichDiff);
    numCorrect_Chao2=zeros(1,numRichDiff);
    numCorrect_ACE=zeros(1,numRichDiff);
    numCorrect_S_aj2=zeros(1,numRichDiff);
    numCorrect_S_ij2=zeros(1,numRichDiff);
    numCorrect_unbiased_apx=zeros(1,numRichDiff);
    numCorrect_unbiased_T=zeros(1,numRichDiff);
    numCorrect_unbiased_0=zeros(1,numRichDiff);
    numCorrect_unbiased_apxC=zeros(1,numRichDiff);
    numCorrect_unbiased_TC=zeros(1,numRichDiff);
    numCorrect_unbiased_0C=zeros(1,numRichDiff);
    for test=1:numTests
        [~,temppos]=sort(all_Species);
        NonZeroPos=temppos(all_DsTrue(temppos)>0); %find indices of communities with detected species, ordered from low to high abundance
        CommunityOrder=randi(sum(all_DsTrue>0),1,2); %draw only from communities with detected species
        CommunityPair=NonZeroPos(CommunityOrder); %match picked community orders to indices
        rich1=all_Species(CommunityPair(1)); %richness of first random community
        for testvar=1:numRichDiff
            pos2=find(all_Species(NonZeroPos)<rich1+maxRichDiff_all(testvar) & all_Species(NonZeroPos)>rich1-maxRichDiff_all(testvar) & all_Species(NonZeroPos)~=rich1); %pick a second community with richness within the defined range when compared to the first community, but does not have identical richness
            if testvar==1 && isempty(pos2)
                test=test-1; %redo test set
                testvar=numRichDiff; %trigger exit from current test loop
            else
                CommunityPair(2)=NonZeroPos(pos2(randperm(length(pos2),1)));
                CommDiff_true=diff(all_Species(CommunityPair)); %difference (richness2-richness1)
                CommDiff_detected=diff(all_DsTrue(CommunityPair));
                CommDiff_Chao1=diff(all_Chao1(CommunityPair));
                CommDiff_GP=diff(all_GP(CommunityPair));
                CommDiff_Chao2=diff(all_Chao2(CommunityPair));
                CommDiff_ACE=diff(all_ACE(CommunityPair));
                CommDiff_S_aj2=diff(all_S_aj2(CommunityPair));
                CommDiff_S_ij2=diff(all_S_ij2(CommunityPair));
                CommDiff_unbiased_apx=diff(all_Ds_unbiased_apx(CommunityPair));
                CommDiff_unbiased_T=diff(all_Ds_unbiased_T(CommunityPair));
                CommDiff_unbiased_0=diff(all_Ds_unbiased_0(CommunityPair));
                CommDiff_unbiased_apxC=diff(all_Ds_unbiased_apxC(CommunityPair));
                CommDiff_unbiased_TC=diff(all_Ds_unbiased_TC(CommunityPair));
                CommDiff_unbiased_0C=diff(all_Ds_unbiased_0C(CommunityPair));
                numCorrect_detected(testvar)=numCorrect_detected(testvar)+((0^CommDiff_true)==(0^CommDiff_detected)); %add 1 if correct
                numCorrect_Chao1(testvar)=numCorrect_Chao1(testvar)+((0^CommDiff_true)==(0^CommDiff_Chao1)); %add 1 if correct
                numCorrect_GP(testvar)=numCorrect_GP(testvar)+((0^CommDiff_true)==(0^CommDiff_GP)); %add 1 if correct
                numCorrect_Chao2(testvar)=numCorrect_Chao2(testvar)+((0^CommDiff_true)==(0^CommDiff_Chao2)); %add 1 if correct
                numCorrect_ACE(testvar)=numCorrect_ACE(testvar)+((0^CommDiff_true)==(0^CommDiff_ACE)); %add 1 if correct
                numCorrect_S_aj2(testvar)=numCorrect_S_aj2(testvar)+((0^CommDiff_true)==(0^CommDiff_S_aj2)); %add 1 if correct
                numCorrect_S_ij2(testvar)=numCorrect_S_ij2(testvar)+((0^CommDiff_true)==(0^CommDiff_S_ij2)); %add 1 if correct
                numCorrect_unbiased_apx(testvar)=numCorrect_unbiased_apx(testvar)+((0^CommDiff_true)==(0^CommDiff_unbiased_apx)); %add 1 if correct
                numCorrect_unbiased_T(testvar)=numCorrect_unbiased_T(testvar)+((0^CommDiff_true)==(0^CommDiff_unbiased_T)); %add 1 if correct
                numCorrect_unbiased_0(testvar)=numCorrect_unbiased_0(testvar)+((0^CommDiff_true)==(0^CommDiff_unbiased_0)); %add 1 if correct
                numCorrect_unbiased_apxC(testvar)=numCorrect_unbiased_apxC(testvar)+((0^CommDiff_true)==(0^CommDiff_unbiased_apxC)); %add 1 if correct
                numCorrect_unbiased_TC(testvar)=numCorrect_unbiased_TC(testvar)+((0^CommDiff_true)==(0^CommDiff_unbiased_TC)); %add 1 if correct
                numCorrect_unbiased_0C(testvar)=numCorrect_unbiased_0C(testvar)+((0^CommDiff_true)==(0^CommDiff_unbiased_0C)); %add 1 if correct
                
                
            end
        end
    end
    
    scenario
    mean_n=nanmean(all_mean_n)
    mean_m=nanmean(all_mean_n_m./all_mean_n)
    %mean_C=nanmean(all_mean_C)
    mean_P=nanmean(all_mean_P)
    
    %regress estimated richness to true richness:
    all_DsTrue(all_DsTrue==0)=NaN;
    mdl_detected=fitlm(all_Species,all_DsTrue);
    mdl_Chao1=fitlm(all_Species,all_Chao1);
    mdl_GP=fitlm(all_Species,all_GP);
    mdl_Chao2=fitlm(all_Species,all_Chao2);
    mdl_ACE=fitlm(all_Species,all_ACE);
    mdl_S_aj2=fitlm(all_Species,all_S_aj2);
    mdl_S_ij2=fitlm(all_Species,all_S_ij2);
    mdl_correction_Apx=fitlm(all_Species,all_Ds_unbiased_apx);
    mdl_correction_T=fitlm(all_Species,all_Ds_unbiased_T);
    mdl_correction_0=fitlm(all_Species,all_Ds_unbiased_0);
    mdl_correction_ApxC=fitlm(all_Species,all_Ds_unbiased_apxC);
    mdl_correction_TC=fitlm(all_Species,all_Ds_unbiased_TC);
    mdl_correction_0C=fitlm(all_Species,all_Ds_unbiased_0C);
    
    RSS_detected=nansum((all_Species-all_DsTrue).^2);
    RSS_Chao1=nansum((all_Species-all_Chao1).^2);
    RSS_GP=nansum((all_Species-all_GP).^2);
    RSS_Chao2=nansum((all_Species-all_Chao2).^2);
    RSS_ACE=nansum((all_Species-all_ACE).^2);
    RSS_S_aj2=nansum((all_Species-all_S_aj2).^2);
    RSS_S_ij2=nansum((all_Species-all_S_ij2).^2);
    RSS_correction_Apx=nansum((all_Species-all_Ds_unbiased_apx).^2);
    RSS_correction_T=nansum((all_Species-all_Ds_unbiased_T).^2);
    RSS_correction_0=nansum((all_Species-all_Ds_unbiased_0).^2);
    RSS_correction_ApxC=nansum((all_Species-all_Ds_unbiased_apxC).^2);
    RSS_correction_TC=nansum((all_Species-all_Ds_unbiased_TC).^2);
    RSS_correction_0C=nansum((all_Species-all_Ds_unbiased_0C).^2);
    TSS=nansum((all_Species-nanmean(all_Species)).^2);
    
    [~,coeffSEs_detected,coeffs_detected]=hac(mdl_detected);
    [~,coeffSEs_Chao1,coeffs_Chao1]=hac(mdl_Chao1);
    [~,coeffSEs_GP,coeffs_GP]=hac(mdl_GP);
    [~,coeffSEs_Chao2,coeffs_Chao2]=hac(mdl_Chao2);
    [~,coeffSEs_ACE,coeffs_ACE]=hac(mdl_ACE);
    [~,coeffSEs_S_aj2,coeffs_S_aj2]=hac(mdl_S_aj2);
    [~,coeffSEs_S_ij2,coeffs_S_ij2]=hac(mdl_S_ij2);
    [~,coeffSEs_correction_Apx,coeffs_correction_Apx]=hac(mdl_correction_Apx);
    [~,coeffSEs_correction_T,coeffs_correction_T]=hac(mdl_correction_T);
    [~,coeffSEs_correction_0,coeffs_correction_0]=hac(mdl_correction_0);
    [~,coeffSEs_correction_ApxC,coeffs_correction_ApxC]=hac(mdl_correction_ApxC);
    [~,coeffSEs_correction_TC,coeffs_correction_TC]=hac(mdl_correction_TC);
    [~,coeffSEs_correction_0C,coeffs_correction_0C]=hac(mdl_correction_0C);
    
    slope_detected=coeffs_detected(2);
    slope_Chao1=coeffs_Chao1(2);
    slope_GP=coeffs_GP(2);
    slope_Chao2=coeffs_Chao2(2);
    slope_ACE=coeffs_ACE(2);
    slope_S_aj2=coeffs_S_aj2(2);
    slope_S_ij2=coeffs_S_ij2(2);
    slope_correction_Apx=coeffs_correction_Apx(2);
    slope_correction_T=coeffs_correction_T(2);
    slope_correction_0=coeffs_correction_0(2);
    slope_correction_ApxC=coeffs_correction_ApxC(2);
    slope_correction_TC=coeffs_correction_TC(2);
    slope_correction_0C=coeffs_correction_0C(2);
    slopeSD_detected=coeffSEs_detected(2)*sqrt(simReps);
    slopeSD_Chao1=coeffSEs_Chao1(2)*sqrt(simReps);
    slopeSD_GP=coeffSEs_GP(2)*sqrt(simReps);
    slopeSD_Chao2=coeffSEs_Chao2(2)*sqrt(simReps);
    slopeSD_ACE=coeffSEs_ACE(2)*sqrt(simReps);
    slopeSD_S_aj2=coeffSEs_S_aj2(2)*sqrt(simReps);
    slopeSD_S_ij2=coeffSEs_S_ij2(2)*sqrt(simReps);
    slopeSD_correction_Apx=coeffSEs_correction_Apx(2)*sqrt(simReps);
    slopeSD_correction_T=coeffSEs_correction_T(2)*sqrt(simReps);
    slopeSD_correction_0=coeffSEs_correction_0(2)*sqrt(simReps);
    slopeSD_correction_ApxC=coeffSEs_correction_ApxC(2)*sqrt(simReps);
    slopeSD_correction_TC=coeffSEs_correction_TC(2)*sqrt(simReps);
    slopeSD_correction_0C=coeffSEs_correction_0C(2)*sqrt(simReps);
    
    %     slope_detected=mdl_detected.Coefficients.Estimate(2);
    %     slope_Chao1=mdl_Chao1.Coefficients.Estimate(2);
    %     slope_GP=mdl_GP.Coefficients.Estimate(2);
    %     slope_Chao2=mdl_Chao2.Coefficients.Estimate(2);
    %     slope_ACE=mdl_ACE.Coefficients.Estimate(2);
    %     slope_S_aj2=mdl_S_aj2.Coefficients.Estimate(2);
    %     slope_S_ij2=mdl_S_ij2.Coefficients.Estimate(2);
    %     slope_correction_Apx=mdl_correction_Apx.Coefficients.Estimate(2);
    %     slope_correction_T=mdl_correction_T.Coefficients.Estimate(2);
    %     slope_correction_0=mdl_correction_0.Coefficients.Estimate(2);
    %     slope_correction_ApxC=mdl_correction_ApxC.Coefficients.Estimate(2);
    %     slope_correction_TC=mdl_correction_TC.Coefficients.Estimate(2);
    %     slope_correction_0C=mdl_correction_0C.Coefficients.Estimate(2);
    %     slopeSD_detected=mdl_detected.Coefficients.SE(2)*sqrt(simReps);
    %     slopeSD_Chao1=mdl_Chao1.Coefficients.SE(2)*sqrt(simReps);
    %     slopeSD_GP=mdl_GP.Coefficients.SE(2)*sqrt(simReps);
    %     slopeSD_Chao2=mdl_Chao2.Coefficients.SE(2)*sqrt(simReps);
    %     slopeSD_ACE=mdl_ACE.Coefficients.SE(2)*sqrt(simReps);
    %     slopeSD_S_aj2=mdl_S_aj2.Coefficients.SE(2)*sqrt(simReps);
    %     slopeSD_S_ij2=mdl_S_ij2.Coefficients.SE(2)*sqrt(simReps);
    %     slopeSD_correction_Apx=mdl_correction_Apx.Coefficients.SE(2)*sqrt(simReps);
    %     slopeSD_correction_T=mdl_correction_T.Coefficients.SE(2)*sqrt(simReps);
    %     slopeSD_correction_0=mdl_correction_0.Coefficients.SE(2)*sqrt(simReps);
    %     slopeSD_correction_ApxC=mdl_correction_ApxC.Coefficients.SE(2)*sqrt(simReps);
    %     slopeSD_correction_TC=mdl_correction_TC.Coefficients.SE(2)*sqrt(simReps);
    %     slopeSD_correction_0C=mdl_correction_0C.Coefficients.SE(2)*sqrt(simReps);
    R2_detected=1-RSS_detected/TSS;
    R2_Chao1=1-RSS_Chao1/TSS;
    R2_GP=1-RSS_GP/TSS;
    R2_Chao2=1-RSS_Chao2/TSS;
    R2_ACE=1-RSS_ACE/TSS;
    R2_S_aj2=1-RSS_S_aj2/TSS;
    R2_S_ij2=1-RSS_S_ij2/TSS;
    R2_correction_Apx=1-RSS_correction_Apx/TSS;
    R2_correction_T=1-RSS_correction_T/TSS;
    R2_correction_0=1-RSS_correction_0/TSS;
    R2_correction_ApxC=1-RSS_correction_ApxC/TSS;
    R2_correction_TC=1-RSS_correction_TC/TSS;
    R2_correction_0C=1-RSS_correction_0C/TSS;
    
    %store all evaluation outcomes
    ScenarioEvals(1,:,scenario)=[slope_detected,slopeSD_detected,R2_detected,numCorrect_detected(1)/numTests,numCorrect_detected(2)/numTests,numCorrect_detected(3)/numTests];
    ScenarioEvals(2,:,scenario)=[slope_Chao1,slopeSD_Chao1,R2_Chao1,numCorrect_Chao1(1)/numTests,numCorrect_Chao1(2)/numTests,numCorrect_Chao1(3)/numTests];
    ScenarioEvals(3,:,scenario)=[slope_GP,slopeSD_GP,R2_GP,numCorrect_GP(1)/numTests,numCorrect_GP(2)/numTests,numCorrect_GP(3)/numTests];
    ScenarioEvals(4,:,scenario)=[slope_Chao2,slopeSD_Chao2,R2_Chao2,numCorrect_Chao2(1)/numTests,numCorrect_Chao2(2)/numTests,numCorrect_Chao2(3)/numTests];
    ScenarioEvals(5,:,scenario)=[slope_ACE,slopeSD_ACE,R2_ACE,numCorrect_ACE(1)/numTests,numCorrect_ACE(2)/numTests,numCorrect_ACE(3)/numTests];
    ScenarioEvals(6,:,scenario)=[slope_S_aj2,slopeSD_S_aj2,R2_S_aj2,numCorrect_S_aj2(1)/numTests,numCorrect_S_aj2(2)/numTests,numCorrect_S_aj2(3)/numTests];
    ScenarioEvals(7,:,scenario)=[slope_S_ij2,slopeSD_S_ij2,R2_S_ij2,numCorrect_S_ij2(1)/numTests,numCorrect_S_ij2(2)/numTests,numCorrect_S_ij2(3)/numTests];
    ScenarioEvals(8,:,scenario)=[slope_correction_Apx,slopeSD_correction_Apx,R2_correction_Apx,numCorrect_unbiased_apx(1)/numTests,numCorrect_unbiased_apx(2)/numTests,numCorrect_unbiased_apx(3)/numTests];
    ScenarioEvals(9,:,scenario)=[slope_correction_T,slopeSD_correction_T,R2_correction_T,numCorrect_unbiased_T(1)/numTests,numCorrect_unbiased_T(2)/numTests,numCorrect_unbiased_T(3)/numTests];
    ScenarioEvals(10,:,scenario)=[slope_correction_0,slopeSD_correction_0,R2_correction_0,numCorrect_unbiased_0(1)/numTests,numCorrect_unbiased_0(2)/numTests,numCorrect_unbiased_0(3)/numTests];
    ScenarioEvals(11,:,scenario)=[slope_correction_ApxC,slopeSD_correction_ApxC,R2_correction_ApxC,numCorrect_unbiased_apxC(1)/numTests,numCorrect_unbiased_apxC(2)/numTests,numCorrect_unbiased_apxC(3)/numTests];
    ScenarioEvals(12,:,scenario)=[slope_correction_TC,slopeSD_correction_TC,R2_correction_TC,numCorrect_unbiased_TC(1)/numTests,numCorrect_unbiased_TC(2)/numTests,numCorrect_unbiased_TC(3)/numTests];
    ScenarioEvals(13,:,scenario)=[slope_correction_0C,slopeSD_correction_0C,R2_correction_0C,numCorrect_unbiased_0C(1)/numTests,numCorrect_unbiased_0C(2)/numTests,numCorrect_unbiased_0C(3)/numTests];
    
    ScenarioEvalsRel(:,:,scenario)=ScenarioEvals(:,:,scenario);
    ScenarioEvalsRel(:,1,scenario)=abs(ScenarioEvals(:,1,scenario)-1);
    ScenarioEvalsRel(:,1:2,scenario)=1-ScenarioEvalsRel(:,1:2,scenario);
    %[~,ScenarioEvalRankIndex(:,3:end,scenario)]=sort(ScenarioEvalsRel(:,3:end,scenario),1,'descend');
    [~,ScenarioEvalRankIndex(:,:,scenario)]=sort(ScenarioEvalsRel(:,:,scenario),1,'descend');
    [~,ScenarioEvalRanks(:,1:end-1,scenario)]=sort(ScenarioEvalRankIndex(:,:,scenario));
    ScenarioEvalRanks(:,end,scenario)=mean(ScenarioEvalRanks(:,1:end-1,scenario),2);
    
    ScenarioEvalsRelSD=std(ScenarioEvalsRel(:,:,scenario),0,1);
    ScenarioEvalsRelDiff(:,1:end-1,scenario)=(ScenarioEvalsRel(:,:,scenario)-mean(ScenarioEvalsRel(:,:,scenario),1))./ScenarioEvalsRelSD;
    ScenarioEvalsRelDiff(:,end,scenario)=mean(ScenarioEvalsRelDiff(:,1:end-1,scenario),2);
    
    %subplot(length(meanAbundance_var)/2,2,scenario)
    subplot(2,2,scenario)
    title(['E[m]=' num2str(mean_m,2) ', E[n]=' num2str(mean_n,2) ', E[P]=' num2str(mean_P,2) ', E[K]=' num2str(k_mean,1)],'fontweight','normal')
    hold on
%     scatter(all_Species,all_DsTrue,8,'k','filled')
%     scatter(all_Species,all_Chao1,8,'b')
%     scatter(all_Species,all_Ds_unbiased_apx,8,'r')
    xvals=[1:100];
    window=10;
    sdBounds_obs=NaN(2,100);
    sdBounds_Chao1=NaN(2,100);
    sdBounds_Apx=NaN(2,100);
    sdBoot_obs=NaN(1,100);
    sdBoot_Chao1=NaN(1,100);
    sdBoot_Apx=NaN(1,100);
    for richCtr=1:100
        %get standard deviation bounds from estimates across communities
        indices=find(all_Species>=richCtr-(window-1)/2 & all_Species<=richCtr+(window-1)/2);
        meanEst=nanmean(all_DsTrue(indices));
        stdEst=nanstd(all_DsTrue(indices));
        sdBounds_obs(:,richCtr)=[meanEst,stdEst];
        meanEst=nanmean(all_Chao1(indices));
        stdEst=nanstd(all_Chao1(indices));
        sdBounds_Chao1(:,richCtr)=[meanEst,stdEst];
        meanEst=nanmean(all_Ds_unbiased_apx(indices));
        stdEst=nanstd(all_Ds_unbiased_apx(indices));
        sdBounds_Apx(:,richCtr)=[meanEst,stdEst];
        
        %get mean standard deviations from bootstrapped estimates within
        %communities
        sdBoot_obs(richCtr)=nanmean(all_bootStd_raw(indices));
        sdBoot_Chao1(richCtr)=nanmean(all_bootStd_Chao1(indices));
        sdBoot_Apx(richCtr)=nanmean(all_bootStd_Apx(indices));
    end
%     CV_obs=sdBounds_obs(2,:)./sdBounds_obs(1,:);
%     CV_Chao1=sdBounds_Chao1(2,:)./sdBounds_Chao1(1,:);
%     CV_Apx=sdBounds_Apx(2,:)./sdBounds_Apx(1,:);
%     CV_trend_obs=fitlm(xvals,CV_obs)
%     CV_trend_Chao1=fitlm(xvals,CV_Chao1)
%     CV_trend_Apx=fitlm(xvals,CV_Apx)
    
%     b1=boundedline(xvals',sdBounds_obs(1,:)',sdBounds_obs(2,:)','k','alpha','transparency', 0.1);
%     b2=boundedline(xvals',sdBounds_Chao1(1,:)',sdBounds_Chao1(2,:)','b','alpha','transparency', 0.1);
%     b3=boundedline(xvals',sdBounds_Apx(1,:)',sdBounds_Apx(2,:)','r','alpha','transparency', 0.1);
%     set(b1,'LineStyle',':')
%     set(b2,'LineStyle',':')
%     set(b3,'LineStyle',':')
%     plot(xvals,sdBounds_obs(1,:)+sdBounds_obs(2,:),'k','LineWidth',2);
%     plot(xvals,sdBounds_obs(1,:)-sdBounds_obs(2,:),'k','LineWidth',2);
    plot(xvals,sdBounds_Chao1(1,:)+sdBounds_Chao1(2,:),'b','LineWidth',2);
    plot(xvals,sdBounds_Chao1(1,:)-sdBounds_Chao1(2,:),'b','LineWidth',2);
    plot(xvals,sdBounds_Apx(1,:)+sdBounds_Apx(2,:),'r','LineWidth',2);
    plot(xvals,sdBounds_Apx(1,:)-sdBounds_Apx(2,:),'r','LineWidth',2);
% 
%     plot(xvals,sdBounds_obs(1,:)+sdBoot_obs,':k','LineWidth',2);
%     plot(xvals,sdBounds_obs(1,:)-sdBoot_obs,':k','LineWidth',2);
    plot(xvals,sdBounds_Chao1(1,:)+sdBoot_Chao1,':b','LineWidth',2);
    plot(xvals,sdBounds_Chao1(1,:)-sdBoot_Chao1,':b','LineWidth',2);
    plot(xvals,sdBounds_Apx(1,:)+sdBoot_Apx,':r','LineWidth',2);
    plot(xvals,sdBounds_Apx(1,:)-sdBoot_Apx,':r','LineWidth',2);
    
    
%     hline1=refline(1,0);
%     hline1.Color = 'c';
%     hline1.LineWidth = 4;
    xlabel 'true species richness'
    ylabel 'estimated species richness'
    
    %ylim([0 180])
    ylims=ylim;
    %     text(minSpecies+5*(maxSpecies-minSpecies)/100,ylims(2)*0.95,['Raw slope=' num2str(slope_detected,2) '\pm' num2str(slopeSD_detected,2) ', R^2*=' num2str(R2_detected,2) ', pairwise_{' num2str(maxRichDiff_all(1)) ',' num2str(maxRichDiff_all(2)) ',' num2str(maxRichDiff_all(3)) '}=' num2str(numCorrect_detected(1)/numTests,2) ',' num2str(numCorrect_detected(2)/numTests,2) ',' num2str(numCorrect_detected(3)/numTests,2)])
    %     text(minSpecies+5*(maxSpecies-minSpecies)/100,ylims(2)*0.9,['Chao1 slope=' num2str(slope_Chao1,2) '\pm' num2str(slopeSD_Chao1,2) ', R^2*=' num2str(R2_Chao1,2) ', pairwise_{' num2str(maxRichDiff_all(1)) ',' num2str(maxRichDiff_all(2)) ',' num2str(maxRichDiff_all(3)) '}=' num2str(numCorrect_Chao1(1)/numTests,2) ',' num2str(numCorrect_Chao1(2)/numTests,2) ',' num2str(numCorrect_Chao1(3)/numTests,2)],'Color','b')
    %     text(minSpecies+5*(maxSpecies-minSpecies)/100,ylims(2)*0.85,['GP slope=' num2str(slope_GP,2) '\pm' num2str(slopeSD_GP,2) ', R^2*=' num2str(R2_GP,2) ', pairwise_{' num2str(maxRichDiff_all(1)) ',' num2str(maxRichDiff_all(2)) ',' num2str(maxRichDiff_all(3)) '}=' num2str(numCorrect_GP(1)/numTests,2) ',' num2str(numCorrect_GP(2)/numTests,2) ',' num2str(numCorrect_GP(3)/numTests,2)],'Color',[0.5 0.5 0.5])
    %     text(minSpecies+5*(maxSpecies-minSpecies)/100,ylims(2)*0.8,['Chao2 slope=' num2str(slope_Chao2,2) '\pm' num2str(slopeSD_Chao2,2) ', R^2*=' num2str(R2_Chao2,2) ', pairwise_{' num2str(maxRichDiff_all(1)) ',' num2str(maxRichDiff_all(2)) ',' num2str(maxRichDiff_all(3)) '}=' num2str(numCorrect_Chao2(1)/numTests,2) ',' num2str(numCorrect_Chao2(2)/numTests,2) ',' num2str(numCorrect_Chao2(3)/numTests,2)],'Color',[0.5 0.5 0.5])
    %     text(minSpecies+5*(maxSpecies-minSpecies)/100,ylims(2)*0.75,['ACE slope=' num2str(slope_ACE,2) '\pm' num2str(slopeSD_ACE,2) ', R^2*=' num2str(R2_ACE,2) ', pairwise_{' num2str(maxRichDiff_all(1)) ',' num2str(maxRichDiff_all(2)) ',' num2str(maxRichDiff_all(3)) '}=' num2str(numCorrect_ACE(1)/numTests,2) ',' num2str(numCorrect_ACE(2)/numTests,2)  ',' num2str(numCorrect_ACE(3)/numTests,2)],'Color',[0.5 0.5 0.5])
    %     text(minSpecies+5*(maxSpecies-minSpecies)/100,ylims(2)*0.7,['JK_a slope=' num2str(slope_S_aj2,2) '\pm' num2str(slopeSD_S_aj2,2) ', R^2*=' num2str(R2_S_aj2,2) ', pairwise_{' num2str(maxRichDiff_all(1)) ',' num2str(maxRichDiff_all(2)) ',' num2str(maxRichDiff_all(3)) '}=' num2str(numCorrect_S_aj2(1)/numTests,2) ',' num2str(numCorrect_S_aj2(2)/numTests,2) ',' num2str(numCorrect_S_aj2(3)/numTests,2)],'Color',[0.5 0.5 0.5])
    %     text(minSpecies+5*(maxSpecies-minSpecies)/100,ylims(2)*0.65,['JK_i slope=' num2str(slope_S_ij2,2) '\pm' num2str(slopeSD_S_ij2,2) ', R^2*=' num2str(R2_S_ij2,2) ', pairwise_{' num2str(maxRichDiff_all(1)) ',' num2str(maxRichDiff_all(2)) ',' num2str(maxRichDiff_all(3)) '}=' num2str(numCorrect_S_ij2(1)/numTests,2) ',' num2str(numCorrect_S_ij2(2)/numTests,2) ',' num2str(numCorrect_S_ij2(3)/numTests,2)],'Color',[0.5 0.5 0.5])
    %     text(minSpecies+5*(maxSpecies-minSpecies)/100,ylims(2)*0.6,['\Omega slope=' num2str(slope_correction_Apx,2) '\pm' num2str(slopeSD_correction_Apx,2) ', R^2*=' num2str(R2_correction_Apx,2) ', pairwise_{' num2str(maxRichDiff_all(1)) ',' num2str(maxRichDiff_all(2)) ',' num2str(maxRichDiff_all(3)) '}=' num2str(numCorrect_unbiased_apx(1)/numTests,2) ',' num2str(numCorrect_unbiased_apx(2)/numTests,2) ',' num2str(numCorrect_unbiased_apx(3)/numTests,2)],'Color','r')
    %     text(minSpecies+5*(maxSpecies-minSpecies)/100,ylims(2)*0.55,['\Omega_T slope=' num2str(slope_correction_T,2) '\pm' num2str(slopeSD_correction_T,2) ', R^2*=' num2str(R2_correction_T,2) ', pairwise_{' num2str(maxRichDiff_all(1)) ',' num2str(maxRichDiff_all(2)) ',' num2str(maxRichDiff_all(3)) '}=' num2str(numCorrect_unbiased_T(1)/numTests,2) ',' num2str(numCorrect_unbiased_T(2)/numTests,2) ',' num2str(numCorrect_unbiased_T(3)/numTests,2)],'Color',[0.5 0.5 0.5])
    %     text(minSpecies+5*(maxSpecies-minSpecies)/100,ylims(2)*0.5,['\Omega_0 slope=' num2str(slope_correction_0,2) '\pm' num2str(slopeSD_correction_0,2) ', R^2*=' num2str(R2_correction_0,2) ', pairwise_{' num2str(maxRichDiff_all(1)) ',' num2str(maxRichDiff_all(2)) ',' num2str(maxRichDiff_all(3)) '}=' num2str(numCorrect_unbiased_0(1)/numTests,2) ',' num2str(numCorrect_unbiased_0(2)/numTests,2) ',' num2str(numCorrect_unbiased_0(3)/numTests,2)],'Color',[0.5 0.5 0.5])
    %     text(minSpecies+5*(maxSpecies-minSpecies)/100,ylims(2)*0.45,['\Omega_C slope=' num2str(slope_correction_ApxC,2) '\pm' num2str(slopeSD_correction_ApxC,2) ', R^2*=' num2str(R2_correction_ApxC,2) ', pairwise_{' num2str(maxRichDiff_all(1)) ',' num2str(maxRichDiff_all(2)) ',' num2str(maxRichDiff_all(3)) '}=' num2str(numCorrect_unbiased_apxC(1)/numTests,2) ',' num2str(numCorrect_unbiased_apxC(2)/numTests,2) ',' num2str(numCorrect_unbiased_apxC(3)/numTests,2)],'Color',[0.5 0.5 0.5])
    %     text(minSpecies+5*(maxSpecies-minSpecies)/100,ylims(2)*0.4,['\Omega_{TC} slope=' num2str(slope_correction_TC,2) '\pm' num2str(slopeSD_correction_TC,2) ', R^2*=' num2str(R2_correction_TC,2) ', pairwise_{' num2str(maxRichDiff_all(1)) ',' num2str(maxRichDiff_all(2)) ',' num2str(maxRichDiff_all(3)) '}=' num2str(numCorrect_unbiased_TC(1)/numTests,2) ',' num2str(numCorrect_unbiased_TC(2)/numTests,2) ',' num2str(numCorrect_unbiased_TC(3)/numTests,2)],'Color',[0.5 0.5 0.5])
    %     text(minSpecies+5*(maxSpecies-minSpecies)/100,ylims(2)*0.35,['\Omega_{0C} slope=' num2str(slope_correction_0C,2) '\pm' num2str(slopeSD_correction_0C,2) ', R^2*=' num2str(R2_correction_0C,2) ', pairwise_{' num2str(maxRichDiff_all(1)) ',' num2str(maxRichDiff_all(2)) ',' num2str(maxRichDiff_all(3)) '}=' num2str(numCorrect_unbiased_0C(1)/numTests,2) ',' num2str(numCorrect_unbiased_0C(2)/numTests,2) ',' num2str(numCorrect_unbiased_0C(3)/numTests,2)],'Color',[0.5 0.5 0.5])
    %
    lgd=legend('observed','Chao1','\Omega');
    lgd.AutoUpdate='off';
    xlims=xlim;
    text(xlims(1)-diff(xlims)*0.18,ylims(2)+diff(ylims)*0.12,char(64+splot),'Fontsize',20)
    splot=splot+2;
    
%     subplot(2,2,scenario+2)
%     hold on
%     xVal=repmat(1:size(all_Apx_T_detectP_terms,1)+1,1,simReps); %mean x-axis positions for correction terms
%     xPos1=xVal+0.4*(rand(size(xVal))-1); %add jitter to x-axis positions (left)
%     xPos2=xVal+0.4*(rand(size(xVal))); %add jitter to x-axis positions (right)
%     scatter(xPos1,reshape([sum(all_Apx_T_detectP_terms); all_Apx_T_detectP_terms],1,[]),5,'o','MarkerEdgeColor','none','MarkerFaceColor',[0.72,0.27,1.00],'MarkerFaceAlpha',0.3);
%     scatter(xPos2,reshape([sum(all_Apx_TC_detectP_terms); all_Apx_TC_detectP_terms],1,[]),5,'o','MarkerEdgeColor','none','MarkerFaceColor',[0.93,.69,.13],'MarkerFaceAlpha',0.3);
%     hline2=refline(0,0);
%     hline2.Color = 'k';
%     xline(1.5)
%     xlabel 'correction factor'
%     ylabel 'observation probability'
%     xticks(1:size(all_Apx_T_detectP_terms,1)+1)
%     xaxisproperties= get(gca, 'XAxis');
%     xaxisproperties.TickLabelInterpreter = 'latex'; % latex for x-axis
%     xticklabels({'sum','D$$_i$$(E$$_i$$[$$\hat{\phi}]$$)','c$$_1$$var$$_i$$($$\hat{mn}$$)','c$$_2$$var$$_i$$($$\hat{P}$$)','c$$_3$$cov$$_i$$($$\hat{mn,P}$$)'});
%     xtickangle(15)
%     xlim([0.5 size(all_Apx_T_detectP_terms,1)+1.5])
%     xlims=xlim;
%     ylims=ylim;
%     lgd=legend('\Omega_T','\Omega_{TC}');
%     lgd.AutoUpdate='off';
%     text(xlims(1)-diff(xlims)*0.18,ylims(2)+diff(ylims)*0.12,char(64+splot),'Fontsize',20)
%     splot=splot-1;
    
    
%     %diagnostic plot for species abundance distributions in scenario:
%     figure('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3)/3.5 scrsz(4)/6]);
%     subplot(1,2,1)
%     hold on
%     histogram(all_mean_n,80,'facecolor','w');
%     xlim([0,8])
%     xlabel 'mn across communities'
%     ylabel 'count'
%     
%     subplot(1,2,2)
%     hold on
%     histogram(mean(nk(:,sum(nk,1)>0),1),'facecolor','w');
%     xlabel 'mn within a community'
%     ylabel 'count'
end