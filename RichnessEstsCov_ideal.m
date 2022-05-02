function [Richness_apx_T,Apx_T_detectP_terms] = RichnessEstsCov_ideal(TransectAbundance,k,Richness_raw)

%RichnessEstsCov.m
%Ed Tekwa Feb 8, 2022 - Apr 11, 2022
%function returns Chao1 and Taylor2 Apx for richness means based on
%the spatial TransectAbundance data: rows=transects, columns=species,
%values=individual counts

TransectAbundance=TransectAbundance(:,sum(TransectAbundance,1)>0); %take out empty species columns

%write observational process model
syms nm C P
Dix=(1-exp(-C*nm))*P; %local detection probability of species i at sampled site x
Dis=1-(1-Dix)^k; %detection probability of species i across all sampled sites x in community s

%second-order partial derivatives for 2nd order Taylor expansion of overall detection probability of any species in the community

%     d2Dis_dnm2=diff(Dis,'nm',2); %partial derivative of Dis w.r.t. nm (1st or 2nd order). eg. second order: - C^2*P*k*exp(-C*nm)*(P*(exp(-C*nm) - 1) + 1)^(k - 1) - C^2*P^2*k*exp(-2*C*nm)*(k - 1)*(P*(exp(-C*nm) - 1) + 1)^(k - 2)
%     d2Dis_dC2=diff(Dis,'C',2); %partial derivative of Dis w.r.t. C (1st or 2nd order). eg. second order: - P*k*nm^2*exp(-C*nm)*(P*(exp(-C*nm) - 1) + 1)^(k - 1) - P^2*k*nm^2*exp(-2*C*nm)*(k - 1)*(P*(exp(-C*nm) - 1) + 1)^(k - 2)
%     d2Dis_dP2=diff(Dis,'P',2); %partial derivative of Dis w.r.t. P (1st or 2nd order). eg. second order: -k*(k - 1)*(exp(-C*nm) - 1)^2*(P*(exp(-C*nm) - 1) + 1)^(k - 2)
%     d2Dis_dnmC=diff(diff(Dis,'nm'),'C');
%     d2Dis_dnmP=diff(diff(Dis,'nm'),'P');
%     d2Dis_dCP=diff(diff(Dis,'C'),'P');

d2Dis_dnm2=- C^2*P*k*exp(-C*nm)*(P*(exp(-C*nm) - 1) + 1)^(k - 1) - C^2*P^2*k*exp(-2*C*nm)*(k - 1)*(P*(exp(-C*nm) - 1) + 1)^(k - 2); %partial derivative of Dis w.r.t. nm (1st or 2nd order). eg. second order: - C^2*P*k*exp(-C*nm)*(P*(exp(-C*nm) - 1) + 1)^(k - 1) - C^2*P^2*k*exp(-2*C*nm)*(k - 1)*(P*(exp(-C*nm) - 1) + 1)^(k - 2)
d2Dis_dC2=- P*k*nm^2*exp(-C*nm)*(P*(exp(-C*nm) - 1) + 1)^(k - 1) - P^2*k*nm^2*exp(-2*C*nm)*(k - 1)*(P*(exp(-C*nm) - 1) + 1)^(k - 2); %partial derivative of Dis w.r.t. C (1st or 2nd order). eg. second order: - P*k*nm^2*exp(-C*nm)*(P*(exp(-C*nm) - 1) + 1)^(k - 1) - P^2*k*nm^2*exp(-2*C*nm)*(k - 1)*(P*(exp(-C*nm) - 1) + 1)^(k - 2)
d2Dis_dP2=-k*(k - 1)*(exp(-C*nm) - 1)^2*(P*(exp(-C*nm) - 1) + 1)^(k - 2); %partial derivative of Dis w.r.t. P (1st or 2nd order). eg. second order: -k*(k - 1)*(exp(-C*nm) - 1)^2*(P*(exp(-C*nm) - 1) + 1)^(k - 2)
d2Dis_dnmC=P*k*exp(-C*nm)*(P*(exp(-C*nm) - 1) + 1)^(k - 1) - C*P*k*nm*exp(-C*nm)*(P*(exp(-C*nm) - 1) + 1)^(k - 1) - C*P^2*k*nm*exp(-2*C*nm)*(k - 1)*(P*(exp(-C*nm) - 1) + 1)^(k - 2);
d2Dis_dnmP=C*k*exp(-C*nm)*(P*(exp(-C*nm) - 1) + 1)^(k - 1) + C*P*k*exp(-C*nm)*(k - 1)*(exp(-C*nm) - 1)*(P*(exp(-C*nm) - 1) + 1)^(k - 2);
d2Dis_dCP=k*nm*exp(-C*nm)*(P*(exp(-C*nm) - 1) + 1)^(k - 1) + P*k*nm*exp(-C*nm)*(k - 1)*(exp(-C*nm) - 1)*(P*(exp(-C*nm) - 1) + 1)^(k - 2);


numTrans=size(TransectAbundance,1); %get number of transects
Richness_true=sum(sum(TransectAbundance,1)>0); %get raw richness
C_true=zeros(1,Richness_true); %array to record clustering for each species
P_true=zeros(1,Richness_true); %array to record occupancy for each species
for species=1:Richness_true
    if numTrans==1
        C_true(species)=1+1./mean(TransectAbundance(:,species)); %if only one sampled site, estimate spatial variance as poisson mean
    else
        C_true(species)=1+var(TransectAbundance(:,species),1)./(mean(TransectAbundance(:,species)).^2); %spatial variance normalized by N (not the normal N-1)
    end
    P_true(species)=sum(TransectAbundance(:,species)>0)/numTrans; %occupancy as number of transects occupied divided by number of transects
end
C_true(C_true==Inf)=NaN;
P_true(P_true==0)=NaN;
%compute on full dataset:
mean_C_true=nanmean(C_true);
var_C_true=nanvar(C_true);
mean_P_true=nanmean(P_true);
var_P_true=nanvar(P_true);
mean_n_m_true=mean(mean(TransectAbundance(:,sum(TransectAbundance)>0)));
var_n_m_true=var(mean(TransectAbundance(:,sum(TransectAbundance)>0)));
cov_nm_C_true=nancov(mean(TransectAbundance(:,sum(TransectAbundance)>0)),C_true);
cov_nm_P_true=nancov(mean(TransectAbundance(:,sum(TransectAbundance)>0)),P_true);
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


%compute correction terms for proposed approximation method
%k=numTrans;
C=mean_C_true;
nm=mean_n_m_true;
P=mean_P_true;
Apx_T_detectP_terms=[eval(Dis)
    eval(d2Dis_dnm2)*var_n_m_true/2
    eval(d2Dis_dC2)*var_C_true/2
    eval(d2Dis_dP2)*var_P_true/2
    eval(d2Dis_dnmC)*cov_nm_C_true
    eval(d2Dis_dnmP)*cov_nm_P_true
    eval(d2Dis_dCP)*cov_C_P_true]; %effects of 7 terms on estimated detection probability

%Ds_apx=eval(Dis)+eval(d2Dis_dnm2)*var_n_m_true/2+eval(d2Dis_dC2)*var_C_true/2+eval(d2Dis_dP2)*var_P_true/2; %Approximated detection probability in community
%Ds_apx=eval(Dis)+eval(d2Dis_dnm2)*var_n_m_true/2+eval(d2Dis_dC2)*var_C_true/2+eval(d2Dis_dP2)*var_P_true/2+eval(d2Dis_dnmC)*cov_nm_C_true+eval(d2Dis_dnmP)*cov_nm_P_true+eval(d2Dis_dCP)*cov_C_P_true; %Approximated detection probability in community
%if sum(abs(Apx_T_detectP_terms)>1)<1 %if no correction term has magnitude over 1
if sum(Apx_T_detectP_terms)>0.1 %if sum of correction terms is positive and greater than a threshold
    Ds_apx=sum(Apx_T_detectP_terms); %use full correction
else
    Ds_apx=Apx_T_detectP_terms(1); %else, use 0th order correction
end
if Richness_raw>0
    Richness_apx_T=Richness_raw/Ds_apx;
else
    Richness_apx_T=NaN;
end