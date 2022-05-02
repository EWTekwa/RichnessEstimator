function [Richness_raw,Chao1,Chao2,ACE,S_aj2,S_ij2,Richness_apx,Apx_detectP_terms] = RichnessEstsCov(TransectAbundance)

%RichnessEstsCov.m
%Ed Tekwa Feb 8, 2022 - Apr 11, 2022
%function returns Chao1 and Taylor2 Apx for richness means based on
%the spatial TransectAbundance data: rows=transects, columns=species,
%values=individual counts

TransectAbundance=TransectAbundance(:,sum(TransectAbundance,1)>0); %take out empty species columns

%write observational process model
syms nm C k P
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
Richness_raw=sum(sum(TransectAbundance,1)>0); %get raw richness
C_detected=zeros(1,Richness_raw); %array to record clustering for each species
P_detected=zeros(1,Richness_raw); %array to record occupancy for each species
for species=1:Richness_raw
    if numTrans==1
        C_detected(species)=1+1./mean(TransectAbundance(:,species)); %if only one sampled site, estimate spatial variance as poisson mean
    else
        C_detected(species)=1+var(TransectAbundance(:,species),1)./(mean(TransectAbundance(:,species)).^2); %spatial variance normalized by N (not the normal N-1)
    end
    P_detected(species)=sum(TransectAbundance(:,species)>0)/numTrans; %occupancy as number of transects occupied divided by number of transects
end
C_detected(C_detected==Inf)=NaN;
P_detected(P_detected==0)=NaN;
%compute on full dataset:
mean_C_detected=nanmean(C_detected);
var_C_detected=nanvar(C_detected);
mean_P_detected=nanmean(P_detected);
var_P_detected=nanvar(P_detected);
mean_n_m_detected=mean(mean(TransectAbundance(:,sum(TransectAbundance)>0)));
var_n_m_detected=var(mean(TransectAbundance(:,sum(TransectAbundance)>0)));
cov_nm_C_detected=nancov(mean(TransectAbundance(:,sum(TransectAbundance)>0)),C_detected);
cov_nm_P_detected=nancov(mean(TransectAbundance(:,sum(TransectAbundance)>0)),P_detected);
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

%compute Chao1 estimate
f1=sum(sum(TransectAbundance,1)==1); %number of singleton species
f2=sum(sum(TransectAbundance,1)==2); %number of doubleton species
Chao1=Richness_raw+f1*(f1-1)/(2*(f2+1)); %Chao1 richness estimator

%compute Chao2 estimate
q1=sum(sum(TransectAbundance>1)==1); %number of species occurring in one transect only
q2=sum(sum(TransectAbundance>1)==2); %number of species occurring in two transect only
m=sum(sum(TransectAbundance>1));
Chao2=Richness_raw+((m-1)/m)*q1*(q1-1)/(2*(q2+1)); %Chao2 richness estimator

%compute abundance-based coverage estimator (ACE)
S_rare=sum(sum(TransectAbundance,1)<=10); %number of rare species (<=10 individuals)
S_abund=sum(sum(TransectAbundance,1)>10); %number of rare species (<=10 individuals)
n_rare=sum(TransectAbundance(:,sum(TransectAbundance,1)<=10),'all'); %total number of individuals in rare species
C_ACE=1-f1/n_rare; %sample coverage estimate
wsumfa=0;
for a=1:10
    wsumfa=wsumfa+a*(a-1)*sum(sum(TransectAbundance,1)==a);
end
V2=max(((S_rare/C_ACE)*wsumfa/(n_rare*(n_rare-1))-1),0); %coefficient of variation
if C_ACE>0
    ACE=S_abund+S_rare/C_ACE+(f1/C_ACE)*V2;
else
    ACE=Chao1;
end

%compute jackknife abundance estimator
S_aj2=Richness_raw+2*f1-f2;

%compute jackknife incidence estimator
S_ij2=Richness_raw+(q1*(2*m-3)/m-q2*((m-2)^2)/(m*(m-1)));

%compute correction terms for proposed approximation method
k=numTrans;
C=mean_C_detected;
nm=mean_n_m_detected;
P=mean_P_detected;
Apx_detectP_terms=[eval(Dis)
    eval(d2Dis_dnm2)*var_n_m_detected/2
    eval(d2Dis_dC2)*var_C_detected/2
    eval(d2Dis_dP2)*var_P_detected/2
    eval(d2Dis_dnmC)*cov_nm_C_detected
    eval(d2Dis_dnmP)*cov_nm_P_detected
    eval(d2Dis_dCP)*cov_C_P_detected]; %effects of 7 terms on estimated detection probability

if sum(Apx_detectP_terms)>0.1 %if sum of correction terms is positive and greater than a threshold
    Ds_apx=sum(Apx_detectP_terms); %use full correction
else
    Ds_apx=Apx_detectP_terms(1); %else, use 0th order correction
end
%Ds_apx=eval(Dis)+eval(d2Dis_dnm2)*var_n_m_detected/2+eval(d2Dis_dC2)*var_C_detected/2+eval(d2Dis_dP2)*var_P_detected/2; %Approximated detection probability in community
%Ds_apx=eval(Dis)+eval(d2Dis_dnm2)*var_n_m_detected/2+eval(d2Dis_dC2)*var_C_detected/2+eval(d2Dis_dP2)*var_P_detected/2+eval(d2Dis_dnmC)*cov_nm_C_detected+eval(d2Dis_dnmP)*cov_nm_P_detected+eval(d2Dis_dCP)*cov_C_P_detected; %Approximated detection probability in community
Richness_apx=Richness_raw/Ds_apx;