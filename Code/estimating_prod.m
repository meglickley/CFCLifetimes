% Production priors

clear all
HomeDir = '/Users/meganlickley/Dropbox (MIT)/Research/Lifetimes/Drafts/NatComms_Resubmission/Code';
cd(HomeDir)

y1 = 1955; % Year Start Date
yswitch = 1988; % Year switching between AFEAS and UNEP
y2 = 2018; % Year End Date
Sim_yrs = y1:y2;
nSamps = 5*10^5;

reportedProduction = NaN(size(Sim_yrs)); % Production is in Tonnes
load('Input/CFC11/WMO2002.mat') % Use until 1988 - WMO data is the adjusted AFEAS data
ytmp1 = find(year==y1); 
ytmp2 = find(year==yswitch);
reportedProduction(1:yswitch+1-y1) = Production(ytmp1:ytmp2);

load('Input/CFC11/a5_na5.mat') % Use Article 5 and non A5 prod total from 1989 onwards
ytmp1 = find(a5_na5(:,1) == yswitch+1);
ytmp2 = find(Sim_yrs == a5_na5(end,1));
reportedProduction(1989-y1+1:ytmp2) = sum(a5_na5(ytmp1:end,[2:3]),2);
reportedProduction(ytmp2:end) = sum(a5_na5(end,[2:3]),2);

%% creating production distribution
nYears = length(y1:yswitch);
rho(1,1,:) = 0.5 + 0.5*betarnd(2,2,nSamps,1);
rho_tmp = repmat(rho,nYears,nYears,1);
exp_val = repmat(abs(repmat([1:nYears],nYears,1)-repmat([1:nYears]',1,nYears)),1,1,nSamps);
rhoMatrix = rho_tmp.^exp_val;
clear rho_tmp exp_val
Zm = mvnrnd(zeros(nSamps,nYears), rhoMatrix, nSamps);
clear Rhom
Um = normcdf(Zm,0,1);
clear Zm
productionSamps(:,1:nYears) = (repmat(0.2*reportedProduction(1:nYears),nSamps,1)).*logninv(Um,zeros(nSamps,nYears),0.5*ones(nSamps,nYears))+0.95*repmat(reportedProduction(1:nYears),nSamps,1);
clear Um

nYears2 = length(Sim_yrs) - nYears; 
rho_tmp = repmat(rho, nYears2, nYears2, 1);

exp_val = repmat(abs(repmat([1:nYears2],nYears2,1)-repmat([1:nYears2]',1,nYears2)),1,1,nSamps);
rhoMatrix = rho_tmp.^exp_val;
clear rho_tmp exp_val
Zm = mvnrnd(zeros(nSamps, nYears2), rhoMatrix, nSamps);
clear Rhom
Um = normcdf(Zm,0,1);
clear Zm
productionSamps(:,nYears+1:length(Sim_yrs)) = (repmat(0.1*reportedProduction(nYears+1:end),nSamps,1)).*logninv(Um,zeros(nSamps,nYears2),0.5*ones(nSamps,nYears2))+0.95*repmat(reportedProduction(nYears+1:end),nSamps,1);
clear Um
save('Input/CFC11/Production_samps.mat','rho','productionSamps','reportedProduction');


%% CFC11 Fugitive Production Scenario

Nyrs3 = length(2000:y2);
rho2(1,1,:) = 0.5+0.5*betarnd(2,2,nSamps,1);
rho_tmp = repmat(rho2,Nyrs3,Nyrs3,1);
exp_val = repmat(abs(repmat([1:Nyrs3],Nyrs3,1)-repmat([1:Nyrs3]',1,Nyrs3)),1,1,nSamps);
rhoMatrix = rho_tmp.^exp_val;
clear rho_tmp exp_val
Zm = mvnrnd(zeros(nSamps,Nyrs3), rhoMatrix, nSamps);
clear Rhom
Um = normcdf(Zm,0,1);
clear Zm

excess_prod = 61*10^3; % Montzka 
AddtoUB = [linspace(0,excess_prod,13),repmat(excess_prod,1,6)];
unexpectedProduction = reportedProduction;
unexpectedProduction(end-Nyrs3+1:end) = reportedProduction(end-Nyrs3+1:end)+AddtoUB;

productionSamps(:,end-Nyrs3+1:end) = repmat(reportedProduction(end-Nyrs3+1:end),nSamps,1)+repmat(unexpectedProduction(end-Nyrs3+1:end)-reportedProduction(end-Nyrs3+1:end),nSamps,1).*betainv(Um,2*ones(nSamps,Nyrs3),2*ones(nSamps,Nyrs3));

clear Um
save('Input/CFC11/unexpectedProduction_samps.mat','rho','productionSamps','reportedProduction','unexpectedProduction');


%%%%%%%%%%%%%%%%   CFC 12 %%%%%%%%%%%%%%%%%%%%
%% 

clearvars -except HomeDir
cd(HomeDir)

y1 = 1955; % Year Start Date
yswitch = 1988;
y2 = 2018; % Year End Date
Sim_yrs = y1:y2;
nSamps = 5*10^5;

reportedProduction = NaN(size(Sim_yrs)); % Production is in Tonnes
load('Input/CFC12/WMO2002.mat') % Use until 1988 - WMO data is the adjusted AFEAS data
ytmp1 = find(year==y1); 
ytmp2 = find(year==yswitch);
reportedProduction(1:yswitch+1-y1) = Production(ytmp1:ytmp2);

load('Input/CFC12/a5_na5.mat') % Use Article 5 and non A5 prod total from 1989 onwards
ytmp1 = find(a5_na5(:,1) == yswitch+1);
ytmp2 = find(Sim_yrs == a5_na5(end,1));
reportedProduction(1989-y1+1:ytmp2) = sum(a5_na5(ytmp1:end,[2:3]),2);
reportedProduction(ytmp2:end) = sum(a5_na5(end,[2:3]),2);

% creating production distribution
nYears = length(y1:yswitch);
rho(1,1,:) = 0.5+0.5*betarnd(2,2,nSamps,1);
rho_tmp = repmat(rho,nYears,nYears,1);
exp_val = repmat(abs(repmat([1:nYears],nYears,1)-repmat([1:nYears]',1,nYears)),1,1,nSamps);
rhoMatrix = rho_tmp.^exp_val;
clear rho_tmp exp_val
Zm = mvnrnd(zeros(nSamps,nYears), rhoMatrix, nSamps);
clear Rhom
Um = normcdf(Zm,0,1);
clear Zm
productionSamps(:,1:nYears) = (repmat(0.2*reportedProduction(1:nYears),nSamps,1)).*logninv(Um,zeros(nSamps,nYears),0.5*ones(nSamps,nYears))+0.95*repmat(reportedProduction(1:nYears),nSamps,1);
clear Um

nYears2 = length(Sim_yrs) - nYears; 
rho_tmp = repmat(rho,nYears2,nYears2,1);

exp_val = repmat(abs(repmat([1:nYears2],nYears2,1)-repmat([1:nYears2]',1,nYears2)),1,1,nSamps);
rhoMatrix = rho_tmp.^exp_val;
clear rho_tmp exp_val
Zm = mvnrnd(zeros(nSamps,nYears2), rhoMatrix, nSamps);
clear Rhom
Um = normcdf(Zm,0,1);
clear Zm
productionSamps(:,nYears+1:length(Sim_yrs)) = (repmat(0.1*reportedProduction(nYears+1:end),nSamps,1)).*logninv(Um,zeros(nSamps,nYears2),0.5*ones(nSamps,nYears2))+0.95*repmat(reportedProduction(nYears+1:end),nSamps,1);
clear Um
save('Input/CFC12/Production_samps.mat','rho','productionSamps','reportedProduction');

%%%%%%%%%%%%%%%%   CFC 113 %%%%%%%%%%%%%%%%%%%%
%% 

clearvars -except HomeDir
cd(HomeDir)

y1 = 1955; % Year Start Date
yswitch = 1988;
y2 = 2018; % Year End Date
Sim_yrs = y1:y2;
nSamps = 5*10^5;

reportedProduction = NaN(size(Sim_yrs)); % Production is in Tonnes
load('Input/CFC113/WMO2002.mat') % Use until 1988 - WMO data is the adjusted AFEAS data
ytmp1 = find(year==y1); 
ytmp2 = find(year==yswitch);
reportedProduction(1:yswitch+1-y1) = Production(ytmp1:ytmp2);

load('Input/CFC113/a5_na5.mat') % Use Article 5 and non A5 prod total from 1989 onwards
ytmp1 = find(a5_na5(:,1) == yswitch+1);
ytmp2 = find(Sim_yrs == a5_na5(end,1));
reportedProduction(1989-y1+1:ytmp2) = sum(a5_na5(ytmp1:end,[2:3]),2);
reportedProduction(ytmp2:end) = sum(a5_na5(end,[2:3]),2);

% creating production distribution
nYears = length(y1:yswitch);
rho(1,1,:) = 0.3+0.7*betarnd(2,2,nSamps,1);
rho_tmp = repmat(rho,nYears,nYears,1);
exp_val = repmat(abs(repmat([1:nYears],nYears,1)-repmat([1:nYears]',1,nYears)),1,1,nSamps);
rhoMatrix = rho_tmp.^exp_val;
clear rho_tmp exp_val
Zm = mvnrnd(zeros(nSamps,nYears), rhoMatrix, nSamps);
clear Rhom
Um = normcdf(Zm,0,1);
clear Zm
productionSamps(:,1:nYears) = (repmat(0.4*reportedProduction(1:nYears),nSamps,1)).*logninv(Um,zeros(nSamps,nYears),0.5*ones(nSamps,nYears))+0.70*repmat(reportedProduction(1:nYears),nSamps,1);
clear Um

nYears2 = length(Sim_yrs) - nYears; 
rho_tmp = repmat(rho,nYears2,nYears2,1);

exp_val = repmat(abs(repmat([1:nYears2],nYears2,1)-repmat([1:nYears2]',1,nYears2)),1,1,nSamps);
rhoMatrix = rho_tmp.^exp_val;
clear rho_tmp exp_val
Zm = mvnrnd(zeros(nSamps,nYears2), rhoMatrix, nSamps);
clear Rhom
Um = normcdf(Zm,0,1);
clear Zm
productionSamps(:,nYears+1:length(Sim_yrs)) = (repmat(0.4*reportedProduction(nYears+1:end),nSamps,1)).*logninv(Um,zeros(nSamps,nYears2),0.5*ones(nSamps,nYears2))+0.70*repmat(reportedProduction(nYears+1:end),nSamps,1);
clear Um
save('Input/CFC113/Production_samps.mat','rho','productionSamps','reportedProduction');

%% CFC113 Fugitive Production Scenario

Nyrs3 = length(2000:y2);
rho2(1,1,:) = 0.3+0.7*betarnd(2,2,nSamps,1);
rho_tmp = repmat(rho2,Nyrs3,Nyrs3,1);
exp_val = repmat(abs(repmat([1:Nyrs3],Nyrs3,1)-repmat([1:Nyrs3]',1,Nyrs3)),1,1,nSamps);
rhoMatrix = rho_tmp.^exp_val;
clear rho_tmp exp_val
Zm = mvnrnd(zeros(nSamps,Nyrs3), rhoMatrix, nSamps);
clear Rhom
Um = normcdf(Zm,0,1);
clear Zm

excess_prod = 15*10^3; % Our upper bound
AddtoUB = repmat(excess_prod,1,Nyrs3);
unexpectedProduction = reportedProduction;
unexpectedProduction(end-Nyrs3+1:end) = reportedProduction(end-Nyrs3+1:end)+AddtoUB;

productionSamps(:,end-Nyrs3+1:end) = repmat(reportedProduction(end-Nyrs3+1:end),nSamps,1)+repmat(unexpectedProduction(end-Nyrs3+1:end)-reportedProduction(end-Nyrs3+1:end),nSamps,1).*betainv(Um,2*ones(nSamps,Nyrs3),2*ones(nSamps,Nyrs3));

clear Um
save('Input/CFC113/unexpectedProduction_samps.mat','rho','productionSamps','reportedProduction','unexpectedProduction');

