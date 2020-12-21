%   Main BPE Script for Lickley et al. (in Revisions) 
%   Title: Joint inference of CFC lifetimes and banks suggests previously 
%   unidentified emissions
%
%   Authors: Megan Lickley, Sarah Fletcher, Matt Rigby and Susan Solomon
%   Corresponding Author: Megan Lickley; mlickley@mit.edu
%  
%   Last udpated Dec 17, 2020 by Megan Lickley.
%   This script runs the BPE model for joint inference of CFC-11, 12 and
%   113 lifetimes, emissions, and banks. 
%   First run the submodelRFandDE.m script. 
%   Then run estimating_prod.m script. 

% 	Run for three scenarios: 
%     - unreportedProd = 1, jointInference = 1; 
%	  - unreportedProd = 0, jointInference = 1; 	
%	  - unreportedProd = 1, jointInference = 0; 
%   Then run MakeFigures.m script to generate figures and table values  


close all
clear all
y1 = 1955; % Year Start Date
y2 = 2010; % Year End Date
simulationYears = y1:y2;
nYears = length(simulationYears);
nSamps = 10^5;  % Default =  10^6
nResamps = 3*10^4; % Default =  3*10^5
nResamps2 = 10^4;  % Default =  10^5

HomeDir = '/Users/meganlickley/Dropbox (MIT)/Research/Lifetimes/Drafts/NatComms_Resubmission/Code'; % Change HomeDirectory name

OutputFileName = 'Run1';  % Can change if running multiple iterations for sensitivity test

% Scenarios:
unreportedProd = 1;     % Either 1 or 0
jointInference = 1;     % Either 1 or 0

%% Output File Names

if unreportedProd
    OutputFileName = strcat(OutputFileName,'_unreportedProd');
end

if jointInference
    OutputFileName = strcat(OutputFileName,'_jointInference');
else
    OutputFileName = strcat(OutputFileName,'_independent');    
end

%% Input parameters
CFC11 = 1;
CFC12 = 2; 
CFC113 = 3;

CFC.ppt_to_tonnes(CFC11) = 22602.38457;
CFC.ppt_to_tonnes(CFC12) = 19895.36;
CFC.ppt_to_tonnes(CFC113) = 30834.68811;

CFC.InputFileName{CFC11} = 'CFC11';
CFC.InputFileName{CFC12} = 'CFC12';
CFC.InputFileName{CFC113} = 'CFC113';

CFC.Sparc_ID(CFC11) = 1;
CFC.Sparc_ID(CFC12) = 2;
CFC.Sparc_ID(CFC113) = 3;

CFC.InitialBankSize(CFC11) = 5893.9; %From WMO 2002
CFC.InitialBankSize(CFC12) = 65198;
CFC.InitialBankSize(CFC113) = 1872;

% Starting values for MF 
CFC.MF1955(CFC11) = 3.29;
CFC.MF1955(CFC12) = 14.29;
CFC.MF1955(CFC113) = 1.25;
%% Reading in DE and RF priors - from bottom-up accounting model
if unreportedProd == 0
    load('Input/CFC11/InferredRFandDE.mat')
elseif  unreportedProd == 1 
    load('Input/CFC11/InferredRFandDE_unreportedProd.mat')
end

ytmp1 = find(yearRFDE==y1);
ytmp2 = find(yearRFDE==y2);

tmp = datasample([1:size(RF,2)],nSamps);
RFDE11sampleIndex = tmp;
RF11_prior = RF(ytmp1:ytmp2,tmp);
DE11_prior = DE(ytmp1:ytmp2,tmp);
RF11_prior(RF11_prior<0) = 0;

load('Input/CFC12/InferredRFandDE.mat')
ytmp1 = find(yearRFDE==y1);
ytmp2 = find(yearRFDE==y2);

tmp = datasample([1:size(RF,2)],nSamps);
RFDE12sampleIndex = tmp;
RF12_prior = RF(ytmp1:ytmp2,tmp);
DE12_prior = DE(ytmp1:ytmp2,tmp);
RF12_prior(RF12_prior<0) = 0;


if unreportedProd == 0
    load('Input/CFC113/InferredRFandDE.mat')
else 
    load('Input/CFC113/InferredRFandDE_unreportedProd.mat')
end

 
tmp = datasample([1:size(RF,2)],nSamps);
RFDE113sampleIndex = tmp;
ytmp1 = find(yearRFDE==y1);
ytmp2 = find(yearRFDE==y2);
RF113_prior = RF(ytmp1:ytmp2,tmp);
DE113_prior = DE(ytmp1:ytmp2,tmp);
RF113_prior(RF113_prior<0) = 0;

clear RF DE

%% Reading in Production priors
if unreportedProd == 0
    load('Input/CFC11/Production_samps.mat')
else 
    load('Input/CFC11/unexpectedProduction_samps.mat')
end

Prod11_prior = datasample(productionSamps,nSamps,1);

load('Input/CFC12/Production_samps.mat')
Prod12_prior = datasample(productionSamps,nSamps,1);

if unreportedProd == 0
    load('Input/CFC113/Production_samps.mat')
else
    load('Input/CFC113/unexpectedProduction_samps.mat')
end
    
Prod113_prior = datasample(productionSamps,nSamps,1);

%% Creating Lifetime Priors
SP_yrs = 1960:2010;
load('Input/SPARC_lifetimes.mat')

tmp = squeeze(LT_sparc(:,:,CFC11)); 
InvLT{1} = 1./tmp;

tmp = squeeze(LT_sparc(:,:,CFC12)); 
InvLT{2} = 1./tmp;

tmp = squeeze(LT_sparc(:,:,CFC113)); 
InvLT{3} = 1./tmp;

% Computing inverse LT 10 year moving average, var name: InvLTma
for mol = 1:3
    for ii = 1:42
        InvLTma{mol}(ii,:) = mean(InvLT{mol}(ii:ii+9,:));
    end
    for mod_ii = 1:7
        indx = isnan(InvLTma{mol}(:,mod_ii));
        InvLTma{mol}(indx,mod_ii) = nanmean(InvLT{mol}(sum(1-indx):end,mod_ii));
    end
end

tmp = [InvLTma{1}; InvLTma{2}; InvLTma{3}];

cov_LT = nancov(tmp');
mean_LT = nanmean(tmp');

% Sample from mean and covariance of 10-yr moving average from SPARC models
LT_priors_short = mvnrnd(mean_LT, cov_LT, nSamps);

% Extend LT priors at the beginning and end of SPARC simualations
for mol = 1:3
    LT_priors{mol}(:,1:SP_yrs(1)-y1+5) = repmat(LT_priors_short(:,1+42*(mol-1)),1,SP_yrs(1)-y1+5); 
    LT_priors{mol}(:,SP_yrs(1)-y1+6:SP_yrs(1)-y1+47) = LT_priors_short(:,1+42*(mol-1):42*mol);
    LT_priors{mol}(:,52:nYears) = repmat(LT_priors_short(:,42*mol),1,nYears-52+1);    
end

LT_prior11 = LT_priors{1};
LT_prior12 = LT_priors{2};
LT_prior113 = LT_priors{3};

%% Reading in Observed Mole Fractions

load('Input/CFC11/wmo2018.mat');
ytmp1 = find(wmo_yr == y1);
ytmp2 = find(wmo_yr == y2);
molefractions{CFC11} = wmo_conc(ytmp1:ytmp2+1);

load('Input/CFC12/wmo2018.mat');
ytmp1 = find(wmo_yr == y1);
ytmp2 = find(wmo_yr == y2);
molefractions{CFC12} = wmo_conc(ytmp1:ytmp2+1);

load('Input/CFC113/wmo2018.mat');
ytmp1 = find(wmo_yr == y1);
ytmp2 = find(wmo_yr == y2);
molefractions{CFC113} = wmo_conc(ytmp1:ytmp2+1);

%% Computing observationally Derived emissions - used in MakeFigures.m file

for mol_ii = 1:3
    MF = repmat(molefractions{mol_ii}',nSamps,1);
    A = CFC.ppt_to_tonnes(mol_ii); % A from eq 1. 
    invLT = LT_priors{mol_ii};
	obsDerivedEmiss{mol_ii} = A*(MF(:,2:end)-MF(:,1:end-1).*exp(-invLT ));
end
clear MF A invLT

if jointInference
    %% Run the Simulation Model for CFC-11 and 12
    [Bank11, Emiss11, MF11,indx1_11, indx2_11] = simulation_model(Prod11_prior, DE11_prior, RF11_prior, CFC.InitialBankSize(CFC11), nSamps, nYears,CFC.MF1955(CFC11),LT_prior11,CFC.ppt_to_tonnes(CFC11) );
    [Bank12, Emiss12, MF12,indx1_12, indx2_12] = simulation_model(Prod12_prior, DE12_prior, RF12_prior, CFC.InitialBankSize(CFC12), nSamps, nYears,CFC.MF1955(CFC12),LT_prior12,CFC.ppt_to_tonnes(CFC12));

    %% Compute Likelihood for each sample 
    ytmp1 = find(simulationYears == 1980);
    ytmp2 = find(simulationYears == y2);
    nRows = length(ytmp1:ytmp2);

    MF_sigma_prior{1} = 0.03*molefractions{1}(ytmp1:ytmp2)';
    MF_sigma_prior{2} = 0.03*molefractions{2}(ytmp1:ytmp2)';
    MF_sigma_prior{3} = 0.055*molefractions{3}(ytmp1:ytmp2)';

    rho_err{1} = 0.99*ones(nSamps,1);
    rho_err{2} = 0.99*ones(nSamps,1);
    rho_err{3} = 0.98*ones(nSamps,1);

    exp_val = abs(repmat([1:nRows],nRows,1)-repmat([1:nRows]',1,nRows));


    parfor ii = 1:nSamps
        rho_tmp = rho_err{1}(ii)*ones(nRows,nRows);
        Rhom11 = rho_tmp.^exp_val;

        Cov_matrix1 = Rhom11.*(MF_sigma_prior{1}'*MF_sigma_prior{1});
        Cov_matrix2 = Rhom11.*(MF_sigma_prior{2}'*MF_sigma_prior{2});

        % likelihood of data (molefractions) given modeled MF (MF11)
        likelihood1(ii) = mvnpdf(MF11(ii,ytmp1:ytmp2)',molefractions{1}(ytmp1:ytmp2),Cov_matrix1);
        likelihood2(ii) = mvnpdf(MF12(ii,ytmp1:ytmp2)',molefractions{2}(ytmp1:ytmp2),Cov_matrix2);
    end

    likelihood1 = (1/sum(likelihood1))*likelihood1;
    likelihood2 = (1/sum(likelihood2))*likelihood2;
    likelihood = likelihood1.*likelihood2; 

    % Resample from the priors based on the relative likelihood
    ResampleIndex1 = nan(nResamps,1);

    % Updating priors based on SIR/likelihood 
    importanceRatioVector = 1/sum(likelihood)*cumsum(likelihood);

    parfor ii = 1:nResamps
        ResampleIndex1(ii) = find(importanceRatioVector>rand,1);
    end


    %% Posteriors from CFC-11,12|Data_11,12 become priors
    % Only lifetime priors for CFC113 are conditioned on data
    LT_prior113 = LT_priors{3}(ResampleIndex1,:);

    % Sampling from priors that are not conditioned on Data 
    Prod113_prior = Prod113_prior(1:nResamps,:);
    DE113_prior = DE113_prior(:,1:nResamps);
    RF113_prior = RF113_prior(:,1:nResamps);

    [Bank113, Emiss113, MF113,indx1_113, indx2_113] = simulation_model(Prod113_prior, DE113_prior, RF113_prior, CFC.InitialBankSize(CFC113), nResamps, nYears,CFC.MF1955(CFC113),LT_prior113,CFC.ppt_to_tonnes(CFC113));

    parfor ii = 1:nResamps

        rho_tmp = rho_err{3}(ii)*ones(nRows,nRows);
        Rhom113 = rho_tmp.^exp_val;
        Cov_matrix3 = Rhom113.*(MF_sigma_prior{3}'*MF_sigma_prior{3});

        % likelihood of data (molefractions) given modeled MF (MF11)
        likelihood3(ii) = mvnpdf(MF113(ii,ytmp1:ytmp2)',molefractions{3}(ytmp1:ytmp2),Cov_matrix3);
    end

    ResampleIndex3 = nan(nResamps2,1);
    importanceRatioVector = 1/sum(likelihood3)*cumsum(likelihood3);

    parfor ii = 1:nResamps2
        ResampleIndex3(ii) = find(importanceRatioVector>rand,1);
    end

    %% Posterior sample indices
    % ResampleIndex3 are based on priors used in simulation model for CFC-113
    % given posteriors from CFC-11 and 12.  Lifetime posteriors are for all
    % molecules sample from priors for CFC-11 and 12.  
    ResampleIndex1 = ResampleIndex1(ResampleIndex3);
    ResampleIndex2 = ResampleIndex1;
    ResampleIndex_LT = ResampleIndex1;
else
    %% Run the Simulation Model
    [Bank11, Emiss11, MF11,indx1_11, indx2_11] = simulation_model(Prod11_prior, DE11_prior, RF11_prior, CFC.InitialBankSize(CFC11), nSamps, nYears,CFC.MF1955(CFC11),LT_prior11,CFC.ppt_to_tonnes(CFC11) );
    [Bank12, Emiss12, MF12,indx1_12, indx2_12] = simulation_model(Prod12_prior, DE12_prior, RF12_prior, CFC.InitialBankSize(CFC12), nSamps, nYears,CFC.MF1955(CFC12),LT_prior12,CFC.ppt_to_tonnes(CFC12));
    [Bank113, Emiss113, MF113,indx1_113, indx2_113] = simulation_model(Prod113_prior, DE113_prior, RF113_prior, CFC.InitialBankSize(CFC113), nSamps, nYears,CFC.MF1955(CFC113),LT_prior113,CFC.ppt_to_tonnes(CFC113));

    %% Compute Likelihood for each sample 
    ytmp1 = find(simulationYears == 1980);
    ytmp2 = find(simulationYears == y2);
    nRows = length(ytmp1:ytmp2);

    MF_sigma_prior{1} = (0.02*ones(nSamps,1)+0.015*betarnd(5,5,nSamps,1))*molefractions{1}(ytmp1:ytmp2)';
    MF_sigma_prior{2} = (0.02*ones(nSamps,1)+0.015*betarnd(5,5,nSamps,1))*molefractions{2}(ytmp1:ytmp2)';
    MF_sigma_prior{3} = (0.03*ones(nSamps,1)+0.03*betarnd(5,5,nSamps,1))*molefractions{3}(ytmp1:ytmp2)';

    rho_err{1} = 0.99*ones(nSamps,1);
    rho_err{2} = 0.99*ones(nSamps,1);
    rho_err{3} = 0.98*ones(nSamps,1);

    exp_val = abs(repmat([1:nRows],nRows,1)-repmat([1:nRows]',1,nRows));


    parfor ii = 1:nSamps
        rho_tmp = rho_err{1}(ii)*ones(nRows,nRows);
        Rhom11 = rho_tmp.^exp_val;

        rho_tmp = rho_err{3}(ii)*ones(nRows,nRows);
        Rhom113 = rho_tmp.^exp_val;
        Cov_matrix1 = Rhom11.*(MF_sigma_prior{1}(ii,:)'*MF_sigma_prior{1}(ii,:));
        Cov_matrix2 = Rhom11.*(MF_sigma_prior{2}(ii,:)'*MF_sigma_prior{2}(ii,:));
        Cov_matrix3 = Rhom113.*(MF_sigma_prior{3}(ii,:)'*MF_sigma_prior{3}(ii,:));

        % likelihood of data (molefractions) given modeled MF (MF11)
        likelihood1(ii) = mvnpdf(MF11(ii,ytmp1:ytmp2)',molefractions{1}(ytmp1:ytmp2),Cov_matrix1);
        likelihood2(ii) = mvnpdf(MF12(ii,ytmp1:ytmp2)',molefractions{2}(ytmp1:ytmp2),Cov_matrix2);
        likelihood3(ii) = mvnpdf(MF113(ii,ytmp1:ytmp2)',molefractions{3}(ytmp1:ytmp2),Cov_matrix3);
    end

    ResampleIndex1 = nan(nResamps2,1);
    ResampleIndex2 = nan(nResamps2,1);
    ResampleIndex3 = nan(nResamps2,1);
    
   
    importanceRatioVector1 = 1/sum(likelihood1)*cumsum(likelihood1);
    importanceRatioVector2 = 1/sum(likelihood2)*cumsum(likelihood2);
    importanceRatioVector3 = 1/sum(likelihood3)*cumsum(likelihood3);

    parfor ii = 1:nResamps2
        ResampleIndex1(ii) = find(importanceRatioVector1>rand,1);
        ResampleIndex2(ii) = find(importanceRatioVector2>rand,1);
        ResampleIndex3(ii) = find(importanceRatioVector3>rand,1);
    end
    ResampleIndex_LT = nan(nResamps2,1);
    
end
%% figures to verify model:

% Plot of priors and posteriors of mole fractions against observations. 
figure; 
for mol_ii = 1:3
    
    if mol_ii == 1;       var = MF11;     tmpIndex = ResampleIndex1; 
    elseif mol_ii == 2;   var = MF12;     tmpIndex = ResampleIndex2; 
    elseif mol_ii == 3;   var = MF113;    tmpIndex = ResampleIndex3; 
    end
    
    subplot(1,3,mol_ii); 
    
    MED = prctile(var,50);  LB = prctile(var,2.5);   UB = prctile(var,97.5);
    p1 = boundedline([y1:y2+1],MED',[MED'-LB',UB'-MED'],'alpha','cmap',[0.3,0.3,0.3]);
    hold on; 
    p2 = plot(y1:y2,molefractions{mol_ii}(1:end-1),'--b','LineWidth',2);
    
    var1 = var(tmpIndex,:); 
    MED = prctile(var1,50); LB = prctile(var1,2.5); UB = prctile(var1,97.5);
    hold on; 
    p3 = boundedline([y1:y2+1],MED',[MED'-LB',UB'-MED'],'alpha','cmap',[0.7,0,0]); box on;
    title(CFC.InputFileName{mol_ii});
    ylabel('molefractions [pmol mol^{-1}]'); xlim([1955,2010]);
end
lgd = legend([p1, p3, p2], 'priors', 'posteriors', 'observed')
lgd.Location = 'southeast'; 


% Figure of priors and posteriors of atmospheric lifetimes
figure; 
for mol_ii = 1:3
    if jointInference
        indexLT = ResampleIndex_LT;
    else
        if mol_ii == 1;     indexLT = ResampleIndex1;
        elseif mol_ii == 2; indexLT = ResampleIndex2;
        elseif mol_ii == 3; indexLT = ResampleIndex3;
        end
    end
        
    subplot(3,1,mol_ii)
    MED = prctile(LT_priors{mol_ii},50);
    LB = prctile(LT_priors{mol_ii},5);
    UB = prctile(LT_priors{mol_ii},95);
    boundedline([y1:y2],1./MED',[1./MED'-1./LB',1./UB'-1./MED'],'alpha','cmap',[0.3,0.3,0.3]);

    hold on;
    MED = prctile(LT_priors{mol_ii}(indexLT,:),50);
    LB = prctile(LT_priors{mol_ii}(indexLT,:),5);
    UB = prctile(LT_priors{mol_ii}(indexLT,:),95);
    boundedline([y1:y2],1./MED',[1./MED'-1./LB',1./UB'-1./MED'],'alpha','cmap',[0.7,0,0]);
    str = strcat(CFC.InputFileName{mol_ii},': lifetime Priors and posteriors'); title(str);
    ylabel('Lifetime [yrs]'); xlim([1955,2010]);
end


OutputFolder_FileName = strcat(HomeDir,'/Output/CFC11/',OutputFileName); 
LT11 = LT_priors{1};
DEsamps11 = DE11_prior(:,indx1_11)';
RFsamps11 = RF11_prior(:,indx1_11)';
Prodsamps11 = Prod11_prior(indx2_11,:);
ObsEmiss11 = obsDerivedEmiss{1};
Sigma11 = MF_sigma_prior{1};
save(OutputFolder_FileName,'DEsamps11','RFsamps11','Sigma11','Prodsamps11','Bank11', 'Emiss11', 'MF11','LT11','ObsEmiss11','molefractions','ResampleIndex1','ResampleIndex_LT','RFDE11sampleIndex')

OutputFolder_FileName = strcat(HomeDir, '/Output/CFC12/',OutputFileName);
LT12 = LT_priors{2};
DEsamps12 = DE12_prior(:,indx1_12)';
RFsamps12 = RF12_prior(:,indx1_12)';
Prodsamps12 = Prod12_prior(indx2_12,:);
ObsEmiss12 = obsDerivedEmiss{2};
Sigma12 = MF_sigma_prior{2};
save(OutputFolder_FileName,'DEsamps12','RFsamps12','Prodsamps12','Bank12', 'Emiss12', 'MF12','LT12','ObsEmiss12','molefractions','ResampleIndex2','ResampleIndex_LT','RFDE12sampleIndex')

OutputFolder_FileName = strcat(HomeDir, '/Output/CFC113/',OutputFileName); 
LT113 = LT_priors{3};
DEsamps113 = DE113_prior(:,indx1_113)';
RFsamps113 = RF113_prior(:,indx1_113)';
Prodsamps113 = Prod113_prior(indx2_113,:);
ObsEmiss113 = obsDerivedEmiss{3};
Sigma113 = MF_sigma_prior{3};
save(OutputFolder_FileName,'DEsamps113','RFsamps113','Sigma113','Prodsamps113','Bank113', 'Emiss113', 'MF113','LT113','ObsEmiss113','molefractions','ResampleIndex3','ResampleIndex_LT','RFDE113sampleIndex')
