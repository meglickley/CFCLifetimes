%estimating Relese Fractions using AFEAS data and Ashford release fractions

clear all
HomeDir = '/Users/meganlickley/Dropbox (MIT)/Research/Lifetimes/Drafts/NatComms_Resubmission/Code';
str = strcat(HomeDir,'/Input/CFC11/AFEAS_cfc11production.mat');
load(str,'closecell','nonhermetic','open_aero','yr')

Yend = 2018;
nYears = length(yr(end)+1:Yend);
nSamps = 5*10^5;

% Production Distribution
% Uncertainty for aerosols and nonhermetic refrigeration
productionDist = 0.2*lognrnd(0,0.5,length(yr)+nYears,nSamps,2)+0.95;

% Uncertainty for closed cell foam
productionDist(:,:,3) = 0.1*lognrnd(0,0.5,length(yr)+nYears,nSamps,1)+0.95;

% Release Fraction for areosols.  All aerosols are emitted in year 1 or 2
DE11.aerosol = betarnd(2,2,1,nSamps);
RF11.aerosol = 1 - DE11.aerosol;

% Release Fraction for closed cell foam. 
m = 0.0366; v = (0.5*m)^2;
mu = log(m^2/sqrt(v+m^2));
sigma = sqrt(log(v/m^2+1));
RF11.closedcell = lognrnd(mu, sigma, 1, nSamps);
DE11.closedcell = 0.03*betarnd(5, 5, 1, nSamps);

% Direct Emissions for nonhermetic refrigeration
m = 0.07; v = (0.5*m)^2;
mu = log(m^2/sqrt(v+m^2));
sigma = sqrt(log(v/m^2+1));
DE11.nonhermetic = lognrnd(mu,sigma,1,nSamps);

% Lifetime for nonhermic refrigrations
m = 10; v = (0.20*m)^2;
mu = log(m^2/sqrt(v+m^2));
sigma = sqrt(log(v/m^2+1));
RF11.nonhermeticLT = lognrnd(mu,sigma,1,nSamps);
bank(1,:) = zeros(1,nSamps);

% Year 1 production distribution
nonhermetic1 = repmat(nonhermetic(1), 1, nSamps).*productionDist(1, :, 1);
open_aero1 = repmat(open_aero(1), 1, nSamps).*productionDist(1, :, 2);
closecell1 = repmat(closecell(1), 1, nSamps).*productionDist(1, :, 3);

annualProduction(1,:) = closecell1 + nonhermetic1 + open_aero1;

annualDE(1,:) = DE11.nonhermetic.*nonhermetic1 + ...
    DE11.aerosol.*open_aero1 + DE11.closedcell.*closecell1;

nonhermeticBank(1, :) = (1 - DE11.nonhermetic).*nonhermetic1;
closedcellBank(1, :)  = (1 - DE11.closedcell).*closecell1;

bankIn(1,:) = annualProduction(1,:) - annualDE(1,:);

for y = 2: length(yr)
    
    nonhermeticAnnualProd(y,:) = repmat(nonhermetic(y) - nonhermetic(y-1), 1, nSamps).*productionDist(y, :, 1);
    opencellAeroAnnualProd(y,:) = repmat(open_aero(y) - open_aero(y-1), 1, nSamps).*productionDist(y, :, 2);
    closedcellAnnualProd(y,:) = repmat(closecell(y) - closecell(y-1), 1, nSamps).*productionDist(y, :, 3);

    annualDE(y,:) = DE11.nonhermetic.*nonhermeticAnnualProd(y,:) + ...
        DE11.aerosol.*opencellAeroAnnualProd(y,:) + ...
        DE11.closedcell.*closedcellAnnualProd(y,:);
    annualProduction(y,:) = nonhermeticAnnualProd(y,:) + ...
        opencellAeroAnnualProd(y,:) + closedcellAnnualProd(y,:);

    bankIn(y,:) = annualProduction(y,:) - annualDE(y,:);
    
    nonhermeticBank(y,:) = exp(-1./RF11.nonhermeticLT).*nonhermeticBank(y-1, :) + ...
        (1 - DE11.nonhermetic).*nonhermeticAnnualProd(y,:);
    
    closedcellBank(y,:) = (1 - RF11.closedcell).*closedcellBank(y-1,:) + ...
        (1 - DE11.closedcell).*closedcellAnnualProd(y,:);
    
    bankEmissions(y,:) = opencellAeroAnnualProd(y-1, :).*RF11.aerosol + ...
        nonhermeticBank(y-1, :).*(1 - exp(-1./RF11.nonhermeticLT)) + ...
        closedcellBank(y-1, :).*RF11.closedcell;
    
    bank(y,:) = bank(y-1,:) + bankIn(y,:) - bankEmissions(y,:);
end

% For years after AFEAS production ends, assume production distribution
% equal to last year of AFEAS data

nonhermeticProd = nonhermetic(end)-nonhermetic(end-1);
openaeroProd = open_aero(end)-open_aero(end-1);
closedcellProd = closecell(end)-closecell(end-1);

for y = length(yr)+1:length(yr)+nYears
    
    nonhermeticAnnualProd(y,:) = repmat(nonhermeticProd, 1, nSamps).*productionDist(y, :, 1);
    opencellAeroAnnualProd(y,:) = repmat(openaeroProd, 1, nSamps).*productionDist(y, :, 2);
    closedcellAnnualProd(y,:) = repmat(closedcellProd, 1, nSamps).*productionDist(y, :, 3);

    annualDE(y,:) = DE11.nonhermetic.*nonhermeticAnnualProd(y,:) + DE11.aerosol.*opencellAeroAnnualProd(y, :) + ...
        DE11.closedcell.*closedcellAnnualProd(y,:);
    
    annualProduction(y,:) = nonhermeticAnnualProd(y, :) + opencellAeroAnnualProd(y, :) + closedcellAnnualProd(y, :);
 
    bankIn(y, :) = annualProduction(y,:) - annualDE(y,:);
    
    nonhermeticBank(y, :) = exp(-1./RF11.nonhermeticLT).*nonhermeticBank(y-1, :) + ...
        (1 - DE11.nonhermetic).*nonhermeticAnnualProd(y, :);
    
    closedcellBank(y, :) = (1 - RF11.closedcell).*closedcellBank(y-1, :) + ...
        (1 - DE11.closedcell).*closedcellAnnualProd(y, :);
    
    bankEmissions(y, :) = opencellAeroAnnualProd(y-1, :).*RF11.aerosol + ...
        nonhermeticBank(y-1, :).*(1 - exp(-1./RF11.nonhermeticLT)) + ...
        closedcellBank(y-1, :).*RF11.closedcell;
    
    bank(y, :) = bank(y-1, :) + bankIn(y, :) - bankEmissions(y, :);
end

RF = bankEmissions(2:end, :)./bank(1:end-1, :);
RF = [NaN(1,nSamps); RF];
DE = annualDE./annualProduction;

yearRFDE = yr(1):Yend;
str = strcat(HomeDir,'/Input/CFC11/InferredRFandDE.mat');

save(str,'RF','yearRFDE', 'DE')

%%  Expected emission scenario - Assuming unreported Production
% Re-evaluate production from 2000 onwards. 

% partition illicit production into equipment types with equal probability
tmp = rand(nSamps,3);
nonhermeticUnreported = tmp(:,1)./sum(tmp,2);
openaeroUnreported = tmp(:,2)./sum(tmp,2);
closedcellUnreported = tmp(:,3)./sum(tmp,2);

Prod_2012 = 71; % From Montzka et al. 2018 estimate
FugitiveADD = linspace(0,Prod_2012,13);
FugitiveADD(end+1:18) = Prod_2012;

newProductionSamps = rand(nSamps,1);

nonhermeticProd = nonhermetic(end)-nonhermetic(end-1);
openaeroProd = open_aero(end)-open_aero(end-1);
closedcellProd = closecell(end)-closecell(end-1);

for y = length(yr)+1:length(yr)+nYears
    
    nonhermeticAnnualProd(y,:) = repmat(nonhermeticProd, 1, nSamps).*productionDist(y, :, 1) + ...
        FugitiveADD(y - length(yr))*nonhermeticUnreported'.*newProductionSamps';
    opencellAeroAnnualProd(y,:) = repmat(openaeroProd, 1, nSamps).*productionDist(y,:,2) + ...
        FugitiveADD(y - length(yr))*openaeroUnreported'.*newProductionSamps';
    closedcellAnnualProd(y,:) = repmat(closedcellProd, 1, nSamps).*productionDist(y,:,3) + ...
        FugitiveADD(y - length(yr))*closedcellUnreported'.*newProductionSamps';

    annualDE(y,:) = DE11.nonhermetic.*nonhermeticAnnualProd(y,:) + ...
        DE11.aerosol.*opencellAeroAnnualProd(y, :) + ...
        DE11.closedcell.*closedcellAnnualProd(y, :);
    
    annualProduction(y,:) = nonhermeticAnnualProd(y, :) + ...
        opencellAeroAnnualProd(y, :) + closedcellAnnualProd(y, :);

    bankIn(y,:) = annualProduction(y, :) - annualDE(y, :);
    
    nonhermeticBank(y,:) = exp(-1./RF11.nonhermeticLT).*nonhermeticBank(y-1, :) + ...
        (1 - DE11.nonhermetic).*nonhermeticAnnualProd(y,:);
    
    closedcellBank(y,:) = (1 - RF11.closedcell).*closedcellBank(y-1,:) + ...
        (1-DE11.closedcell).*closedcellAnnualProd(y,:);
    
    bankEmissions(y,:) = opencellAeroAnnualProd(y-1,:).*RF11.aerosol + ...
        nonhermeticBank(y-1, :).*(1 - exp(-1./RF11.nonhermeticLT)) + ...
        closedcellBank(y-1, :).*RF11.closedcell;
    
    bank(y,:) = bank(y-1,:) + bankIn(y,:) - bankEmissions(y,:);
end

RF = bankEmissions(2:end,:)./bank(1:end-1,:);
RF = [NaN(1,nSamps); RF];
DE = annualDE./annualProduction;

yearRFDE = yr(1):Yend;
str = strcat(HomeDir,'/Input/CFC11/InferredRFandDE_unreportedProd.mat');

save(str,'RF','yearRFDE', 'DE')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CFC-12 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars -except HomeDir  nSamps Yend

str = strcat(HomeDir,'/Input/CFC12/AFEAS_cfc12production.mat');
load(str,'closecell','nonhermetic','open_aero','refrigeration','yr')

nYears = length(yr(end)+1 : Yend);

productionDist = 0.2*lognrnd(0, 0.5, length(yr)+nYears, nSamps, 4) + 0.95;

% model aerosol and opencell foam together bc of AFEAS grouping
DE12.aerosol = betarnd(2, 2, 1, nSamps);
RF12.aerosol = 1 - DE12.aerosol;

DE12.closedcell = betarnd(2, 2, 1, nSamps);
RF12.closedcell = 1 - DE12.closedcell;

m = 0.07; v = (0.5*m)^2;
mu = log(m^2/sqrt(v+m^2));
sigma = sqrt(log(v/m^2+1));
DE12.nonhermetic = lognrnd(mu,sigma,1,nSamps);%0.07*DEu(5,:);

m = 10; v = (0.20*m)^2;
mu = log(m^2/sqrt(v+m^2));
sigma = sqrt(log(v/m^2+1));
RF12.nonhermeticLT = lognrnd(mu,sigma,1,nSamps);

m = 0.02; v = (0.5*m)^2;
mu = log(m^2/sqrt(v+m^2));
sigma = sqrt(log(v/m^2+1));
DE12.hermetic = lognrnd(mu,sigma,1,nSamps); %0.02*DEu(7,:);

m = 20; v = (0.20*m)^2;
mu = log(m^2/sqrt(v+m^2));
sigma = sqrt(log(v/m^2+1));
RF12.hermeticLT = lognrnd(mu,sigma,1,nSamps);

closecell1 = repmat(closecell(1),1,nSamps).*productionDist(1,:,1);
nonhermetic1 = repmat(nonhermetic(1),1,nSamps).*productionDist(1,:,2);
open_aero1 = repmat(open_aero(1),1,nSamps).*productionDist(1,:,3);
refrigeration1 =  repmat(refrigeration(1),1,nSamps).*productionDist(1,:,4);

annualProduction(1,:) = closecell1 + nonhermetic1 + open_aero1 + refrigeration1;

annualDE(1,:)         = DE12.nonhermetic.*nonhermetic1 + ...
    DE12.aerosol.*open_aero1 + DE12.hermetic.*refrigeration1 + ...
    DE12.closedcell.*closecell1;

nonhermeticBank(1, :) = (1 - DE12.nonhermetic).*nonhermetic1;
hermeticBank(1, :)    = (1 - DE12.hermetic).*refrigeration1;

bankIn(1,:) = annualProduction(1,:)-annualDE(1,:);
bank(1,:) = bankIn(1,:);
for y = 2: length(yr)
    nonhermeticAnnualProd(y,:) = repmat((nonhermetic(y)-nonhermetic(y-1)),1,nSamps).*productionDist(y,:,1);
    hermeticRefrigeration(y,:) = repmat((refrigeration(y) - refrigeration(y-1)),1,nSamps).*productionDist(y,:,2);
    opencellAeroAnnualProd(y,:) = repmat((open_aero(y)-open_aero(y-1)),1,nSamps).*productionDist(y,:,3);
    closedcellAnnualProd(y,:) = repmat((closecell(y)-closecell(y-1)),1,nSamps).*productionDist(y,:,4);
    
    annualDE(y,:) = DE12.nonhermetic.*nonhermeticAnnualProd(y,:) + ...
        DE12.aerosol.*opencellAeroAnnualProd(y,:) + ...
        DE12.hermetic.*hermeticRefrigeration(y,:) + ...
        DE12.closedcell.*closedcellAnnualProd(y,:);
    
    annualProduction(y, :) = nonhermeticAnnualProd(y, :) + ...
        opencellAeroAnnualProd(y, :) + closedcellAnnualProd(y, :) + ...
        hermeticRefrigeration(y, :);
    
    bankIn(y,:) = annualProduction(y,:) - annualDE(y,:);
    
    nonhermeticBank(y,:) = exp(-1./RF12.nonhermeticLT).*nonhermeticBank(y-1,:) + ...
        (1 - DE12.nonhermetic).*nonhermeticAnnualProd(y,:);
    
    hermeticBank(y,:) = exp(-1./RF12.hermeticLT).*hermeticBank(y-1,:) + ...
        (1-DE12.hermetic).*hermeticRefrigeration(y,:);
    
    bankEmissions(y,:) = opencellAeroAnnualProd(y-1,:).*RF12.aerosol + ...
        closedcellAnnualProd(y-1, :).*RF12.closedcell + ...
        nonhermeticBank(y-1, :).*(1 - exp(-1./RF12.nonhermeticLT)) + ...
        hermeticBank(y-1, :).*(1 - exp(-1./RF12.hermeticLT));
    
    bank(y,:) = bank(y-1,:) + bankIn(y,:) - bankEmissions(y,:);
end

nonhermeticProd   = nonhermetic(end) - nonhermetic(end-1);
refrigerationProd = refrigeration(end) - refrigeration(end-1);
openaeroProd      = open_aero(end) - open_aero(end-1);
closedcellProd    = closecell(end) - closecell(end-1);

for y = length(yr)+1 : length(yr)+nYears
    
    nonhermeticAnnualProd(y, :)  = repmat(nonhermeticProd, 1, nSamps).*productionDist(y,:,1);
    hermeticRefrigeration(y, :)  = repmat(refrigerationProd, 1, nSamps).*productionDist(y,:,2);
    opencellAeroAnnualProd(y, :) = repmat(openaeroProd, 1, nSamps).*productionDist(y,:,3);
    closedcellAnnualProd(y, :)   = repmat(closedcellProd, 1, nSamps).*productionDist(y,:,4);
    
    annualDE(y,:) = DE12.nonhermetic.*nonhermeticAnnualProd(y, :) + ...
        DE12.aerosol.*opencellAeroAnnualProd(y, :) + ...
        DE12.hermetic.*hermeticRefrigeration(y, :) + ...
        DE12.closedcell.*closedcellAnnualProd(y, :);
    
    annualProduction(y,:) = nonhermeticAnnualProd(y, :) + ...
        opencellAeroAnnualProd(y,:) + closedcellAnnualProd(y,:) + ...
        hermeticRefrigeration(y,:);
    
    bankIn(y,:) = annualProduction(y,:) - annualDE(y,:); 
    
    nonhermeticBank(y,:) = exp(-1./RF12.nonhermeticLT).*nonhermeticBank(y-1, :) + ...
        (1 - DE12.nonhermetic).*nonhermeticAnnualProd(y, :);
    hermeticBank(y,:) = exp(-1./RF12.hermeticLT).*hermeticBank(y-1, :) + ...
        (1 - DE12.hermetic).*hermeticRefrigeration(y, :);
    
    bankEmissions(y, :) = opencellAeroAnnualProd(y-1, :).*RF12.aerosol + ...
        closedcellAnnualProd(y-1, :).*RF12.closedcell + ...        
        nonhermeticBank(y-1, :).*(1 - exp(-1./RF12.nonhermeticLT)) + ...
        hermeticBank(y-1, :).*(1 - exp(-1./RF12.hermeticLT));
    
    bank(y,:) = bank(y-1,:) + bankIn(y,:) - bankEmissions(y,:);
end

RF = bankEmissions(2:end,:)./bank(1:end-1,:);
RF = [NaN(1,nSamps);RF];

DE = annualDE./annualProduction;

yearRFDE = yr(1):Yend;
str = strcat(HomeDir,'/Input/CFC12/InferredRFandDE.mat');

save(str,'RF','yearRFDE', 'DE')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CFC-113 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars -except HomeDir nSamps Yend 

str = strcat(HomeDir,'/Input/CFC113/AFEAS_cfc113production.mat');
load(str,'longbank','shortbank','yr')

nYears = length(yr(end)+1:Yend);

productionDist = 0.27*lognrnd(0, 0.5, length(yr)+nYears, nSamps, 4) + 0.7;

DE113.shortbank = betarnd(2,2,1,nSamps);
RF113.shortbank = 1 - DE113.shortbank;

m = 0.02; v = (m)^2;
mu = log(m^2/sqrt(v+m^2));
sigma = sqrt(log(v/m^2+1));
DE113.longbank = lognrnd(mu,sigma,1,nSamps); 

m = 20; v = (0.3*m)^2;
mu = log(m^2/sqrt(v+m^2));
sigma = sqrt(log(v/m^2+1));
RF113.longbankLT = lognrnd(mu,sigma,1,nSamps); 

shortbank1 = repmat(shortbank(1), 1, nSamps).*productionDist(1,:,1);
longbank1 = repmat(longbank(1), 1, nSamps).*productionDist(1,:,2);

annualProduction(1,:) = shortbank1 + longbank1;
annualDE(1,:) = DE113.shortbank.*shortbank1 + DE113.longbank.*longbank1;

longbankBank(1,:) = (1 - DE113.longbank).*longbank1;
bankIn(1,:) = annualProduction(1,:) - annualDE(1,:);
bank(1,:) = bankIn(1,:);

for y = 2: length(yr)
    
    shortbankAnnualProd(y, :) = repmat(shortbank(y) - shortbank(y-1),1,nSamps).*productionDist(y,:,1);
    longbankAnnualProd(y, :) = repmat(longbank(y) - longbank(y-1),1,nSamps).*productionDist(y,:,2);
    
    annualDE(y, :) = DE113.shortbank.*shortbankAnnualProd(y,:) + ...
        DE113.longbank.*longbankAnnualProd(y,:);
    
    annualProduction(y, :) = shortbankAnnualProd(y,:) + longbankAnnualProd(y,:);
    
    bankIn(y, :) = annualProduction(y,:) - annualDE(y,:);
    
    longbankBank(y, :) = exp(-1./RF113.longbankLT).*longbankBank(y-1,:) + ...
        (1 - DE113.longbank).*longbankAnnualProd(y,:);
    
    bankEmissions(y, :) = RF113.shortbank.*shortbankAnnualProd(y-1,:) + ...
        longbankBank(y-1, :).*(1 - exp(-1./RF113.longbankLT));
    
    bank(y,:) = bank(y-1,:) + bankIn(y,:) - bankEmissions(y,:);
end

shortbankProd = shortbank(end) - shortbank(end-1);
longbankProd = longbank(end) - longbank(end-1);

for y = length(yr)+1: length(yr)+nYears
    
    shortbankAnnualProd(y,:) = repmat(shortbankProd, 1, nSamps).*productionDist(y, :, 1);
    longbankAnnualProd(y,:) = repmat(longbankProd, 1, nSamps).*productionDist(y, :, 2);
    
    annualDE(y,:) = DE113.shortbank.*shortbankAnnualProd(y, :) + ...
        DE113.longbank.*longbankAnnualProd(y, :);
    annualProduction(y,:) = shortbankAnnualProd(y, :) + longbankAnnualProd(y, :);
    
    bankIn(y,:) = annualProduction(y, :) - annualDE(y, :);
    
    longbankBank(y,:) = exp(-1./RF113.longbankLT).*longbankBank(y-1, :) + ...
        (1 - DE113.longbank).*longbankAnnualProd(y,:);
    
    bankEmissions(y,:) = RF113.shortbank.*shortbankAnnualProd(y-1,:) + ...
        longbankBank(y-1,:).*(1-exp(-1./RF113.longbankLT));
    bank(y,:) = bank(y-1,:) + bankIn(y,:) - bankEmissions(y,:);
end

RF = bankEmissions(2:end,:)./bank(1:end-1,:);
RF = [NaN(1,nSamps);RF];

DE = annualDE./annualProduction;

yearRFDE = yr(1):Yend;
str = strcat(HomeDir,'/Input/CFC113/InferredRFandDE.mat');
save(str,'RF','yearRFDE', 'DE')

%%  Fugitive  CFC-113, using upper bound based on Lickley et al. 2020
% Assume all of the emissions come from a short bank for feedstock use. 

Fugitive_prod = 15; %Gg/yr as upper limit 
Fugitive_add = Fugitive_prod*betarnd(2,2,nSamps,1);

shortbankProd = shortbank(end) - shortbank(end-1);
longbankProd = longbank(end) - longbank(end-1);
for y = length(yr)+1: length(yr)+nYears
    
    shortbankAnnualProd(y,:) = repmat(shortbankProd, 1, nSamps).*productionDist(y, :, 1) + Fugitive_add';
    longbankAnnualProd(y,:) = repmat(longbankProd, 1, nSamps).*productionDist(y,:,2);
    
    annualDE(y,:) = DE113.shortbank.*shortbankAnnualProd(y,:) + ...
        DE113.longbank.*longbankAnnualProd(y,:);
    annualProduction(y,:) = shortbankAnnualProd(y,:) + longbankAnnualProd(y,:);
    
    bankIn(y,:) = annualProduction(y,:)  - annualDE(y,:);
    
    longbankBank(y,:) = exp(-1./RF113.longbankLT).*longbankBank(y-1,:) + ...
        (1 - DE113.longbank).*longbankAnnualProd(y,:);
    
    bankEmissions(y,:) = RF113.shortbank.*shortbankAnnualProd(y-1,:) + ...
        longbankBank(y-1,:).*(1-exp(-1./RF113.longbankLT));
    bank(y,:) = bank(y-1,:) + bankIn(y,:) - bankEmissions(y,:);
end

RF = bankEmissions(2:end,:)./bank(1:end-1,:);
RF = [NaN(1,nSamps);RF];

annualDE = annualDE./annualProduction;

yearRFDE = yr(1):Yend;
str = strcat(HomeDir,'/Input/CFC113/InferredRFandDE_unreportedProd.mat');
save(str,'RF','yearRFDE', 'DE')


