% Main BPE Script for Lickley et al. (in Revisions) 
% Title: Joint inference of CFC lifetimes and banks suggests previously 
% unidentified emissions
%
% Authors: Megan Lickley, Sarah Fletcher, Matt Rigby and Susan Solomon
% Corresponding Author: Megan Lickley; mlickley@mit.edu
%
% Last udpated Dec 18, 2020 by Megan Lickley.
% This script creates the figures for the Main text and supplement

close all
clear all

HomeDir = '/Users/meganlickley/Dropbox (MIT)/Research/Lifetimes/Drafts/NatComms_Resubmission/Code/';
FolderName = strcat(HomeDir,'Output/');
FigureFolderName = strcat(HomeDir, 'FigureFolder/');

CFC11_FileName{1}  = strcat(FolderName, 'CFC11/Run1_unreportedProd_jointInference.mat');
CFC12_FileName{1}  = strcat(FolderName, 'CFC12/Run1_unreportedProd_jointInference.mat');
CFC113_FileName{1} = strcat(FolderName, 'CFC113/Run1_unreportedProd_jointInference.mat');

CFC11_FileName{2}   = strcat(FolderName, 'CFC11/Run1_jointInference.mat');
CFC12_FileName{2}   = strcat(FolderName, 'CFC12/Run1_jointInference.mat');
CFC113_FileName{2}  = strcat(FolderName, 'CFC113/Run1_jointInference.mat');

CFC11_FileName{3}   = strcat(FolderName, 'CFC11/Run1_unreportedProd_independent.mat');
CFC12_FileName{3}   = strcat(FolderName, 'CFC12/Run1_unreportedProd_independent.mat');
CFC113_FileName{3}  = strcat(FolderName, 'CFC113/Run1_unreportedProd_independent.mat');


%% Fig 1: Assumptions skewing our interpretation of results
% Indices 1,3,5,7 correspond to SPARC timevarying MMM lifetime
% Indices 2,4,6,8 correspond to constant lifetime of 52 yrs
% Indices 1,2,5,6 correspond to reported production
% Indices 3,4,7,8 correspond to  1.1 x reported production
% Indices 1,2,3,4 correspond to median RF value from prior distribution
% Indices 5,6,7,8 corresponds to 16th percentile of RF value 
y1 = 1955; % Year Start Date
yswitch = 1988;
y2 = 2018; % Year End Date
simulationyrs = y1:y2;

load('Input/CFC11/wmo2018.mat');
ytmp1 = find(wmo_yr == 1955);
ytmp2 = find(wmo_yr == 2017);
MF = wmo_conc(ytmp1:ytmp2+1);

yearsSPARC = 1960:2010;
load('Input/SPARC_lifetimes.mat')

tmp = squeeze(LT_sparc(:,:,1)); 
InvLT = 1./tmp;

% Computing 10 year moving average
for ii = 1:42
    InvLTma(ii,:) = nanmean(InvLT(ii:ii+9,:));
end

InvLT = nanmean(InvLTma,2);
LTp(1:yearsSPARC(1)-y1+5) = InvLT(1); 
LTp(yearsSPARC(1)-y1+6:yearsSPARC(1)-y1+47) = InvLT; 
LTp(52:64) = InvLT(end);

InitialBankSize = 5893.9; %From WMO 2002
A = 22602.38457; % Conversion factor of mole fractions to emissions
% Top-down emissions are either from time varying lifetime or constant
% lifetime
topdownEmissions = A*(MF(2:end)'-MF(1:end-1)'.*exp(-LTp(1:end-1)));
topdownEmissions(2,:) = A*(MF(2:end)'-MF(1:end-1)'.*exp(-1./(52*ones(size(LTp(1:end-1))))));

reportedProduction = NaN(size(simulationyrs)); % Production is in Tonnes
% Use until 1988 - WMO data is the adjusted AFEAS data
load('Input/CFC11/WMO2002.mat') 
ytmp1 = find(year==y1); 
ytmp2 = find(year==yswitch);
reportedProduction(1:yswitch+1-y1) = Production(ytmp1:ytmp2);

% Use Article 5 and non A5 prod total from 1989 onwards
load('Input/CFC11/a5_na5.mat') 
ytmp1 = find(a5_na5(:,1) == yswitch+1);
ytmp2 = find(simulationyrs == a5_na5(end,1));
reportedProduction(1989-y1+1:ytmp2) = sum(a5_na5(ytmp1:end,[2:3]),2);
reportedProduction(ytmp2:end) = sum(a5_na5(end,[2:3]),2);

% top-down bank estimates given different assumptions
clear Bank

for LT_ii = 1:2
    Bank(LT_ii,1) = InitialBankSize;
    for yy = 2:length(topdownEmissions)
        Bank(LT_ii,yy) = sum(reportedProduction(1:yy)) - ...
            sum(topdownEmissions(LT_ii,1:yy)) + InitialBankSize;
    end
end

for LT_ii = 1:2
    for yy = 2:length(topdownEmissions)
        Bank(LT_ii+2,yy) = sum(1.1*reportedProduction(1:yy)) - ...
            sum(topdownEmissions(LT_ii,1:yy)) + InitialBankSize;
    end
end

load(CFC11_FileName{1},'RFsamps11');
RF_opts = median(RFsamps11);
RF_opts(2,:) = prctile(RFsamps11,16);
RF_opts(:,57:64) = repmat(RF_opts(:,56),1,8);

bankEmiss(1:4,:) = repmat(RF_opts(1,2:end),4,1).*Bank;
bankEmiss(5:8,:) = repmat(RF_opts(2,2:end),4,1).*Bank;

DEmiss = repmat(topdownEmissions,4,1) - bankEmiss;


FigHandle = figure(1)
subplot(2,2,1); 
yrs = y1:y2-1;
plot(yrs, 0.001*topdownEmissions(1,:),   'k', 'LineWidth', 1);   hold on; 
plot(yrs, 0.001*topdownEmissions(2,:), '--k', 'LineWidth', 1); 
xlim([1955,2010]); ylabel('[Gg/yr]'); set(gca, 'FontSize', 12);
title('CFC-11 Emissions');

subplot(2,2,2); 
plot(yrs, 0.001*Bank(1,:),   'k',  'LineWidth', 3);    hold on; 
plot(yrs, 0.001*Bank(2,:), '--k',  'LineWidth', 3);    hold on; 
plot(yrs, 0.001*Bank(3,:),   'k',  'LineWidth', 1);    hold on; 
plot(yrs, 0.001*Bank(4,:),  '--k', 'LineWidth', 1); 
xlim([1955,2010]); ylabel('[Gg]'); set(gca, 'FontSize', 12);
title('CFC-11 Banks'); 

subplot(2,2,3); 
colRF1 = [0, 0.447, 0.741];
colRF2 = [0.85, 0.325, 0.098];

plot(yrs, 0.001*bankEmiss(1,:),     'Color',colRF1,'LineWidth', 3); hold on; 
plot(yrs, 0.001*bankEmiss(2,:),'--','Color',colRF1,'LineWidth', 3); hold on; 
plot(yrs, 0.001*bankEmiss(3,:),     'Color',colRF1,'LineWidth', 1); hold on;
plot(yrs, 0.001*bankEmiss(4,:),'--','Color',colRF1,'LineWidth', 1); hold on; 
plot(yrs, 0.001*bankEmiss(5,:),     'Color',colRF2,'LineWidth', 3); hold on; 
plot(yrs, 0.001*bankEmiss(6,:),'--','Color',colRF2,'LineWidth', 3); hold on; 
plot(yrs, 0.001*bankEmiss(7,:),     'Color',colRF2,'LineWidth', 1); hold on; 
plot(yrs, 0.001*bankEmiss(8,:),'--','Color',colRF2,'LineWidth', 1);
xlim([1955,2010]); ylabel(' [Gg/yr]'); set(gca, 'FontSize', 12);
title('Bank Emissions'); 

subplot(2,2,4); 
plot(yrs, 0.001*DEmiss(1,:),     'Color',colRF1,'LineWidth', 3); hold on; 
plot(yrs, 0.001*DEmiss(2,:),'--','Color',colRF1,'LineWidth', 3); hold on; 
plot(yrs, 0.001*DEmiss(3,:),     'Color',colRF1,'LineWidth', 1); hold on;
plot(yrs, 0.001*DEmiss(4,:),'--','Color',colRF1,'LineWidth', 1); hold on; 
plot(yrs, 0.001*DEmiss(5,:),     'Color',colRF2,'LineWidth', 3); hold on; 
plot(yrs, 0.001*DEmiss(6,:),'--','Color',colRF2,'LineWidth', 3); hold on; 
plot(yrs, 0.001*DEmiss(7,:),     'Color',colRF2,'LineWidth', 1); hold on; 
plot(yrs, 0.001*DEmiss(8,:),'--','Color',colRF2,'LineWidth', 1);
xlim([1955,2010]);  ylabel(' [Gg/yr]'); set(gca, 'FontSize', 12);
title('Total Direct Emissions'); 

figure_width = 14; % in inches
figure_height = 8; % in inches
screen_ppi = 72; % inherent property of Matlab. CHANGING THIS WILL MAKE FIGURE LOOK INCONSISTENT BETWEEN SCREEN AND EXPORT!
set(FigHandle, 'Position', [100, 100, round(figure_width*screen_ppi), round(figure_height*screen_ppi)]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [figure_width figure_height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 figure_width figure_height]);

str = strcat(FigureFolderName,'Figure1.pdf');
print(gcf, '-dpdf', str);


%%  Fig 2: Mole Fraction prior, observations and posterior
clearvars -except CFC11_FileName CFC12_FileName CFC113_FileName ...
    FigureFolderName FigureFolderName HomeDir

load(CFC11_FileName{1}, 'MF11', 'ResampleIndex1','molefractions');
load(CFC12_FileName{1}, 'MF12', 'ResampleIndex2');
load(CFC113_FileName{1},'MF113','ResampleIndex3');

y1 = 1955; % Year Start Date
y2 = 2010; % Year End Date
nYears = length(y1:y2);

FigHandle = figure(2); 

for mol_ii = 1:3

    if mol_ii == 1;     var = MF11;  indx = ResampleIndex1; str = 'CFC11';
    elseif mol_ii == 2; var = MF12;  indx = ResampleIndex2; str = 'CFC12';
    elseif mol_ii == 3; var = MF113; indx = ResampleIndex3; str = 'CFC113';
    end
    
    subplot(1,3,mol_ii); 
    MED = prctile(var,50); LB = prctile(var,2.5); UB = prctile(var,97.5);
    p1 = boundedline([y1:y2+1], MED',[MED'-LB',UB'-MED'],'alpha','cmap',[0.3,0.3,0.3]); 
    
    hold on; 
    var1 = var(indx,:); 
    MED = prctile(var1,50); LB = prctile(var1,2.5); UB = prctile(var1,97.5);
    p3 = boundedline([y1:y2+1],MED',[MED'-LB',UB'-MED'],'alpha','cmap',[0.7,0,0]);  
    
    hold on; 
    p2 = plot(y1:y2+1,molefractions{mol_ii}(1:end),'--b','LineWidth',2);
    
    ylabel('Mole Fractions [pmol mol^{-1}]'); 
    xlabel('Year'); xlim([y1,2010]);
    title(str); box on; set(gca, 'FontSize', 12); 
    lgd = legend([p1 p2 p3],'prior','observed','posterior')
    lgd.Location = 'northwest';
end

figure_width = 18; % in inches
figure_height = 5; % in inches
screen_ppi = 72; % inherent property of Matlab. CHANGING THIS WILL MAKE FIGURE LOOK INCONSISTENT BETWEEN SCREEN AND EXPORT!
set(FigHandle, 'Position', [100, 100, round(figure_width*screen_ppi), round(figure_height*screen_ppi)]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [figure_width figure_height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 figure_width figure_height]);

str = strcat(FigureFolderName,'Figure2.pdf');
print(gcf, '-dpdf', str);

%%  Fig 3 (top): Emissions priors for top-down and bottom up.  bottom (posteriors)

clearvars -except CFC11_FileName CFC12_FileName CFC113_FileName ...
    FigureFolderName FigureFolderName HomeDir

load(CFC11_FileName{1}, 'Emiss11',  'ResampleIndex1', 'ObsEmiss11', 'ResampleIndex_LT');
load(CFC12_FileName{1}, 'Emiss12',  'ResampleIndex2', 'ObsEmiss12');
load(CFC113_FileName{1},'Emiss113', 'ResampleIndex3', 'ObsEmiss113');

load('Figure_Data/WMO2018Emissions')
y1 = 1955; % Year Start Date
y2 = 2010; % Year End Date
yrs = y1:y2;
nYears = length(yrs);

FigHandle = figure(3);

for mol_ii = 1:3
    if mol_ii == 1; 
        var1 = Emiss11;  var2 = ObsEmiss11;  ind = ResampleIndex1; ymax = 450;
    elseif mol_ii == 2; 
        var1 = Emiss12;  var2 = ObsEmiss12;  ind = ResampleIndex2; ymax = 600;
    elseif mol_ii == 3; 
        var1 = Emiss113; var2 = ObsEmiss113; ind = ResampleIndex3; ymax = 400;
    end
    
    wmoyrs = WMO2018Emissions(:,1);
    wmoObs = WMO2018Emissions(:,1+mol_ii);
    
    
    subplot(2,3,mol_ii)
    var1 = 0.001*var1;  var2 = 0.001*var2;
    MED = prctile(var1,50); 
    LB = prctile(var1,2.5); UB = prctile(var1,97.5);
    p1 = boundedline(yrs, MED',[MED'-LB',UB'-MED'],'alpha','cmap',[0.3,0.3,0.3]);
    
    hold on;
    MED = prctile(var2,50); 
    LB = prctile(var2,2.5); UB = prctile(var2,97.5);
    p2 = boundedline(yrs, MED',[MED'-LB',UB'-MED'],'alpha','cmap',[0.7,0,0]);
   
    hold on;
    p3 = plot(wmoyrs, 0.001*wmoObs,'--b','LineWidth',2);
    xlim([y1,2010]); xlabel('Year'); ylabel('[Gg yr^{-1}]'); 
    title('CFC 11 Prior Emissions'); 
    box on; set(gca, 'FontSize', 12); ylim([0,ymax]);
    
    
    subplot(2,3,mol_ii+3)
    var1 = var1(ind,:);  var2 = var2(ResampleIndex_LT,:);
    MED = prctile(var1, 50); 
    LB = prctile(var1, 2.5); UB = prctile(var1, 97.5);
    p1 = boundedline(yrs, MED',[MED'-LB',UB'-MED'],'alpha','cmap',[0.3,0.3,0.3]);
    
    hold on;
    MED = prctile(var2, 50); 
    LB =  prctile(var2, 2.5); UB = prctile(var2, 97.5);
    p2 = boundedline([y1:y2],MED',[MED'-LB',UB'-MED'],'alpha','cmap',[0.7,0,0]);
    
    hold on;
    p3 = plot(wmoyrs, 0.001*wmoObs, '--b','LineWidth',2);
    ylabel('[Gg yr^{-1}]'); xlabel('Year'); xlim([1955,2010]); ylim([0,ymax]);
    title('CFC 11 Posterior Emissions'); 
    set(gca, 'FontSize', 12); box on; 
end

figure_width = 14; % in inches
figure_height = 10; % in inches
screen_ppi = 72; % inherent property of Matlab. CHANGING THIS WILL MAKE FIGURE LOOK INCONSISTENT BETWEEN SCREEN AND EXPORT!
set(FigHandle, 'Position', [100, 100, round(figure_width*screen_ppi), round(figure_height*screen_ppi)]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [figure_width figure_height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 figure_width figure_height]);

str = strcat(FigureFolderName,'Figure3.pdf');
print(gcf, '-dpdf', str);

%% Fig 4: Lifetime prior and posteriors

clearvars -except CFC11_FileName CFC12_FileName CFC113_FileName ...
    FigureFolderName FigureFolderName HomeDir

% scenario = 1 : Figure 4 and table 1 in main text
% scenario = 2 : table S1; 
% scenario = 3 : Figure S2 and table S2

scenario = 1;
load(CFC11_FileName{scenario},'LT11','ResampleIndex1','ResampleIndex_LT');
load(CFC12_FileName{scenario},'LT12','ResampleIndex2');
load(CFC113_FileName{scenario},'LT113','ResampleIndex3');

y1 = 1955; % Year Start Date
y2 = 2010; % Year End Date
yrs = y1:y2;
nYears = length(yrs);

yearsSPARC = 1960:2010;
load('Input/SPARC_lifetimes.mat')

% Computing 10 year moving average
for mol = 1:3
    LT = squeeze(LT_sparc(:,:,mol)); 
    
    for ii = 1:42
        LTma{mol}(ii,:) = mean(LT(ii:ii+9,:));
    end
    for mod_ii = 1:7
        indx = isnan(LTma{mol}(:,mod_ii));
        LTma{mol}(indx,mod_ii) = nanmean(LT(sum(1-indx):end,mod_ii));
    end
    
    LTp{mol}(1:yearsSPARC(1)-y1+5) = mean(LTma{mol}(1,:)); 
    LTp{mol}(yearsSPARC(1)-y1+6:yearsSPARC(1)-y1+47) = mean(LTma{mol},2); 
    LTp{mol}(52:nYears) = mean(LTma{mol}(end,:));
end

FigHandle = figure(4); 

for mol_ii = 1:3
    if scenario == 1 || scenario == 2
        FigName = 'Figure4.pdf';
        if mol_ii == 1;   var1 = LT11;  str = strcat('CFC-11 Lifetimes'); 
            indx = ResampleIndex_LT; 
            indx1 = indx; indx2 = indx; indx3 = indx;
        elseif mol_ii == 2; var1 = LT12;  str = strcat('CFC-12 Lifetimes');
        elseif mol_ii == 3; var1 = LT113; str = strcat('CFC-113 Lifetimes');
        end
    else
        FigName = 'FigureS2.pdf';
        if mol_ii == 1;     var1 = LT11;  str = strcat('CFC-11 Lifetimes');
           indx = ResampleIndex1;  indx1 = indx;
        elseif mol_ii == 2; var1 = LT12;  str = strcat('CFC-12 Lifetimes');
           indx = ResampleIndex2;  indx2 = indx;
        elseif mol_ii == 3; var1 = LT113; str = strcat('CFC-113 Lifetimes');
           indx = ResampleIndex3;  indx3 = indx;
        end
    end
    
    subplot(1,3,mol_ii)

    var2 = var1(indx,:);

    MED = prctile(var1,50); med_tmp = 1./MED(1:nYears)';
    LB = prctile(var1,2.5); UB = prctile(var1,97.5);
    lb_tmp = 1./MED(1:nYears)'-1./LB(1:nYears)';
    ub_tmp = 1./UB(1:nYears)'-1./MED(1:nYears)';
    p1 = boundedline(yrs, med_tmp, [lb_tmp, ub_tmp],'alpha','cmap',[0.3,0.3,0.3]);
    hold on; plot(yrs,med_tmp,'Color',[0.3,0.3,0.3],'LineWidth',2);
    
    hold on;
    MED = prctile(var2,50); med_tmp = 1./MED(1:nYears)';
    LB = prctile(var2,2.5); UB = prctile(var2,97.5);
    lb_tmp = 1./MED(1:nYears)'-1./LB(1:nYears)';
    ub_tmp = 1./UB(1:nYears)'-1./MED(1:nYears)';
    p2 = boundedline(yrs, med_tmp, [lb_tmp, ub_tmp],'alpha','cmap',[0.7,0,0]);
    hold on; plot(yrs, med_tmp,'Color',[0.7,0,0],'LineWidth',2);

    hold on; 
    p3 = plot(y1:y2,LTp{mol_ii},'--b','LineWidth',2);
    title(str); ylabel('[years]'); xlabel('Year');
    box on; set(gca, 'FontSize', 12); xlim([y1,2010]);
    xlim([yearsSPARC(5),yearsSPARC(end-5)]);
end

figure_width = 18; % in inches
figure_height = 4; % in inches
screen_ppi = 72; % inherent property of Matlab. CHANGING THIS WILL MAKE FIGURE LOOK INCONSISTENT BETWEEN SCREEN AND EXPORT!
set(FigHandle, 'Position', [100, 100, round(figure_width*screen_ppi), round(figure_height*screen_ppi)]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [figure_width figure_height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 figure_width figure_height]);
str = strcat(FigureFolderName,FigName);
print(gcf, '-dpdf', str);


tmp = LT11(indx1, end);
tmp11end = [1./median(tmp),  1./prctile(tmp,2.5),  1./prctile(tmp,97.5)];

tmp = mean(LT11(indx1, :)');
tmp11_mean = [1./median(tmp),  1./prctile(tmp,2.5),  1./prctile(tmp,97.5)];

tmp = LT12(indx2, end);
tmp12end = [1./median(tmp),  1./prctile(tmp,2.5),  1./prctile(tmp,97.5)];

tmp = mean(LT12(indx2, :)');
tmp12_mean = [1./median(tmp),  1./prctile(tmp,2.5),  1./prctile(tmp,97.5)];

tmp = LT113(indx3, end);
tmp113end = [1./median(tmp),  1./prctile(tmp,2.5),  1./prctile(tmp,97.5)];

tmp = mean(LT113(indx3, :)');
tmp113_mean = [1./median(tmp),  1./prctile(tmp,2.5),  1./prctile(tmp,97.5)];

Data_end = [11 tmp11end
            12 tmp12end
            113 tmp113end];
VarNames = {'CFC', 'med_2010LT','prctile2p5','prctile97p5'};
T = table(Data_end(:,1),Data_end(:,2),Data_end(:,3),Data_end(:,4), 'VariableNames',VarNames)

Data_all = [11 tmp11_mean
            12 tmp12_mean
            113 tmp113_mean];
VarNames = {'CFC', 'med_meanLT','prctile2p5','prctile97p5'};
T = table(Data_all(:,1),Data_all(:,2),Data_all(:,3),Data_all(:,4), 'VariableNames',VarNames)

%% Figure 5 Emissions projections going forward from 2011 onwards
clearvars -except CFC11_FileName CFC12_FileName CFC113_FileName ...
     FigureFolderName HomeDir

y1 = 1955; % Year Start Date
y2 = 2010; % Year End Date
yrs = y1:y2;
nYears = length(yrs);

% scenario = 1 : Figure 4 and table 1 in main text
% scenario = 2 : table S3; Figure S7

scenario = 1; 
load(CFC11_FileName{scenario},'LT11','ResampleIndex1','Bank11','RFsamps11','Emiss11','ObsEmiss11','ResampleIndex_LT');
load(CFC12_FileName{scenario},'LT12','ResampleIndex2','Bank12','RFsamps12','Emiss12','ObsEmiss12');
load(CFC113_FileName{scenario},'LT113','ResampleIndex3','Bank113','RFsamps113','Emiss113','ObsEmiss113');
A(1) = 22602.38457;
A(2) = 19895.36;
A(3) = 30834.68811;

load('Input/CFC11/wmo2018.mat');
ytmp1 = find(wmo_yr == 2010);
ytmp2 = find(wmo_yr == 2016);
obsMF{1} = wmo_conc(ytmp1:ytmp2+1);

load('Input/CFC12/wmo2018.mat');
ytmp1 = find(wmo_yr == 2010);
ytmp2 = find(wmo_yr == 2016);
obsMF{2} = wmo_conc(ytmp1:ytmp2+1);

load('Input/CFC113/wmo2018.mat');
ytmp1 = find(wmo_yr == 2010);
ytmp2 = find(wmo_yr == 2016);
obsMF{3} = wmo_conc(ytmp1:ytmp2+1);

FigHandle = figure(5); 

for mol_ii = 1:3
    
    if mol_ii == 1; 
        varLT = LT11(ResampleIndex_LT,end);     RFsamps = RFsamps11; 
        indx = ResampleIndex1;                  Banksamps = Bank11; 
        ObsEmissSamps = ObsEmiss11;             EmissSamps = Emiss11;
        ymax = 400;    str = 'CFC 11 Posterior Emissions';
    elseif mol_ii == 2; 
        varLT = LT12(ResampleIndex_LT,end);     RFsamps = RFsamps12; 
        indx = ResampleIndex2;                  Banksamps = Bank12; 
        ObsEmissSamps = ObsEmiss12;             EmissSamps = Emiss12;
        ymax = 600;    str = 'CFC 12 Posterior Emissions';
    elseif mol_ii == 3; 
        varLT = LT113(ResampleIndex_LT,end);    RFsamps = RFsamps113; 
        indx = ResampleIndex3;                  Banksamps = Bank113; 
        ObsEmissSamps = ObsEmiss113;            EmissSamps = Emiss113;
        ymax = 300;    str = 'CFC 12 Posterior Emissions';
    end
    
    clear Banktmp
    nSamps = size(varLT);
    BankEmiss{mol_ii} = Banksamps(indx,1:end-1).*RFsamps(indx,2:end);
    Banktmp(:,1) = Banksamps(indx,55);   

    topDownEmissProj{mol_ii} = A(mol_ii)*(repmat(obsMF{mol_ii}(2:end),1,nSamps) - ...
        repmat(obsMF{mol_ii}(1:end-1),1,nSamps).*exp(-repmat(varLT',7,1)));

    for ii = 56:62
        BankEmiss{mol_ii}(:,ii) = Banktmp(:,ii-55).*RFsamps(indx,end);
        Banktmp(:,ii-54) = Banktmp(:,ii-55).*(1-RFsamps(indx,end));
    end

    tmp_obsemiss{mol_ii} = cat(2,ObsEmissSamps(ResampleIndex_LT,:),topDownEmissProj{mol_ii}(2:end,:)')';

    subplot(1,3,mol_ii)
    clear vartmp
    var1 = 0.001*EmissSamps(indx,:);
    var2 = 0.001*ObsEmissSamps(ResampleIndex_LT,:);
    var3 = 0.001*topDownEmissProj{mol_ii}';
    var4 = 0.001*BankEmiss{mol_ii};

    MED = prctile(var1,50); 
    LB = MED - prctile(var1,2.5);    UB = prctile(var1,97.5) - MED;
    p1 = boundedline(y1:y2, MED', [LB',UB'],'alpha','cmap',[0.3,0.3,0.3]);
    
    hold on;
    MED = prctile(var2,50); 
    LB = MED - prctile(var2,2.5);    UB = prctile(var2,97.5) - MED;
    p2 = boundedline(y1:y2, MED',[LB', UB'],'alpha','cmap',[0.7,0,0]);
    
    hold on; 
    MED = prctile(var3,50); 
    LB = MED - prctile(var3,2.5);    UB = prctile(var3,97.5) - MED;
    p3 = boundedline(y2:2016, MED',[LB', UB'],'alpha','cmap',[0,0,0.7]);
   
    hold on; 
    MED = prctile(var4,50); 
    p4 = plot(y1:2016, MED', 'k', 'LineWidth', 2);
    
    xlim([1955,2016]); ylim([0,ymax]);box on; 
    title(str); ylabel('Emissions [Gg yr^{-1}]'); xlabel('Year');
    box on; set(gca, 'FontSize', 12); 
end

lgd = legend([p1 p2 p3 p4],'bottom-up','top-down','projected','Bank Emiss')
lgd.Location = 'northwest';

figure_width = 18; % in inches
figure_height = 4; % in inches
screen_ppi = 72; % inherent property of Matlab. CHANGING THIS WILL MAKE FIGURE LOOK INCONSISTENT BETWEEN SCREEN AND EXPORT!
set(FigHandle, 'Position', [100, 100, round(figure_width*screen_ppi), round(figure_height*screen_ppi)]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [figure_width figure_height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 figure_width figure_height]);
str = strcat(FigureFolderName, 'Figure5.pdf');
print(gcf, '-dpdf', str);


prctile_opts = [97.5, 84, 50, 16, 2.5];
for mol_ii = 1:3
    DE14to16tmp = 0.001*(tmp_obsemiss{mol_ii}(end-2:end,:)' - ...
        BankEmiss{mol_ii}(:,end-2:end));
    DE02to12 = 0.001*(tmp_obsemiss{mol_ii}(end-14:end-4,:)' - ...
        BankEmiss{mol_ii}(:,end-14:end-4));
    DiffEmisstmp = 0.001*(mean(tmp_obsemiss{mol_ii}(end-2:end,:),1) - ...
        mean(tmp_obsemiss{mol_ii}(end-14:end-4,:),1));
    val1 = mean(tmp_obsemiss{mol_ii}(end-2:end,:) - ...
        BankEmiss{mol_ii}(:,end-2:end)',1);
    val2 = mean(tmp_obsemiss{mol_ii}(end-14:end-4,:) - ...
        BankEmiss{mol_ii}(:,end-14:end-4)',1);
    DiffTotalDEtmp = 0.001*(val1 - val2);
   
    for ii = 1:5
        tmp = prctile_opts(ii);
        BankEmiss2014to2016(mol_ii,ii)  = prctile(mean(0.001*BankEmiss{mol_ii}(:,end-2:end),2),tmp);
        BankEmiss2002to2012(mol_ii,ii)  = prctile(mean(0.001*BankEmiss{mol_ii}(:,end-14:end-4),2),tmp);
        TotalEmiss2014to2016(mol_ii,ii) = prctile(mean(0.001*tmp_obsemiss{mol_ii}(end-2:end,:),1),tmp);
        TotalEmiss2002to2012(mol_ii,ii) = prctile(mean(0.001*tmp_obsemiss{mol_ii}(end-14:end-4,:),1),tmp);
        TotalDE2014to2016(mol_ii,ii)    = prctile(mean(DE14to16tmp,2),tmp);
        TotalDE2002to2012(mol_ii,ii)    = prctile(mean(DE02to12,2),tmp);
        DiffTotalEmiss(mol_ii,ii)       = prctile(DiffEmisstmp,tmp);
        DiffDirTotalEmiss(mol_ii,ii)    = prctile(DiffTotalDEtmp,tmp);
    end
end

% Table Creation
RowNames = {'BankEmiss2014to2016', 'TotalEmiss2014to2016', ...
    'DirectEmiss2014to2016', 'BankEmiss2002to2012', ...
    'TotalEmiss2002to2012','DirectEmiss2002to2012'}
VarNames{1} = {'CFC11pctile97p5',  'ptile84', 'median', 'ptile16', 'ptile2p5'};
VarNames{2} = {'CFC12pctile97p5',  'ptile84', 'median', 'ptile16', 'ptile2p5'};
VarNames{3} = {'CFC113pctile97p5', 'ptile84', 'median', 'ptile16', 'ptile2p5'};

for mol_ii = 1:3;
    Data_tab = [BankEmiss2014to2016(mol_ii,:)
               TotalEmiss2014to2016(mol_ii,:)
               TotalDE2014to2016(mol_ii,:)
               BankEmiss2002to2012(mol_ii,:)
               TotalEmiss2002to2012(mol_ii,:)
               TotalDE2002to2012(mol_ii,:)];

    T = table(Data_tab(:,1),Data_tab(:,2),Data_tab(:,3), Data_tab(:,4), ...
        Data_tab(:,5),'VariableNames',VarNames{mol_ii},'RowNames',RowNames)
end

% Table Creation
RowNames2 = {'97.5th','84th','50th','16th','2.5th'}
VarNames2 = {'CFC11','CFC12','CFC113'};

T = table(DiffTotalEmiss(1,:)', DiffTotalEmiss(2,:)', DiffTotalEmiss(3,:)', ...
    'VariableNames',VarNames2,'RowNames',RowNames2)

T = table(DiffDirTotalEmiss(1,:)', DiffDirTotalEmiss(2,:)', DiffDirTotalEmiss(3,:)', ...
    'VariableNames',VarNames2,'RowNames',RowNames2)


%%

FigHandle = figure(6)
for mol_ii = 1:3
    
    if mol_ii == 1;      ObsEmissSamps = ObsEmiss11;
        str = 'CFC-11 Direct Total Emissions';  ymax = 240; 
    elseif mol_ii == 2;  ObsEmissSamps = ObsEmiss12;
        str = 'CFC-12 Direct Total Emissions';  ymax = 350;
    elseif mol_ii == 3;  ObsEmissSamps = ObsEmiss113;
        str = 'CFC-113 Direct Total Emissions'; ymax = 220;
    end
    
    subplot(1, 3, mol_ii)
    
    var1 = 0.001*BankEmiss{mol_ii};
    var2 = 0.001*ObsEmissSamps(ResampleIndex_LT,:);
    var3 = 0.001*topDownEmissProj{mol_ii}(2:end,:)';
    vartmp1  = cat(2,var2,var3);
    vartmp = vartmp1 - var1;
    MED = prctile(vartmp,50); 
    LB  = MED - prctile(vartmp,2.5);      UB = prctile(vartmp,97.5) - MED;
    p1 = boundedline(1955:2016, MED', [LB', UB'],'alpha','cmap',[0.3,0.3,0.3]);
    
    ylim([0, ymax]);  xlim([1955,2016]); 
    title(str); box on; set(gca, 'FontSize', 12);  ylabel('[Gg/yr]');

end

figure_width = 18; % in inches
figure_height = 4; % in inches
screen_ppi = 72; % inherent property of Matlab. CHANGING THIS WILL MAKE FIGURE LOOK INCONSISTENT BETWEEN SCREEN AND EXPORT!
set(FigHandle, 'Position', [100, 100, round(figure_width*screen_ppi), round(figure_height*screen_ppi)]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [figure_width figure_height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 figure_width figure_height]);

FigFoldName = '/Users/meganlickley/Dropbox (MIT)/Research/Lifetimes/Figures/May18/';
str = strcat(FigFoldName, 'Figure6.pdf');
print(gcf, '-dpdf', str);

%%  Inset for 2000-2016

FigHandle = figure(6)
subplot(1,3,1)
ylim([-20,40]); xlim([2010,2016]); 
hold on; plot([2010,2016],[0,0],'k'); 
box on; set(gca, 'FontSize', 12);

subplot(1,3,2)
ylim([-15,50]); xlim([2010,2016]); 
hold on; plot([2010,2016],[0,0],'k')
box on; set(gca, 'FontSize', 12); 

subplot(1,3,3)
ylim([-5,15]); xlim([2010,2016]); 
box on; set(gca, 'FontSize', 12); ylabel('[Gg/yr]'); 
hold on; plot([2010,2016],[0,0],'k')

figure_width = 12; % in inches
figure_height = 3; % in inches
screen_ppi = 72; % inherent property of Matlab. CHANGING THIS WILL MAKE FIGURE LOOK INCONSISTENT BETWEEN SCREEN AND EXPORT!
set(FigHandle, 'Position', [100, 100, round(figure_width*screen_ppi), round(figure_height*screen_ppi)]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [figure_width figure_height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 figure_width figure_height]);
str = strcat(FigFoldName,'Figure6inset.pdf');
print(gcf, '-dpdf', str);


%%  Same as Fig  6 but accounting for 113a going from 0.5 ppt in 2012 to 0.7 ppt in 2017 and removing this from the concentrations. 

MF113a = [0.5, 0.5, 0.5, 0.54, 0.58, 0.62, 0.66, 0.7];
load('Input/CFC113/wmo2018.mat');
ytmp1 = find(wmo_yr == 2010);
ytmp2 = find(wmo_yr == 2016);
obsMF{3} = wmo_conc(ytmp1:ytmp2+1) - MF113a';
topDownEmissProj{3} = A(3)*(repmat(obsMF{3}(2:end),1,nSamps)-repmat(obsMF{3}(1:end-1),1,nSamps).*exp(-repmat(varLT',7,1)));
  
BankEmiss{3} = Bank113(ResampleIndex3,1:end-1).*RFsamps113(ResampleIndex3,2:end);
Bank113tmp(:,1) = Bank113(ResampleIndex3,55);

for ii = 56:62  
    BankEmiss{3}(:,ii) = Bank113tmp(:,ii-55).*(RFsamps113(ResampleIndex3,end));
    Bank113tmp(:,ii-54)= Bank113tmp(:,ii-55).*(1-RFsamps113(ResampleIndex3,end));
end

FigHandle = figure(13); 

var1 = 0.001*BankEmiss{3};
var2 = 0.001*ObsEmiss113(ResampleIndex_LT,:);
var3 = 0.001*topDownEmissProj{3}(2:end,:)';

vartmp1  = cat(2,var2,var3);
vartmp = vartmp1 - var1;

MED113 = MED;   LB113 = LB;   UB113 = UB; 
p2 = boundedline([1955:2016],MED113',[LB113',UB113'],'alpha','cmap',[0,0,0.7]);

hold on; 
MED = prctile(vartmp,50); 
LB = MED - prctile(vartmp,2.5);    UB = prctile(vartmp,97.5) - MED;
p1 = boundedline([1955:2016],MED',[LB',UB'],'alpha','cmap',[0.3,0.3,0.3]);

xlim([2010,2016]);
title('CFC-113 Direct Total Emissions')
box on; set(gca, 'FontSize', 12); ylabel('[Gg/yr]');
legend([p1, p2],'113 emissions only','113 + 113a emissions')

figure_width = 6; % in inches
figure_height = 4; % in inches
screen_ppi = 72; % inherent property of Matlab. CHANGING THIS WILL MAKE FIGURE LOOK INCONSISTENT BETWEEN SCREEN AND EXPORT!
set(FigHandle, 'Position', [100, 100, round(figure_width*screen_ppi), round(figure_height*screen_ppi)]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [figure_width figure_height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 figure_width figure_height]);

str = strcat(FigureFolderName,'FigureS11.pdf');
print(gcf, '-dpdf', str);


%% Supplementary Figure 1 Residuals

load(CFC11_FileName{1},'MF11','ResampleIndex1','molefractions');
load(CFC12_FileName{1},'MF12','ResampleIndex2');
load(CFC113_FileName{1},'MF113','ResampleIndex3');

FigHandle = figure(7); 

for mol_ii = 1:3
    if mol_ii == 1;         str1 = 'CFC-11 residuals,';
        MFsamps = MF11;  indx = ResampleIndex1; 
    elseif mol_ii == 2;     str1 = 'CFC-12 residuals,';
        MFsamps = MF12;  indx = ResampleIndex2; 
    elseif mol_ii == 3;     str1 = 'CFC-113 residuals,';
        MFsamps = MF113; indx = ResampleIndex3; 
    end
    
    obs_tmp = molefractions{mol_ii}(1:end);
    mod_tmp = MFsamps(indx,:);
    resids = mod_tmp - repmat(obs_tmp', size(mod_tmp,1),1);

    for decade_ii = 1:3
        subplot(3, 3, 3*(mol_ii - 1) + decade_ii)
        hist(resids(:,20+10*decade_ii),40); 
        str = strcat(str1, num2str(1975+10*decade_ii));
        title(str);   xlabel('mod - obs MF [ppt]');
        set(gca, 'FontSize', 12);
    end
end

figure_width = 14; % in inches
figure_height = 10; % in inches
screen_ppi = 72; % inherent property of Matlab. CHANGING THIS WILL MAKE FIGURE LOOK INCONSISTENT BETWEEN SCREEN AND EXPORT!
set(FigHandle, 'Position', [100, 100, round(figure_width*screen_ppi), round(figure_height*screen_ppi)]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [figure_width figure_height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 figure_width figure_height]);

str = strcat(FigureFolderName,'FigureS1.pdf');
print(gcf, '-dpdf', str);

%%  Fig S3: Bank estimate posterior against Nat Comms estimate
clearvars -except CFC11_FileName CFC12_FileName CFC113_FileName ...
    FigureFolderName HomeDir

y1 = 1955; % Year Start Date
y2 = 2010; % Year End Date
nYears = length(y1:y2);

load(CFC11_FileName{1},'Bank11','ResampleIndex1');
load(CFC12_FileName{1},'Bank12','ResampleIndex2');
load(CFC113_FileName{1},'Bank113','ResampleIndex3');
load('Figure_Data/Lickleyetal2020.mat')

FigHandle = figure(8)

subplot(1,3,1)
var = 0.001*Bank11(ResampleIndex1,:);
MED = prctile(var, 50); LB = prctile(var,2.5); UB = prctile(var,97.5);
p1 = boundedline([y1:y2], MED', [MED' - LB',UB'-MED'],'alpha','cmap',[0.7,0,0]);
hold on; 
p2 = boundedline([y1:2016], med11', [med11' - lb11', ub11' - med11'], 'alpha','cmap',[0,0,1]);
title('CFC-11 Bank size'); xlabel('Year'); ylabel('Bank [Gg]');
box on; set(gca, 'FontSize', 12); xlim([y1,2010]);

subplot(1,3,2)
var = 0.001*Bank12(ResampleIndex2,:);
MED = prctile(var,50); LB = prctile(var,2.5); UB = prctile(var,97.5);
p1 = boundedline([y1:y2],MED',[MED'-LB',UB'-MED'],'alpha','cmap',[0.7,0,0]);
hold on; 
p2 = boundedline([y1:2016], med12',[med12' - lb12', ub12' - med12'],'alpha','cmap',[0,0,1]);
title('CFC-12 Bank size'); xlabel('Year'); ylabel('Bank [Gg]');
box on; set(gca, 'FontSize', 12); xlim([y1,2010]);

subplot(1,3,3)
var = 0.001*Bank113(ResampleIndex3,:);
MED = prctile(var,50); LB = prctile(var,2.5); UB = prctile(var,97.5);
p1 = boundedline([y1:y2],MED',[MED'-LB',UB'-MED'],'alpha','cmap',[0.7,0,0]);
hold on; 
p2 = boundedline([y1:2016], med113', [med113' - lb113', ub113' - med113'],'alpha','cmap',[0,0,1]);
title('CFC-113 Bank'); xlabel('Year'); ylabel('Bank [Gg]');
box on; set(gca, 'FontSize', 12); xlim([y1,2010]);

lgd = legend([p1 p2],'Inferred LTs','Nat Comms')
lgd.Location = 'northwest';


figure_width = 18; % in inches
figure_height = 5; % in inches
screen_ppi = 72; % inherent property of Matlab. CHANGING THIS WILL MAKE FIGURE LOOK INCONSISTENT BETWEEN SCREEN AND EXPORT!
set(FigHandle, 'Position', [100, 100, round(figure_width*screen_ppi), round(figure_height*screen_ppi)]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [figure_width figure_height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 figure_width figure_height]);

str = strcat(FigureFolderName,'FigureS3.pdf');
print(gcf, '-dpdf', str);


%% Supplementary Figure 3 - Bank comparison for different BPE assumptions

FigHandle = figure(9)
load('Figure_Data/Lickleyetal2020.mat')
load('Figure_Data/FigS4_data.mat')

% Lickley et al., 2020
p1 = boundedline([1955:2016],med11',[med11' - lb11',ub11'-med11'],'alpha','cmap',[0.3,0.3,0.3]);

% New LT model, old RF
hold on;
p2 = boundedline([1955:2010],0.001*MED_schem2',0.001*[LB_schem2, UB_schem2],'alpha','cmap',[0.7,0,0]);

% Prescribed RF, prescribed LT
hold on;
p4 = boundedline([1955:2010],0.001*MED_schem1',0.001*[LB_schem1, UB_schem1],'alpha','cmap',[0.2,0.4,0]);

% New Model
hold on; 
load(CFC11_FileName{1},'Bank11','ResampleIndex1');
Bank_BPE = Bank11(ResampleIndex1,:);
MED_schem3 = prctile(Bank_BPE,50);
UB_schem3 = (prctile(Bank_BPE,97.5)-MED_schem3)';
LB_schem3 = (MED_schem3 - prctile(Bank_BPE,2.5))';
p3 = boundedline([1955:2010],0.001*MED_schem3',0.001*[LB_schem3 , UB_schem3],'alpha','cmap',[0,0,0.7]);


lgd = legend([p1 p4 p2 p3],'Lickley et al. 2020','New Model, Schem1','New Model, Schem2', 'New Model, Schem3')
lgd.Location = 'northwest';
box on; xlim([1955,2016]);  xlabel('Year'); ylabel('Bank [Gg]'); 
title('CFC-11 Bank Comparison')
set(gca, 'FontSize', 12)

figure_width = 7; % in inches
figure_height = 5; % in inches
screen_ppi = 72; 
set(FigHandle, 'Position', [100, 100, round(figure_width*screen_ppi), round(figure_height*screen_ppi)]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [figure_width figure_height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 figure_width figure_height]);

str = strcat(FigureFolderName,'FigureS4.pdf');
print(gcf, '-dpdf', str);


%% Figure S5 - WACCM modeled lifeitmes
load('Figure_Data/CFC11Lifetimes.mat')
FigHandle = figure(10); 
subplot(2,1,1)
p1 = plot(years,yearLifetimes(:,1:end),'Color',[0.7,0.7,0.7],'LineWidth',2);
hold on; 
p2 = plot(years,mean(yearLifetimes'),'k', 'LineWidth',4);
legend([p1(1), p2], 'Individual Member','Ensemble Mean');
set(gca, 'FontSize', 12); ylabel('years'); title('WACCM CFC-11 Lifetime');
xlim([min(years),max(years)]); ylim([46,60])

subplot(2,1,2)
for yr = 1:26
    runningmean(yr,:) = mean(yearLifetimes(yr:yr+4,:),1);
end
p1 = plot(years(3:end-2),runningmean(:,1:end),'Color',[0.7,0.7,0.7],'LineWidth',2);
hold on; 
p2 = plot(years(3:end-2),mean(runningmean'),'k', 'LineWidth',4);
legend([p1(1), p2], 'Individual Member','Ensemble Mean');
xlim([min(years(3:end-2)),max(years(3:end-2))]);
set(gca, 'FontSize', 12); ylabel('years'); title('WACCM CFC-11 Lifetime (5-year running mean)');

figure_width = 10; % in inches
figure_height = 6; % in inches
screen_ppi = 72; % inherent property of Matlab. CHANGING THIS WILL MAKE FIGURE LOOK INCONSISTENT BETWEEN SCREEN AND EXPORT!
set(FigHandle, 'Position', [100, 100, round(figure_width*screen_ppi), round(figure_height*screen_ppi)]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [figure_width figure_height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 figure_width figure_height]);

str = strcat(FigureFolderName,'FigureS5.pdf');
print(gcf, '-dpdf', str);


%% Do Same For Emissions for comparison

load('Figure_Data/FigureS8_data.mat')

FigHandle = figure(11)

MED = MED_lickleyetal;   UB = UB_lickleyetal;   LB = LB_lickleyetal;

subplot(1,3,1)
p1a = boundedline([1955:2016],0.001*MED',0.001*[LB,UB],'alpha','cmap',[0.3,0.3,0.3]);
subplot(1,3,2)
p1b = boundedline([1955:2016],0.001*MED',0.001*[LB,UB],'alpha','cmap',[0.3,0.3,0.3]);
subplot(1,3,3)
p1c = boundedline([1955:2016],0.001*MED',0.001*[LB,UB],'alpha','cmap',[0.3,0.3,0.3]);


subplot(1,3,1)
hold on;
MED = MED_schem1;   UB = UB_schem1;   LB = LB_schem1;
p3 = boundedline([1955:2010],0.001*MED',0.001*[LB,UB],'alpha','cmap',[0.2,0.4,0]);
lgd = legend([p1a p3],'NatComms','New Model, Schem1')
lgd.Location = 'northwest';
box on; xlim([1955,2016]);  xlabel('Year'); ylabel('Bank [Gg]'); 
title('CFC-11 Bank Comparison');  set(gca, 'FontSize', 12)


subplot(1,3,2)
hold on;
MED = MED_schem2;   UB = UB_schem2;   LB = LB_schem2;
p2 = boundedline([1955:2010],0.001*MED',0.001*[LB,UB],'alpha','cmap',[0.7,0,0]);
lgd = legend([p1b p2],'NatComms','New Model, Schem2')
lgd.Location = 'northwest';
box on; xlim([1955,2016]);  xlabel('Year'); ylabel('Bank [Gg]'); 
title('CFC-11 Bank Comparison');  set(gca, 'FontSize', 12)


subplot(1,3,3)
MED = MED_schem3;   UB = UB_schem3;   LB = LB_schem3;
hold on;
p2a = boundedline([1955:2010],0.001*MED',0.001*[LB,UB],'alpha','cmap',[0,0,0.7]);
lgd = legend([p1c p2a],'NatComms','New Model, Schem 3')
lgd.Location = 'northwest';
box on; xlim([1955,2016]);  xlabel('Year'); ylabel('Bank [Gg]'); 
title('CFC-11 Bank Comparison');  set(gca, 'FontSize', 12)


figure_width = 18; % in inches
figure_height = 6; % in inches
screen_ppi = 72; 
set(FigHandle, 'Position', [100, 100, round(figure_width*screen_ppi), round(figure_height*screen_ppi)]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [figure_width figure_height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 figure_width figure_height]);

str = strcat(FigureFolderName,'FigureS8.pdf');
print(gcf, '-dpdf', str);


%%
%Posterior RF for each gas (same as above but add CFC-113)
FigHandle = figure(12)

for mol_ii = 1:3
    if mol_ii == 1
        load(CFC11_FileName{1},'RFDE11sampleIndex','ResampleIndex1');
        load(strcat(HomeDir, 'Input/CFC11/InferredRFandDE_unreportedProd.mat'));
        indx1 = RFDE11sampleIndex;  indx2 = ResampleIndex1; str = 'CFC-11';
    elseif mol_ii == 2
        load(CFC12_FileName{1},'RFDE12sampleIndex','ResampleIndex2');
        load(strcat(HomeDir, 'Input/CFC12/InferredRFandDE.mat'));
        indx1 = RFDE12sampleIndex;  indx2 = ResampleIndex2; str = 'CFC-12';
    elseif mol_ii == 3
        load(CFC113_FileName{1},'RFDE113sampleIndex','ResampleIndex3');
        load(strcat(HomeDir, 'Input/CFC113/InferredRFandDE.mat'));
        indx1 = RFDE113sampleIndex;  indx2 = ResampleIndex3; str = 'CFC-113';
    end
    
    subplot(1, 3, mol_ii)
    rf = RF(:, indx1);
    rf = rf(:, indx2);
    plot(yearRFDE, prctile(rf',50),'k','LineWidth',2); hold on; 
    plot(yearRFDE, prctile(rf',97.5),'--k'); hold on; 
    plot(yearRFDE, prctile(rf',2.5),'--k'); 
    xlim([1955, 2018]); ylim([0,1]); ylabel('Posterior Release Fraction'); 
    legend('median','95% CI'); title(str);  set(gca, 'FontSize', 12);
end

figure_width = 18; % in inches
figure_height = 5; % in inches
screen_ppi = 72; % inherent property of Matlab. CHANGING THIS WILL MAKE FIGURE LOOK INCONSISTENT BETWEEN SCREEN AND EXPORT!
set(FigHandle, 'Position', [100, 100, round(figure_width*screen_ppi), round(figure_height*screen_ppi)]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [figure_width figure_height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 figure_width figure_height]);

str = strcat(FigureFolderName,'FigureS5.pdf');
print(gcf, '-dpdf', str);


%% Supplementary Figure 3 or more, how different lifetime assumptions result in different posterior estimates for F11
clearvars -except CFC11_FileName CFC12_FileName CFC113_FileName ...
    FigureFolderName HomeDir

FolderName = '/Volumes/G-DRIVE mobile USB/Lifetimes/';
load('Figure_Data/FigureS9andS10.mat')
%  Scenario1 = CFC11_FileName{1}; % Fugitive emissions allowed for CFC11 and 113
%  Scenario2 = CFC11_FileName{2}; % No fugitive
%  Scenario3 = assuming no LT uncertainty and LT = SPARC MMM
%  Scenario4 = assuming no LT uncertainty and LT = 52
%  Scenario5 = assuming no LT uncertainty and LT = 57
%  Scenario6 = assuming no LT uncertainty and LT = 62
%  All scenarios are contained here: 'Figure_Data/FigureS9andS10.mat'  
%  fReplace the first two scenarios with the new simulation runs:

for Scen_ii = 1:2
    FileName_tmp = strcat(CFC11_FileName{1});
    load(FileName_tmp)
    
    topdownEmissions(Scen_ii,:) = median(ObsEmiss11(ResampleIndex1,:));
    Emiss_BU(Scen_ii,:) = median(Emiss11(ResampleIndex1,:));
    MF(Scen_ii,:) = median(MF11(ResampleIndex1,:));
    MF_obs(Scen_ii,:)= molefractions{1};
    RF(Scen_ii,:)= median(RFsamps11(ResampleIndex1,:));
    Prod(Scen_ii,:) = median(Prodsamps11(ResampleIndex1,1:56));
    DE(Scen_ii,:)= median(Prodsamps11(ResampleIndex1,1:56).*DEsamps11(ResampleIndex1,:));
    Bank(Scen_ii,:) = median(Bank11(ResampleIndex1,:));
    Bankemiss(Scen_ii,:)= median(Bank11(ResampleIndex1,1:end-1).*RFsamps11(ResampleIndex1,2:end));

end
    

FigHandle = figure(13)
subplot(1,3,1)
plot([1955:2010],RF([2:end],:),'LineWidth',2);
hold on; plot([1955:2010],RF(1,:),'--','LineWidth',2);
title('Release Fraction');xlim([1955,2010]); ylabel('Release Fraction');
set(gca, 'FontSize', 14);

subplot(1,3,2)
plot([1955:2010],0.001*Prod([2:end],:),'LineWidth',2);
hold on; plot([1955:2010],0.001*Prod(1,:),'--','LineWidth',2); 
title('Production');xlim([1955,2010]); ylabel('[Gg/yr]');
set(gca, 'FontSize', 14);

subplot(1,3,3)
plot([1955:2010],0.001*Bank([2:end],:),'LineWidth',2);
hold on; plot([1955:2010],0.001*Bank(1,:),'--','LineWidth',2);
title('Banks');xlim([1955,2010]); ylabel('[Gg]');
legend('inferred LT No fugitive','SPARC MMM','52 yrs','57 yrs','62 yrs','inferred LT fugitive')
set(gca, 'FontSize', 14);

figure_width = 18; % in inches
figure_height = 5; % in inches
screen_ppi = 72; 

screen_figure_width = round(figure_width*screen_ppi); % in pixels
screen_figure_height = round(figure_height*screen_ppi); % in pixels
set(FigHandle, 'Position', [100, 100, round(figure_width*screen_ppi), round(figure_height*screen_ppi)]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [figure_width figure_height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 figure_width figure_height]);

str = strcat(FigureFolderName,'FigureS9.pdf');
print(gcf, '-dpdf', str);
   

FigHandle = figure(14); 

subplot(2,2,1)
plot([1955:2010],0.001*topdownEmissions([2:end],:),'LineWidth',2);
hold on; plot([1955:2010],0.001*topdownEmissions(1,:),'--','LineWidth',2); 
title('TD Emissions');xlim([1955,2010]);  ylabel('[Gg/yr]');
set(gca, 'FontSize', 14);

subplot(2,2,2)
plot([1955:2010],0.001*Emiss_BU([2:end],:),'LineWidth',2);
hold on; plot([1955:2010],0.001*Emiss_BU(1,:),'--','LineWidth',2); 
title('BU Emissions');xlim([1955,2010]); ylabel('[Gg/yr]');
set(gca, 'FontSize', 14);

subplot(2,2,3)
plot([1955:2010],0.001*DE([2:end],:),'LineWidth',2);
hold on; plot([1955:2010],0.001*DE(1,:),'--','LineWidth',2); 
title('Direct Emissions');xlim([1955,2010]); ylabel('[Gg/yr]');
set(gca, 'FontSize', 14);

subplot(2,2,4)
plot([1956:2010],0.001*Bankemiss([2:end],:),'LineWidth',2);
hold on; plot([1956:2010],0.001*Bankemiss(1,:),'--','LineWidth',2); 
title('Bank Emissions');xlim([1955,2010]); ylabel('[Gg/yr]');
legend('inferred LT No fugitive','SPARC MMM','52 yrs','57 yrs','62 yrs','inferred LT fugitive')
set(gca, 'FontSize', 14);


figure_width = 14; % in inches
figure_height = 12; % in inches
screen_ppi = 72; 
set(FigHandle, 'Position', [100, 100, round(figure_width*screen_ppi), round(figure_height*screen_ppi)]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [figure_width figure_height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 figure_width figure_height]);

str = strcat(FigureFolderName,'FigureS10.pdf');
print(gcf, '-dpdf', str);

