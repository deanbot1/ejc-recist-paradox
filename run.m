% run.m is a master run file to load, filter, fit, and process clinical
% trial data for every cancer type
% calls loadDataPar, filterDataPar, and fitDataPar, then outputs this file
% as an excel fitted data cancer type file

clear
close all
clc

%% Run fit on all indications and save the output 
indications = {'colon','ovarian1','ovarian2','prostate','head neck'};
for i = 1:length(indications)
    patfit = [];
    list = [];
    stats = [];
cancer_type = indications{i};
%cancer_type = 'colon';


%load std deviation in imaging measurement from Schwartz paper
avg_sd = load('../out/avg_sd.mat');
avg_sd = struct2cell(avg_sd);
avg_sd = cell2mat(avg_sd);
[pat] = loadDataPar(cancer_type); % input:cancer type, output: cleaned patient data for specific cancer
[patf, list] = filterDataPar(pat); % input: patient data 
%output: list of patients to be fit, patf cleaned patient data with
%censored interval, tumor volume, and time

% fit each patient data set
[patfit, stats] = fitDatakCase2(patf, list, avg_sd); % input: filtered data set
% output: patfit containing patient parameter sets and stats for each
% indication

switch cancer_type

    case 'colon'
        save('../out/patfitcol.mat','patfit');
        save('../out/listcol.mat','list');
        save('../out/statscol.mat','stats');
        
    case 'ovarian1'
            save('../out/patfitova1.mat','patfit');
            save('../out/listova1.mat','list');
            save('../out/statsova1.mat','stats');

    case 'ovarian2'
            save('../out/patfitova2.mat','patfit');
            save('../out/listova2.mat','list');
            save('../out/statsova2.mat','stats');

   case 'prostate'
        save('../out/patfitpros.mat','patfit');
        save('../out/listpros.mat','list');
        save('../out/statspros.mat','stats');
   
    case 'head neck'
        save('../out/patfithn.mat','patfit');
        save('../out/listhn.mat','list');
        save('../out/statshn.mat','stats');
        end








end

