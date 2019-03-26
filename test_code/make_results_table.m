%% results table
clear all; close all
indications = {'colon','ovarian1','ovarian2','prostate','head neck'};
filenames = {'statscol.mat','statsova1.mat','statsova2.mat','statspros.mat','statshn.mat'};

clear S

for j = 1:length(indications)
    foo = load(['../out/' filenames{j}]);
    S(j) = table2struct(dataset2table(foo.stats));
end


T = struct2table(S,'rownames',indications)
    

%% plot the table

nums = T(:,2:end-1);
for j = 1:size(nums,1)
    nums{j,:} = round(T.N_total(j)*nums{j,:}/100);
end
colors = {'y','c','m','g','r','b'};
columns = {'<2 time points','1 exponential, increasing','1 exponential, decreasing','biexponential:TTP consistent with k','biexponential:TTP INCONSISTENT with k','poor fits'};

figure
hh = bar(table2array(nums),.5,'stacked');
for j = 1:length(hh)
    hh(j).FaceColor = colors{j};
end

%set(gca,'Ydir','reverse');
set(gca,'XTickLabel',T.Row);
ylabel('number of patients');
legend(hh,columns,'location','north');
title('EJC figure 2');

%% manuscript table

Tman = T(:,{'N_total','pct_par'});
Tman.N_par = Tman.N_total.*Tman.pct_par*0.01;
Tman
