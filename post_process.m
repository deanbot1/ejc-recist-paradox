%post-processing of fitted patient data to single and double exponential
%models
% just need to load patfit.xls file and/or matlab file
close all
clear all 
clc
green = [0 .75 0];
fitflag = 0; % whether to refit every individual and plot 
alp = 0.05; % alpha value, ie, plot the 100-2*alp CI's


%% Load patient data: parameters and simulations
cancer_type = 'colon';
trial_num = 1; % default to 1 when only one trial for that cancer type
switch cancer_type

    case 'colon'
    S = load('../out/patfitcol.mat');
    pat = S.patfit;
    list = load('../out/listcol.mat','list');    list = struct2cell(list);
    list = cell2mat(list);
    stats = load('../out/statscol.mat');
    
    case 'prostate'
    S = load('../out/patfitpros.mat');
    pat = S.patfit;
    list = load('../out/listpros.mat','list');
    list = struct2cell(list);
    list = cell2mat(list);
    stats = load('../out/statspros.mat');
    
    case 'head neck'
    S = load('../out/patfithn.mat');
    pat = S.patfit;
    list = load('../out/listhn.mat','list');
    list = struct2cell(list);
    list = cell2mat(list);
    stats = load('../out/statshn.mat');


    case 'ovarian'
        switch trial_num
        case 1
        S = load('../out/patfitova1.mat');
        pat = S.patfit;
        list = load('../out/listova1.mat','list');
        list = struct2cell(list);
        list = cell2mat(list);
        stats = load('../out/statsova1.mat');
        case 2
        S = load('../out/patfitova2.mat');
        pat = S.patfit;
        list = load('../out/listova2.mat','list');
        list = struct2cell(list);
        list = cell2mat(list);
        stats = load('../out/statsova2.mat');
    end
end

avg_sd = load('../out/avg_sd.mat');
avg_sd = struct2cell(avg_sd);
avg_sd = cell2mat(avg_sd);
sigma = .01* avg_sd;

%% Spaghetti plot of raw SOD data
% need to filter these by only including patients contained in list, but
% for now leave all patients


% plot filtered data
flist=[];
figure(1)
for i=1:length(list)
    j=list(i); % only plot patients whos numbers in list
      
    if length(unique(pat(j).time))>=3% get rid of weird numbers
       flist=[flist j];
        plot(pat(j).tmeas,pat(j).V,'color',[0 0 .7]);    
        hold on            
    end
    xlabel('time (days)');
    ylabel('sum of the longest dimensions');
end
hold off
title(['Tumor "Volume" plot (N=',num2str(length(list)),')'])
% set(gcf, 'Position', get(0, 'Screensize'));

disp([num2str(length(flist)), ' patients have at least 3 time points measured.']);
%% Spaghetti plot normalized by SOD, colored by paradox non-paradox
% plot filtered data
dlist=[];

for i=1:length(list)
    j=list(i); % only plot patients whos numbers in list
      
    if pat(j).double == 1
       dlist=[dlist j];
    end
end
figure();
 for i = 1:length(dlist)
      j = dlist(i);
      t = pat(j).tmeas;
      vol = pat(j).V./pat(j).V(1);
      %V = smooth(t,vol);
      V = vol;
    if pat(j).par == 0
    plot(t, V, 'g-');    
        hold on            
    end
    if pat(j).par == 1
        plot(t, V, 'r-');
    end
    xlabel('time (days)');
    ylabel('sum of the longest dimensions(t)/baseline');
end
hold off
title(['Tumor SLD(t) raw data plot (N=',num2str(length(dlist)),')'])
xlim ([ -10 700])
% set(gcf, 'Position', get(0, 'Screensize'));

%% Individual patient fits to winning model
fitflag = 1; % change to 1 to run this code
if fitflag
%figure;
% need to plot censored interval separately, color censored points
% differently, and label each point with its contribution to the
% loglikelihood function
T = [0:.1:1000]';  


    for i = 5
    figure;
    % first define cases for model fit
    %pause 
    %drawnow
    txt = [];
    dim = [.2 .6 .3 .3];
    hold off
    j=flist(i); % only plot patients whos numbers in list
     if pat(j).badfit ==0 
        icens = find(pat(j).cens_int ~= 0);
        igood = find(pat(j).cens_int == 0); % uncensored data
        errh = sigma.* pat(j).V;
        errh(icens) = sigma.* pat(j).V(icens) + pat(j).cens_int(icens);
        errl = sigma.* pat(j).V;
        
        errorbar(pat(j).tmeas(igood),pat(j).V(igood), errl(igood), errh(igood),'ro');
        hold on
        if length(icens)~= 0
        errorbar(pat(j).tmeas(icens),pat(j).V(icens), errl(icens), errh(icens),'mo', 'linewidth', 2);  
        end
        xlabel('time (days)');
        ylabel('sum of the longest dimensions');
        xlim ([ pat(j).tmeas(1), (pat(j).tmeas(end)+10)])
        title(['Patient ',num2str(j)])
        
        %Calculate Contribution of loglikeihood function of each point
        lw1 = 1;
        lw2 = 1;
        lwb = 1;
        lwg = 1;
        if pat(j).single == 1 && pat(j).single_inc == 0
            lw1 = 3;
             txt2 = ['V_0 = ', num2str(round(pat(j).phat1(1),2)),', k = ', num2str(round(pat(j).phat1(2),5))];
             plot(pat(j).tseries,pat(j).yseries1,'b','linewidth',lw1);
             if length(icens)~=0
            legend('patient data uncensored', 'patient data censored', 'single exp shrinking')
            end
            if length(icens) == 0
                legend('patient data ', 'single exp shrinking')
            end
        end
        if pat(j).twotps == 1 && pat(j).twotps_single_inc == 0
            lw1 = 3;
            txt2 = ['V_0 = ', num2str(round(pat(j).phat1(1)), 2),', k = ', num2str(round(pat(j).phat1(2)), 5)];
            plot(pat(j).tseries,pat(j).yseries1,'b','linewidth',lw1);
            if length(icens)~=0
            legend('patient data uncensored', 'patient data censored', 'single exp shrinking')
            end
            if length(icens) == 0
                legend('patient data ', 'single exp shrinking')
            end
           
        end
        if pat(j).single == 1 && pat(j).single_inc == 1
            lwg = 3;
            txt2 = ['V_0 = ', num2str(round(pat(j).phatg(1), 2)),', g = ', num2str(round(pat(j).phatg(2), 5))];
                  plot(pat(j).tseries, pat(j).yseriesg,'m', 'linewidth', lwg);
            if length(icens)~=0
            legend('patient data uncensored', 'patient data censored', 'single exp growing')
            end
            if length(icens) == 0
                legend('patient data ', 'single exp growing')
            end     
        end
        if pat(j).twotps == 1 && pat(j).twotps_single_inc == 1
            lwg = 3;
            txt2 = ['V_0 = ', num2str(round(pat(j).phatg(1)), 2),', g = ', num2str(round(pat(j).phatg(2)), 5)];
            plot(pat(j).tseries, pat(j).yseriesg,'m', 'linewidth', lwg);
            if length(icens)~=0
                legend('patient data uncensored', 'patient data censored', 'single exp growing')
            end
            if length(icens) == 0
                legend('patient data ', 'single exp growing')
            end

        end
        if pat(j).double == 1
            lw2 = 3;
            txt2 = ['V_0 = ', num2str(round(pat(j).phat2(1), 1)),', phi = ', num2str((pat(j).phat2(2))), ', k = ', num2str(round(pat(j).phat2(3),3)),', g = ', num2str(round(pat(j).phat2(4),3))];
                it = find(pat(j).tseries == round(pat(j).TTP));
                itk = find(pat(j).tseries == round(pat(j).TTPplusdk));
                 itime = find(T==max(pat(j).tseries));
                itk = find(ismember(T, round(pat(j).TTPplusdk)));
                 plot(pat(j).tseries,pat(j).yseries2,'b','linewidth',lw2);
                plot(T(1:itime), pat(j).Vplusdk(1:itime), 'color',green, 'linewidth', 2);
                if ~isnan( pat(j).TTP)
                    if pat(j).tseries(end) > pat(j).TTP % for now only plot within range of measured data
                        plot(pat(j).TTP, pat(j).yseries2(it), 'k*', 'linewidth',5);
                        plot(pat(j).TTPplusdk, pat(j).Vplusdk(itk), 'r*', 'linewidth', 5);
                    end
                end
            if length(icens)~=0
                legend('patient data uncensored', 'patient data censored', 'double exp', 'double exp k+dk',...
                    ['TTP20=', num2str(round(pat(j).TTP))], ['TTP20(k +dk)= ', num2str(round(pat(j).TTPplusdk))])
            end
            if length(icens) == 0
                legend('patient data ', 'double exp', 'double exp k+dk',...
                    ['TTP20=', num2str(round(pat(j).TTP))], ['TTP20(k +dk)= ', num2str(round(pat(j).TTPplusdk))])
            end
        end
     end
        
end
end

%% FIGURE 5A Patient example to illustrate TTB20

figure;
% need to plot censored interval separately, color censored points
% differently, and label each point with its contribution to the
% loglikelihood function
T = [0:.1:1000]';  


i=23%1:length(list)
    % first define cases for model fit
    dim = [.2 .6 .3 .3];

    j=list(i); % only plot patients whos numbers in list
%      if pat(j).badfit ==0 
        icens = find(pat(j).cens_int ~= 0);
        igood = find(pat(j).cens_int == 0); % uncensored data
        errh = sigma.* pat(j).V;
        errh(icens) = sigma.* pat(j).V(icens) + pat(j).cens_int(icens);
        errl = sigma.* pat(j).V;
        
        %Calculate Contribution of loglikeihood function of each point
      
        lw1 = 1;
        lw2 = 1;
        lwb = 1;
        lwg = 1;
        if pat(j).single == 1 && pat(j).single_inc == 0
            lw1 = 3;
             txt2 = ['V_0 = ', num2str(round(pat(j).phat1(1),2)),', k = ', num2str(round(pat(j).phat1(2),5))];
        end
        if pat(j).twotps == 1 && pat(j).twotps_single_inc == 0
            lw1 = 3;
            txt2 = ['V_0 = ', num2str(round(pat(j).phat1(1)), 2),', k = ', num2str(round(pat(j).phat1(2)), 5)];
        end
        if pat(j).single == 1 && pat(j).single_inc == 1
            lwg = 3;
            txt2 = ['V_0 = ', num2str(round(pat(j).phatg(1), 2)),', g = ', num2str(round(pat(j).phatg(2), 5))];
        end
        if pat(j).twotps == 1 && pat(j).twotps_single_inc == 1
            lwg = 3;
            txt2 = ['V_0 = ', num2str(round(pat(j).phatg(1)), 2),', g = ', num2str(round(pat(j).phatg(2)), 5)];
        end
        if pat(j).double == 1
            lw2 = 3;
            txt2 = ['V_0 = ', num2str(round(pat(j).phat2(1), 1)),', phi = ', num2str((pat(j).phat2(2))), ', k = ', num2str(round(pat(j).phat2(3),3)),', g = ', num2str(round(pat(j).phat2(4),3))];
        end

        
        
        errorbar(pat(j).tmeas(igood),pat(j).V(igood), errl(igood), errh(igood),'ro');  
        xlabel('time (days)');
        ylabel('sum of the longest dimensions');
        xlim ([ pat(j).tmeas(1), (pat(j).tmeas(end)+10)])
        ylim ([ 0, 40])
        title(['FIGURE 5A: Patient ',num2str(j)])

        hold on
        if length(icens)~= 0
        errorbar(pat(j).tmeas(icens),pat(j).V(icens), errl(icens), errh(icens),'mo', 'linewidth', 2);  
        end

        for i = 1:length(pat(j).tmeas)
            if pat(j).single == 1
            txt = [num2str(pat(j).LL1(i))];
            end
            if pat(j).double == 1
            txt= [num2str(pat(j).LL2(i))];
            end
            if pat(j).single == 0 && pat(j). double == 0
                txt = 'No model parameters';
                txt2 = 'No model parameters describe kinetics';
            end
            %text(pat(j).tmeas(i), pat(j).V(i), txt, 'HorizontalAlignment','left');
        end
        imeas = find(pat(j).tmeas(end) == pat(j).tseries);
        plot(pat(j).tseries(1:imeas),pat(j).yseries2(1:imeas),'b','linewidth',lw2);
        hold on
        plot(pat(j).tseries(imeas:end), pat(j).yseries2(imeas:end), '--b', 'linewidth',lw2);
%         text(100, 100, txt2)
       annotation('textbox',dim,'String',txt2,'FitBoxToText','on');


        itime = find(T==max(pat(j).tseries));

        xlabel('time (days)');
        ylabel('SOD(t)');
        xlim ([ pat(j).tmeas(1), 350])
        ylim ([ 0, 120])
        title(['FIGURE 5A: Patient ',num2str(j)])
        hold off
        
%% Make Waterfall plot of dTTP/dk values and best response for FIGYRE 3A

bestresp = [];
bestrespdk = [];
bestrespdkh = [];
bestrespdoublek = [];
dTTPdk = [];
for j = 1:length(pat)
   if ~isempty(pat(j).dTTPdk) && ~isnan(pat(j).dTTPdk)
       both = horzcat(pat(j).dTTPdk, pat(j).BORFLAG);
       dTTPdk = vertcat(dTTPdk, both);
   end
   if ~isempty(pat(j).bestrespw) && ~isnan(pat(j).bestrespw) %&& pat(j).double == 1
        both = horzcat(pat(j).bestrespw, pat(j).BORFLAG, pat(j).double);
       bestresp = vertcat(bestresp, both);
   end
   if ~isempty(pat(j).bestrespdkw) && ~isnan(pat(j).bestrespdkw) %&& pat(j). double == 1
        both = horzcat(pat(j).bestrespdkw, pat(j).BORFLAG);
       bestrespdk = vertcat(bestrespdk, both);
   end
   if ~isempty(pat(j).bestrespdkhw) && ~isnan(pat(j).bestrespdkhw) %&& pat(j). double == 1
        both = horzcat(pat(j).bestrespdkhw, pat(j).BORFLAG);
       bestrespdkh = vertcat(bestrespdkh, both);
   end
   if ~isempty(pat(j).bestrespdoublekw) && ~isnan(pat(j).bestrespdoublekw) %&& pat(j). double == 1
        both = horzcat(pat(j).bestrespdoublekw, pat(j).BORFLAG, pat(j).double);
       bestrespdoublek = vertcat(bestrespdoublek, both);
   end
end

figure;
hold on
[Y, i]= sort(dTTPdk(:,1), 'descend');
dTTPdk_val = dTTPdk(i,:);
for i = 1:length(dTTPdk)
    if dTTPdk(i,2) ==0
    bar(i, dTTPdk_val(i,1), 'b');
    end
    if dTTPdk(i,2) ==1
    bar(i,dTTPdk_val(i,1),'r');
    end
end
xlim([ 0 length(dTTPdk)])
title('Waterfall plot of \partial TTP20/\partial k')
ylabel ('\partial TTP20/\partial k')
xlabel('patients')
set(gcf,'color','w');

ind2 = bestresp(:,3) == 1;
ir = bestresp(:,1)<-0.3;
ORR = length(bestresp(ir,1))./length(bestresp);
bestrespdk = sort(bestrespdk, 'descend');
irdk = bestrespdk(:,1)<-0.3;
ORRdk = length(bestrespdk(irdk,1))./length(bestrespdk);
irdkh = bestrespdkh(:,1)<-0.3;
ORRdkh = length(bestrespdkh(irdkh,1))./length(bestrespdkh);
ird = bestrespdoublek(:,1)<-0.3;
ORRd = length(bestrespdoublek(ird,1))./length(bestrespdoublek);

bestresp2 = (bestresp(ind2,1));
bestrespdoublek2 = bestrespdoublek(ind2,1);
ir2 = bestresp2(:,1)<-0.3;
ird2 = bestrespdoublek2(:,1) <-0.3;

ORR2 = length(bestresp2(ir2,1))./length(bestresp2);
ORRd2 = length(bestrespdoublek2(ird2,1))./length(bestrespdoublek2);


% draw lines at 20% and -30 %

ct_zero = 0;
ct_zerodk = 0;
ct_zerod = 0;
ct_zerod2 = 0;
ct_zerodkh = 0;
ct_zero_2 = 0;
for j = 1:length(pat)
    if pat(j).single == 1 && pat(j).BORFLAG == 0
        if min(pat(j).yseries1) < max(1.5*pat(j).num_meas, 7.5)
        ct_zero = ct_zero + 1;
        end
    end
    if pat(j).double == 1 && pat(j).BORFLAG == 0
        if min(pat(j).yseries2) < max(1.5*pat(j).num_meas, 7.5)
        ct_zero = ct_zero + 1;
        ct_zero_2 = ct_zero_2 + 1;
        end
        if min(pat(j).Vdoublek) < max(1.5*pat(j).num_meas, 7.5)
        ct_zerod2= ct_zerod2 + 1;
        end
    end
    

    if min(pat(j).Vplusdk) < max(1.5*pat(j).num_meas, 7.5)
        ct_zerodk = ct_zerodk + 1;
    end
    if min(pat(j).Vplusdkh) < max(1.5*pat(j).num_meas, 7.5)
        ct_zerodkh = ct_zerodkh + 1;
    end
    if min(pat(j).Vdoublek) < max(1.5*pat(j).num_meas, 7.5)
        ct_zerod= ct_zerod + 1;
    end
end

CR = ct_zero./length(bestresp);
CRdk = ct_zerodk./length(bestrespdk);
CRdkh = ct_zerodkh./length(bestrespdkh);
CRd = ct_zerod./length(bestrespdoublek);
CR2 = ct_zero_2./length(bestresp2);
CRd2 = ct_zerod2./length(bestrespdoublek2);

upper_lim_resp = ones([1 length(bestresp)]).*20;
lower_lim_resp = ones([1 length(bestresp)]).*(-30);
 %FIGURE 3A
figure;
hold on
subplot(2,1,1)
bar(100*sort(bestresp(:,1), 'descend'), 'b')
% text(70, 45, ['ORR = ',num2str(round(100*ORR,1)),'%'])
% text(70, 80, ['CRR = ',num2str(round(100*CR,1)),'%'])
hold on
plot(1:1:length(bestresp), upper_lim_resp, 'r-')
plot(1:1:length(bestresp), lower_lim_resp, 'r-')
xlim([ 0 length(bestresp)])
ylim([ -100 100])
ylabel ('best % reponse')
xlabel('patients')
title('FIGURE 3A')

subplot(2,1,2)
bar(100*sort(bestrespdoublek(:,1), 'descend'),'facecolor',green,'edgecolor','none')
hold on
% text(70, 45, ['ORR = ',num2str(round(100*ORRd, 1)),'%'])
% text(70, 80, ['CRR = ',num2str(round(100*CRd, 1)),'%'])
plot(1:1:length(bestresp), upper_lim_resp, 'r-')
plot(1:1:length(bestresp), lower_lim_resp, 'r-')
xlim([ 0 length(bestresp)])
ylim([ -100 100])
%title('Waterfall plot of best response k=10k')
ylabel ('best % response')
xlabel('patients')


set(gcf,'color','w');
% plot only patients who respond and relapse



%% Make Kaplan Meyer Curve, k histogram, and OS curve

% initialize vectors and matrices of TTP values for each patient

TTP = [];
TTPplusdk = [];
TTPplusdkh = [];
TTPdoublek = [];
TTPallk = []; % each row will be a patient, each column a k value
TTB20 = [];
TTB20doublek = [];
TTB20allk = []; % each row will be a patient, each column a k value
OSobs = [];
list_TTP = [];
list_TTPplusdk = [];
list_TTPdoublek = [];
list_TTPplusdkh = [];
list_OS = [];
k = [];
dk = 0.005;
dkh = 0.05;
dkm = 10;


% assume that 
for j = 1:length(pat)
    if ~isempty(pat(j).TTP)
        both = horzcat(pat(j).TTP, pat(j).BORFLAG, pat(j).double);
        TTP = vertcat(TTP, pat(j).TTP);
    
        % find k values
        if pat(j).double == 1
            k = vertcat(k, pat(j).phat2(4));
        elseif pat(j).BORFLAG == 1
            k = vertcat(k, pat(j).phat1(2));
        else 
            k = vertcat(k,0);
        end
    end
    
    if~isempty(pat(j).TTPplusdk)
        both = horzcat(pat(j).TTPplusdk, pat(j).BORFLAG);
        TTPplusdk = vertcat(TTPplusdk, pat(j).TTPplusdk);
    end  
    if~isempty(pat(j).TTB20)
        both = horzcat(pat(j).TTB20, pat(j).BORFLAG);
        TTB20 = vertcat(TTB20, pat(j).TTB20);
    end  
    if~isempty(pat(j).TTB20doublek)
        both = horzcat(pat(j).TTB20doublek, pat(j).BORFLAG);
        TTB20doublek = vertcat(TTB20doublek, pat(j).TTB20doublek);
    end  
    if~isempty(pat(j).TTPplusdkh)
        both = horzcat(pat(j).TTPplusdkh, pat(j).BORFLAG);
        TTPplusdkh = vertcat(TTPplusdkh, pat(j).TTPplusdkh);
        
    end  
    if~isempty(pat(j).TTPdoublek)
        both = horzcat(pat(j).TTPdoublek(end), pat(j).BORFLAG, pat(j).double);
        TTPdoublek = vertcat(TTPdoublek, pat(j).TTPdoublek(end));
       
    end  
    if~isempty(pat(j).TTP)
        both = horzcat(pat(j).OSobs, pat(j).BORFLAG, pat(j).double);
        OSobs = vertcat(OSobs, pat(j).OSobs);
    end
    
    if ~isempty(pat(j).TTPallk)
        TTPallk = vertcat(TTPallk, pat(j).TTPallk);
    end
    
    if ~isempty(pat(j).TTB20allk)
        TTB20allk = vertcat(TTB20allk, pat(j).TTB20allk);
    end
    
end
% get censored intervals, then change TTP = NaN to last day measured

TTP_tot = horzcat(TTP(:,1), bestresp(:,2), bestresp(:,3));
TTPdoublek_tot = horzcat(TTPdoublek(:,1),bestresp(:,2), bestresp(:,3));
TTB20_tot = horzcat(TTB20(:,1), bestresp(:,2), bestresp(:,3));
TTB20doublek_tot = horzcat(TTB20doublek(:,1),bestresp(:,2), bestresp(:,3));
OSobs_tot = horzcat(OSobs, bestresp(:,2), bestresp(:,3));

ind = TTP_tot(:,2) == 1; % patient not responding i.e. BORFLAG == 1
for i = 1:length(TTP_tot)
    ind2(i) = TTP_tot(i,2) == 0 && TTP_tot(i,3) == 1; % patient responding, and double exponential
    ind3(i) = TTP_tot(i,2) == 0 && TTP_tot(i,3) == 0; % patient responding, but single exponential 
end
TTP_nr = TTP_tot(ind,1);
TTP_resp = TTP_tot(~ind, 1);
TTP_2exp = TTP_tot(ind2,1);
TTP_1expr = TTP_tot(ind3,1);

TTPdoublek_nr = TTPdoublek_tot(ind,1);
TTPdoublek_resp = TTPdoublek_tot(~ind,1);
TTPdoublek_2exp = TTPdoublek_tot(ind2,1);
TTPdoublek_1expr = TTPdoublek_tot(ind3,1);

TTB20_nr = TTB20_tot(ind,1);
TTB20_resp = TTB20_tot(~ind, 1);
TTB20_2exp = TTB20_tot(ind2,1);
TTB20_1expr = TTB20_tot(ind3,1);

TTB20doublek_nr = TTB20doublek_tot(ind,1);
TTB20doublek_resp = TTB20doublek_tot(~ind,1);
TTB20doublek_2exp = TTB20doublek_tot(ind2,1);
TTB20doublek_1expr = TTB20doublek_tot(ind3,1);


% get censored intervals, then change TTP = NaN to last day measured
cens_br = isnan(TTB20_resp); % should be 0 for observations fully observed and 1 for right-censored
cens_bnr = isnan(TTB20_nr); % these always progress, no censoring needed
cens_bdkmr = isnan(TTB20doublek_resp); % should be 0 for observations fully observed and 1 for right-censored
cens_bdkmnr = isnan(TTB20doublek_nr);

cens_r = isnan(TTP_resp); % should be 0 for observations fully observed and 1 for right-censored
cens_nr = isnan(TTP_nr); % these always progress, no censoring needed
cens_dkmr = isnan(TTPdoublek_resp); % should be 0 for observations fully observed and 1 for right-censored
cens_dkmnr = isnan(TTPdoublek_nr);

cens_r1 = isnan(TTP_1expr);
cens_r2 = isnan(TTP_2exp);

censkm_r1 = isnan(TTPdoublek_1expr);
censkm_r2 = isnan(TTPdoublek_2exp);

cens_br1 = isnan(TTB20_1expr);
cens_br2 = isnan(TTB20_2exp);

censkm_br1 = isnan(TTB20doublek_1expr);
censkm_br2 = isnan(TTB20doublek_2exp);

cens = isnan(TTP); % if no TTP, single exponential with positive kill rate wins, i.e. no TTP
censb = isnan(TTB20);
ctnans = sum(cens)
censk = isnan(TTPplusdk); % if no TTP, single exponential with positive kill rate wins, i.e. no TTP
ctnansk = sum(censk)
censkm = isnan(TTPdoublek); % if no TTP, single exponential with positive kill rate wins, i.e. no TTP
censbkm = isnan(TTB20doublek);
ctnanskm = sum(censkm)
censkh = isnan(TTPplusdkh);
ctnanskh = sum(censkh)

censOS = isnan(OSobs);
% Change TTP = NaN to last day measured
% make lists of last day measured
for j = 1:length(pat)
    if ~isempty(pat(j).TTP)
        list_TTP = vertcat(list_TTP, pat(j).tmeas(end));
    end
    if ~isempty(pat(j).TTPplusdk)
        list_TTPplusdk = vertcat(list_TTPplusdk, pat(j).tmeas(end));
    end
    if ~isempty(pat(j).TTPplusdkh)
        list_TTPplusdkh = vertcat(list_TTPplusdkh, pat(j).tmeas(end));
    end
    if ~isempty(pat(j).TTPdoublek)
        list_TTPdoublek = vertcat(list_TTPdoublek, pat(j).tmeas(end));
    end
 
end
% if TTP = NaN, let TTP + last day measured
for i = 1:length(TTP)
    if isnan(TTP(i,1))
        TTP(i,1) = list_TTP(i);
    end
    if isnan(TTPplusdk(i,1))
        TTPplusdk(i,1) = list_TTPplusdk(i);
    end
    if isnan(TTPplusdkh(i,1))
        TTPplusdkh(i,1) = list_TTPplusdkh(i);
    end
    if isnan(TTPdoublek(i,1))
        TTPdoublek(i,1) = list_TTPdoublek(i);
    end
    if isnan(OSobs(i,1))
        OSobs(i,1) = list_TTP(i);
    end
    if isnan(TTB20(i,1))
        TTB20(i,1) = list_TTP(i);
    end
    if isnan(TTB20doublek(i,1))
        TTB20doublek(i,1) = list_TTP(i);
    end
    if isnan(TTPallk(i,1))
        TTPallk(i,:) = list_TTP(i);
    end
    if isnan(TTB20allk(i,1))
        TTB20allk(i,:) = list_TTP(i);
    end
end

% Make k distributions

for i = 1:length(k)
    if k(i) >0 % add this if
    kdk(i) = k(i) + dk;
    kdkh(i) = k(i) + dkh;
    km(i) = k(i).*dkm;
    end
end
% Make Kaplan-Meier cuves with TTP and censored

[d,td, dlo, dup] = ecdf(TTP,'censoring',cens);
[dk,tdk, dlok, dupk] = ecdf(TTPplusdk,'censoring',censk);
[dkh,tdkh, dlokh, dupkh] = ecdf(TTPplusdkh,'censoring',censkh);
[dkm,tdkm, dlokm, dupkm] = ecdf(TTPdoublek,'censoring',censkm);

[do,tdo, dloo, dupo] = ecdf(OSobs,'censoring',censOS);

ireal = find(d>0.5,1, 'first')
isim = find(dkm >0.5,1, 'first')

medpfsreal = td(ireal)
medpfssim = tdkm(isim)

all_TTP = vertcat(TTP, TTPdoublek);
lab = vertcat(zeros(length(TTP),1), ones(length(TTPdoublek),1));
censall = logical(vertcat(cens, censkm));

[b, logl, H, stats] = coxphfit(lab, all_TTP,  'Censoring', censall)
hr = exp(b)
dev = icdf('Normal',1-alp,0,1);
hrl = exp(b-dev*stats.se)
hrh = exp(b+dev*stats.se)
figure;
set(gcf,'color','w');
stairs(td, 1-d, 'b','linewidth', 1.5)
hold on
% stairs(td, 1-dlo, 'b','linewidth', .5)
% stairs(td, 1-dup, 'b','linewidth', .5)

stairs(tdkm, 1-dkm, 'color',green,'linewidth', 1.5)
% stairs(tdkm, 1-dlokm, 'r','linewidth', .5)
% stairs(tdkm, 1-dupkm, 'r','linewidth', .5)

ylim([0 1])
xlim([0 800])
ylabel ('RECIST PFS of 1.2SODmin')
xlabel('Time (days)')
title('FIGURE 3B: Simulated effect of increasing k on RECIST PFS')
%  legend('true PFS','true upper bound', 'true lower bound','k = 10k', 'simulated upper bound', 'simulated lower bound')

 figure;
set(gcf,'color','w');
stairs(tdo, 1-do, 'b','linewidth', 1.5)
hold on
stairs(tdo, 1-dloo, 'b','linewidth', .5)
stairs(tdo, 1-dupo, 'b','linewidth', .5)
ylim([0 1])
ylabel ('OS observed')
xlabel('Time (days)')
title('True OS')
 legend('true OS','true upper bound', 'true lower bound')

k = sort(k);
kdk = sort(kdk);

edges = [-5:.25:1];
figure;
set(gcf,'color','w');
subplot(2,1, 1)
histogram(log10(k),edges)
hold on
xlabel ('log k')
xlim( [-4 1])
ylabel('frequency')
% title( 'Distribution of k-values from observed data')
% subplot(4,1,2)
% histogram(log10(kdk),edges)
% xlabel ('log k (negative growth rate)')
% xlim( [-5 1])
% ylabel('frequency')
% title( 'Distribution of k-values, k+dk = 0.005')
subplot(2,1,2)
histogram(log10(km),edges)
xlabel ('log k')
xlim( [-4 1])
ylabel('frequency')
% title( 'Distribution of k-values, k= 10k')
% subplot(4,1,4)
% histogram(log10(kdkh),edges)
% xlabel ('log k (negative growth rate)')
% xlim( [-5 1])
% ylabel('frequency')
% title( 'Distribution of k-values, k+dk = 0.05')
% plot(log(0.005), 0 ,'y*','linewidth',5)
% plot(log(0.05), 0, 'r*', 'linewidth',5)
%% Make Kaplan Meyer curve with "new progression" TTB20
[db,tdb, dlob, dupb] = ecdf(TTB20,'censoring',censb);

[dkmb,tdkmb, dlokmb, dupkmb] = ecdf(TTB20allk(:,end),'censoring',censbkm);

ireal = find(db>0.5,1, 'first');
isim = find(dkmb >0.5,1, 'first');
mednpfsreal = tdb(ireal)
mednpfssim = tdkmb(isim)

all_TTB20 = vertcat(TTB20, TTB20allk(:,end));
lab = vertcat(zeros(length(TTB20),1), ones(length(TTB20doublek),1));
censall = logical(vertcat(censb, censbkm));

[b, logl, H, stats] = coxphfit(lab, all_TTB20, 'Censoring', censall)
hr = exp(b)
dev = icdf('Normal',1-alp,0,1);
hrh = exp(b+dev*stats.se)
hrl = exp(b-dev*stats.se)


figure;
set(gcf,'color','w');
stairs(tdb, 1-db, 'b','linewidth', 1.5)
hold on
% stairs(td, 1-dlo, 'b','linewidth', .5)
% stairs(td, 1-dup, 'b','linewidth', .5)
stairs(tdkmb, 1-dkmb, 'm','linewidth', 1.5)
% stairs(tdkm, 1-dlokm, 'r','linewidth', .5)
% stairs(tdkm, 1-dupkm, 'r','linewidth', .5)
ylim([0 1])
xlim([0 800])
ylabel ('"New PFS" (patients not at 1.2V_0)')
xlabel('Time (days)')
title('FIGURE 6: Simulated effect of increasing k on " new PFS" curve')
%  legend('true new PFS','k = 10k new PFS')



%% Plot Hazard ratio as a function of k for TTP20 and TTB20

% censoring -- do not uncomment this block!
% cens = isnan(TTP); % if no TTP, single exponential with positive kill rate wins, i.e. no TTP
% censb = isnan(TTB20);
% ctnans = sum(cens)
% censk = isnan(TTPplusdk); % if no TTP, single exponential with positive kill rate wins, i.e. no TTP
% ctnansk = sum(censk)
% censkm = isnan(TTPdoublek); % if no TTP, single exponential with positive kill rate wins, i.e. no TTP
% censbkm = isnan(TTB20doublek);
% ctnanskm = sum(censkm)
% censkh = isnan(TTPplusdkh);
% ctnanskh = sum(censkh)

thinby = 1;
kfactor = logspace(-1,1,100/thinby);
foo = find(cens ~= censb);  % find weird patients who have TTB120 censored and TTP20 not censored!
jgood = find(~ismember(1:length(TTP),foo));
labels = vertcat(zeros(length(TTP(jgood)),1), ones(length(TTP(jgood)),1));
censall = [cens(jgood);cens(jgood)];
censallb = [censb(jgood);censb(jgood)];

Nboot = 2500*0;
clear b bnew
for j = 1:length(kfactor)
    
    if Nboot > 0
            
    b(:,j) = bootstrp(Nboot,@coxphfit,labels,vertcat(TTP(jgood), TTPallk(jgood,(j-1)*thinby + 1)));
    bnew(:,j) = bootstrp(Nboot,@coxphfit,labels,vertcat(TTB20(jgood), TTB20allk(jgood,(j-1)*thinby + 1)));

    hr(j) = exp(median(b))';
    hrlo(j) = exp(quantile(b,alp))';
    hrhi(j)= exp(quantile(b,1-alp))';
    hrnew(j) = exp(median(bnew))';
    hrlonew(j) = exp(quantile(bnew,alp))';
    hrhinew(j) = exp(quantile(bnew,1-alp))';
    save('../out/daboot.mat','b','bnew','kfactor','Nboot');

    else
    [b,~,~,stats] = coxphfit(labels,vertcat(TTP(jgood), TTPallk(jgood,(j-1)*thinby + 1)),'censoring',censall);
    hr(j) = exp(b);
    dev = icdf('Normal',1-alp,0,1);
    hrlo(j) = exp(b-dev*stats.se);
    hrhi(j) = exp(b+dev*stats.se);
    [b,~,~,stats] = coxphfit(labels, vertcat(TTB20(jgood), TTB20allk(jgood,(j-1)*thinby + 1)),'censoring',censallb);
        hrnew(j) = exp(b); 
    hrlonew(j) = exp(b-dev*stats.se);
    hrhinew(j) = exp(b+dev*stats.se);

    end
    disp(j);
    
end


%beep; beep; beep;

% plot the results (which take about 30+ minutes to generate for 10x1000)
    

figure;
hold on;
plot(log10(kfactor([1 end])),[1 1],'k-');
plot([0 0],[0.7 1.3],'k-'); 
plot(log10(kfactor), hr, 'b-','LineWidth',3);
text(log10(kfactor(end)),hr(end),' TTP20','Color','b');
fill(log10(kfactor([1:end end:-1:1])),[hrlo';hrhi(end:-1:1)'],'b','edgecolor','none');
plot(log10(kfactor), hrnew, 'm-','LineWidth',3)
text(log10(kfactor(end)),hrnew(end),' TTB120','Color','m');
fill(log10(kfactor([1:end end:-1:1])),[hrlonew';hrhinew(end:-1:1)'],'m','edgecolor','none');
alpha(0.15)



ylim([0.7 1.3])
xlabel(' \leftarrow active less effective     k_{active}/k_{control}      active more effective \rightarrow')
ylabel( '\leftarrow favors active HR favors control \rightarrow')
title('FIGURE 7: Simulated Hazard Ratio by Endpoint')
set(gcf,'color','w');
xlim([log10(kfactor(1)),log10(kfactor(end))])

xtix = [0.1 0.25 0.5 1.0 2.5 5.0 10];
set(gca,'Xtick',log10(xtix),'XtickLabel',num2str(xtix','%2.1f'));

%%
%stop
%% Make pat table and fit predictors
kdbl = [];

% first find and plot k -distribution for all double exponential patients
for j = 1:length(pat)
    if pat(j). double == 1
        kdbl = vertcat(kdbl, pat(j).phat2(4));
    end
end

edges = [-5:.25:1];
edges2 = [0:0.004:.1];
figure;
subplot(2,1,1)
set(gcf,'color','w');
histogram(log10(kdbl),edges)
hold on
xlabel ('log k (negative growth rate)')
xlim( [-5 1])
ylabel('frequency')
title( 'Distribution of k-values for double exponentials')
subplot(2,1,2)
set(gcf,'color','w');
histogram(kdbl, edges2)
hold on
xlabel (' k (negative growth rate)')
xlim( [0 .1])
ylabel('frequency')
title( 'Distribution of k-values for double exponentials')

kmed = median(kdbl);
% Now separate patients into trt group (k> kmed) and control group (k<kmed)
ct_trt = 0;
ct_ctrl = 0;
for j = 1:length(pat)
    if pat(j).double == 1
        if pat(j).phat2(4) > kmed
            pat(j).trt = 1;
            ct_trt = ct_trt +1;
        end
        if pat(j).phat2(4) <kmed
            pat(j).trt = 0;
            ct_ctrl = ct_ctrl + 1;
        end
    end
end

% Make a separate structure of patient parameters, and treatment/control
% status
% i = 1:1:length(kdbl);
dbls = [];
for j = 1:length(pat)
    if pat(j).double == 1
        dbls = [dbls j];
    end
end

for i = 1:length(dbls)
    j = dbls(i);
       pat_match(i).ID = pat(j).ID;
       pat_match(i).V_o = pat(j).phat2(1);
       pat_match(i).phi = pat(j).phat2(2);
       pat_match(i).gr = pat(j).phat2(3);
       pat_match(i).k = pat(j).phat2(4);
       pat_match(i).Vmin = pat(j).bestresp;
       pat_match(i).TTP = pat(j).TTP;
       pat_match(i).trt = pat(j).trt;
       pat_match(i).TTPobs = pat(j).TTPobs;
       pat_match(i).OSobs = pat(j).OSobs;
       pat_match(i).Tlast = pat(j).tmeas(end);
       pat_match(i).Vminobs = pat(j).Vminobs;
       pat_match(i).PFS = pat(j).PFS;
end
% make treatment and control the same length 
for j = 1:length(pat_match)
    ia(j) = isempty(pat_match(j).trt);
end

pat_match = pat_match(~ia);


% Fit the predictors (phi, g, k, Vo) to logistic regression model of TTP
pat_table = struct2table(pat_match);
l = height(pat_table);
phi = log(pat_table.phi./(1-pat_table.phi));
gr = log(pat_table.gr);
k = log(pat_table.k);
V_o = log(pat_table.V_o);

% phi = pat_table.phi;
% gr = pat_table.gr;
% k = pat_table.k;
% V_o = pat_table.V_o;


col_ones = ones( l, 1);
TTP = pat_table.TTP;
 



pred = horzcat(col_ones, k, phi, gr, V_o);
[b, bint,r, rint, stats] = regress(log(TTP), pred);
[bcox, logl, H, stats] = coxphfit(pred(:, 2:end), TTP);
HR = exp(bcox(1).*(k-mean(k)) + bcox(2).*(phi-mean(phi)) + bcox(3).*(gr-mean(gr)) + bcox(4).*(V_o-mean(V_o)));
yfit = pred*b;

figure;
plot(yfit, log(TTP), 'b*')
hold on
plot(yfit,yfit, 'g-')
ylabel ('TTP20 (days)')
xlim([ 2.5 8])
ylim([2 9])
xlabel( [num2str(round(b(1),2)),' + ',num2str(round(b(2),2)),'k + ', num2str(round(b(3),2)), '\phi + ', num2str(round(b(4),2)),'g + ', num2str(round(b(5),2)),'SLD_{0}' ])
title ('FIGURE 4C: Multiple linear regression of parameters on TTP')


%% Stratify high and linear combo predictor (yfit) and simulate Kaplan-Meyer curves
% record average k and Vmin for each curve
pat_table.predictor = yfit;
for j =1:length(yfit)
    if yfit(j) < median(yfit)
        predstrat(j) = 0;
    end
    if yfit(j) > median(yfit)
        predstrat(j) = 1;
    end
end
    
pat_table.predstrat = predstrat';

% Now divide up into matrices that give TTP for each group
% let 0 = lower on predictor variab;e
% let 1 = higher on predictor variable
TTP0 = [];
TTP1 = [];
cens0 = [];
cens1 = [];
k0 = [];
k1 = [];

for j = 1:height(pat_table)
    if pat_table.predstrat(j) == 0
        if ~isnan(pat_table.TTP(j))
            TTP0 = vertcat(TTP0, pat_table.TTP(j));
            k0 = vertcat(k0, pat_table.k(j));
            cens0 = vertcat(cens0, 0);
        end
        if isnan(pat_table.TTP(j))
            TTP0 = vertcat(TTP0, pat_table.Tlast(j));
            k0 = vertcat(k0, pat_table.k(j));
            cens0 = vertcat(cens0, 1);
        end
            
    end
    if pat_table.predstrat(j) == 1
        if ~isnan(pat_table.TTP(j))
            TTP1 = vertcat(TTP1, pat_table.TTP(j));
            k1 = vertcat(k1, pat_table.k(j));
            cens1 = vertcat(cens1,0);
        end
         if isnan(pat_table.TTP(j))
             TTP1 = vertcat(TTP1, pat_table.Tlast(j));
             k1 = vertcat(k1, pat_table.k(j));
             cens1 = vertcat(cens1,1);
         end
    end
end

% make censory indicators logicals
cens0 = logical(cens0);
cens1 = logical(cens1);

[d0,td0, dlo0, dup0] = ecdf(TTP0,'censoring',cens0);
[d1,td1, dlo1, dup1] = ecdf(TTP1,'censoring',cens1);

figure;
set(gcf,'color','w');
stairs(td0, 1-d0, 'b','linewidth', 1.5)
hold on
stairs(td0, 1-dlo0, 'b','linewidth', .5)
stairs(td0, 1-dup0, 'b','linewidth', .5)
stairs(td1, 1-d1, 'r','linewidth', 1.5)
stairs(td1, 1-dlo1, 'r','linewidth', .5)
stairs(td1, 1-dup1, 'r','linewidth', .5)
xlabel('time (days)')
ylabel('"PFS" (patients not at 1.2Vmin)')
xlim([0 800])
% text(td0(5), 1-d0(5), [' median k  = ', num2str(median(k0))],'horizontalalignment','right','verticalalignment','bottom')
% text(td1(5), 1-d1(5), [' median k  = ', num2str(median(k1))],'horizontalalignment','left','verticalalignment','bottom')
%  

i0 = find(d0>0.5,1, 'first');
i1 = find(d1 >0.5,1, 'first');

medpfs0 = td0(i0)
medpfs1 = td1(i1)
%% Plot TTP vs Vmin stratified by phi
% Fit each stratum with smoothing spline w/ CIs
% also make boxplots of TTP for x -axis lumped into quartiles
% also digitize data of true PFS vs Vmin observed from He et al paper
pat_table_phi_sort = sortrows(pat_table, {'phi'});
qsize = round(height(pat_table)/4);
yy1 = horzcat(pat_table_phi_sort.Vmin(1:33), pat_table_phi_sort.TTP(1:33), pat_table_phi_sort.phi(1:33));
yy2 = horzcat(pat_table_phi_sort.Vmin(34:67), pat_table_phi_sort.TTP(34:67), pat_table_phi_sort.phi(34:67));
yy3 = horzcat(pat_table_phi_sort.Vmin(68:100), pat_table_phi_sort.TTP(68:100), pat_table_phi_sort.phi(68:100));
yy4 = horzcat(pat_table_phi_sort.Vmin(101:134), pat_table_phi_sort.TTP(101:134), pat_table_phi_sort.phi(101:134));
figure;
% plot( pat_table_phi_sort.Vmin(1:33), pat_table_phi_sort.TTP(1:33),'r*','Markersize', 8)
hold on
plot( yy1(:,1),yy1(:,2), 'r.','Markersize', 12)
plot(yy2(:,1),yy2(:,2), 'm.','Markersize', 12)
plot(yy3(:,1),yy3(:,2), 'g.','Markersize', 12)
plot(yy4(:,1),yy4(:,2), 'b.','Markersize', 12)
% plot(pat_table_phi_sort.Vmin(34:68), pat_table_phi_sort.TTP(34:68), 'm.','Markersize', 8)
% plot(pat_table_phi_sort.Vmin(69:101), pat_table_phi_sort.TTP(69:101), 'g.', 'Markersize', 8)
% plot(pat_table_phi_sort.Vmin(102:134), pat_table_phi_sort.TTP(102:134), 'b.', 'Markersize', 8)
ylabel ('TTP')
xlabel ('Vmin/Vo')
ylim ([0 600])

phi1 = mean(yy1(:,3))
phi2 = mean(yy2(:,3))
phi3 = mean(yy3(:,3))
phi4 = mean(yy4(:,3))
figure;
% plot( pat_table_phi_sort.Vmin(1:33), pat_table_phi_sort.TTP(1:33),'r*','Markersize', 8)
hold on
plot(pat_table_phi_sort.Vmin, pat_table_phi_sort.TTP, 'k.','Markersize', 12)
ylabel ('TTP')
xlabel ('Vmin/Vo')
ylim ([0 600])

