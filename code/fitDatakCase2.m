 function [pat, stat_table] = fitDatakCase2( pat, list, avg_sd)
% this function fits the filtered patient data to the single and double
% exponential models. It determines which model wins for each patient, and
% if the two population evolutionary model wins, it calculates dTTPdk to
% determine if that patient falls within the paradox region

% Case 2: This time when we perturb k, we assume that a change in k has no
% effect on the phi = 1 population (all resistant)

% it outputs the patient data, along with the parameters for each model
% fit, the AIC values, the flags for which model wins, and a stats table of
% the patient fit

% write settings for parameter and data transforms

% transforms continuous variables into lognormal space and limited values
% % % into logit space
% S = load('../out/patf.mat');
% pat = S.patf;

% single exponential
pfxform1 = @(pval)[1 1].*log(pval); %'forward' parameter transform into Reals
pbxform1 = @(phat)[1 1].*exp(phat);  %'backward' parameter transform into model space
yfxform = @(y)log(y); % 'forward' transform for data and model output
ybxform = @(yhat)exp(yhat); % 'inverse' transform for data and model output

% double exponential
pfxform2 = @(pval)[1 0 1 1].*log(pval)+[0 1 0 0].*log(pval./(1-pval)); %'forward' parameter transform into Reals
pbxform2 = @(phat)[1 0 1 1].*exp(phat)+[0 1 0 0].*(exp(phat)./(1+exp(phat)));  %'backward' parameter transform into model space

 

% Bayesian fit


sigma = .01* avg_sd; % std_dev from Schwartz reproducibility paper
dk = 0.005; 
dkh = 0.05;
dkm = 10; % stands for dk multiplier use this to simulate effect of 2 x as effective kil rate

ct_pats = 0;
ct_single = 0;
ct_single_inc = 0;
ct_double = 0;
ct_par = 0;
ct_badfit = 0;
ct_lessthan2 = 0;
ct_greaterthan2 = 0;
ct_2tps = 0;
ct_2tps_inc = 0;
ct_BORFLAG_single = 0;
ct_BORFLAG_double = 0;
ct_small = 0;
kfactor = logspace(-1,1,100);


for i=1:length(list)
    ct_pats = ct_pats +1;
    j = list(i); % only fit patients in list 
     % set all markers = 0
       pat(j).lessthan2 = 0; %
       pat(j).twotps = 0;
       pat(j).single = 0;
       pat(j).double = 0;
       pat(j).par = 0;
       pat(j).badfit = 0;
       pat(j).BORFLAG = 0;
       pat(j).twotps_single_inc = 0;
       pat(j).small = 0;
    if length(pat(j).tmeas) ==1
       ct_lessthan2 = ct_lessthan2 + 1;
       pat(j).lessthan2 = 1; %
       % these do not get fit to a model
    end
if length(pat(j).tmeas) >=2
    ydata = pat(j).V';
    ytime = pat(j).tmeas;
    % make toy data set to test fitting algorithm
%     ytime = -20:40:220;

    icens = find(pat(j).cens_int ~= 0); % left censored time points
    igood = find(pat(j).cens_int == 0); % uncensored data
    
        
     modelfungood1 = @(p)simmodel1(p, ytime(igood)); % single exponential model with death
     modelfuncens1 = @(p)simmodel1(p, ytime(icens)); 
     modelfungoodg = @(p)simmodelg(p, ytime(igood)); % single exponential model with growth
     modelfuncensg = @(p)simmodelg(p, ytime(icens));
     modelfungood2 = @(p)simmodel2(p,ytime(igood)); % double exponential model function with ytime
     modelfuncens2 = @(p)simmodel2(p,ytime(icens));
     
    % INITIAL GUESSES BASED ON DATA
    V_0guess = ydata(1);
    if ydata(2)~= 0
        kguess = -(yfxform(ydata(2)) - yfxform(ydata(1)))/(ytime(2));
    end
    if ydata(2) == 0
        kguess = -(yfxform(1e-5) - yfxform(ydata(1)))/(ytime(2));
    end
    if kguess <= 0
        kguess = 1e-5;
    end
    gguess = (yfxform(ydata(end))-yfxform(ydata(end-1)))/(ytime(end)-ytime(end-1)); 
    % alter initial guesses to prevent NaNs and zeros
    if isnan(gguess)
        gguess = 1e-5;
    end
    if isinf(gguess)
        gguess = 0.8;
    end
    if gguess <= 0
        gguess = 1e-5;
    end
    phiguess = (ydata(end)/ydata(1)).*exp(-gguess*ytime(end));
    if phiguess <= 0
       phiguess = 1e-5;
    end
    % Initial guess matrices
    theta1 = [V_0guess, kguess]; % V_0 and k
    thetag = [V_0guess, gguess];
    theta2 = [V_0guess, phiguess, gguess, kguess];

    
    % for censored data, calculate upper and lower bound of interval to
    % integrate over
    ydata_high = ydata(icens) + pat(j).cens_int(icens)';
    ydata_low = ydata(icens);
 

    % need to know contribution of each point to the log likelihood
    % function... right now separately calculate contributions
    if length(icens) == 0 % any case where there are no zeros
    %single exponential with death
    loglikelihood1 = @(phat)sum(log(normpdf(yfxform(ydata),yfxform(modelfungood1(pbxform1(phat))), sigma)));
    LL1contribg= @(phat)(log(normpdf(yfxform(ydata),yfxform(modelfungood1(pbxform1(phat))), sigma)));
    %single exponential with growth
    loglikelihoodg = @(phat)sum(log(normpdf(yfxform(ydata),yfxform(modelfungoodg(pbxform1(phat))), sigma)));
    LLgcontribg= @(phat)(log(normpdf(yfxform(ydata),yfxform(modelfungoodg(pbxform1(phat))), sigma)));
    % double exponential
    loglikelihood2 = @(phat)sum(log(normpdf(yfxform(ydata),yfxform(modelfungood2(pbxform2(phat))), sigma)));
    LL2contribg = @(phat)(log(normpdf(yfxform(ydata),yfxform(modelfungood2(pbxform2(phat))), sigma)));
  
    else
    % computes log likelihood over censored interval 
    % single exponential with death
    loglikelihood1 = @(phat)sum(log(normpdf(yfxform(ydata(igood)),yfxform(modelfungood1(pbxform1(phat))), sigma)))+...
        sum(log(normcdf(yfxform(ydata_high), yfxform(modelfuncens1(pbxform1(phat))),sigma)-...
        (normcdf(yfxform(ydata_low), yfxform(modelfuncens1(pbxform1(phat))), sigma))));
    LL1contribg = @(phat)(log(normpdf(yfxform(ydata(igood)),yfxform(modelfungood1(pbxform1(phat))), sigma)));
    LL1contribcens = @(phat)(log(normcdf(yfxform(ydata_high), yfxform(modelfuncens1(pbxform1(phat))),sigma)-...
        normcdf(yfxform(ydata_low), yfxform(modelfuncens1(pbxform1(phat))), sigma)));
     % single exponential with growth
    loglikelihoodg = @(phat)sum(log(normpdf(yfxform(ydata(igood)),yfxform(modelfungoodg(pbxform1(phat))), sigma)))+...
        sum(log(normcdf(yfxform(ydata_high), yfxform(modelfuncensg(pbxform1(phat))),sigma)-...
        (normcdf(yfxform(ydata_low), yfxform(modelfuncensg(pbxform1(phat))), sigma))));
    LLgcontribg = @(phat)(log(normpdf(yfxform(ydata(igood)),yfxform(modelfungoodg(pbxform1(phat))), sigma)));
    LLgcontribcens = @(phat)(log(normcdf(yfxform(ydata_high), yfxform(modelfuncensg(pbxform1(phat))),sigma)-...
        normcdf(yfxform(ydata_low), yfxform(modelfuncensg(pbxform1(phat))), sigma)));
    % double exponential
     loglikelihood2 = @(phat)sum(log(normpdf(yfxform(ydata(igood)),yfxform(modelfungood2(pbxform2(phat))), sigma)))+...
        sum(log(normcdf(yfxform(ydata_high), yfxform(modelfuncens2(pbxform2(phat))),sigma)-...
        (normcdf(yfxform(ydata_low), yfxform(modelfuncens2(pbxform2(phat))), sigma))));
    LL2contribg = @(phat)(log(normpdf(yfxform(ydata(igood)),yfxform(modelfungood2(pbxform2(phat))), sigma)));
    LL2contribcens = @(phat)(log(normcdf(yfxform(ydata_high), yfxform(modelfuncens2(pbxform2(phat))),sigma)-...
        normcdf(yfxform(ydata_low), yfxform(modelfuncens2(pbxform2(phat))), sigma)));
    end
    
    % minimize objective function for each structural model
    objfun1 = @(phat)-loglikelihood1(phat);
    objfung = @(phat)-loglikelihoodg(phat);
    objfun2 = @(phat)-loglikelihood2(phat);
    phatbest1 = fminsearch(objfun1, pfxform1(theta1)); % find best fitting parameters
    phatbestg = fminsearch(objfung, pfxform1(thetag));
    phatbest2 = fminsearch(objfun2, pfxform2(theta2));

    
    % save parameters and model into patient structure
    pat(j).phat1 = pbxform1(phatbest1);
    %pat(j).tseries = ytime(1):1:ytime(end); % limits TTP to be within
    %measured time
    pat(j).tseries = ytime(1):1:1000;
    pat(j).yseries1 = simmodel1(pbxform1(phatbest1), pat(j).tseries); % model fit for plotting
    pat(j).LL1(igood) = LL1contribg(phatbest1);
    if length(icens)~=0
    pat(j).LL1(icens) = LL1contribcens(phatbest1);
    end
    % indices on model data where measurement occurred
    isamp = find(ismember(pat(j).tseries,ytime));
    pat(j).yhat1 = pat(j).yseries1(isamp);
    pat(j).residuals1 = pat(j).yhat1 - pat(j).V';
    pat(j).ybar = mean(pat(j).V); % average value of volumes over all measured times
    pat(j).Rsq1 = 1- (sum((pat(j).residuals1).^2)./(sum((pat(j).ybar-pat(j).V).^2))); 
    % this will be used to determine goodness of fit to count in our
    % analysis
    
    % single exponential with growth
    pat(j).phatg = pbxform1(phatbestg);
    pat(j).yseriesg = simmodelg(pbxform1(phatbestg), pat(j).tseries); % model fit for plotting
    pat(j). LLg(igood) = LLgcontribg(phatbestg);
    if length(icens)~= 0
    pat(j).LLg(icens) = LLgcontribcens(phatbestg);
    end
    pat(j).yhatg = pat(j).yseriesg(isamp);
    pat(j).residualsg = pat(j).yhatg - pat(j).V';
    pat(j).Rsqg = 1- (sum((pat(j).residualsg).^2)./(sum((pat(j).ybar-pat(j).V).^2)));
    
    % double exponential
    pat(j).phat2 = pbxform2(phatbest2);
    pat(j).yseries2 = simmodel2(pbxform2(phatbest2), pat(j).tseries); % model fit for plotting
    pat(j). LL2(igood) = LL2contribg(phatbest2);
    if length(icens)~= 0
    pat(j).LL2(icens) = LL2contribcens(phatbest2);
    end
    pat(j).yhat2 = pat(j).yseries2(isamp);
    pat(j).residuals2 = pat(j).yhat2 - pat(j).V';
    pat(j).Rsq2 = 1- (sum((pat(j).residuals2).^2)./(sum((pat(j).ybar-pat(j).V).^2)));
    
    
    % AIC for single exponential model
    num_params1 = 2;
    n = length(pat(j).tmeas);
    AIC1 = -2*loglikelihood1(phatbest1) + 2*num_params1;
    AICc1 = -2*loglikelihood1(phatbest1) + (2*num_params1*(num_params1+1))./(n-num_params1-1);
    pat(j).AIC1 = AIC1;
    pat(j).AICc1 = AICc1;
    
    num_paramsg = 2;
    n = length(pat(j).tmeas);
    AICg = -2*loglikelihoodg(phatbestg) + 2*num_paramsg;
    AICcg = -2*loglikelihoodg(phatbestg) + (2*num_paramsg*(num_paramsg+1))./(n-num_paramsg-1);
    pat(j).AICg = AICg;
    pat(j).AICcg = AICcg;

    
    % AIC for double exponential model
    num_params2 = 4;
    AIC2 = -2*loglikelihood2(phatbest2) + 2*num_params2;
    AICc2 = -2*loglikelihood2(phatbest2) + (2*num_params2*(num_params2+1))./(n-num_params2-1);
    pat(j).AIC2 = AIC2;
    pat(j).AICc2 = AICc2;
    
   
   % Count number that  have less than 2 time points, and fit 1 or 2 exponential
   if all(pat(j).size(1:pat(j).num_meas)< 2.5)
       pat(j).small = 1;
       ct_small = ct_small + 1;
   elseif (~isreal(pat(j).Rsq1) || pat(j).Rsq1<0.6 || isnan(pat(j).Rsq1)) && (~isreal(pat(j).Rsqg) || pat(j).Rsqg<0.6 || isnan(pat(j).Rsqg))&& (~isreal(pat(j).Rsq2) ||pat(j).Rsq2<0.6 || isnan(pat(j).Rsq2))
        pat(j).badfit = 1;
       ct_badfit = ct_badfit+1;
   elseif ~isfinite(AIC1) && ~isfinite(AIC2) && ~isfinite(AICg)
        pat(j).badfit = 1;
       ct_badfit = ct_badfit+1;
   elseif length(pat(j).tmeas) == 2
       ct_2tps = ct_2tps + 1;
       pat(j).twotps = 1;
       if AICg< AIC1
           pat(j).BORFLAG =1;
           pat(j).twotps_single_inc = 1;
           ct_BORFLAG_single = ct_BORFLAG_single + 1;
           ct_2tps_inc = ct_2tps_inc + 1;
       % these default to single exponential analyses
       end
   elseif length(pat(j).tmeas) >2
       ct_greaterthan2 = ct_greaterthan2 +1;
   
        if AIC1 < AIC2 || AICg < AIC2 || (AICg < AIC1 && isnan(AIC2)) || (AIC1 < AICg && isnan(AIC2))
             pat(j).single = 1;
             ct_single = ct_single + 1; % counts only singles with more than 2 tps
        end 

        if AIC2 < AIC1 && AIC2 < AICg && isreal(pat(j).Rsqg) 
            pat(j).double = 1;
            ct_double = ct_double + 1;
        end
   end
    % allocates all patients to less than 2, twotps, single, double, or bad
    % fit 
%     T = [0:.1:1000]';  
    T = [0:.1:3000]'; % make T longer to find TTB20
    iv= find(ismember(T,56));
    if pat(j).single == 1 || pat(j).twotps == 1 && pat(j).badfit == 0  && pat(j).small == 0
        % Calculate TTP and Vplusdk for all single exponentials
        
        %****** DEATH *******
        if AIC1 < AICg % single exponential death model
            pat(j).BORFLAG = 0;
            if length(pat(j).tmeas) >1
                pat(j).single_inc = 0;
            end
            P1 = num2cell(pat(j).phat1); 
            [V_0, k1] = deal(P1{:}); % our parameters
            V = V_0.*(exp(-T*k1));
            Vmin = min(V);
            Vplusdk = V_0.*(exp(-T*(k1+dk)));
            Vplusdkh = V_0.*(exp(-T*(k1+dkh)));
            Vdoublek = V_0.*(exp(-T*(dkm.*k1)));
            Vvisitdk = Vplusdk(iv);
            Vminusdk= V_0.*(exp(-T*(k1-dk)));
            Vminplusdk = min(Vplusdk);
            Vminplusdkh = min(Vplusdkh);
            Vmindoublek = min(Vdoublek);
            Vminminusdk = min(Vminusdk);
            pat(j).bestresp = (Vmin)./V_0; 
            pat(j).bestrespdk = (Vminplusdk)./V_0; 
            pat(j).bestrespdkh = (Vminplusdkh)./V_0;
            pat(j).bestrespdoublek = (Vmindoublek)./V_0;
            pat(j).bestrespw = (Vmin-V_0)./V_0; 
            pat(j).bestrespdkw = (Vminplusdk-V_0)./V_0; 
            pat(j).bestrespdkhw = (Vminplusdkh-V_0)./V_0;
            pat(j).bestrespdoublekw = (Vmindoublek-V_0)./V_0;
            pat(j).dTTPdk = NaN;
            pat(j).par = NaN;
            pat(j).BORFLAG = 0;
            pat(j).Vplusdk = Vplusdk;
            pat(j).Vplusdkh = Vplusdkh;
            pat(j).Vdoublek = Vdoublek;
            pat(j).Vminusdk = Vminusdk;
            pat(j).TTP = NaN;
            pat(j).TTB20 = NaN;
            pat(j).TTB20doublek = NaN;
            pat(j).TTPplusdk = NaN;
            pat(j). TTPplusdkh = NaN;
            pat(j).TTPdoublek = NaN;
            pat(j).TTPallk = NaN(1, length(kfactor));
            pat(j).TTB20allk = NaN(1,length(kfactor));
        end
        
        %*****GROWTH*****
        if AICg < AIC1 % single exponential growth model
            pat(j).BORFLAG = 1;
            if length(pat(j).tmeas) >1
                pat(j).single_inc = 1;
            end
            ct_BORFLAG_single = ct_BORFLAG_single +1;
            ct_single_inc = ct_single_inc + 1;
            P1 = num2cell(pat(j).phatg); 
            [V_0, g] = deal(P1{:}); % our parameters
            V = V_0.*(exp(T*g));
            Vmin = min(V);
             iv= find(ismember(T,56));
             Vvisit =V(iv);
             imin = find(V==Vmin,1,'first');
             ipro = find(V>1.2*Vmin & T > T(imin) & V-Vmin > 5 ,1,'first'); % require 5mm minimum increase
             inpro = find(V>1.2*V(1) & T> T(imin) & V-V(1) > 5, 1, 'first');% require 5mm minimum increase
            if ~isempty(ipro) % if this exists, TTP is that time
                TTP = T(ipro);
                pat(j).TTP = TTP;
            end
            if isempty(ipro)
                TTP = NaN;
                pat(j).TTP = NaN; % aka tumor does not progress 
            end
            if ~isempty(inpro) % if this exists, TTP is that time
                TTB20 = T(inpro);
                pat(j).TTB20 = TTB20;
            end
            if isempty(inpro)
                TTB20 = NaN;
                pat(j).TTB20 = NaN; % aka tumor does not progress within time frame
            end
            if ~isempty(inpro) % if this exists, TTP is that time
                TTB20doublek = T(inpro);
                pat(j).TTB20doublek = TTB20doublek;
            end
            if isempty(inpro)
                TTB20doublek = NaN;
                pat(j).TTB20doublek = NaN; % aka tumor does not progress within time frame
            end
            pat(j).dTTPdk = 0;
            pat(j).Vplusdk = V;
            pat(j).Vplusdkh = V;
            pat(j).Vdoublek = V;
            pat(j).Vminusdk = V;
            pat(j).bestresp = (Vvisit)./V_0; % took this from Dean's analyze paradox, best resp = V(8wks)
            pat(j).bestrespdk = (Vvisit)./V_0; % Vol at visit if k is increased
            pat(j).bestrespdkh = (Vvisit)./V_0;
            pat(j).bestrespdoublek = (Vvisit)./V_0;
            pat(j).bestrespw = (Vvisit-V_0)./V_0; % took this from Dean's analyze paradox, best resp = V(8wks)
            pat(j).bestrespdkw = (Vvisit-V_0)./V_0; % Vol at visit if k is increased
            pat(j).bestrespdkhw = (Vvisit-V_0)./V_0;
            pat(j).bestrespdoublekw = (Vvisit-V_0)./V_0;
            pat(j).TTPallk = TTP.*ones(1,length(kfactor));
            pat(j).TTB20allk = TTB20.*ones(1,length(kfactor));
            pat(j).TTPdoublek = TTP;
            pat(j).TTPplusdk = TTP;
            pat(j).TTPplusdkh = TTP;
        end
    end
    if pat(j).double == 1  && pat(j).small == 0
        P2 = num2cell(pat(j).phat2); 
        [V_0, phi, g, k2] = deal(P2{:}); % our parameters
        % find TTP
        T = [0:.1:3000]';  
        BORFLAG = 0; % assume V not only increasing
        pat(j).double_inc = 0;
        V = V_0.*(phi*exp(T*g) + (1-phi)*exp(-T*k2));
        Vmin = min(V);
        imin = find(V==Vmin,1,'first');
        iv= find(ismember(T,56));
        Vvisit =V(iv);
        pat(j).bestresp = (Vmin)./V_0;
        pat(j).bestrespw = (Vmin-V_0)./V_0;
   if imin==1 % V is strictly increasing
        BORFLAG = 1;
        pat(j).double_inc = 1;
        ct_BORFLAG_double = ct_BORFLAG_double + 1;
        pat(j).bestresp = (Vvisit)./V_0;
        pat(j).bestrespw = (Vvisit-V_0)./V_0;
   end
     % finds first index after Vmin where V=1.2*Vmin
     ipro = find(V>1.2*Vmin & T > T(imin) & V-Vmin > 5,1,'first');
     inpro = find(V>1.2*V(1) & T>T(imin) & V-V(1) > 5, 1, 'first');
    if ~isempty(ipro) % if this exists, TTP is that time
        TTP = T(ipro);
        pat(j).TTP = TTP;
    end
    if isempty(ipro)
        TTP = NaN;
        pat(j).TTP = NaN; % aka tumor does not progress (with two exponential fit this should never happen)
    end
    if ~isempty(inpro) % if this exists, TTP is that time
        TTB20 = T(inpro);
        pat(j).TTB20 = TTB20;
    end
    if isempty(inpro)
        TTB20 = NaN;
        pat(j).TTB20 = NaN; 
    end
 
   % find Vplusdk, Vminusdk, TTP(k+dk) and TTP(k-dk)
   Vplusdk = V_0.*(phi*exp(T*g) + (1-phi)*exp(-T*(k2+dk)));
   Vplusdkh = V_0.*(phi*exp(T*g) + (1-phi)*exp(-T*(k2+dkh)));
   Vdoublek = V_0.*(phi*exp(T*g) + (1-phi)*exp(-T*(k2.*dkm)));
   Vminusdk= V_0.*(phi*exp(T*g) + (1-phi)*exp(-T*(k2-dk)));
   Vminplusdk = min(Vplusdk);
   Vminplusdkh = min(Vplusdkh);
   Vmindoublek = min(Vdoublek);
   Vminminusdk = min(Vminusdk);
   iminplus = find(Vplusdk==Vminplusdk,1,'first');
   iminplush = find(Vplusdkh==Vminplusdkh,1,'first');
   imindoublek = find(Vdoublek==Vmindoublek,1,'first');
   iminminus = find(Vminusdk==Vminminusdk,1,'first');
   pat(j).bestrespdk = (Vminplusdk)./V_0;
   pat(j).bestrespdkh = (Vminplusdkh)./V_0;
   pat(j).bestrespdoublek = (Vmindoublek)./V_0;
   pat(j).bestrespdkw = (Vminplusdk-V_0)./V_0;
   pat(j).bestrespdkhw = (Vminplusdkh-V_0)./V_0;
   pat(j).bestrespdoublekw = (Vmindoublek-V_0)./V_0;
   % loop through kfactors, calculate V, Vmin, TTP, TTB20 and save as
   % vectors
           for i = 1:length(kfactor)
               Vk = V_0.*(phi*exp(T*g) + (1-phi)*exp(-T*(k2.*kfactor(i))));
               Vmink = min(Vk);
               imink = find(Vk==Vmink,1,'first');
               ivk= find(ismember(T,56));
               Vvisitk =V(ivk);
               iprok = find(Vk>1.2*Vmink & T > T(imink) & Vk-Vmink > 5,1,'first');
               inprok = find(Vk>1.2*Vk(1) & T>T(imink) & Vk-Vk(1) > 5, 1, 'first');
                if ~isempty(iprok), TTPallk(1,i) = T(iprok); end
                if isempty(iprok),TTPallk(1,i) = NaN; end
                if ~isempty(inprok), TTB20allk(1,i) = T(inprok);end
                if isempty(inprok),TTB20allk(1,i) = NaN; end
           end
     pat(j).TTPallk = TTPallk;
     pat(j).TTB20allk = TTB20allk;
     iproplus = find(Vplusdk>1.2*Vminplusdk & T > T(iminplus) & Vplusdk-Vminplusdk > 5,1,'first');
     iproplush = find(Vplusdkh>1.2*Vminplusdkh & T > T(iminplush) & Vplusdkh-Vminplusdkh > 5,1,'first');
     iprodoublek = find(Vdoublek>1.2*Vmindoublek & T > T(imindoublek) & Vdoublek-Vmindoublek > 5,1,'first');
     inprodoublek = find(Vdoublek>1.2*Vdoublek(1) & T>T(imindoublek) & Vdoublek-Vdoublek(1) > 5,1,'first');
     iprominus = find(Vminusdk>1.2*Vminminusdk & T > T(iminminus) & Vminusdk-Vminminusdk > 5,1,'first');
    if ~isempty(iproplus)
        TTPplusdk = T(iproplus);
        TTPplusdkh = T(iproplush);
        TTPdoublek = T(iprodoublek);
        pat(j).TTPplusdk = TTPplusdk;
        pat(j).TTPplusdkh = TTPplusdkh;
        pat(j).TTPdoublek = TTPdoublek;
    end
    if ~isempty(inprodoublek)
        TTB20doublek = T(inprodoublek);
        pat(j).TTB20doublek= TTB20doublek;
    end
    if isempty(inprodoublek)
        TTB20doublek =NaN;
        pat(j).TTB20doublek = NaN;
    end
    if isempty(iproplus)
        TTPplusdk = NaN; 
        pat(j).TTPplusdk = TTPplusdk;
        TTPplusdkh = NaN; 
        pat(j).TTPplusdkh = TTPplusdkh;
        TTPdoublek = NaN;
        pat(j).TTPdoublek = TTPdoublek;
    end
    if ~isempty(iprominus)
        TTPminusdk = T(iprominus);
    end
    if isempty(iprominus)
        TTPminusdk = Inf; 
    end
    dttpdk = (TTPplusdk - TTPminusdk)/(2*dk); % approx partial derivative of TTP wrt k
    pat(j).BORFLAG = BORFLAG;
    pat(j).dTTPdk = dttpdk;
    pat(j).Vplusdk = Vplusdk;
    pat(j).Vplusdkh = Vplusdkh;
    pat(j).Vdoublek = Vdoublek;
    pat(j).Vminusdk = Vminusdk;

if dttpdk <0 && BORFLAG == 0
    pat(j).par = 1;
    ct_par = ct_par +1;
end
    end

end
end

pct_lessthan2 = 100 * ct_lessthan2/ct_pats;
pct_2tps = 100 *ct_2tps/ct_pats;
pct_single = 100* ct_single/ct_pats;
pct_double = 100 * ct_double/ct_pats;
pct_par_tot = 100 * ct_par/ct_pats;
pct_Ok = 100 * (ct_double - ct_par)./ct_pats;
pct_par_double = 100* ct_par/ct_double;
pct_badfit = 100 * ct_badfit/ct_pats;
pct_BORFLAG = 100 * (ct_BORFLAG_single + ct_BORFLAG_double)./ct_pats;
pct_2tps_inc = 100*(ct_2tps_inc)/ct_pats;
pct_2tps_dec = 100*(ct_2tps-ct_2tps_inc)/ct_pats;
pct_single_inc = 100 * (ct_single_inc)/ct_pats;
pct_single_dec = 100 * (ct_single - ct_single_inc)/ct_pats;
pct_small = 100 * ct_small/ct_pats;

pct_single_inc = 100 *( ct_single_inc )/ct_pats;
pct_single_dec = 100 * (ct_single + ct_2tps-ct_single_inc)/ct_pats;

stats = horzcat( ct_pats, pct_lessthan2, pct_single_inc, pct_single_dec, pct_Ok, pct_par_tot, pct_badfit, pct_small);

stat_table = dataset({stats, 'N_total', 'pct_less_than_2_TPs', 'pct_single_inc','pct_single_dec','pct_OK', 'pct_par','pct_badfit', 'pct_small'});


 end
