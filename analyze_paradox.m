%% this code explores the evolutionary kill rate/TTP paradox
% specifically, it has been noted (see MLN8237/task_09_SNP_TUK project)
% that certain parameter ranges in the two exponential evolutionary model
% can lead to time to RECIST progression varying *inversely* with kill rate, which
% is counterintuitive to many.
%
% The model looks like this
%
% $$ V(t) = V_0 [ \phi e^{gt} + (1-\phi) e^{-kt} ] $$
% 
% where:
%
% * $V_0$ is the initial volume
% * $\phi$ is the initial resistant fraction
% * g>0 is the resistant growth rate
% * k>0 is the kill rate (actually -net "growth rate" on treatment) on sensitive cells
%
% If we define TTP to be the time > $t_{min}$ at which $V(t) = 1.2*V_{min}$,
% then we are interested in finding regions in g,k,phi parameter space for
% which
% 
% $$ \frac{\partial TTP}{\partial k} < 0 $$
% 

clear all; close all

%% first we brute force it
% Although we can find a closed form expression for time to nadir, I know
% of no analytic way to get time to progression ie time when V(t) = 1.2*Vmin. 

gam = [-1:.02:3]; % gamma the ratio of k/g 
phi = [0:.001:.1]; % phi the resistant fraction at treatment start
phi = [0:.01:1]; % try with larger range of phi
g = .01;    % growth rate (fixed)

[GAM,PHI] = meshgrid(gam,phi); % big ol grid of parameters

GAMflat = reshape(GAM,1,[]);
PHIflat = reshape(PHI,1,[]); % resistant fraction vector
Gflat = g*ones(size(GAMflat)); % growth rate vector
Kflat = Gflat.*GAMflat; % kill rate vector

T = [0:.1:1000]';  OST = ones(size(T)); 
samplefreq = 8*7; % day of first follow up/sampling frequency
Tsamp = [0:samplefreq:max(T)];
isamp = find(ismember(T,Tsamp));

%V = OST*PHIflat.*exp(T*Gflat) + OST*(1-PHIflat).*exp(-T*Kflat);
dk = 0.005;
%Vk = OST*PHIflat.*exp(T*Gflat) + OST*(1-PHIflat).*exp(-T*(Kflat+dk));

TTP = NaN*ones(size(PHIflat));  YBEST = TTP;
BORFLAG = zeros(size(TTP)); % whether a particular trajectory is boring, ie, doesn't have a minimum for t>0 
for j = 1:length(PHIflat);
	V = PHIflat(j)*exp(T*Gflat(j)) + (1-PHIflat(j))*exp(-T*Kflat(j));
   Vmin = min(V);
   imin = find(V==Vmin,1,'first');
   if imin==1 % V is strictly increasing
	   BORFLAG(j) = 1;
	   ibest = isamp(min(find(V(isamp)>1.2*V(1))));
	   ibest = isamp(2); % keep it simple...
	   if ~isempty(ibest)
		   YBEST(j) = V(ibest); % best response by waterfall definition is really worst best first bad response
	   else
		   YBEST(j) = V(1);
	   end
   end
   ipro = find(V>1.2*Vmin & T > T(imin),1,'first');
   if ~isempty(ipro)
        TTP(j) = T(ipro);
		if BORFLAG(j)==0
		YBEST(j) = V(imin); % best observed? response
		end
   end
end

TTPk = NaN*ones(size(PHIflat));
for j = 1:length(PHIflat);
	V = PHIflat(j)*exp(T*Gflat(j)) + (1-PHIflat(j))*exp(-T*(Kflat(j)+dk));
   Vmin = min(V);
   imin = find(V==Vmin,1,'first');
   ipro = find(V>1.2*Vmin & T > T(imin),1,'first');
   if ~isempty(ipro)
        TTPk(j) = T(ipro);
   end
end


dttpdk = (TTPk-TTP)/dk; % approx partial derivative of TTP wrt k 

DTTPDK = reshape(dttpdk,size(GAM));
YBEST = reshape(YBEST,size(GAM));

%% Fit patient data to model
cancer_type = 'colon';

if strcmp(cancer_type,'ovarian')
S = load('../out/pat_params_ova.mat'); % loads patient data with fitting params
pat_params = S.pat_params;
end

if strcmp(cancer_type,'colon')
    S = load('../out/pat_params_col.mat'); % loads patient data with fitting params
    pat_params = S.pat_params;
    %pats_dbl = load('../out/pat_table_col.mat');
    list = load('../out/listcol.mat','list');
    list = struct2cell(list);
    list = cell2mat(list);
end


%% plot the results

% put fitted patient parameter fit in figure


figure; 
minmax = @(x)([min(x) max(x)]);
imagesc(minmax(1*(g*gam)),minmax(phi),DTTPDK);
%[C,h]=contourf(GAM,PHI,DTTPDK); clabel(C,h); colorbar
caxis(10000*[-1 1]);
cmap = [[1-[0:0.8/31:.8]';zeros(32,1)] [zeros(32,1);[.2:.8/31:1]']  zeros(64,1)]; % green to red colormap
colormap(cmap); colorbar;
xlabel('k (net tumor shrinkage rate)');
ylabel('resistant fraction \phi');
title('\partial TTP20/\partial k');
set(gca,'Ydir','normal'); hold on;
% plot patient parameters on search of parameter space
% for i = 1:length(pat_params)
%  % if pat_params(i).double == 1
%         plot(pat_params(i).k,pat_params(i).phi,'k*');    
%         hold on  
%    %end
% end
% plot(pat_params(10).k,pat_params(10).phi,'k*');  
[C,h] = contour(g*GAM,PHI,DTTPDK,[-100000 0 100000],'m','linewidth',2); clabel(C,h,'Color','m');


% overlay boring indicators
BORFLAG = reshape(BORFLAG,size(GAM));
hold on; 
plot(1*(g*GAM(find(BORFLAG))),PHI(find(BORFLAG)),'b.');
text(mean(mean(1*(g*GAM(find(BORFLAG))))),mean(mean(PHI(find(BORFLAG)))),'t_{min} \leq 0','color','w','FontSize',14,...
	'fontweight','bold');
%h2 = plot(gam,1./(1+gam),'b-','LineWidth',2);

%% plot just plain ol' TTP as a function of the parameters
% something is a little off here, in particular the DTTP/DK = 0 level curve in
% magenta should be going thru the "peaks" of the black level curves of
% TTP.
% This might be an artifact of DTTP/DK being estimated as a right
% difference rather than a center difference...
% this is possibly the most unintuitive graph I've ever made

TTP = reshape(TTP,size(GAM));

% plot the results
figure; 
minmax = @(x)([min(x) max(x)]);
imagesc(minmax(1*(g*gam)),minmax(phi),TTP);
%[C,h]=contourf(GAM,PHI,DTTPDK); clabel(C,h); colorbar
caxis([00 300]);
cmap = [[1-[0:0.8/31:.8]';zeros(32,1)] [zeros(32,1);[.2:.8/31:1]']  zeros(64,1)]; % green to red colormap
colormap('pink'); colorbar;
xlabel('k = net tumor shrinkage rate');
ylabel('resistant fraction \phi');
title('TTP20');
set(gca,'Ydir','normal'); hold on;
[C,h] = contour(g*GAM,PHI,TTP,[0:50:300],'k','linewidth',1); clabel(C,h,'Color','k'); hold on
[C,h] = contour(g*GAM,PHI,DTTPDK,[-100000 0 100000],'m','linewidth',2); %clabel(C,h,'Color','w');

% overlay boring indicators
BORFLAG = reshape(BORFLAG,size(GAM));
hold on; 
plot(1*(g*GAM(find(BORFLAG))),PHI(find(BORFLAG)),'b.');
text(mean(mean(g*GAM(find(BORFLAG)))),mean(mean(PHI(find(BORFLAG)))),'t_{min} \leq 0','color','w','FontSize',14,...
	'fontweight','bold');
%h2 = plot(gam,1./(1+gam),'b-','LineWidth',2);

%% plot just plain ol' SCALED TTP as a function of the parameters

TTP = reshape(TTP,size(GAM));

% plot the results
figure; 
minmax = @(x)([min(x) max(x)]);
imagesc(minmax(1*(g*gam)),minmax(phi),TTP);
mx = max(TTP,[],2);
TTPflag = TTP == mx*ones(size(gam));

%[C,h]=contourf(GAM,PHI,DTTPDK); clabel(C,h); colorbar
caxis([00 300]);
cmap = [[1-[0:0.8/31:.8]';zeros(32,1)] [zeros(32,1);[.2:.8/31:1]']  zeros(64,1)]; % green to red colormap
colormap('pink'); colorbar;
xlabel('k = net tumor shrinkage rate');
ylabel('resistant fraction \phi');
title('TTP20');
set(gca,'Ydir','normal'); hold on;
% [C,h] = contour(g*GAM,PHI,TTP,[0:50:300],'k','linewidth',1); clabel(C,h,'Color','k'); hold on
% [C,h] = contour(g*GAM,PHI,DTTPDK,[-100000 0 100000],'m','linewidth',2); %clabel(C,h,'Color','w');

% overlay boring indicators
BORFLAG = reshape(BORFLAG,size(GAM));
hold on; 
plot(1*(g*GAM(find(BORFLAG))),PHI(find(BORFLAG)),'b.');
text(mean(mean(g*GAM(find(BORFLAG)))),mean(mean(PHI(find(BORFLAG)))),'t_{min} \leq 0','color','w','FontSize',14,...
	'fontweight','bold');
%h2 = plot(gam,1./(1+gam),'b-','LineWidth',2);

%% plot TTP as a function of k
figure; clear h
ymax = 400;
jloop = find(ismember(phi,[0 0.01 0.05 0.1 0.25 0.5 0.75 1]));
kvec = 1*(g*gam);
for j = jloop
	
	plot(kvec,TTP(j,:),'k-','LineWidth',2); hold on
	iboring = find(BORFLAG(j,:));
	if ~isempty(iboring), h(1)=plot(kvec(iboring),TTP(j,iboring),'b-','LineWIdth',2); end
	ikgmax = find(TTP(j,:)==max(TTP(j,:)),1,'last');
	iinc = iboring(end):ikgmax-1;
	idec = ikgmax+1:length(TTP(j,:));
	if ~isempty(iinc), h(2)=plot(kvec(iinc),TTP(j,iinc),'g-','LineWidth',2); end
	if ~isempty(idec), h(3)=plot(kvec(idec),TTP(j,idec),'r-','LineWidth',2); end
	if ikgmax < length(kvec)
	plot(kvec(ikgmax),TTP(j,ikgmax),'ko','MarkerFaceColor','w');
	end
	phistr = num2str(phi(j),'%3.2f');
	if TTP(j,end) <= ymax
		text(kvec(end),TTP(j,end),[' \phi = ' phistr]);
	else
		ilast = find(TTP(j,:)<=ymax,1,'last');
		text(kvec(ilast),ymax,['\phi = ' phistr],'HorizontalAlignment','center','VerticalAlignment','bottom');
	end
end
hold on
% throw on patient parameters
% for i = 1:length(pat_params)
%     j = list(i);
%     plot(pat_params(j).k, pat_params(j).TTP, 'k*')
% end
set(gca,'Ylim',[0 ymax]);
xlim ([0 kvec(end)])
ylabel('TTP20 (days)');
xlabel('k (net kill rate)');
legend(h,{'t_{min} < 0','\uparrow k \rightarrow \uparrow TTP20','\uparrow k \rightarrow \downarrow TTP20'});
%% plot Vmin as a function of k
% Analytical functional form for Vmin(k):

% Vmin(k) = V_0[phi*(k/g*(1-phi)/phi)^(g/(k+g)) + (1-phi)*(k/g*(1-phi)/phi)^(-k/(k+g))

phigloop = [ 0.5];
philoop = [ 0.005 0.01 0.02 0.03 0.04];
Vmin = [];
figure; clear h


for j = 1:length(philoop)
    for m = 1:length(phigloop)
     phi = philoop(j);
     g = phi./phigloop;
     phistr = num2str(phi,'%3.2f');
     gstr = num2str(g,'%3.2f');
     phigstr = num2str(phi./g, '%3.2f');
     V_0 = 100;
     k = 0:.001:0.1;
     for i = 1:length(k)
     Vmin(i, j) = V_0.*(phi*(((k(i)./g).*(1-phi)/phi))^(g./(k(i)+g)) + (1-phi).*((k(i)./g).*((1-phi)/phi)).^(-k(i)./(k(i)+g)));
     end
     plot(k, Vmin(:, j), 'linewidth', 2);
     hold on
     text(k(end-12),Vmin(end-20,j),[ ' \phi = ' phistr, ' g = ' gstr],'horizontalalignment','center','verticalalignment','top');
     xlabel ('k')
     ylabel('Minimum Volume')
     title(['Vmin vs. k for constant \phi/g = ' phigstr])
     hold on
    end  
end
 hold on
for i = 1:length(pat_params)
    j = list(i);
    if pat_params(j).double == 1
    plot(pat_params(j).k,pat_params(j).Vmin*100, 'k*')
    end
    xlim( [ k(1) k(end)])
end

%% plot TTP as a function of YBEST and perform univariate linear regression
figure; clear h
ymax = 400;
phi = [0:.01:1]; % go back to original phi
jloop = find(ismember(phi,[0 0.005 0.01 0.05 0.1]));
% jloop = find(ismember(phi,[0, 0.01 0.14 0.45]));
%kvec = 1*(g*gam);
% stratify real patients
% pat_dbl = pats_dbl.pat_table;
% pat_table_phi_sort = sortrows(pat_dbl, {'phi'});
% qsize = round(height(pat_table_phi_sort)/4);

% throw on patient parameters
hold on
% plot(yy1(:,1),yy1(:,2), 'r.','Markersize', 12)
% plot(yy2(:,1),yy2(:,2), 'm.','Markersize', 12)
% plot(yy3(:,1),yy3(:,2), 'g.','Markersize', 12)
% plot(yy4(:,1),yy4(:,2), 'b.','Markersize', 12)

% *****Uncomment to plot overlay of patient parameters****
%evomod = @(t,v0,phi,k,g)v0*(phi*exp(g*t) + (1-phi)*exp(-k*t));
%t_{min} = \frac{ln(\frac{k}{g})-ln(\frac{\phi}{1-\phi})}{k+g}
model = @(t,phi,k,g)(phi*exp(g*t)+(1-phi)*exp(-k*t)); % normalized model
tminfun = @(p,k)(log(k/p.gr) - log(p.phi/(1-p.phi)))/(k+p.gr);
vminfun = @(p,k)model(tminfun(p,k),p.phi,p.k,p.gr);
ttpfun = @(p,k)fzero(@(t)model(t,p.phi,k,p.gr) - 1.2*vminfun(p,k),...
    tminfun(p,k) + [0 log(100)/p.gr]);
dk = 0.0001;
deriv = @(fun,p)(feval(fun,p,p.k+dk)-feval(fun,p,p.k-dk))/(2*dk);
XX = []; YY = []; UU = []; VV = [];

Vmins = cat(1,pat_params.Vmin);
TTPs = cat(1,pat_params.TTP);


plot(Vmins,TTPs,'k*');

set(gca,'Ylim',[0 ymax]);
xlim( [ 0 2])
ylabel('TTP20 (days)');
xlabel('SLD_{min}/SLD_0 (or SLD_{8wks}/SLD_0 if growing)');

[b,bint,r,rint,stats] = regress(TTPs,[ones(size(Vmins)) Vmins]);

plot([0 2],b(1) + b(2)*[0 2],'m-','linewidth',2);
title(['FIGURE 4A: R^2 = ' num2str(stats(1)) ', p = ' num2str(stats(3))])


% mags  = sqrt(UU.^2 + VV.^2);
% quiver(XX,YY,UU./mags,VV./mags,.001);
%legend(h,{'t_{min} < 0','\downarrow V_{min} \rightarrow \uparrow TTP20','\downarrow V_{min} \rightarrow \downarrow TTP20'});
% legend( ['mean \phi = ', num2str(mean(yy1(:,3)))])
% meanphi1 = mean(yy1(:,3))
% meanphi2 = mean(yy2(:,3))
% meanphi3 = mean(yy3(:,3))
% meanphi4 = mean(yy4(:,3))

%% plot TTP as a function of YBEST (best depth of response)
figure; clear h
ymax = 400;
phi = [0:.01:1]; % go back to original phi
jloop = find(ismember(phi,[0 0.005 0.01 0.05 0.1]));
% jloop = find(ismember(phi,[0, 0.01 0.14 0.45]));
%kvec = 1*(g*gam);
% stratify real patients
% pat_dbl = pats_dbl.pat_table;
% pat_table_phi_sort = sortrows(pat_dbl, {'phi'});
% qsize = round(height(pat_table_phi_sort)/4);
% yy1 = horzcat(pat_table_phi_sort.Vmin(1:33), pat_table_phi_sort.TTP(1:33), pat_table_phi_sort.phi(1:33));
% yy2 = horzcat(pat_table_phi_sort.Vmin(34:67), pat_table_phi_sort.TTP(34:67), pat_table_phi_sort.phi(34:67));
% yy3 = horzcat(pat_table_phi_sort.Vmin(68:100), pat_table_phi_sort.TTP(68:100), pat_table_phi_sort.phi(68:100));
% yy4 = horzcat(pat_table_phi_sort.Vmin(101:134), pat_table_phi_sort.TTP(101:134), pat_table_phi_sort.phi(101:134));
for j = jloop
	kvec = YBEST(j,:);
	plot(kvec,TTP(j,:),'k-','LineWidth',2); hold on
	iboring = find(BORFLAG(j,:));
	if ~isempty(iboring), h(1)=plot(kvec(iboring),TTP(j,iboring),'b-','LineWIdth',2); end
	ikgmax = find(TTP(j,:)==max(TTP(j,:)),1,'last');
	iinc = iboring(end):ikgmax-1;
	idec = ikgmax+1:length(TTP(j,:));
	if ~isempty(iinc), h(2)=plot(kvec(iinc),TTP(j,iinc),'g-','LineWidth',2); end
	if ~isempty(idec), h(3)=plot(kvec(idec),TTP(j,idec),'r-','LineWidth',2); end
	if ikgmax < length(kvec)
	plot(kvec(ikgmax),TTP(j,ikgmax),'ko','MarkerFaceColor','w');
	end
	phistr = num2str(phi(j),'%3.2f');
	if TTP(j,end) <= ymax
		text(kvec(end),TTP(j,end),[' \phi = ' phistr],'horizontalalignment','center','verticalalignment','top');
	else
		ilast = find(TTP(j,:)<=ymax,1,'last');
		text(kvec(ilast),ymax,['\phi = ' phistr],'HorizontalAlignment','center','VerticalAlignment','bottom');
	end
end
% throw on patient parameters
hold on
% plot(yy1(:,1),yy1(:,2), 'r.','Markersize', 12)
% plot(yy2(:,1),yy2(:,2), 'm.','Markersize', 12)
% plot(yy3(:,1),yy3(:,2), 'g.','Markersize', 12)
% plot(yy4(:,1),yy4(:,2), 'b.','Markersize', 12)

% *****Uncomment to plot overlay of patient parameters****
%evomod = @(t,v0,phi,k,g)v0*(phi*exp(g*t) + (1-phi)*exp(-k*t));
%t_{min} = \frac{ln(\frac{k}{g})-ln(\frac{\phi}{1-\phi})}{k+g}
model = @(t,phi,k,g)(phi*exp(g*t)+(1-phi)*exp(-k*t)); % normalized model
tminfun = @(p,k)(log(k/p.gr) - log(p.phi/(1-p.phi)))/(k+p.gr);
vminfun = @(p,k)model(tminfun(p,k),p.phi,p.k,p.gr);
ttpfun = @(p,k)fzero(@(t)model(t,p.phi,k,p.gr) - 1.2*vminfun(p,k),...
    tminfun(p,k) + [0 log(100)/p.gr]);
dk = 0.0001;
deriv = @(fun,p)(feval(fun,p,p.k+dk)-feval(fun,p,p.k-dk))/(2*dk);
XX = []; YY = []; UU = []; VV = [];
for i = 1:length(pat_params)
    j = list(i);
    plot(pat_params(j).Vmin, pat_params(j).TTP, 'k*');  
    if ~isempty(pat_params(j).phi)
        if abs(tminfun(pat_params(j),pat_params(j).k)) < Inf
        dttpdk = deriv(ttpfun,pat_params(j));   
        dvmindk = deriv(vminfun,pat_params(j));
        XX(end+1) = vminfun(pat_params(j),pat_params(j).k);
        YY(end+1) = ttpfun(pat_params(j),pat_params(j).k);
        UU(end+1) = dvmindk;
        VV(end+1) = dttpdk;
        %quiver(XX(end),YY(end),UU(end),VV(end));
        end
    end
end
set(gca,'Ylim',[0 ymax]);
xlim( [ 0 2])
ylabel('TTP20 (days)');
xlabel('SLD_{min}/SLD_0 (or SLD_{8wks}/SLD_0 if growing)');

% mags  = sqrt(UU.^2 + VV.^2);
% quiver(XX,YY,UU./mags,VV./mags,.001);
legend(h,{'t_{min} < 0','\downarrow SLD_{min} \rightarrow \uparrow TTP20','\downarrow SLD_{min} \rightarrow \downarrow TTP20'});
% legend( ['mean \phi = ', num2str(mean(yy1(:,3)))])
% meanphi1 = mean(yy1(:,3))
% meanphi2 = mean(yy2(:,3))
% meanphi3 = mean(yy3(:,3))
% meanphi4 = mean(yy4(:,3))
title('FIGURE 4B                           ')
%% make a phi-stratified plot of pat-params
figure
hold on;
Vmins = cat(1,pat_params.Vmin);
TTP20 = cat(1,pat_params.TTP);
phifit = cat(1,pat_params.phi);
Nq = 4;
quants = quantile(phifit(phifit>eps & phifit < 1),Nq);
quants = [0 eps quants 1 1+eps]; Nq = length(quants)-1;
qcolors = {'r','g','b','k','m','c','y'};
for j = 1:Nq
    igood = find(phifit >= quants(j) & phifit < quants(j+1) & ~isnan(TTP20));
% switch j
%     case 1
%         igood = find(Vmins < 1);
%     case 2
%         igood = find(Vmins >= 1);
% end
    plot(Vmins(igood),TTP20(igood),'*','color',qcolors{j}); hold on;
   % disp(length(igood));
   b = robustfit(Vmins(igood),TTP20(igood));
   x = [min(Vmins(igood)) max(Vmins(igood))];
   plot(x,b(1)+b(2)*x,'-','color',qcolors{j});
end




%% plot some key trajectories to illustrate the point!
T = [0:.1:1000]';
gam = [-1:.02:3]; % gamma the ratio of k/g 
phi = [0:.001:.1]; % phi the resistant fraction at treatment start
phi = [0:.01:1]; % try with larger range of phi
g = .01;    % growth rate (fixed)
k2 = 0.03;
k1 = 0.003;
kk = [k1;k2];
phitest = 0.01;
model = @(t,phi,k,g)(phi*exp(g*t)+(1-phi)*exp(-k*t)); % normalized model

figure('Name','EJC GRAPHICAL ABSTRACT FIGURE');
colors = {[0 0 1];[0 .75 0]};
for j = 1:length(kk)
	Ymod = @(t)(model(t,phitest,kk(j),g));
	Y = Ymod(T);
	plot(T,Y,'-','LineWidth',3,'color',colors{j});
	hold on;
    imin = find(Y == min(Y))
	ittp = find(abs(g*gam-kk(j))<100*eps);
	jttp = find(phi==phitest);
	TTPj = TTP(jttp,ittp);
	plot(TTPj,Ymod(TTPj),'s','color',colors{j},'markerfacecolor',colors{j},'MarkerSize',10);
	text(TTPj,Ymod(TTPj),['  TTP20=' num2str(round(TTPj)) 'd'],'color',colors{j},'rotation',90);
 	text(25,Ymod(25),['   k=' num2str(kk(j))],'color',colors{j},'rotation',-35*j,'verticalalignment','bottom');
% 	plot([0:TTPj],Ymod([0:TTPj]),'-','LineWidth',3,'color',colors(j));
%     text(T(imin), min(Y), ['  Vmin=' num2str(round(min(Y),2)) ],'verticalalignment', 'top','horizontalalignment', 'right','color','k');
%     plot(T(imin),min(Y),'s','color','k','markerfacecolor','k','MarkerSize',10);
end

set(gca,'Xlim',[0 max(T(Y<=.8))],'Ylim',[0 1.1]);
xlabel('days on treatment');
ylabel('tumor size/baseline');
title(sprintf('g = %3.3g, \\phi = %3.3g',g,phitest));

%% plot bunch of key trajectories to illustrate the point!
model = @(t,phi,k,g)(phi*exp(g*t)+(1-phi)*exp(-k*t)); % normalized model
T = [0:.1:400]'; 
gam = [-1:.02:3]; % gamma the ratio of k/g 
phi = [0:.001:.1]; % phi the resistant fraction at treatment start
phi = [0:.01:1]; % try with larger range of phi
g = .01;    % growth rate (fixed)
k2 = 0.03;
k1 = 0.003;
kk = [0 10.^[-4:.2:-1]];
tt = zeros(size(kk));
yy = tt;
%kk = [0:.01:.1];
phitest = 0.01;
modelf = @(t,phi,k,g)(phi*exp(g*t)+(1-phi)*exp(-k*t)); % normalized model
vminf = @(phi,k,g)phi*(k/g*(1-phi)/phi)^(g/(k+g)) + (1-phi)*(k/g*(1-phi)/phi)^(-k/(k+g))
tminf = @(phi,k,g)(log(k/g) - log(phi/(1-phi)))/(k+g);
TTB120f = @(phi,k,g)fzero(@(t)(model(t,phi,k,g)-1.2),1000);
TTP20f = @(phi,k,g)fzero(@(t)(model(t,phi,k,g)-1.2*vminf(phi,k,g)),[tminf(phi,k,g) TTB120f(phi,k,g)]);


figure('Position',[307   212   389   288],'Name','EJC Figure 1 inset');

for j = 1:length(kk)
    plot(T,modelf(T,phitest,kk(j),g),'-','LineWidth',2,'COlor',.95*[0 0 1]); hold on
    if kk(j) > 0
        tt(j)  = TTP20f(phitest,kk(j),g);
    else
        tt(j) = log((phitest+.2)/phitest)/g;
    end
    yy(j)=modelf(tt(j),phitest,kk(j),g);
    %plot(tt(j),yy(j),'ro','MarkerFaceColor','r');
end

plot(tt,yy,'ro','MarkerFaceColor','r','MarkerSize',5);

ylabel('normalized tumor size');
xlabel('days on treatment');
grid on
legend(h,'TTP20','Location','NorthWest');
set(gca,'Xlim',[-10 400],'Ylim',[-.1 1.4]);
title(sprintf('g = %3.3g, \\phi = %3.3g',g,phitest));
http20 = scaledarrows(tt,yy,'color','r','headstyle','vback3');
text(tt(end),yy(end),'TTP20','color','r','verticalalignment','top','horizontalalignment','center');

figure('Name','EJC Figure 1 main');
plot(kk,tt,'r-o','MarkerFaceColor','r'); hold on; 
text(kk(end),tt(end),'  TTP20','color','r');

set(gca,'Xscale','log');
xlabel('cell kill rate k');
ylabel('time to event');
title(sprintf('g = %3.3g, \\phi = %3.3g',g,phitest));

%% plot bunch of key trajectories to illustrate the point AND ADD TTB120
model = @(t,phi,k,g)(phi*exp(g*t)+(1-phi)*exp(-k*t)); % normalized model
T = [0:.1:800]'; 
gam = [-1:.02:3]; % gamma the ratio of k/g 
phi = [0:.001:.1]; % phi the resistant fraction at treatment start
phi = [0:.01:1]; % try with larger range of phi
g = .01;    % growth rate (fixed)
k2 = 0.03;
k1 = 0.003;
kk = [0 10.^[-4:.2:-1]];
tt = zeros(size(kk));
yy = tt;
%kk = [0:.01:.1];
phitest = 0.01;
modelf = @(t,phi,k,g)(phi*exp(g*t)+(1-phi)*exp(-k*t)); % normalized model
vminf = @(phi,k,g)phi*(k/g*(1-phi)/phi)^(g/(k+g)) + (1-phi)*(k/g*(1-phi)/phi)^(-k/(k+g))
tminf = @(phi,k,g)(log(k/g) - log(phi/(1-phi)))/(k+g);
TTB120f = @(phi,k,g)fzero(@(t)(model(t,phi,k,g)-1.2),1000);
TTP20f = @(phi,k,g)fzero(@(t)(model(t,phi,k,g)-1.2*vminf(phi,k,g)),[tminf(phi,k,g) TTB120f(phi,k,g)]);

figure('Position',[307   212   389   288],'Name','EJC Figure 5 inset');
for j = 1:length(kk)
    plot(T,modelf(T,phitest,kk(j),g),'-','LineWidth',2,'COlor',.95*[0 0 1]); hold on
    if kk(j) > 0
        tt(j)  = TTP20f(phitest,kk(j),g);
    else
        tt(j) = log((phitest+.2)/phitest)/g;
    end
    ttb(j) = TTB120f(phitest,kk(j),g);
    yyb(j) = modelf(ttb(j),phitest,kk(j),g);
    yy(j)=modelf(tt(j),phitest,kk(j),g);
    %plot(tt(j),yy(j),'ro','MarkerFaceColor','r');
end

http = plot(tt,yy,'ro','markerfacecolor','r','markersize',5); hold on;
green = [0 .75 0];
httb = plot(ttb,yyb,'o','markerfacecolor',green,'color',green,'markersize',5); hold on;



ylabel('normalized tumor size');
xlabel('days on treatment');
grid on
%legend([h;hb],
set(gca,'Xlim',[-10 500],'Ylim',[-.1 1.4]);
title(sprintf('g = %3.3g, \\phi = %3.3g',g,phitest));

scaledarrows(tt,yy,'color','r','headstyle','vback3');
%text(tt(end),yy(end),'TTP20','color','r','verticalalignment','top','horizontalalignment','center');
scaledarrows(ttb,yyb,'color',green,'headstyle','vback3');
%text(ttb(end),yyb(end),' TTB120','color',green);

legend([http;httb],{'TTP20','TTB120'},'Location','NorthWest');
figure('Name','EJC Figure 5B Main');
plot(kk,tt,'r-o','MarkerFaceColor','r'); hold on; text(kk(end),tt(end),'  TTP20','color','r');
plot(kk,ttb,'-o','MarkerfaceColor',green,'color',green);text(kk(end),ttb(end),'  TTB120','color',green);

set(gca,'Xscale','log');
xlabel('cell kill rate k');
ylabel('time to event');
title(sprintf('g = %3.3g, \\phi = %3.3g',g,phitest));

%% another version of TTP field plot
figure; 
INC = zeros(size(BORFLAG)); 
DEC = zeros(size(BORFLAG));
MAX = zeros(size(BORFLAG));
for j = 1:length(phi)
	iboring = find(BORFLAG(j,:));
	ikgmax = find(TTP(j,:)==max(TTP(j,:)),1,'last');
	iinc = iboring(end):ikgmax-1;
	idec = ikgmax+1:length(TTP(j,:));
	INC(j,iinc)=1;
	DEC(j,idec)=1;
	MAX(j,ikgmax)=1;
end

TTPmod = zeros(size(TTP));
TTPmod(DEC==1) = 1;
TTPmod(INC==1) = 2;
TTPmod(BORFLAG==1) = 3;
TTPmod(MAX==1)=4;
cmap = [eye(3);[1 1 1]];
imagesc(minmax(kvec),minmax(phi),TTPmod); hold on
colormap(cmap);
caxis([1 4]);
ylabel('\phi');
xlabel('k (net tumor shrinkage rate)');
set(gca,'Ydir','normal')

plot([0 0],[0 1],'m:','Linewidth',2);
text(0,1,'|','color','m','verticalalignment','bottom','horizontalalignment','center');
text(0,1,'\leftarrow initial growth ','color','m','verticalalignment','bottom','horizontalalignment','right');
text(0,1,' initial regression \rightarrow','color','m','verticalalignment','bottom','horizontalalignment','left');

text(1*(1+0.5),0.75,'t_{min}<0','color','w','fontsize',14,'fontweight','bold');
text(1*(1+1),.4,'\uparrow k \rightarrow \uparrow TTP','color','w','fontsize',14,'fontweight','bold','Rotation',30);
text(1*(1+1.5),.2,'\uparrow k \rightarrow \downarrow TTP','color','w','fontsize',14,'fontweight','bold');

%% interactive visualization
% next step, would be nice to either port this to R-shiny or using GINPUT,
% let user click around in parameter space 
% and time course plot of TUK with those parameters gets plotted, with some
% curves with slightly higher and lower k values... 

%% theoretical approximation
% a closely related and easier to analyse question is where in parameter
% space does $\partial t_{min} / \partial k < 0 $ ?
% if you derive the expression for $t_{min}$ (by setting V'(t)=0), you get
%
% $$ t_{min} = \frac{ln(\frac{k}{g})-ln(\frac{\phi}{1-\phi})}{k+g} $$
%
% differentiating this with respect to k and doing some algebra we get that the derivative < 0 whenever:
%
% $$ ln(\frac{\phi}{1-\phi}) < -1-\frac{g}{k}-ln(\frac{g}{k}) $$
%
% we plot this separatrix over the last plot as a red curve. It's not
% exactly the same but it gives you an idea...

% alpha = -1 - gri - log(gri);
% phicrit = exp(alpha)./(1+exp(alpha));
% hold on;
% h1 = plot(gri,phicrit,'m:','linewidth',2);

% theoretical boringness separatrix

%legend([h1;h2],{'\partial t_{min}/\partial{k} = 0','t_{min} = 0'});

%% closing thoughts
% Can we show analytically that given any initial distribution
% of x=ic50's $\rho(x)$ and at a constant exposure C, and step-like pharmacology, no drift 
% ie 
% $$ g_{net}(x) = g, x>C, =-k, x<C $$
% that our
% 'generalized' model
%
% $$ V(t) = \int_0^\infty \rho(x) e^{g_{net}(x) t} dx $$
%
% does or doesn't have the artifact problem, ie, 
% $\partial t_{min}/\partial C < 0$ in this case
%
% I think we showed somewhere (see Guanyu's slides) that in this case the
% model reduces to our two-exponential model but with $\phi$ depending on C:
%
% $$ \phi(c) = \int_C^\infty \rho(x) dx $$