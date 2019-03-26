function [Y] = simmodel1(p,T)
% SIMMODEL returns Y values corresponding to times T for parameter vector p
% for single exponential model

% Single Exponential
%$$ V(t) = V_0*e^{-kt} ] $$
% where:
% * $V_0$ is the initial volume
% * k>0 is the kill rate (actually -net "growth rate" on treatment) on all cells

P = num2cell(p); 
[V_0, k] = deal(P{:}); % our parameters

TT = T(1):1:T(end)';
Ntimes = length(TT);

% initialize solutions
Y(1,1) = V_0.*exp(-TT(1)*k);

for j = 1:length(TT)
    if TT(j)<=0
        Y(j,1) = V_0;
    end
    if TT(j) > 0 
    Y(j,1) = V_0 .* exp(-k*TT(j));
    % need to write something to account for when Y(j,1) <= 0...
    end
    
  
end

% return only predictions at time points originally passed in T
ikeep = find(ismember(TT,T));
Y = Y(ikeep,1);
end