function [pat, list] = filterDataPar(pat)
% load, filter, and sum tumor lesion data

%S = load('../out/pat.mat');
%pat = S.pat;

% add a cell for each patient to record
% V = Sum of the longest dimensions at each time
% cens_int = number of 0s in tumor lesions sizes * 0.15 cm 
%            this gives interval to integrate over in Bayes fit
% tmeas = time of measurement

list = [];
LLOQ = 1.5; % note LLOQ is based on minimum size of lesion detected on CT scan
ct_zeroV = 0;
for j = 1:length(pat)

       if length(unique(pat(j).time))>=0
           list = [list j];
       end
       
       % find unique time values
        pat(j).tmeas = unique(pat(j).time);
        tbl = tabulate(pat(j).time);
        pat(j).num_meas = max(tbl(:,2));
        % at each unique time value, calculate sum of size
        % if volume is 0, want to calculate the censored interval
        % cens_int = num-zeros * LLOQ
        for i = 1:length(pat(j).tmeas)
            iV=find((pat(j).time == pat(j).tmeas(i)));
            pat(j).V(i) = sum(pat(j).size(iV));
            pat(j).cens_int(i) = LLOQ.*sum(pat(j).V == 0);
            % if want to censor all zeros
%             pat(j).cens_int(i) = LLOQ.*sum(pat(j).size(iV) == 0);
        end
            if sum(pat(j).V==0)~=0
                ct_zeroV = ct_zeroV +1;
            end
        pat(j).Vminobs = min(pat(j).V)./pat(j).V(1);
        pat(j).bestrespobs = (min(pat(j).V) - pat(j).V(1))./pat(j).V(1);
end

end
