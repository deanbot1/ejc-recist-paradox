function  [ pat ] = loadDataPar( cancer_type )
% load data and create structured array
switch cancer_type

    case 'colon'
        
        disp('Loading data from first-line treatment of metastatic colorectal carcinoma patients in control arm of Phase 3 study.');
        disp('Treatment: Ox, Lu, 5FU every 2 weeks')

    %load data (not needed for paper results)
    % N: numerical data; T: text fields
    [NO, TO] = xlsread('../data/OS_colon_data.xls','','','basic');
    [NP, TP] = xlsread('../data/PFS_colon_data.xls','','','basic');

    colP = @(s)(find(strcmp(s,TP(1,:))));
    rP = @(str)(find(strcmp(str,TP(2:end,colP('ENDRSON')))));
    rowsP=rP('DISEASE PROGRESSION/LACK OF EFFICACY');

    NP=NP(rowsP,:);

     N =readtable('../data/Les_colon_data.xls','TreatAsEmpty', {'NA', 'Na', 'na', 'none'});
   % only use rows that have 'PLACEBO' under TRTNAME
   ikeep = ismember(N.TRTNAME, 'PLACEBO');
    N1 = N(ikeep, :);


% Remove rows that don't have a time point or a tumor lesion measurement
    torv1=isnan(N1.MSDY);
    torv2 = isnan(N1.MSDIM1);
    torvall = torv1 | torv2;
    N1(torvall,:) = [];
    ID = N1.RSUBJID;
    uID = unique(ID);
    tname = N1.TRTNAME;
    tlesID = N1.TLESID;
    day = N1.MSDY;
    tlessize = N1. MSDIM1;

    num_pats = length(uID);
    for j = 1:num_pats
            igood = find(ismember(ID, uID(j)));
            pat(j).ID = uID(j);

            % to create new spreadsheet
            pat(j).IDlist = ID(igood);
            pat(j).drug = tname(igood);
            pat(j).time = day(igood);
            pat(j).size = tlessize(igood);

    end
    % Loop through OS and tack on if matching uID, else censor

    for j=1:num_pats
        % note 10 is specific to this data set (hard coded in)
        % this is location of RSUBJID
        ind = find(ismember(NO(:,10), pat(j).ID));
        if isempty(ind)
            pat(j).OSobs = NaN;
        else
            pat(j).OSobs = NO(ind, 12)*7;
        end
    end

    % Loop through PFS and tack on if matching uID, else censor
    for j=1:num_pats
        % note 15 is specific to this data set (hard coded in)
        % this is location of RSUBJID
        ind = find(ismember(NP(:,15), pat(j).ID));
        if isempty(ind)
            pat(j).TTPobs = NaN;
            pat(j).PFS = pat(j).OSobs;
        else
            pat(j).TTPobs = NP(ind, 11);
            pat(j).PFS = NP(ind,11);
        end
    end


    
    case 'head neck'
           disp('Loading data from treatment of advanced head and neck cancer in control arm of Phase 3 study.');
            disp('Treatment: Radiation of 70 Gy in 35 fractions over 7 weeks plus 100 mg/m2 cisplatin once on weeks 1, 4, and 7')

         N =readtable('../data/headnecktrial_clinexa.xls','TreatAsEmpty', {'NA', 'Na', 'na', 'none'});
       % only use rows that have 'PLACEBO' under TRTNAME
       ikeep = ismember(N.ROWS, 'Y');
        N1 = N(ikeep, :);

    % Remove rows that don't have a time point or a tumor lesion measurement
        torv1=isnan(N1.CLDY);
        torv2 = isnan(N1.CLDIM);
        torvall = torv1 | torv2;
        N1(torvall,:) = [];
        ID = N1.RSUBJID;
        uID = unique(ID);
        tlesID = N1.TLESID;
        day = N1.CLDY;
        tlessize = N1. CLDIM;

        num_pats = length(uID);

        for j = 1:num_pats
             igood = find(ismember(ID, uID(j)));

            pat(j).ID = uID(j);

            % to create new spreadsheet
            pat(j).IDlist = ID(igood);
            pat(j).time = day(igood);
            pat(j).size = tlessize(igood);
        end

    case 'prostate'
        disp('Loading data from progressive metastatic castration-resistant prostate cancer patients after failure of a docetaxel-based chemotherapy regimen in control arm of Phase 3 study.');
        disp('Treatment: Prednisone orally daily at 5 mg PO BID')
           N =readtable('../data/prostatetrial_tmm_p.xls','TreatAsEmpty', {'NA', 'Na', 'na', 'none'});
       % only use rows that have 'PLACEBO' under TRTNAME
       ikeep = ismember(N.ROWS, 'Y');
        N1 = N(ikeep, :);

    % Remove rows that don't have a time point or a tumor lesion measurement
        torv1=isnan(N1.EFDAY);
        torv2 = isnan(N1.TMMDIA);
        torvall = torv1 | torv2;
        N1(torvall,:) = [];
        ID = N1.PID_A;
        uID = unique(ID);
        day = N1.EFDAY;
        tlessize = N1.TMMDIA;

        num_pats = length(uID);
        for j = 1:num_pats
                igood = find(ismember(ID, uID(j)));
                pat(j).ID = uID(j);

                % to create new spreadsheet
                pat(j).IDlist = ID(igood);
                pat(j).time = day(igood);
                pat(j).size = tlessize(igood);

        end       
        
  case 'ovarian1'
            N =readtable('../data/ovariantrial_aomd1.xls','TreatAsEmpty', {'NA', 'Na', 'na', 'none'});
        disp('Loading data from relapsed epithelial ovarian carcinoma patients in control arm of Phase 3 study.');
        disp ('Treatment:Paclitaxel HCl every 3 weeks');

     % only use rows that have 'PLACEBO' under TRTNAME
       ikeep = ismember(N.ROWS, 'Y');
        N1 = N(ikeep, :);

    % Remove rows that don't have a time point or a tumor lesion measurement
        torv1=isnan(N1.AOMDDY);
        torv2 = isnan(N1.LSNAREAN);
        torvall = torv1 | torv2;
        N1(torvall,:) = [];
        ID = N1.DSUBJID;
        uID = unique(ID);
        day = N1.AOMDDY;
        tlesarea = N1.LSNAREAN;
        tlessize = sqrt(tlesarea);
        num_pats = length(uID);
        for j = 1:num_pats
                igood = find(ismember(ID, uID(j)));
                pat(j).ID = uID(j);

                % to create new spreadsheet
                pat(j).IDlist = ID(igood);
                pat(j).time = day(igood);
                pat(j).size = tlessize(igood);

        end       

  case 'ovarian2'
      N =readtable('../data/ovariantrial_aomd2.xls','TreatAsEmpty', {'NA', 'Na', 'na', 'none'});
            disp('Loading data from relapsed epithelial ovarian carcinoma patients in control arm of Phase 3 study.');
            disp ('Treatment:Topotecan HCl every 3 weeks');
            % only use rows that have 'PLACEBO' under TRTNAME
            ikeep = ismember(N.ROWS, 'Y');
            N1 = N(ikeep, :);
            
            % Remove rows that don't have a time point or a tumor lesion measurement
            torv1=isnan(N1.AOMDDY);
            torv2 = isnan(N1.LSNAREAN);
            torvall = torv1 | torv2;
            N1(torvall,:) = [];
            ID = N1.DSUBJID;
            uID = unique(ID);
            day = N1.AOMDDY;
            tlesarea = N1.LSNAREAN;
            tlessize = sqrt(tlesarea);
            
            num_pats = length(uID);
            for j = 1:num_pats
                igood = find(ismember(ID, uID(j)));
                pat(j).ID = uID(j);
                
                % to create new spreadsheet
                pat(j).IDlist = ID(igood);
                pat(j).time = day(igood);
                pat(j).size = tlessize(igood);
                
            end
            
    otherwise
        disp('We do not have data for that cancer type')
end

end
