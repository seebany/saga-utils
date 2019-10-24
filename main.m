function [] = main(yearin, doyin, lronlyflag, prnlist_in, init_time_in, xtime_in)
%% usage
% yearin, doyin : required; i.e. [2015:2017], [1:365]
% lronlyflag [1/0] : required; 1, do low-rate only processing; 0, both hr and lr
% prnlist_in [a vector] : optional; specify a list of PRN to override detected event satellites
% init_time_in [datevec], xtime_in [s] : optional; i.e. [2015, 10, 7, 3, 6, 2, 0, 0], [200, 1200]
%% Initialization
close all;
dbstop if error;
init();
% maxNumCompThreads(4);

%file separator "/" in linux
global sep;
sep = filesep;

%where cases data are stored
global cases_folder;
cases_folder = '/data1/public/Data/cases/pfrr/';
%case_folder = '/data1/public/Data/cases/calg/';

%where output data and plots are stored
global home_dir;
global mat_dir;
[~, dummy] = system('echo $HOME');
home_dir = [dummy(1:end-1), sep];
 
if strcmp(cases_folder(end-4:end-1), 'pfrr')
    %path for Poker Flat data
    mat_dir = ['PFRR_Data', sep];
    mat_dir = ['matfiles', sep];%, num2str(yearin), sep, num2str(doyin), sep, '181011', sep];%datestr(date, 'yymmdd'),sep];
else
    %folder_path for 2013 Calgary data
    mat_dir = ['Calgary_Data', sep];
end

[~, op_path] = ver_chk();

%create parent directories
comm = strjoin({'mkdir -p', [home_dir, sep, mat_dir, sep]});
system(comm);

%find yesterday in doy
comm = 'date -u -d "a day ago" +%j';
[~, yesterday] = system(comm);
yesterday = str2num(yesterday);

comm1 = 'date -u -d "a day ago" +%Y';
[~, year] = system(comm1);

switch nargin
    % run for yesterday
    case {0}
        yearlist = str2num(year);
        doylist = yesterday;
        doyin = yesterday;
        lronlyflag = 1;
    % run for year and doy, with lronlyflag
    case {3}
        yearlist = yearin;
        doylist = doyin;
    % run for year and doy and one single specific prn w/o specific period
    case {4, 6}
        yearlist = yearin;
        doylist = doyin;
        % looks at high rate data, so lronlyflag should always be 0
        lronlyflag = 0;
        prnlist = prnlist_in;
    otherwise       
        error('Not enough inputs');
end

SD = [];
MEGA_MSP = [];
MSP_days = [];
SCINTEVENTS = [];
%specify signal: L1C1 = 0; L2CL = 2.
for signal_type = 0 %[0, 2]
    for yearnum = yearlist
        year = num2str(yearnum, '%04i')
        % for doy = yesterday
        for doynum = unique(doylist)
            %day of year in string
            doy = num2str(doynum, '%03i')
            year;
            
            %% Process low-rate data
            %first sigmaphi threshold in [rad]
            spmask = 0;
            s4mask = 0;
            
            %check if low rate data has already been processed
            matfile = dir([op_path, 'lrdata_', num2str(signal_type), '_', year, '_', doy, '.mat'])
            if ~isempty(matfile) && matfile.datenum >= datenum([2015, 5, 30, 0, 0, 0])
                disp([op_path, 'lrdata_', num2str(signal_type), '_', year, '_', doy, '.mat already exists']);
            else
                %filter low rate sigma phi data for all operational receivers
                %with sigmaphi threshold and save the results into a large array
                [MSP, MS4, HITDATA, CST, rcvr_op, tlim, splim, s4lim, signal] = scint_el_stackplot(cases_folder, home_dir, signal_type, doy, year, spmask, s4mask);
                lrdata = {MSP, HITDATA, CST, rcvr_op, tlim, splim, signal, spmask, MS4, s4mask, s4lim};
                if ~isempty(MSP) && ~isempty(MS4)
                    save([op_path, 'lrdata_', num2str(signal_type), '_', year, '_', doy, '.mat'], 'lrdata');
                else
                    disp('No data for this day, exiting....');
                    continue;
                end
            end
            
            %% Generate low-rate sigmaphi stackplots
            %load low rate results
            lrfile = [op_path, 'lrdata_', num2str(signal_type), '_', year, '_', doy, '.mat'];
            load(lrfile);
            tlim_vec = datevec(lrdata{5})
            hr_times = [];
            hrtimes = [];
            MSP = lrdata{1};
            HITDATA = lrdata{2};
            CST = lrdata{3};
            rcvr_op = lrdata{4};
            tlim = lrdata{5};
            splim = lrdata{6};
            signal = lrdata{7};
            spmask = lrdata{8};
            MS4 = lrdata{9};
            s4mask = lrdata{10};
            s4lim = lrdata{11};
            span = 60;
            
            disp(['There are totally ', num2str(size(MSP, 1)), ...
                ' points of sigma_phi data']);
            disp(['There are totally ', num2str(size(MS4, 1)), ...
                ' points of s4 data']);
            
            % receiver structure including all poker flat receivers
            rcvr_struct = rx_dirs(cases_folder, year, doy);
            rcvr_struct = ['grid108'; 'grid154'; 'grid160'; 'grid161'; 'grid162'; 'grid163'; 'ASTRArx'];
            % rcvr_struct = ['grid108';'grid154';'grid160';'grid161';'grid162';'grid163'];
            
            if ~isempty(MSP) && lronlyflag
                lr_rx_stackplot(lrdata, year, doy, rcvr_struct, op_path, 'sp');
                lr_prn_stackplot(lrdata, year, doy, rcvr_struct, op_path, 'sp');
            else % no enough scintillation data to conitue, considered a quiet day
                disp(['doy:', doy, ' of ', year, ' is a quiet day']);
            end
            
            if ~isempty(MS4) && lronlyflag
                lr_rx_stackplot(lrdata, year, doy, rcvr_struct, op_path, 's4');
                lr_prn_stackplot(lrdata, year, doy, rcvr_struct, op_path, 's4');
            else % no enough scintillation data to conitue, considered a quiet day
                disp(['doy:', doy, ' of ', year, ' is a quiet day']);
            end
            
            %         exit;
            disp('Quick look plots generated.');
            pause(2);
            if lronlyflag
                continue;
            end
            
            
            %save MSP for each day into MEGA_MSP
            % datevec(MSP([1 end],1))
            if strcmp(year, '2014') == 1 || strcmp(year, '2013') == 1
                MSP = MSP(MSP(:, 1) <= datenum([2015, 1, 1, 0, 0, 0]) & MSP(:, 1) >= datenum([2013, 12, 1, 0, 0, 0]), :);
            end
            
            MSP_hours = MSP(:, 1) - str2num(doy) + 1;
            % MSP_hours = MSP(:,1);
            % datevec(MSP_hours([1 end],:))
            MEGA_MSP = [MEGA_MSP; [MSP_hours, MSP(:, 2:end)]];
            % datevec(MEGA_MSP([1 end],1))
            
            %% Discard scintillation values below certain thresholds
            %daily mean
            meansp_doy = mean(MSP(:, 2));
            stdsp_doy = std(MSP(:, 2));
            means4_doy = mean(MS4(:, 2));
            stds4_doy = std(MS4(:, 2));
            %fixed threshold
            sp_fixed = 0.6; % Yang had as 0.6; I'm temporarily lowering to 0 just to get Vaisnavi's requested HR plot.
            s4_fixed = 0.2; % Yang had as 0.4; Vaishnavi lowered to 0.2
            
            RCVR_OP = [];
            MSP_NUM = [];
            MSP_NUM0 = [];
            for rr = 1:size(rcvr_op, 1)
                if ~isempty(MSP(MSP(:, 4) == rr, :))
                    MSP_NUM = [MSP_NUM; size(MSP(MSP(:, 4) == rr, :), 1)];
                    RCVR_OP = [RCVR_OP; rcvr_op(rr, :)];
                end
                MSP_NUM0 = [MSP_NUM0; size(MSP(MSP(:, 4) == rr, :), 1)];
            end
            
            RCVR_OP
            rcvr_op
            % There are two thresholds used to compute the weighted
            % scintillation number WSN.  The two thresholds are spth_hr0
            % and spth_hr1. Here the first one is set to a fixed value and
            % the second is set to the mean value for that day.
            spth_hr0 = sp_fixed;
            MSP_hr0 = MSP(MSP(:, 2) >= spth_hr0, :);
            MSP_hr0 = sortrows(MSP_hr0, 1);
            spth_hr1 = meansp_doy;
            MSP_hr1 = MSP(MSP(:, 2) >= spth_hr1, :);
            MSP_hr1 = sortrows(MSP_hr1, 1);
            
            MSP_hr = floor([size(MSP_hr0, 1); size(MSP_hr1, 1)]/size(RCVR_OP, 1));
            spth_hr = [spth_hr0; spth_hr1];
            WSN = floor((MSP_hr' * spth_hr)/sum(spth_hr));
            
            disp(['There are totally ', num2str(size(MSP, 1)), ' points for this day']);
            disp([num2str(size(MSP_hr0, 1)), ...
                ' points of sigma_phi data exceeding the fixed threshold ', num2str(spth_hr0)]);
            disp([num2str(size(MSP_hr1, 1)), ...
                ' points of sigma_phi data exceeding the daily average threshold ', num2str(spth_hr1)]);
            disp(['Weighted Scintillation Number is ', num2str(WSN)]);
            
            sd_num = [yearnum, doynum, ...
                size(rcvr_op, 1), floor(mean(MSP_NUM0)), size(RCVR_OP, 1), floor(mean(MSP_NUM)), ...
                spth_hr(1), MSP_hr(1), spth_hr(2), MSP_hr(2), WSN]
            SD = [SD; sd_num];
            ONES = ones(size(MSP, 1), 1);
            MSP_day = [doynum * ONES, length(unique(MSP(:, 4))) * ONES, MSP(:, 2:end)];
            MSP_days = [MSP_days; MSP_day];
            fid1 = fopen([op_path, 'sd.txt'], 'a+');
            fprintf(fid1, '%d %03d %d %d %d %d %.3f %d %.3f %d %d\n', sd_num');
            % continue;
            
            %% Generate a list of time intervals used for downloading high-rate data
            dt = 3600 / 24 / 3600 / 2;
            if spth_hr1 <= spth_hr0
                MSP_hr = MSP_hr0;
                spth_hr = spth_hr0
            else
                MSP_hr = MSP_hr1;
                spth_hr = spth_hr1
            end
            
            %spth_hr is always the larger one of the two
            
            lrtimesfile = [op_path, 'lrtimes', '_', year, '_', doy, '.mat']
            if isempty(dir(lrtimesfile))
                disp([lrtimesfile, 'doesn not exist']);
                [mega_t, TSP_hr0, TSP_hrv0] = find_general_times(MSP, rcvr_op, spth_hr);
                TSP_hrv0
                save(lrtimesfile, 'MSP', 'rcvr_op', 'spth_hr', 'TSP_hr0', 'TSP_hrv0');
            else
                load(lrtimesfile);
                TSP_hrv0
            end
            
            %         keyboard;
            e_common = findhrtimes(year, doy);
            t0 = datevec(e_common(:, 2));
            tf = datevec(e_common(:, 3));
            e_common0 = [e_common(:, 1), t0(:, 4), t0(:, 5), ...
                tf(:, 4), tf(:, 5), e_common(:, end)]
            %         keyboard;
            TSP_hr0_short = TSP_hr0(1:min(10,size(TSP_hr0,1)), :);
            TSP_hrv0_short = TSP_hrv0(1:min(10,size(TSP_hr0,1)), :);
            TSP_hr = TSP_hr0;
            
            %         keyboard;
            % continue;
            
            %%
            %receiver structure of high rate data, different from that of low rate as
            %receivers are in site ID order instead
            rcvr_struct_hr = ['grid108'; ...
                'grid163'; ...
                'grid160'; ...
                'grid162'; ...
                'grid161'; ...
                'grid154'];
            rcvr_op_hr = [];
            for rr = 1:size(rcvr_struct_hr, 1)
                rcvr_name = rcvr_struct_hr(rr, :);
                if ismember(rcvr_name, RCVR_OP, 'rows')
                    [~, rr_op] = ismember(rcvr_name, RCVR_OP, 'rows');
                    rcvr_op_hr = [rcvr_op_hr; RCVR_OP(rr_op, :)];
                end
            end
            
            if yearnum == 2013 && doynum == 342
                rcvr_op_hr = ['grid108'; ...
                    'grid163'; ...
                    'grid162'; ...
                    'ASTRArx'; ...
                    'grid161';]
            end
            
            rcvr_op_hr;
            %         continue;
            
            %% High-rate processing w/{w/o} specified PRNs or time intervals
            flag = 'single';
            %         flag = 'multiple';
            if ~exist('prnlist', 'var') 
                if strcmp(flag, 'single')
                    prnlist = unique(TSP_hr(:, 1), 'stable');
                    prnlist = prnlist(1);
                else
                    prnlist = e_common(:, 1);
                    %             t_common = unique(e_common(:,2:3),'rows','stable');
                    t_common = e_common(:, 2:3);
                end
            end
            %         keyboard;
            switch doy
                %             case '342'
                %                 prnlist = [23,10,13];
                %     case '050'
                %         prnlist = [17];
%                case '051'
%                    prnlist = [29];
                 case '076'
                    % prnlist = 27;%18;%[19, 18, 22, 27];
                    %         prnlist = [27];
                    %     case '077'
                    %         prnlist = [25 29 31];
                    %         prnlist = 1;                    
%                 case '280'
%                     prnlist = [3];
% %                 case '326'
% %                     prnlist = [14];
% %                 case '049'
% %                     prnlist = [5];
            end
            length(prnlist)
            for kk = 1:min(3, length(prnlist))
                %  for kk = 1
                % for kk = 1:min(length(prnlist),5)
                prn = prnlist(kk)
                if strcmp(flag, 'single')
                    trow = find(TSP_hr(:, 1) == prn);
                else
                    trow = find(prnlist == prn);
                end
                length(trow)
                % only take the first ranked event for a certain satellite,
                % uncomment to get all events
                for irow = 1 %:length(trow)
                    if strcmp(flag, 'single')
                        tt = TSP_hr(trow(irow), 2:3)';
                        duration = TSP_hr(trow(irow), 4);
                        sp_median = TSP_hr(trow(irow), 5);
                    else
                        tt = t_common(trow(irow), :)';
                        duration = diff(tt) * 24 * 60;
                    end
                    
                    %         datevec(tt)
                    %         if (isempty(t_old) && diff(tt)<=dt)
                    %             tt_new(1,:) = tt(1)-(dt-diff(tt))/2;
                    %             tt_new(2,:) = tt(2)+(dt-diff(tt))/2;
                    %         else
                    %             tt_new(1,:) = tt(1)-300/24/3600;
                    %             tt_new(2,:) = tt(2)+300/24/3600;
                    %         end
                    %         tt = tt_new;
                    datevec(tt)
                    %                                 keyboard;
                    init_time = datevec(tt(1));
                    init_time = datenum([init_time(1:4), 0, 0])
                    xtime = (tt - init_time) * 24 * 3600;
                    
                    if exist('init_time_in', 'var')
                        init_time = datenum(init_time_in);
                    end
                    if exist('xtime_in', 'var')
                        xtime = xtime_in;
                    end
                    %
                    %         specify times of interest used in the case study
                    switch prn
                        %             case {5,25,29}
                        %                 init_time = datenum([2015 3 17 9 0 0]);
                        %                 xtime = [-1800;1800];
                                     case {18,27,22}
                                         init_time = datenum([2015 3 17 13 0 0]);
                                         xtime = [60*8;60*25];
%                         case {19}
%                             init_time = datenum([2015 3 17 15 0 0]);
%                             xtime = [60*0; 60*30];
%                         case {14}
%                             init_time = datenum([2014 11 22 22 0 0]);
%                             xtime = [60*29; 60*49];
%                         case {5}
%                             init_time = datenum([2015 2 18  11 0 0]);
%                             xtime = [60*17; 60*29];
                        %             case {21}
                        %                 init_time = datenum([2015 3 17 11 0 0]);
                        %                 xtime = [0;3600*3];
                        %                                 case {10,13}
                        %                                     init_time = datenum([2013 12 8 3 0 0]);
                        %                                     xtime = [2604.4;4648.9];
                        %             case {31,29,25}
                        %                 init_time = datenum([2015 3 18 7 50 0]);
                        %                 xtime = [0;2400];
                        %             case 1
                        %                 init_time = datenum([2015 3 18 15 2 0]);
                        %                 xtime = [-600;600];
                        case {3}
                            init_time = datenum([2015 10 7 6 0 0]);
                            xtime = [120; 250];
%                         case {6}
                            
                    end
                    
                    %                         keyboard;
                    %         xtime = [3600;2*3600+45*60];
                    %         xtime = [3600+1800;2*3600];
                    tt = init_time + xtime / 24 / 3600;
                    disp(['Time: ', num2str(xtime(1)), '-', num2str(xtime(end)), ...
                        ' seconds after ', datestr(init_time, 'yyyy/mm/dd HH:MM:SS')]);
                    %         lr_hr_sigmaphi_plot(tspan,prn,home_dir,cases_folder,rcvr_op_hr,year,doy,sep,signal,signal_type);
                    %         dt = 3600/24/3600;
                    dt = diff(tt) * 1.1;
                    % test with all 5 receivers including ASTRA
                    %         disp('All 5 receivers including ASTRA');
                    [tstt, tend] = zoomin_hrplot_iq_carrier(irow, tt, dt, signal_type, home_dir, cases_folder, year, doy, prn, rcvr_op_hr);
                    %
                    try
                        SCINTEVENTS = [SCINTEVENTS; ...
                            [yearnum, doynum, prn, tstt(4:5), tend(4:5), duration, veststats']]
                    catch
                        SCINTEVENTS = [SCINTEVENTS; ...
                            [yearnum, doynum, prn, tstt(4:5), tend(4:5), duration]]
                        
                    end
                end
            end
            % plotprnvs(prnlist,year,doy);
        end
        
        length(doyin)
        if length(doyin) >= 365 && ~isempty(SD) && ~isempty(MEGA_MSP) && ~isempty(MSP_days)
            save([op_path, 'sd_', year, '.mat'], 'SD', 'MEGA_MSP', 'MSP_days', 'SCINTEVENTS');
        end
        
        csvwrite([op_path, 'ScintillationEvents_short_', year, '.csv'], SCINTEVENTS);
        
    end
end
% fclose(fid1);
