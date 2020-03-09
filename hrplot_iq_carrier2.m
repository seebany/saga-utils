function [tstt, tend] = hrplot_iq_carrier2(signal_type, home_dir, cases_folder, year, doy, ...
    prn, tspan_d, rcvr_op, zcounter, set_plot, fluct)

close all
sep = filesep;
% try
%     load('../../local.mat');
%     [~, op_path] = ver_chk;
% catch
signal = getSignal(signal_type);

% High-rate data processing and plotting

% tspan_d(1,:) = tspan_d(1,:)-600/24/3600;
% tspan_d(2,:) = tspan_d(2,:)+600/24/3600;

tspan_utc = datevec(tspan_d);
tstt = tspan_utc(1, :);
tend = tspan_utc(2, :);
hour = tstt(4);
init_time = datenum([tspan_utc(1, 1:3), hour, 0, 0]);
tlim = (tspan_d' - init_time') * 24 * 3600;
RCVRNAME = {};
[sitenum_op] = rx2site(rcvr_op);

global mat_dir;

hr_results = [home_dir, sep, mat_dir, sep, ...
    'hrplot_', year, '_', doy, '_PRN', num2str(prn), datestr(tspan_d(1, :), '_HHMMUT'), '_zoom', num2str(zcounter), '.mat']
tic;
% keyboard;
% if isempty(dir(hr_results))
if 1
    for rr = 1:size(rcvr_op, 1)
        rcvr_name = rcvr_op(rr, :);
        sitenum = sitenum_op{rr, :};
        %     if strcmp(case_folder(end-4:end-1),'pfrr')
        %         %folder_path for 2013 Poker Flat data
        %         op_path = strcat([home_dir,'PFRR_Data/',rcvr_name,sep,year,sep,doy,sep]);
        %         in_path = strcat([case_folder,rcvr_name,sep,year,sep,doy,sep]);
        %     else
        %         %folder_path for 2013 Calgary data
        %         op_path = strcat([home_dir,'Calgary_Data/',rcvr_name,sep,doy,sep]);
        %         in_path = strcat([case_folder,rcvr_name,sep,doy,sep]);
        %     end
        [~, op_path] = inoutpath(cases_folder, home_dir, year, doy, rcvr_name);
        %     if strcmp(rcvr_name,'ASTRArx')
        %         keyboard;
        %     end
        
        
        infilename = strcat(op_path, 'hr_prn_files_', signal, sep, ...
            'PRN', num2str(prn), datestr(tstt, '_HHMMUT'), '_zoom', num2str(zcounter), '.mat')
        if isempty(dir(infilename))
            disp('No existing data stored; first time reading data for this particular time period')
            DATAM2 = Fn_ReadHighRate_CASESdata_sdb(prn, op_path, cases_folder, rcvr_name, signal_type, tstt, tend);
            if ~isempty(DATAM2)
                save(infilename, 'DATAM2');
            end
        else
            load(infilename);
        end
        
        %         if rcvr_name == 'ASTRArx'
        %             keyboard;
        %         end
        rcvr_name
        size(DATAM2)
        %     DATAM2([1 end],:)
        %         keyboard;
        if ~isempty(DATAM2)
            RCVRNAME = [RCVRNAME; sitenum];
            if strcmp(set_plot, 'A') == 1
                outfilename = [op_path, 'hr_prn_files_', signal, sep, 'FilteredData_PRN', num2str(prn), ...
                    datestr(tstt, '_HHMMUT'), '_zoom', num2str(zcounter), '.mat']
            else
                outfilename = [op_path, 'hr_prn_files_', signal, sep, 'HR_Scintdata_PRN', num2str(prn), ...
                    datestr(tstt, '_HHMMUT'), '_zoom', num2str(zcounter), '.mat']
            end
            if isempty(dir(outfilename))
                disp('No existing data stored; first time processing data for this particular time period');
                Fn_Plot_HighRate_CASESdata_sdb(prn, tstt, init_time, op_path, signal_type, set_plot, infilename, zcounter);
            end
            load(outfilename);
            data_PRN = data_PRN';
            size(data_PRN);
            data_PRN([1, end], :);
            
            %specify time interval
            if 1 < 0
            elseif strcmp(year, '2017') && strcmp(doy, '233')
                ttt = data_PRN([1, end], 1);
                ttt = [20 * 60; 40 * 60];
                                 %elseif prn == 27 && strcmp(year, '2015') && strcmp(doy, '076')
                %         ttt = data_PRN([1 end],1);
                         %ttt = [693;753];
                %     elseif prn == 22 && strcmp(year,'2015') && strcmp(doy,'076')
                %         ttt = [600;900];
%                      elseif prn == 19 && strcmp(year,'2015') && strcmp(doy,'076')
%                          ttt = [183;214];
%              elseif (prn == 23 || prn == 10 || prn == 13) && strcmp(year,'2013') && strcmp(doy,'342')
%                  ttt = [2615;2660]; %CASE y1  (after 3UT)
%              elseif (prn == 23 || prn == 10 || prn == 13) && strcmp(year,'2013') && strcmp(doy,'342')
%                  ttt = [2625;2660]; %CASE y1 SHORTER  (after 3UT)
%              elseif (prn == 32) && strcmp(year,'2014') && strcmp(doy,'320')
%                 ttt = [1078;1108];  %Pralay case 3b (after 1UT)
%              elseif (prn == 32) && strcmp(year,'2014') && strcmp(doy,'320')
%                  ttt = [1030;1071];  %Pralay case 3A (after 1UT)
%                elseif (prn == 32) && strcmp(year,'2014') && strcmp(doy,'320')
%                    ttt = [975;1011];  % TEST ERROR L,Z(after 1UT)
%              elseif prn == 29 && strcmp(year,'2014') && strcmp(doy,'051')
%                          ttt = [2685;2729];
%                          ttt = [656;697];
%                          ttt = [535;570];
            else
                ttt = data_PRN([1, end], 1);
            end
            % %             override the time limits from low rate detection
            tlim = ttt';
            data_PRN = data_PRN(data_PRN(:, 1) <= ttt(end) & data_PRN(:, 1) >= ttt(1), :);
            if size(data_PRN, 1) == 0
                data_PRN = NaN * ones(1, 4);
            end
            obstime = data_PRN(:, 1);
            %1/12/2015 make time labels in :MM:SS format
            obstime = obstime / 24 / 3600 + init_time;
            piqpowdata = data_PRN(:, 2);
            piqphdata = data_PRN(:, 3);
            
            dt = max(unique(diff(obstime)));
            dt = 0.015 / 24 / 3600;
            [obstime_e, power_e] = discont_proc(obstime, piqpowdata, dt);
            maxpwr(rr) = max(abs(power_e))
            [obstime_e, phase_e] = discont_proc(obstime, piqphdata, dt);
            maxph(rr) = max(abs(phase_e))
            [color] = rx_color(rcvr_name);
            %high-rate s4 and sigmaphi
            if strcmp(set_plot, 'B') == 1
                s4 = power_e;
                sp = phase_e;
                subplot(2, 1, 1)
                plot(obstime_e, s4, 'Color', color, 'Linewidth', 0.5);
                grid on;
                str = strcat('100Hz S_4 and', ...
                    {' \sigma_{\Phi} for '}, signal, ', PRN:', num2str(prn));
                title(str);
                axis([tlim / 24 / 3600 + init_time, 0, 1]);
                ylabel('(a) S_4');
                if diff(tlim) <= 300
                    ticklbl = 'HH:MM:SS';
                else
                    ticklbl = 'HH:MM';
                end
                datetick('x', ticklbl, 'keeplimits');
                hold on;
                subplot(2, 1, 2)
                plot(obstime_e, sp, 'Color', color, 'Linewidth', 0.5);
                grid on
                ylabel('(b) \sigma_{\Phi} [rad]')
                axis([tlim / 24 / 3600 + init_time, 0, 2 * pi]);
                hold on;
                datetick('x', ticklbl, 'keeplimits');
                lg = legend(gca, RCVRNAME, 'Location', ...
                    'north', 'Orientation', 'horizontal');
                lgpos = get(lg, 'Position');
                lg = legend(gca, RCVRNAME, 'Position', ...
                    [lgpos(1), 0.5, lgpos(3:4)], 'Orientation', 'horizontal');
                set(lg, 'FontSize', 8);
                lgpos = get(lg, 'Position');
                xstring = ['Time ', num2str(tlim(1)), '-', num2str(tlim(2)), ...
                    '[s] after ', num2str(hour), ...
                    ':00 UT on: ', datestr(tstt, 'mm/dd/yy')];
                %1/12
                xstring = ['Time [', ticklbl, '] on: ', datestr(tstt, 'mm/dd/yy')];
                %
                xlabel(xstring);
                xdata_PRN{rr} = data_PRN;
            end
            %high-rate filtered power and phase
            if strcmp(set_plot, 'A') == 1
                %1/12
                if diff(tlim) <= 300
                    ticklbl = 'HH:MM:SS';
                    rotang = 0;
                else
                    ticklbl = 'HH:MM';
                    rotang = 25;
                end
                
                subplot(2, 1, 1);
                
                %no db:
                %plot(obstime_e, (power_e), 'Color', color, 'Linewidth', 0.25);
                plot(obstime_e, 10*log10(power_e), 'Color', color, 'Linewidth', 0.25);
                set(gca, 'xticklabelrotation', rotang);
                str = strcat('Detrended Power $P_{f}$ and', ...
                    {' Phase $\Phi_f$ for '}, signal, ', PRN:', num2str(prn));
                title(str);
                %         axis([tlim/24/3600+init_time -log10(max(maxpwr))*10*1.5 log10(max(maxpwr))*10*1.5]);
                axis([tlim / 24 / 3600 + init_time, -10, 5]);
                datetick('x', ticklbl, 'keeplimits');
                %         set(gca,'XTick',(ttt(1):10:ttt(2))/24/3600+init_time);
                ylabel('(a) Power $P_f$ [dB]');
                
                
                hold on;
                
                %         legend(RCVRNAME,'Location','NorthEastOutside');
                if 2 < 1
                    %         if rr == size(rcvr_op,1) && strcmp(doy,'342')
                    %1/12/2015 red line snapshot times
                    tpp_s = [2641.789, 2654.335, 2666.819, 2679.305];
                    for tpp = tpp_s / 24 / 3600 + init_time
                        plot([tpp, tpp], [-log10(max(maxpwr)) * 10 * 1.5, log10(max(maxpwr)) * 10 * 1.5], 'color', [0.5, 0.5, 0.5], 'Linewidth', 0.5);
                        plot([tpp + 1 / 24 / 3600, tpp + 1 / 24 / 3600], [-log10(max(maxpwr)) * 10 * 1.5, log10(max(maxpwr)) * 10 * 1.5], 'color', [0.5, 0.5, 0.5], 'Linewidth', 0.5);
                        hold on;
                    end
                    %             %1/21/2015 green line snapshot times
                    %             for tpp = [2634,2646,2659,2671,2684]/24/3600+init_time
                    %                 plot([tpp,tpp],[-log10(max(maxpwr))*10*1.5 log10(max(maxpwr))*10*1.5],'color',[0.5 0.5 0.5],'Linewidth',0.5);
                    %                 plot([tpp+1/24/3600,tpp+1/24/3600],[-log10(max(maxpwr))*10*1.5 log10(max(maxpwr))*10*1.5],'color',[0.5 0.5 0.5],'Linewidth',0.5);
                    %                 hold on;
                    %             end
                end
                phasesp = subplot(2, 1, 2);
                h(rr) = plot(gca, obstime_e, phase_e, 'Color', color, 'Linewidth', 0.25);
                set(gca, 'xticklabelrotation', rotang);
                ylabel('Phase $\Phi_f$ [rad]');
                hold on;
                %phase peaks
                %                 [pks, locs, width , prominence] = ...
                %                     findpeaks(piqphdata, obstime, 'MinPeakHeight', 1.5, ...
                %                     'Annotate','extents');
                
                %         axis([tlim/24/3600+init_time -max(maxph)*1.5 max(maxph)*1.5]);
                axis([tlim / 24 / 3600 + init_time, -2 * pi, 2 * pi]);
                datetick('x', ticklbl, 'keeplimits');
                %         set(gca,'XTick',(ttt(1):10:ttt(2))/24/3600+init_time);
                
                if rr == size(rcvr_op, 1) && strcmp(doy, '342') && 2 < 1
                    %1/12/2015 red line snapshot times
                    for tpp = tpp_s / 24 / 3600 + init_time
                        plot([tpp, tpp], [-max(maxph) * 1.5, max(maxph) * 1.5], 'color', [0.5, 0.5, 0.5], 'Linewidth', 0.5);
                        plot([tpp + 1 / 24 / 3600, tpp + 1 / 24 / 3600], [-max(maxph) * 1.5, max(maxph) * 1.5], 'color', [0.5, 0.5, 0.5], 'Linewidth', 0.5);
                        text(tpp, max(maxph)*1.5, datestr(round((tpp - init_time)*24*3600)/24/3600+init_time, 'HH:MM:SS'), ...
                            'VerticalAlignment', 'Bottom', 'color', 'k');
                        hold on;
                    end
                    %             %1/21/2015 green line snapshot times
                    %             for tpp = [2634,2646,2659,2671,2684]/24/3600+init_time
                    %                 plot([tpp,tpp],[-max(maxph)*1.5 max(maxph)*1.5],'color',[0.5 0.5 0.5],'Linewidth',0.5);
                    %                 plot([tpp+1/24/3600,tpp+1/24/3600],[-max(maxph)*1.5 max(maxph)*1.5],'color',[0.5 0.5 0.5],'Linewidth',0.5);
                    %                 text(tpp,max(maxph)*1.5,datestr(tpp,':MM:SS'),...
                    %                     'VerticalAlignment','Top','color','g');
                    %                 hold on;
                    %             end
                end
                
                %         blanks = repmat({' '},9,1);
                %         lb1 = [2620;blanks;2630;blanks;
                %             2640;blanks;2650;blanks;
                %             2660;blanks;2670;blanks;
                %             2680;blanks;2690];
                %         set(gca,'XTick',ttt(1):ttt(2),'XTickLabel',lb1);
                %     subplot(3,1,3)
                %         plot(obstime_e,carrier_e,'Color',color,'Linewidth',0.8);
                %         grid on
                %         ylabel('Carrier Phase w/o IQ [r]');
                %         axis([tlim -max(maxcph)*1.5 max(maxcph)*1.5]);
                %         hold on;
                %         legend(RCVRNAME,'Location','NorthEastOutside');
                xstring = ['Time ', num2str(tlim(1)), '-', num2str(tlim(2)), ...
                    '[s] after ', num2str(hour), ...
                    ':00 UT on: ', datestr(tstt, 'mm/dd/yy')];
                %1/12
                xstring = ['Time [', ticklbl, ' UT] on: ', datestr(tstt, 'mm/dd/yy')];
                %
                xlabel(xstring);
                %         set(gca,'XTick',ttt(1):ttt(2),'XTickLabel',lb1);
                
                %     disc = find(diff(obstime)>dt);
                %     if size(obstime,1)<=2
                %         xtdata{rr} = sortrows([obstime(1);obstime(end)]);
                %     else
                %         xtdata{rr} = sogrtrows([obstime(1);obstime(disc);obstime(disc+1);obstime(end)]);
                %     end
                %     xtdata{rr} = xtdata{rr}';
                xdata_PRN{rr} = data_PRN;
            end
        else
            xdata_PRN{rr} = [];
        end
    end
    gcabottom = get(phasesp, 'outerposition');
   lg = legend(h, RCVRNAME, 'Location', ...
       'north', 'Orientation', 'horizontal');
   lgpos = get(lg, 'position');
   lg = legend(h, RCVRNAME, 'Position', ...
       [lgpos(1) + 0.001, gcabottom(2) + gcabottom(4), lgpos(3:4)], 'Orientation', 'horizontal');
    
    %save the plot
    if strcmp(set_plot, 'A') == 1
        prefix = 'IQ_';
    elseif strcmp(set_plot, 'B') == 1
        prefix = 'S4SP_';
    end
    % plot_name = [signal,'_CorrIQ&CorrCarrPh_PRN',num2str(prn),...
    plot_name = [prefix, signal, '_PRN', num2str(prn), ...
        '_', year, '_', doy, ...
        '_zoom', num2str(zcounter), ...
        '_', num2str(tlim(1)), '-', num2str(tlim(2)), ...
        's_after_', datestr(init_time, 'HHMM'), 'UT'];
    %     plotpath = [op_path, plot_name,'.eps'];
    %     plotpath = [op_path, plot_name, '_phonly','.eps'];
    %     saveas(gcf, plotpath, 'epsc2');
    [~, hr_path, ~] = ver_chk;
    % save high_rate data as matfiles
    save([hr_path, plot_name, '.mat'], 'xdata_PRN', 'init_time', 'RCVRNAME', 'rcvr_op');
    
    plotpath = [hr_path, plot_name, '.png'];
    saveas(gcf, plotpath, 'png');
    close;
    
    init_t_utc = datevec(init_time);
    save(hr_results);
else
    load(hr_results);
end
init_t_utc = datevec(init_time);
hrplot_te = toc;

disp(['Finished preprocessing for PRN', num2str(prn)]);
% tic;
continuteflag = input('Proceed to estimation? [y]/n', 's');
if strcmp(continuteflag, 'n')
    return;
end

% Cross-correlation for pairs of receivers
% xcorr_results = [home_dir,'/PFRR_Data/','xcorr_',...
%     year,'_',doy,'_PRN',num2str(prn),'_zoom',num2str(zcounter),'.mat']
xcorr_results = [home_dir, sep, mat_dir, sep, 'xcorr_', ...
    year, '_', doy, '_PRN', num2str(prn), datestr(tspan_d(1, :), '_HHMMUT'), ...
    '_case', num2str(fluct),'_zoom', num2str(zcounter), '_60s.mat']
tic;

xdata_PRN
%make sure no empty set in the data to be cross-correlated
DRX = {};
rcvr_op_xcorr = rcvr_op
for iii = 1:size(xdata_PRN, 2)
    if ~isempty(xdata_PRN{iii})
        DRX = [DRX, xdata_PRN{iii}];
    else
        rcvr_op_xcorr = setdiff(rcvr_op_xcorr, rcvr_op(iii, :), 'rows');
    end
end
rcvr_op = rcvr_op_xcorr;
xdata_PRN = DRX

%check if there is data of only one/no receiver data, if
if isempty(rcvr_op) || size(rcvr_op, 1) == 1
    disp(['Caution! There is data of only 1 or no receiver available, ', ...
        'unable to perform cross-correlation']);
    return;
end

if ~isempty(dir(xcorr_results))
    fprintf('Estimation analysis results exist\n');
    %     plotSAGAvsPFISR(prn, tstt, 'debug',fluct);
    % plotprnvs(prn,year,doy);
    disp(['Finished processing for PRN', num2str(prn)]);
    rerunflag = input('Rerun the analysis? y/[n]', 's');
    if strcmp(rerunflag, 'y')
        renamecomm = strjoin({'mv', xcorr_results, [xcorr_results, '.bak']});
        system(renamecomm);
    else
        %plotSAGAvsPFISR(prn, tstt, 'debug');
        plotSAGAvsPFISR(prn, tstt, 'debug',fluct);
        return;
    end
end
% end

dt = 0.015;
%make all receiver data have equal length
for rr = 1:size(rcvr_op, 1)
    xdata_PRN{rr}(:, 1) = round(xdata_PRN{rr}(:, 1)*100) / 100;
    checkrows = find(diff(xdata_PRN{rr}(:, 1) == 0));
    uniquerows = setdiff(1:numel(xdata_PRN{rr}(:, 1)), checkrows);
    xdata_PRN{rr} = xdata_PRN{rr}(uniquerows, :);
    %the same with timestamps
    xt = xdata_PRN{rr}(:, 1);
    disc = find(diff(xt) >= dt);
    if size(xt, 1) <= 2
        xtdata{rr} = sortrows([xt(1); xt(end)]);
    else
        xtdata{rr} = sortrows([xt(1); xt(disc); xt(disc+1); xt(end)]);
    end
    xtdata{rr} = xtdata{rr}';
end

% return;

% Find continuous segment for all operational receivers
t = find_common_times(xtdata); 
% t([1 end])
% tlim([1 end])

%True combination of receiver pairs
combos_fig = nchoosek(1:size(rcvr_op, 1), 2);

combos = combos_fig;

v_ccmin = [0.6];
v_dtau = 60;
for i_dtau = 1:length(v_dtau)
    dtau = v_dtau(i_dtau);
    %[tslist, telist] = dividet_v2(t, dtau, 10);
    [tslist, telist] = dividet_v1(t, dtau, 10);
     if ((prn == 32) && strcmp(year,'2014') && strcmp(doy,'320'))==1 %Case 3b Pralay
         [tslist, telist] = dividet_v2(t, dtau, 10);
     elseif  isempty(tslist)
         [tslist, telist] = dividet_v2(t, dtau, 10);
     end
    %[tslist, telist] = dividet_v3(t, dtau*3/4, 10);
    [tslist, telist, telist - tslist]
    for tt = 1:length(telist)
        for rr = 1:size(rcvr_op, 1)
            time = xdata_PRN{rr}(:, 1);
            xdata{rr} = xdata_PRN{rr}(time <= telist(tt) & time >= tslist(tt), :);
            [~, ia, ~] = unique(xdata{rr}(:, 1), 'stable');
            %         size(c)
            xdata{rr} = xdata{rr}(ia, :);
            size(xdata{rr});
        end
        xtime = xdata{1}(:, 1);
        
        if fluct==1
             flagfluc='power';
        else
            flagfluc='phase';
        end
        fprintf('Begin estimation for period %i/%i \n', tt, length(telist));
        
       
        %read origin receiver location
        [~, op_path_0] = inoutpath(cases_folder, home_dir, year, doy, rcvr_op(1, :));
        load([op_path_0, 'prn_files_', signal, sep, 'navdata.mat']);
        NAVDATA_O = NAVDATA;
        [~, lla_0, ~] = compute_baselines(NAVDATA, NAVDATA_O, init_time, xtime);
        for rr = 1:size(rcvr_op, 1)
            rcvr_name = rcvr_op(rr, :);
            %         if strcmp(cases_folder(end-4:end-1),'pfrr')
            %             %folder_path for 2013 Poker Flat data
            %             op_path = strcat([home_dir,'PFRR_Data/',rcvr_name,sep,year,sep,doy,sep]);
            %         else
            %             %folder_path for 2013 Calgary data
            %             op_path = strcat([home_dir,'Calgary_Data/',rcvr_name,sep,doy,sep]);
            %         end
            [~, op_path] = inoutpath(cases_folder, home_dir, year, doy, rcvr_name);
            load([op_path, 'prn_files_', signal, sep, 'navdata.mat']);
            load([op_path, 'prn_files_', signal, sep, 'txinfodata.mat']);
            load([op_path, 'prn_files_', signal, sep, 'ionodata.mat']);
            load([op_path, 'prn_files_', signal, sep, 'scintdata.mat']);
            %         %plot STEC
            %         if ~isempty(IONODATA)
            %             ts_iono = (datenum(gps2utc(IONODATA(:,1:2)))-init_time)*24*3600;
            %             ind = ts_iono<=telist(tt) & ts_iono>=telist(tt) & IONODATA(:,5)==prn;
            %             t_iono = datenum(gps2utc(IONODATA(ind,1:2)));
            %             STEC(:,rr) = IONODATA(ind,3);
            %         end
            %         %
            %         if strcmp(rcvr_op(rr,:),'grid108') == 1
            %             NAVDATA_O = NAVDATA;
            %         end
            [enu, lla, xyz] = ...
                compute_baselines(NAVDATA, NAVDATA_O, init_time, xtime);
  %                  [~, DATA_el,~] = read_data(doy,year,in_path,op_path,sep,signal_type,signal);
%                         plot(enu(1), enu(2), 'o', 'Color', rx_color(rcvr_name), 'MarkerFaceColor', rx_color(rcvr_name));
%                         text(enu(1), enu(2), [num2str(lla(1:2)', '%.6g'), ['\circ N'; '\circ E']], 'VerticalAlignment', 'Bottom');
%             hold on;
%             grid on;
            ENU(:, rr) = enu';
            XYZ(:, rr) = xyz';
            LLA(:, rr) = lla';
            if ~isempty(DATA_el)
                time_azel = (datenum(gps2utc(DATA_el(:, 1:2))) - init_time) * 24 * 3600;
                iii = time_azel <= telist(tt) & time_azel >= tslist(tt) & DATA_el(:, 4) == prn;
                t_azel = time_azel(iii, :);
                el = DATA_el(iii, 3) * pi / 180;
                az = DATA_el(iii, 5) * pi / 180;
                EL(tt, rr) = mean(DATA_el(iii, 3)) * pi / 180;
                ZE(tt, rr) = pi / 2 - EL(tt, rr);
                AZ(tt, rr) = mean(DATA_el(iii, 5)) * pi / 180;
            else
                EL(tt, rr) = NaN;
                ZE(tt, rr) = NaN;
                AZ(tt, rr) = NaN;
            end
            if ~isempty(DATA)
                time_scint = (datenum(gps2utc(DATA(:, 1:2))) - init_time) * 24 * 3600;
                iii = time_scint <= telist(tt) & time_scint >= tslist(tt) & DATA(:, 4) == prn;
                %                 SP(tt, rr) = DATA(iii, 3);
            else
                SP(tt, rr) = NaN;
            end
            if rr == 1
                hvec = 1e3 * (250 + 50 * randn(100, 1));
                [ipplat, ipplon] = deal(zeros(length(az), length(hvec)));
                [vge_prn, vgn_prn] = deal(zeros(length(az)-1, length(hvec)));
                for kkk = 1:length(hvec)
                    %mean ipp latitude should be the same for different ccmin
                    [ipplat(:, kkk), ipplon(:, kkk)] = sill(az, el, ...
                        lla_0(1), lla_0(2), lla_0(3), hvec(kkk));
                    ipp_xyz = wgslla2xyz(ipplat(:, kkk), ipplon(:, kkk), ...
                        hvec(kkk)*ones(size(ipplat(:, kkk))));
                    ipp_enu = xyz2enu_new(ipp_xyz, lla_0(1), lla_0(2), lla_0(3));
                    denu_ipp = diff(ipp_enu);
                    %             dist = sqrt(denu_ipp(:,1).^2 + denu_ipp(:,2).^2);
                    vge_prn(:, kkk) = denu_ipp(:, 1) ./ diff(t_azel);
                    vgn_prn(:, kkk) = denu_ipp(:, 2) ./ diff(t_azel);
                end
                vge_prnbar = mean(mean(vge_prn, 2));
                vgn_prnbar = mean(mean(vgn_prn, 2));
                evge_prn = mean(std(vge_prn, 0, 2));
                evgn_prn = mean(std(vgn_prn, 0, 2));
            end
        end
%                 axis([-100, 3500, -1200, 400]);
%                 legend(sitenum_op);
%                 title('Operational Receviers of SAGA at PFRR');
%                 xlabel(['Time after ', datestr(init_time, 'HH:MM'), ...
%                     ' UT ', datestr(init_time, 'mm/dd/yyyy')]);
%                 plotpath = [op_path, 'PRN', num2str(prn), '_PFRR_array', '.eps'];
%                 saveas(gcf, plotpath, 'epsc2');
%                 close;
        %             return;
        
         [peak, tpeak, altpeak] = ...
            dataxcorr_alt(sitenum_op, xdata, combos_fig, flagfluc);
        %save the plot
        title(['Cross-correlation over ', ...
            num2str(tslist(tt), '%.0f'), '-', num2str(telist(tt), '%.0f'), 's after ', ...
            datestr(init_time, 'HHMM'), 'UT']);
        plotname = ['PRN', num2str(prn), '_Lag_plot_', ...
            num2str(tslist(tt), '%.0f'), '-', num2str(telist(tt), '%.0f'), 's_after_', ...
            datestr(init_time, 'HHMM'), 'UT'];
        plotpath = [op_path, plotname, '.png']
        saveas(gcf, plotpath, 'png');
        close;
        Fs=round(100);
               % get time-shifted signals %EACH %ALL
         [ph_truncated, pwr_truncated] = ...
                    plot_shifted(xdata, tpeak, rcvr_op, sitenum_op, combos, Fs, fluct);
%         
         scale=dev(rcvr_op,ph_truncated,pwr_truncated,sitenum_op,fluct);
 %       scale=0.25;
        numsim=10; 
        [tauaarrn, taucarrn, ccvalarrn, ccerrarrn, ...
             tauaarr, taucarr, ccvalarr, ccerrarr, A_nu, phi_nu] = ...
             estimate_obs2(rcvr_op, xdata, combos, flagfluc,scale,numsim); 
        
        % Plot input and output of NOISE 
        plot_noisesignal=0;
        if plot_noisesignal==1
        for rr = 1:size(rcvr_op, 1)
            rcvr_name = rcvr_op(rr, :);
            obstime_e = xdata{rr}(:,1);
            power_e = xdata{rr}(:, 2);
            phase_e = xdata{rr}(:, 3);
            [color] = rx_color(rcvr_name);
            if diff(tlim) <= 300
                ticklbl = 'HH:MM:SS';
                rotang = 0;
            else
                ticklbl = 'HH:MM';
                rotang = 25;
            end
            subplot(2, 1, 1);
            for n=1:numsim
                scint_n = power_e.*exp(1i.*phase_e) + A_nu(:,rr,n).* exp(1i.* phi_nu(:,rr,n));
                plot(obstime_e/ 24 / 3600 + init_time, 10*log10(abs(scint_n)),'*g','MarkerSize',0.5);                
                hold on;
            end
            plot(obstime_e/ 24 / 3600 + init_time, 10*log10(power_e), 'Color', color, 'Linewidth', 0.25);
            set(gca, 'xticklabelrotation', rotang);
            str = strcat('Detrended Power $P_{f}$ and', ...
                {' Phase $\Phi_f$ for '}, signal, ', PRN:', num2str(prn));
            title(str);
            axis([obstime_e(1)/ 24 / 3600 + init_time, obstime_e(end)/ 24 / 3600 + init_time, -10, 5]);
            datetick('x', ticklbl, 'keeplimits');
            ylabel('(a) Power $P_f$ [dB]');
            hold on;
            phasesp = subplot(2, 1, 2);
             for n=1:numsim
                scint_n = power_e.*exp(1i.*phase_e) + A_nu(:,rr,n).* exp(1i.* phi_nu(:,rr,n));
                scint_n2 = angle(scint_n);
                ind = find(((phase_e>pi)&(scint_n2<0))|(phase_e-scint_n2>2));
                scint_n2(ind) = scint_n2(ind) + 2*pi;
                ind2 = find(((phase_e<(-pi))&(scint_n2>0))|(phase_e-scint_n2<-2));
                scint_n2(ind2) = scint_n2(ind2) - 2*pi;
                h2(rr)=plot(gca, obstime_e/ 24 / 3600 + init_time,(scint_n2),'*g', 'MarkerSize',0.5);
                hold on;
             end
            h(rr) = plot(gca, obstime_e/ 24 / 3600 + init_time, phase_e, 'Color', color, 'Linewidth', 0.25);
            set(gca, 'xticklabelrotation', rotang);
            ylabel('Phase $\Phi_f$ [rad]');
            hold on;
            axis([obstime_e(1)/ 24 / 3600 + init_time, obstime_e(end)/ 24 / 3600 + init_time, -2 * pi, 2 * pi]);
            datetick('x', ticklbl, 'keeplimits');
            xstring = ['Time [', ticklbl, ' UT] on: ', datestr(tstt, 'mm/dd/yy')];
            xlabel(xstring);
            gcabottom = get(phasesp, 'outerposition');
            lg = legend([h(rr) h2(rr)], [RCVRNAME(rr) 'With noise'], 'Location', ...
                'north', 'Orientation', 'horizontal');
            lgpos = get(lg, 'position');
            lg = legend([h(rr) h2(rr)], [RCVRNAME(rr) 'With noise'], 'Position', ...
                [lgpos(1) + 0.001, gcabottom(2) + gcabottom(4), lgpos(3:4)], 'Orientation', 'horizontal');
            plot_name = [signal, '_PRN', num2str(prn), 'RCVR', num2str(rr), ...
                num2str(tslist(tt), '%.0f'), '-', num2str(telist(tt), '%.0f'), 's_after_', ...
                datestr(init_time, 'HHMM'), 'UT_withNOISE'];
            [~, hr_path, ~] = ver_chk;
            plotpath = [hr_path, plot_name, '.png'];
            saveas(gcf, plotpath, 'png');
            close;
        end
        end      
        % Solve for drift velocity
        %     combos = nchoosek(1:size(rcvr_op,1), 2);
        [H, YN, Y, CCVALN, CCVAL, CCERR, COMBOS, RXILOC, RHO0] = deal([]);
        for i = 1:size(combos, 1)
            denu(:, i) = ENU(:, combos(i, 1)) - ENU(:, combos(i, 2));
            x_ij = denu(1, i);
            y_ij = denu(2, i);
            %Use SAGA phase measurements to solve system for velocity
            tau_an = tauaarrn{i, :};
            tau_cn = taucarrn{i, :};
            ccvaln = ccvalarrn{i, :};
            ccerr = ccerrarrn{i, :};
            
            tau_a = tauaarr{i, :};
            tau_c = taucarr{i, :};
            ccval = ccvalarr{i, :};
            
            N = size(unique(tau_cn));
            h = [ones(N) * x_ij^2, ...
                ones(N) * 2 * x_ij * y_ij, ...
                ones(N) * y_ij^2, ...
                2 * x_ij * unique(tau_cn), ...
                2 * y_ij * unique(tau_cn)];
            yn = tau_an.^2 - tau_cn.^2;
            y = tau_a.^2 - tau_c.^2;
            H = [H; h]; %D in costa's description Nx5
            YN = [YN; yn]; %T in costa's description Nx1
            Y = [Y; y];
            CCVALN = [CCVALN; ccvaln];
            CCVAL = [CCVAL; ccval];
            CCERR = [CCERR; ccerr]; %Nx1
            COMBOS = [COMBOS; repmat([combos(i, 1), combos(i, 2)], N)];
            RHO0 = [RHO0; repmat([peak(i), altpeak(i)], N)];
            RXILOC = [RXILOC; repmat([ENU(1:2, combos(i, 1))', ENU(1:2, combos(i, 2))'], N)];
        end
        
        %         exit;
        %         if tt == 7
        %             exit;
        %         else
        %             continue;
        %         end
        if isempty(H)
            disp(['H is empty'])
            continue
        end
        [estbar, covest, percentage, estbaro] ...
            = estimate_SAGA(H, YN, Y, CCVALN, CCVAL, CCERR, 0, COMBOS, RHO0, RXILOC, rcvr_op);
      
       %[vmaghat, vanghat, vgehat, vgnhat, arhat, psiahat, vctovhat, argshat]= plotva(H, Y, CCVAL, CCERR, COMBOS, peak, denu)
        estbarvec = num2cell(estbar);
        eestvec = num2cell(diag(sqrt(covest)));
        [~, ~, vge_sc, vgn_sc, ar, psia, vc] = deal(estbarvec{:});
        [~, ~, ~, ~, ear, epsia, evc] = deal(eestvec{:});
        % Mapping matrix from ^scV^g to ^scV^prn
        % ^scV^prn = ^scV^g - ^prnV^g
        M = [1, 0, -1, 0; 0, 1, 0, -1];
        vg = M * [vge_sc; vgn_sc; vge_prnbar; vgn_prnbar];
        vge = vg(1);
        vgn = vg(2);
        covscen = covest(3:4, 3:4);
        covdrift = blkdiag(covscen, evge_prn^2, evgn_prn^2);
        coven = M * covdrift * M';
        evge = sqrt(coven(1, 1));
        evgn = sqrt(coven(2, 2));
        % Jacobian from [vge;vgn] to [vmag;vang]
        vmag = sqrt(vge.^2+vgn.^2);
        vang = atan2(vgn, vge);
        if vang < 0
            vang = 2 * pi + vang;
        end
        J = [vge ./ vmag, vgn ./ vmag; ...
            -vgn ./ vmag.^2, vge ./ vmag.^2];
        vge
        vgn
        vmag
        vang
        ar %axial ratio
        psia %Psi_axialratiopsi
        vc/vmag 
        vmag_mean = vmag;
        vang_mean = vang;
        
        % covariances of ^scV^prn
        covmagang = J * M * covdrift * (J * M)';
        evmag = sqrt(covmagang(1, 1));
        evang = sqrt(covmagang(2, 2));
        
        ESTV(tt, :) = [datenum(tslist(tt)/24/3600+init_time), ...
            datenum(telist(tt)/24/3600+init_time), ...
            vmag, evmag, rad2deg(vang), rad2deg(evang), vge, evge, vgn, evgn, ...
            ar, ear, rad2deg(psia), rad2deg(epsia), vc, evc, percentage];
        
        ESTO(:, tt) = estbaro;
        %save the ellipse
%         plotname = ['Ellipse_PRN', num2str(prn), '_', ...
%             num2str(tslist(tt), '%.0f'), '-', num2str(telist(tt), '%.0f'), 's_after', ...
%             datestr(init_time, '_HHMMUT')];
%         plotpath = [op_path, plotname, '.png'];
%         saveas(gcf, plotpath, 'png');
%         close;
        
        %         save('lz.mat');
        %         return;
        
        %skyplot
%         if tt == length(telist) && i_dtau == 0
%             tt
%             tstt_sky = datevec(init_time+t(1)/24/3600);
%             tend_sky = datevec(init_time+t(end)/24/3600);
%             [AZb, ELb, beamid] = plotPFISR_NeTe(tstt_sky, tend_sky, 'Ne');
%             AZp = [repmat(AZb, [1, length(telist)]); AZ' * 180 / pi];
%             ELp = [repmat(ELb, [1, length(telist)]); EL' * 180 / pi];
%             prnp = [beamid; prn * ones(size(rcvr_op, 1), 1)];
%             beamAZ = repmat(AZb, [1, length(telist)]);
%             beamEL = repmat(ELb, [1, length(telist)]);
%             azprn = AZ' * 180 / pi;
%             elprn = EL' * 180 / pi;
%             prnprn = prn * ones(size(rcvr_op, 1), 1);
%             save([home_dir, sep, mat_dir, sep, 'skyplotdata_PRN', ...
%                 num2str(prn), '_', year, '_', doy, '_zoom', num2str(zcounter), '.mat'], ...
%                 'beamid', 'beamAZ', 'beamEL', 'azprn', 'elprn', 'prnprn');
%                              skyPlot(AZ'*180/pi,EL'*180/pi,prn*ones(size(rcvr_op,1),1));
%             figure;
%             skyPlot(AZp, ELp, prnp);
%             
            %             skyplot_v2(AZp, ELp, scintdata, prnp);
            
%             plotpath = [op_path, 'PRN', num2str(prn), ...
%                 '_', year, '_', doy, ...
%                 '_zoom', num2str(zcounter)', ...
%                 '_', num2str(tlim(1)), '-', num2str(tlim(2)), ...
%                 's_after_', datestr(init_time, 'HHMM'), 'UT', ...
%                 '_skyplot', '.png'];
%             saveas(gcf, plotpath, 'png');
%             close;
%             %                 return;
              %end
        
        
        % Spectral analysis
      [~, op_path, ~] = ver_chk;
        
        % flags used in the analysis
        
        % normalize the mean squared errors?
        normalize = 0;
        if normalize
            normflag = 'normalized';
        else
            normflag = 'non-normalized';
        end
        
        % time shift received signals?
        shifted = 0;
        if shifted
            shiftedflag = 'timeshifted';
        else
            shiftedflag = 'nottimeshifted';
        end
        
        % do the fit for all receivers or individually?
        fitall = 0;
        if fitall
            fitflag = 'fitall';
        else
            fitflag = 'fiteach';
        end
        
        % show estimates from Dr. Bust's implementation?
        debug = 0;
        if debug
            debugflag = 'yesBust';
        else
            debugflag = 'noBust';
        end
        
        zeropadding = 1;
        if zeropadding
            zeropadflag = 'zeropadding';
        else
            zeropadflag = 'nopadding';
        end
        
        % welch method or basic periodogram?
        welch = 1;
        if welch
            spectrumflag = 'welch';
        else
            spectrumflag = 'periodogram';
        end
        
        % if welch, do overlap?
        overlap = 1;
        if overlap
            overlapflag = '50overlap';
        else
            overlapflag = '0overlap';
        end
        
        % filter the signal with windows?
        windowed = 1;
        
        % power or amplitude?
        if fluct==1
           factor = 0.5; 
        elseif fluct==0
            factor= 1;
        end
                
        % iteration grid
        zmin = 125e3;
        zmax = 1000e3;
        Lmin = 25e3;
        step = 25e3;
        
        % figure format, png/pdf/eps2c
        format = 'png';
        suffix = '.png';
        %         format = 'pdf'; suffix = '.pdf';
        %         format = 'epsc2'; suffix = '.eps';
        
        %         keyboard;
        
        %creat folders for figures
        op_path = strjoin({op_path, 'spectral', ...
            normflag, debugflag, zeropadflag, overlapflag, fitflag, filesep}, filesep);
        mkdir = strjoin({'mkdir -p', op_path});
        system(mkdir);
        
        if ~isnan(vmag) && ~isnan(vang) && (vc/vmag<1) && (vc/vmag>0)
            for n = 1:numsim %defined for estimate_obs
                for rr = 1:size(rcvr_op, 1)
                    if strcmp(rcvr_op(rr, :), 'ASTRArx') && strcmp(doy, '342') && strcmp(year, '2013')
                        AZ(tt, rr) = mean(AZ(tt, [1:rr - 1, rr + 1:end]));
                        ZE(tt, rr) = mean(ZE(tt, [1:rr - 1, rr + 1:end]));
                    end        

                    %Amplitude and phase of the receivered signal

                    %original signals
                    pwr_o = xdata{rr}(:, 2); %+ 0.25*randn(l,1);
                    ph_o = xdata{rr}(:, 3); %+ 0.25*randn(l,1);
                    %time-shifted, aligned and truncated signals
                    if shifted  
                        pwr = pwr_truncated{rr}; %+ 0.25*randn(l,1);
                        ph = ph_truncated{rr}; %+ 0.25*randn(l,1);
                    else             
                        pwr =pwr_o; %+ 0.25*randn(l,1);
                        ph = ph_o; %+ 0.25*randn(l,1);
                    end
                    
                    scint = pwr.*exp(1i.*ph);
                    noise = A_nu(:,rr,n).* exp(1i.* phi_nu(:,rr,n));
                    scint_n = scint + noise;
                    pwr_n = abs(scint_n);
                    ph_n = angle(scint_n);
                    
                    clear ind ind2
                    ind = find(((ph>pi)&(ph_n<0))|(ph-ph_n>2));
                    ph_n(ind) = ph_n(ind) + 2*pi;
                    ind2 = (((ph<(-pi))&(ph_n>0))|(ph-ph_n<-2));
                    ph_n(ind2) = ph_n(ind2) - 2*pi;
                    
                    l = length(pwr);
                    if zeropadding
                        NFFT = 2^nextpow2(l);
                    else
                        NFFT = l;
                    end

                    if windowed
                        window_welch = [];
                        window_period = [];
                    else
                        window_welch = ones(NFFT, 1);
                        window_period = ones(NFFT, 1);
                    end

                    if overlap
                        noverlap = [];
                    else
                        noverlap = 0;
                    end
                    
                    if n==1
                        pwr = pwr;
                        ph =ph; %no noise added
                        [Spwr_obs_welch_nn{tt, rr}, ~] = pwelch(factor * log(pwr), window_welch, noverlap, NFFT, Fs, 'psd');
                        [Sph_obs_welch_nn{tt, rr}, f] = pwelch(ph, window_welch, noverlap, NFFT, Fs, 'psd');
                        [Spwr_obs_period_nn{tt, rr}, ~] = periodogram(factor * log(pwr), window_period, NFFT, Fs, 'psd');
                        [Sph_obs_period_nn{tt, rr}, f] = periodogram(ph, window_period, NFFT, Fs, 'psd');
                        R_obs_welch_nn{tt}(:, rr) = Spwr_obs_welch_nn{tt, rr} ./ Sph_obs_welch_nn{tt, rr};
                        R_obs_period_nn{tt}(:, rr) = Spwr_obs_period_nn{tt, rr} ./ Sph_obs_period_nn{tt, rr};
                    end
                    
                    pwr = pwr_n;
                    ph = ph_n;   
                    [Spwr_obs_welch{tt, rr}, ~] = pwelch(factor * log(pwr), window_welch, noverlap, NFFT, Fs, 'psd');
                    [Sph_obs_welch{tt, rr}, f] = pwelch(ph, window_welch, noverlap, NFFT, Fs, 'psd');

                    %                 keyboard;

                    [Spwr_obs_period{tt, rr}, ~] = periodogram(factor * log(pwr), window_period, NFFT, Fs, 'psd');
                    [Sph_obs_period{tt, rr}, f] = periodogram(ph, window_period, NFFT, Fs, 'psd');
                    
                    R_obs_welch{tt}(:, rr) = Spwr_obs_welch{tt, rr} ./ Sph_obs_welch{tt, rr};
                    R_obs_period{tt}(:, rr) = Spwr_obs_period{tt, rr} ./ Sph_obs_period{tt, rr};
                end
                %Welch or Periodogram?
                if welch
                    R_obs = R_obs_welch{tt};
                    R_obs_n{n} = R_obs_welch{tt};
                    if n == 1
                        R_obs_nn = R_obs_welch_nn{tt}; %no noise added
                    end
                else
                    R_obs = R_obs_period{tt};
                    R_obs_n{n} = R_obs_period{tt};
                    if n == 1
                        R_obs_nn = R_obs_period_nn{tt};
                    end
                end

                [Lgrid, zgrid] = meshgrid(Lmin:step:zmax, zmin:step:zmax);

                for rr = 1:size(rcvr_op, 1)
                    vmag = evmag*randn(1) + vmag_mean;
                    vang = evang*randn(1) + vang_mean;
                    ep = NaN(size(Lgrid));
                    for i = 1:size(Lgrid, 1)
                        for ii = 1:size(Lgrid, 2)
                            L = Lgrid(i, ii);
                            z = zgrid(i, ii);
                            if z - L >= 100e3
                                %add a distribution to velocity with its
                                %mean and std value
                               
                                if fitall == 0
                                    [R_rytov, k_par, k_par_index, R_Bust] = Lz(vmag, vang, AZ(tt, rr), ZE(tt, rr), L, z, f);
                                    R_obs_c = R_obs(k_par_index, rr);
                                    if n ==1
                                        R_obs_c_nn = R_obs_nn(k_par_index, rr);
                                    end
                                else
                                    [R_rytov, k_par, k_par_index, R_Bust] = Lz(vmag, vang, AZ(tt, :), ZE(tt, :), L, z, f);
                                    R_obs_c = R_obs(k_par_index, rr);
                                    if n ==1
                                        R_obs_c_nn = R_obs_nn(k_par_index, rr);
                                    end
                                end
                                
                                k_par_index_005 = find(k_par(k_par_index) < 0.05);
                                if isempty(k_par_index_005)==1
                                      k_par_index2 = k_par_index;
                                      R_obs_c2 = R_obs_c;
                                else
                                    R_obs_c_max = max(R_obs_c(k_par_index_005));
                                    k_par_index_R_obs_max = find(R_obs_c== R_obs_c_max);
                                    k_par_R_obs_max = k_par(k_par_index_R_obs_max);
                                    k_fit_limit = 2* k_par_R_obs_max;
                                    k_index_fit_limit = find(k_par == k_fit_limit);
                                    for ki=1:size(k_par_index)
                                        if k_par_index(ki) <= k_index_fit_limit
                                            k_par_index2(ki) = k_par_index(ki); %index to fit R_rytov to R_observed in k<k(R_obs_max*2)
                                        end
                                    end
                                    if exist('k_par_index2','var')==0
                                        k_par_index2=k_par_index;
                                    end
                                    R_obs_c2 = R_obs_c(k_par_index2, :);
                                end
                               
                                R_rytov_c = R_rytov(k_par_index2, :);
                                R_Bust_c = R_Bust(k_par_index2, :);
                                
                                clear k_par_index2
                                % normalized?
                                if normalize == 0
                                    sumsquared = (R_obs_c2 - R_rytov_c).^2;
                                    sumsquared = (log10(R_obs_c2) - log10(R_rytov_c)).^2;
                                else 
                                    sumsquared = ((log10(R_obs_c2) - log10(R_rytov_c)) ./ log10(R_obs_c2)).^2;
                                end  

                                epsqr = mean(sumsquared(:), 'omitnan');
                                ep(i, ii) = epsqr;
                                                                %                             if epsqr > 10
                                %                                 ep(i,j) = NaN;
                                %                                 continue;
                                %                             end
                            end
                        end
                    end
                    
                    ep_min_n(n, rr) = min(min(ep));
                    ep_max_n(n, rr) = max(max(ep));
                    
                    if isnan(ep_min_n(n,rr))
                        disp(['MSError Nan'])
                        continue
                    end

                    L_hat_n(n, rr) = Lgrid(ep == ep_min_n(n, rr));
                    z_hat_n(n, rr) = zgrid(ep == ep_min_n(n, rr));
                    %                 mesh(Lgrid/10^3, zgrid/10^3, ep, ep)
                    flagplotmin = 0;
                    if flagplotmin == 1 && n == numsim
                    figj = figure;
                    pcolor(Lgrid/10^3, zgrid/10^3, ep);
                    %                 shading(gca, 'flat');
                    %                 set(gca, 'layer', 'top');
                    caxis([0, 2]);
                    hold on;
                    %                 plot3(L_hat/10^3,z_hat/10^3,ep_min,'ro');
                    plot(L_hat_n(n, rr)/10^3, z_hat_n(n, rr)/10^3, 'ro');
                    xlabel('Thickness $L$ [km] ');
                    ylabel('Top height $z$ [km]');
                    zlabel('$\epsilon^2$');
                    title(['$\hat{L} = $', num2str(L_hat_n(n, rr)/10^3), ', ', ...
                        '$\hat{z} = $', num2str(z_hat_n(n, rr)/10^3), ', ', ...
                        '$\epsilon^2_{min} = $', num2str(ep_min_n(n, rr)), ', ', ...
                        '$v =$', num2str(vmag), ...
                        'for ', sitenum_op{rr, :}]);
                    %                 set(gca, 'Zscale', 'log');
                    %         zlim([0 10]);
                    %                 view([45,15]);
                    tightfig;
                    cb = colorbar;
                    set(get(cb, 'YLabel'), 'String', '$\epsilon^2$', 'interpreter', 'latex');
                    plotname = [year, '_', doy, '_PRN', num2str(prn), '_', sitenum_op{rr, :}, '_CostFunction_', ...
                        num2str(tslist(tt), '%.0f'), '-', num2str(telist(tt), '%.0f'), 's_after_', ...
                        datestr(init_time, 'HHMM'), 'UT_', num2str(factor), '_', num2str(zmax/1000)];
                    plotpath = [op_path, plotname, suffix];
                    saveas(gcf, plotpath, format);
                    close;
                    end
                end
            end
            for rr = 1:size(rcvr_op, 1)
                if isnan(ep_min_n(n,rr))
                continue
                end
                    L_hat_mean(:, rr) = mean(L_hat_n(:, rr));
                    L_hat_std(:, rr) = std(L_hat_n(:, rr));
                    z_hat_mean(:, rr) = mean(z_hat_n(:, rr));
                    z_hat_std(:,rr) =  std(z_hat_n(:, rr)); 
                    ep_min(:, rr) = mean(ep_min_n(:, rr));
            end
            xl = [-inf, inf]; %xl = [1e-3 1e-1];
            %             return;
            for rr = 1:size(rcvr_op, 1)
                if isnan(ep_min_n(n,rr))
                continue
                end
                L_hat(:, rr) = L_hat_mean(:, rr);
                z_hat(:, rr) = z_hat_mean(:, rr);
                if flagplotmin == 1
                figall = figure;
                if welch
                    loglog(k_par, Spwr_obs_welch{tt, rr}, '*c', ...
                        k_par, Sph_obs_welch{tt, rr}, '*g', ...
                        k_par, R_obs_welch{tt}(:, rr), '*m');
                    hold on
                    loglog(k_par, Spwr_obs_welch_nn{tt, rr}, 'b', ...
                        k_par, Sph_obs_welch_nn{tt, rr}, 'k', ...
                        k_par, R_obs_welch_nn{tt}(:, rr), 'r');
%                     loglog(k_par, R_obs_welch{tt}(:, rr), '*g')                    
%                     hold on
%                     loglog(k_par, Spwr_obs_welch_nn{tt, rr}, 'b', ...
%                         k_par, Sph_obs_welch_nn{tt, rr}, 'k', ...
%                         k_par, R_obs_welch_nn{tt}(:, rr), 'r');
                    
                else
                    % hold on;
                    loglog(k_par, Spwr_obs_period{tt, rr}, '*b', ...
                        k_par, Sph_obs_period{tt, rr}, '*k', ...
                        k_par, R_obs_period{tt}(:, rr), '*r');
                    hold on
                    loglog(k_par, Spwr_obs_period_nn{tt, rr}, 'c', ...
                        k_par, Sph_obs_period_nn{tt, rr}, 'g', ...
                        k_par, R_obs_period_nn{tt}(:, rr), 'm');
                end
                xlim(xl);
                ylim([10^-6, 10^2]);
                %                     set(gca,'YTick',[1e-5 1e-4 1e-3 1e-2 1e-1 1e0 1e1 1e2]);
                title('Observed Log-Amplitude to Phase Power Spectrum Ratio');
                legend({
                    'Noise Ratio', 'Noise Log-Amplitude', 'Noise Phase',...
                    'Noise Free Log-Amplitude', 'Noise Free Phase', 'Noise Free Ratio'}, 'location', 'southwest');
                    
                    
                xlabel(['Wavenumber along Drift Velocity Direction $\kappa_v$ [rad/m], ', sitenum_op{rr, :}]);
                %                     legend({'Phase','Log_{10} Power'},'location','best')
                tightfig;
                plotname = [year, '_', doy, '_PRN', num2str(prn), '_', sitenum_op{rr, :}, '_ObservedRatio_', ...
                    num2str(tslist(tt), '%.0f'), '-', num2str(telist(tt), '%.0f'), 's_after_', ...
                    datestr(init_time, 'HHMM'), 'UT'];
                plotpath = [op_path, plotname, suffix];
                saveas(gcf, plotpath, format);
                close;
                end
                
                % PLOT OF OBSERVED RATIO WITH NOISE
%                 if welch %No noise added analysis
%                     loglog(k_par, R_obs_welch{tt}(:, rr), 'g*', 'MarkerSize',0.5);
%                     hold on;
%                     loglog(k_par, R_obs_welch_nn{tt}(:, rr), 'r');
%                  else
%                     loglog(k_par, R_obs_period{tt}(:, rr), 'g*','MarkerSize',0.5);
%                     hold on
%                     loglog(k_par, R_obs_period_nn{tt}(:, rr), 'r'); %nn no noise added
%                 end
%                 xlim(xl);
%                 ylim([10^-6, 10^2]);
%                 %                        set(gca,'YTick',[1e-5 1e-4 1e-3 1e-2 1e-1 1e0 1e1 1e2]);
%                 title('Observed Log-Amplitude to Phase Power Spectrum Ratio');
%                 legend({'Observed Ratio with noise data','Observed Ratio'}, 'location', 'southwest');
%                 xlabel(['Wavenumber along Drift Velocity Direction $\kappa_v$ [rad/m], ', sitenum_op{rr, :}]);
%                 %                     legend({'Phase','Log_{10} Power'},'location','best')
%                 tightfig;
%                 plotname = [year, '_', doy, '_PRN', num2str(prn), '_', sitenum_op{rr, :}, '_ObservedRatio_', ...
%                     num2str(tslist(tt), '%.0f'), '-', num2str(telist(tt), '%.0f'), 's_after_', ...
%                     datestr(init_time, 'HHMM'), 'UT', '_nonoise'];
%                 plotpath = [op_path, plotname, suffix];
%                 saveas(gcf, plotpath, format);
%                 close;
            end
            fighat = figure;
            if size(rcvr_op, 1) > 2
                set(gcf, 'papersize', [8, 2 * size(rcvr_op, 1)], ...
                    'paperposition', [0, 0, 8, 2 * size(rcvr_op, 1)], ...
                    'paperpositionmode', 'auto', ...
                    'position', [0, 0, 8, 2 * size(rcvr_op, 1)]);
            end
            [sp, ~] = tight_subplot(size(rcvr_op, 1), 1, [0, 0.03], [0.11, 0.05], [0.11, 0.05]);
            for rr = 1:size(rcvr_op, 1)
                 if isnan(ep_min_n(n,rr))
                    disp(['No L found'])
                    continue
                 end
                for n=1:numsim 
                [R_rytov_hat, k_par, k_par_index, R_Bust] = Lz(vmag, vang, ...
                    AZ(tt, rr), ZE(tt, rr), L_hat_n(n, rr), z_hat_n(n, rr), f);
                [R_rytov_fixed, ~, k_par_index_Bust, R_Bust_fixed] = Lz(vmag, vang, ...
                    AZ(tt, rr), ZE(tt, rr), 500e3, 200e3, f);
                k_par_c = k_par(k_par_index, :);
                R_rytov_hat_c{tt}(:, rr) = R_rytov_hat(k_par_index, :);
                %R_obs_c(:, rr) = R_obs(k_par_index, rr);
                R_obs_c(:, rr) = R_obs_n{n}(k_par_index, rr);
                if debug
                    loglog(sp(rr), k_par_c, R_obs_c(:, rr), 'r', ...
                        k_par_c, R_rytov_hat_c{tt}(:, rr), 'c', ...
                        k_par_c, R_Bust(k_par_index), 'k', ...
                        k_par_c, R_rytov_fixed(k_par_index), 'g', ...
                        k_par_c, R_Bust_fixed(k_par_index), 'b');
                    if n==1
                    legend(sp(rr), ['Observed, ', sitenum_op{rr, :}], ...
                        ['R, $\hat{L}=', num2str(L_hat(:, rr)/10^3), '$, $\hat{z}=', num2str(z_hat(:, rr)/10^3), '$'], ...
                        'B, -', ...
                        'R, $\hat{L}=500$, $\hat{z}=200$', ...
                        'B, -', ...
                        'location', 'northwest', 'orientation', 'horizontal');
                    end
                else
                    loglog(sp(rr), k_par_c, R_obs_c(:, rr), 'r', ...
                        k_par_c, R_rytov_hat_c{tt}(:, rr), 'c');
                    hold(sp(rr),'on')
                    if n==1
                    legend(sp(rr), ['Observed, ', sitenum_op{rr, :}], ...
                        'Rytov', 'location', 'southeast'); %, $\hat{L}$ mean=', num2str(L_hat(:, rr)/10^3), ...
                        %'$\hat{z}$ mean=', num2str(z_hat(:, rr)/10^3), '$'], ...
                        
                    end 
                end
                end
                xlim(sp(rr), xl);
                ylim(sp(rr), [10^-5, 10^2.5]);
                                   %  set(gca,'YTick',[1e-5 1e-4 1e-3 1e-2 1e-1 1e0 1e1 1e2]);
                                    % ylim([1e-5*0.99 1e2*1.01]);
                
                
%                 MEGA_LZ(tt, rr, :) = [datenum(tslist(tt)/24/3600+init_time), ...
%                     datenum(telist(tt)/24/3600+init_time), ...
%                     L_hat(:, rr) / 10^3, z_hat(:, rr) / 10^3, ...
%                     ep_min(:, rr)];
                MEGA_LZ(tt, rr, :) = [datenum(tslist(tt)/24/3600+init_time), ...
                    datenum(telist(tt)/24/3600+init_time), ...
                    L_hat(:, rr) / 10^3, z_hat(:, rr) / 10^3, ...
                    L_hat_n(:, rr)' /10^3, z_hat_n(:, rr)' / 10^3, ...
                    ep_min(:, rr)];
            end
%             MEGA_LZ
            title(sp(1), 'Rytov and Observed Log-Amplitude to Phase Power Spectrum Ratio');
            xlabel(sp(rr), 'Wavenumber along Drift Velocity Direction $\kappa_v$ [rad/m]');
            set(sp(1:(rr - 1)), 'xticklabel', []);
            set(sp(2:2:rr), 'yticklabel', []);
            tightfig;
            plotname = [year, '_', doy, '_PRN', num2str(prn), '_', sitenum_op{rr, :}, '_RytovObserved_', ...
                num2str(tslist(tt), '%.0f'), '-', num2str(telist(tt), '%.0f'), 's_after_', ...
                datestr(init_time, 'HHMM'), 'UT_', num2str(factor), '_', num2str(zmax/1000)];
            plotpath = [op_path, plotname, suffix];
            saveas(gcf, plotpath, format);
            %             print(gcf, [plotpath, '_pver.png'], '-dpng', '-r600');
            %                 close;
            
            close;
        else
            for rr = 1:size(rcvr_op, 1)
                MEGA_LZ(tt, rr, :) = [datenum(tslist(tt)/24/3600+init_time), ...
                    datenum(telist(tt)/24/3600+init_time), ...
                    NaN, NaN, NaN(1,numsim), NaN(1,numsim), NaN,];
            end
        end
        %         keyboard;
    end
    if isempty(tslist) && isempty(telist)
        for rr = 1:size(rcvr_op, 1)
            MEGA_LZ(tt, rr, :) = NaN(1, 5);
        end
    end
end

if flagplotmin == 1
L_hat_vector = L_hat_n(:);
z_hat_vector = z_hat_n(:);
fighis3=figure;
sp1 = subplot(2,1,1);
histogram(L_hat_vector/1000);
sp2 = subplot(2,1,2);
histogram(z_hat_vector/1000);
title(sp1,'Histogram');
ylabel(sp1, 'L_hat') 
ylabel(sp2, 'z_hat')
saveas(gcf, [hr_path,'Lz_histogram'], 'png');
close
end

eststruct = struct('prn', prn, 't0', ESTV(:, 1), 'tf', ESTV(:, 2), ...
    'v', ESTV(:, 3), 'ev', ESTV(:, 4), ...
    'theta', ESTV(:, 5), 'etheta', ESTV(:, 6), ...
    've', ESTV(:, 7), 'eve', ESTV(:, 8), 'vn', ESTV(:, 9), 'evn', ESTV(:, 10), ...
    'ar', ESTV(:, 11), 'ear', ESTV(:, 12), ...
    'psia', ESTV(:, 13), 'epsia', ESTV(:, 14), ...
    'vc', ESTV(:, 15), 'evc', ESTV(:, 16), ...
    'percentage', ESTV(:, 17));

save(xcorr_results, 'ESTV', 'eststruct', 'ESTO', ...
    'zcounter', 'tlim', 'init_time', ...
    'v_ccmin', 'v_dtau', 'sitenum_op', 'rcvr_op');

save(xcorr_results, 'MEGA_LZ', '-append');

xcorr_te = toc;

disp(['Plotting lasted ', num2str(hrplot_te), 's'])
disp(['Drift estimation lasted ', num2str(xcorr_te), 's']);
close all;


% return;
%actual interval taken for cross-correlation
ta = datevec(init_time+t([1, end])/24/3600);
%set time axis limit for SAGA
% ts_cc = [tlim(1)-max(v_dtau)/24/3600 tlim(2)+max(v_dtau)/24/3600];
ts_cc = [tlim(1), tlim(2)];
%     %override time invtervals
 ts_cc = ([datenum([2014 2 20 11 0 0]) datenum([2014 2 20 12 0 0])]-init_time)*24*3600;
% save(['/data1/home/ysu27/Dropbox/research/MEGAVEST_',num2str(prn),'.mat']);
save([home_dir,sep,'high_rate.mat'],'RCVRNAME','init_time','rcvr_op','xdata_PRN')
% plotmisc(xcorr_results,year,doy,prn);
if numsim==1
    plotlz(prn,tstt,'debug')
else
    plotlz2(prn, tstt, 'debug',numsim)
end

% plotSAGAvsPFISR(prn,tstt,'ve_vn',fluct);
plotSAGAvsPFISR(prn, tstt, 'debug',fluct);
return;

% plotprnvs(prn,year,doy);
disp(['Finished processing for PRN', num2str(prn)]);

save([home_dir,sep,'high_rate.mat'],'RCVRNAME','init_time','rcvr_op','xdata_PRN')
end
