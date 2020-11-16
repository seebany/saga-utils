function [AZb, ELb, beamid, op_path] = plotPFISR_NeTe(t1, t2, flag, beam, codetype, scintstart, scintend)
% This function plots PFISR electron density Ne or temperature Te,
% specified by string flag 'Ne' or 'Te' or 'Ti' respectively.  The inputs t1 and t2
% are datenum values for the desired start and end times on the x-axis.
%
% Y. Su 2017
% Commented by S. Datta-Barua 29 Nov 2018.
% S. Datta-Barua updating plots 1 Apr 2019. Additional input specifying
% beamID, which is 64016 for vertical and 64157 for up-B beam. Load and use
% both long pulse data (>200 km) and alternating code data (< 200 km).
% 27 Aug 2019 Add codetype as an input string, either 'lp' or 'ac'.
% 7 Sept 2019 Add flag option of Ti to plot ion temperatures.
% 17 Sept 2019 Add the 12/8/13 L,z study boxes from Aurora's data.
% 23 Sept 2019 Add input datenums scintstart, scintend as the times to be
% marked with vertical lines for the scintillation.

hold on

% Select plotting time interval based on input dates.
if floor(t1) == datenum([2015 3 17])
    % For this date, there are multiple experiments so this has to be
    % changed manually and run four times: 2 experiments * 2 pulse types.
%    t0 = [2015, 3, 17, 12, 45, 0];
%    tf = [2015, 3, 17, 13, 0, 0];
%    t0 = [2015, 3, 17, 13, 0, 0];
%    tf = [2015, 3, 17, 14, 15, 0];
    t0 = [2015, 3, 17, 14, 30, 0];
    tf = [2015, 3, 17, 16, 0, 0];
elseif floor(t1) == datenum([2015 10 7])
    % For this date, because there are multiple experiments during the time
    % interval of interest, this has to be changed manually and run four
    % times: 2 experiments * 2 pulse types (lp and ac). 4/3/19 SDB.
    t0 = [2015, 10, 7, 5, 30, 0];%[2015, 10, 7, 6, 0, 0];%[2015, 10, 7, 5, 30, 0];
    tf = [2015, 10, 7, 6, 0, 0];%[2015, 10, 7, 7, 0, 0];%[2015, 10, 7, 6, 0, 0];
%    t0 = [2015, 10, 7, 6, 0, 0];
%    tf = [2015, 10, 7, 7, 0, 0];
else
    t0 = t1;
    tf = t2;
end

%flag = 'Te';
[~, op_path, mat_path] = ver_chk;
%codetype ='lp'; %'ac'
alt_cutoff = 195; % km
mat_path = '~/matfiles/pfisr/';
switch codetype
    case 'lp'
% Load long pulse data.
try
% This was downloaded using Yang's dlPFISR.m which saves in the format of
% 'YEAR,MONTH,DAY,HOUR,MIN,SEC,AZM,ELM,BEAMID,RANGE,NEL,DNEL,TI,DTI,TE,DTE'.
    load(['PFISR_lp_' datestr(t1, 'yymmdd') '.mat'], 'PFISR_data');
catch
    try
    PFISR_data = importdata([datestr(t1, 'mmm-dd,yyyy'), '.txt']);%'Oct-7,2015.txt'); %'Mar-17,2015.txt');
    save(['PFISR_lp_', datestr(t1,'yymmdd'), '.mat'], 'PFISR_data');%151007.mat'], 'PFISR_data');
    catch
        load([mat_path, 'Madrigal_' datestr(t1, 'dd_mm_yyyy'), '.mat'], 'Madrigal');
        % Rename so that sorting columns works below. SDB 4/3/19.
        Madrigalac = Madrigal; clear Madrigal;
    end
end
    case 'ac'
% Load alternating code data, which was downloaded using Vaishnavi's
% script.
        load([mat_path, 'Madrigalac_' datestr(t1, 'dd_mm_yyyy'), '.mat'], 'Madrigalac');%7_10_2015.mat'],'Madrigalac');
end

if exist('Madrigalac','var')
% When downloaded from Madrigal using Vaishnavi's script, these are listed as
%[azm, minute, dayno, gdalt, gdlat, hour, tr, Ne, sec, Ti, Te, day, glon, elm, vO, range, year, ut, month, kp]
% so re-order with columns as
% [year, month, day, hour, minute, sec, azm, elm, beamid(hardcode as 64157), range, Nel (i.e., log10(Ne)), zeros(dNel), Ti, zeros(dTi), Te, zeros(dTe)]
PFISR_data = zeros(size(Madrigalac, 1),16);
cols = [17, 19, 12, 6, 2, 9, 1, 14]; cols2 = 16; cols3 = 8;
        PFISR_data(:,1:8) = Madrigalac(:,cols);%, ...
        
        upBrows = find(Madrigalac(:,1) == -154.3 & Madrigalac(:,14) == 77.5);
        PFISR_data(upBrows,9) = 64157;%repmat(64157,size(Madrigalac,1),1);%, ...
        PFISR_data(:,10) = Madrigalac(:,cols2);%, ...
        PFISR_data(:,11) = log10(Madrigalac(:,cols3));
        PFISR_data(:,13) = Madrigalac(:,10);
        PFISR_data(:,15) = Madrigalac(:,11);
            %];
        clear cols cols2 cols3
end    
% columns
% 1:DATENUM,2:AZM,3:ELM,4:BEAMID,5:RANGE,6:NEL,7:DNEL,8:TI,9:DTI,10:TE,11:DTE
data = [datenum(PFISR_data(:,1:6)), PFISR_data(:,7:end)];

% set up flags
if strcmp(flag, 'Ne')
    cols = 6;
elseif strcmp(flag, 'Te')
    cols = 10;
elseif strcmp(flag, 'Ti');
    cols = 8;
end

delta = 0;
if delta
    deltaflag = 'delta';
else
    deltaflag = 'normal';
end

%create folders for figures
% SDB commenting this out because I don't know which script defines these
% paths and just want to get the plots quickly. 11/29/18.
%op_path = strjoin({op_path, 'spectral', ...
%    flag, deltaflag, filesep}, filesep);
%mkdir = strjoin({'mkdir -p', op_path});
%system(mkdir);
%
% Yang previously widened the time window to be plotted, but to get around
% the multiple-ranges for different times problem, force this to be exact.
% SDB 11/29/18.
%data = data(data(:, 1) <= datenum(tf)+360/24/3600 & data(:, 1) >= datenum(t0)-360/24/3600, :);
data = data(data(:,1) <= datenum(tf) & data(:,1) >= datenum(t0),:);
% Also enforce an altitude limit for lp to be > 200 km and ac < 200 km.
% SDB 4/1/19.
switch codetype
    case 'lp'
data = data(data(:,5) > alt_cutoff,:);
    case 'ac'
        data = data(data(:,5) <= alt_cutoff, :);
end
if ~isempty(data)
    beamid = unique(data(:, 4), 'stable');
    % Select the vertical beam which has ID 64016. Southern beam has beam
    % ID 64157 with az == -154.3 deg and el == 77.5 deg.
    for ibeam = find(beamid == beam);%length(beamid)
        AZb(ibeam, :) = unique(data(data(:, 4) == beamid(ibeam), 2));
        ELb(ibeam, :) = unique(data(data(:, 4) == beamid(ibeam), 3));
        beamstr = strjoin({num2str(beamid(ibeam)), ...
            [num2str(AZb(ibeam, :)), '$^\circ$ az'], [num2str(ELb(ibeam, :)), '$^\circ$ el']}, ', ');
        %         subplot(2,1,2);
        data_beam = data(data(:, 4) == beamid(ibeam), :);
        % Change from using ranges to times to help reshape ac data. 4/2/19.
        ranges = unique(data_beam(:, 5), 'stable');
        % Sometimes ranges are nans.
        nonnanrangerows = find(isfinite(data_beam(:,5)));
        ranges = ranges(isfinite(ranges))
        times = unique(data_beam(:, 1), 'stable');
        switch flag
            case 'Ne'
        ne = 10.^data_beam(nonnanrangerows, cols);
        negrid = reshape(ne, [], length(times));
        %negrid = reshape(ne, length(ranges), []);
            case 'Te'
        te = data_beam(:, cols);
        tegrid = reshape(te, [], length(times));
        %tegrid = reshape(te, length(ranges), []);
            case 'Ti'
                ti = data_beam(:,cols);
                tigrid = reshape(ti, [], length(times));
        end
        
        % For some data sets, the ranges change partway through the time
        % interval so the data are not divisible by number of unique
        % ranges.

        %         continue;
%         difftimenegrid = [diff(negrid, 1, 2), NaN(length(ranges), 1)] ./ ...
%             [diff(times) * 24 * 3600; NaN]';
%         diffrangenegrid = [diff(negrid, 1, 1); NaN(1, length(times))] ./ ...
%             [diff(ranges) * 10e3; NaN];
        
        [timegrid, rangegrid] = meshgrid(times, ranges);
        if delta
            pcolor(timegrid, rangegrid, difftimenegrid);
            caxis([-1.5 * 10 ^ 9, 1.5 * 10 ^ 9]);
            title(['\parbox{6in}{\centering ', beamstr, '\\', ...
                'PFISR Electron Density Variation $\Delta {N_e} / \Delta t, \Delta t \approx 180 s$}'],...
                'interpreter', 'latex');
            cb = colorbar;
            set(get(cb, 'YLabel'), 'String', ...
                '$\frac{\Delta N_e}{\Delta t} [m^{-3} s^{-1}]$', ...
                'interpreter', 'latex');
        else
            switch flag
                case 'Ne'
                    % This if statement was needed before I eliminated
                    % rows whose ranges that were nans.
%                    if size(negrid) == size(timegrid)
            pcolor(timegrid, rangegrid, log10(negrid));
%                    else
%                        pcolor(timegrid, rangegrid, log10(negrid(1:end-1,:)));
%                    end
            caxis([10, 12]);
            title(['\parbox{6in}{\centering ', beamstr, '\\', ...
                'PFISR Electron Density $N_e$}'],'interpreter','latex');
            cb = colorbar;
            set(get(cb, 'YLabel'), 'String', ...
                '$\log_{10}{N_e} [m^{-3}]$', ...
                'interpreter', 'latex');
                case 'Te'
                    pcolor(timegrid, rangegrid, tegrid);
                    caxis([0, 3000]);
                    title(['\parbox{6in}{\centering ', beamstr, '\\', ...
                        'PFISR Electron Temperature $T_e$}'], ...
                        'interpreter', 'latex');
                    cb = colorbar;
                    set(get(cb, 'YLabel'), 'String', ...
                        '$T_e$ [K]', ...
                        'interpreter', 'latex');
                case 'Ti'
                    pcolor(timegrid, rangegrid, tigrid);
                    caxis([0, 3000]);
                    title(['\parbox{6in}{\centering ', beamstr, '\\', ...
                        'PFISR Ion Temperature $T_i$}'], ...
                        'interpreter', 'latex');
                    cb = colorbar;
                    set(get(cb, 'YLabel'), 'String', ...
                        '$T_i$ [K]', ...
                        'interpreter', 'latex');
            end
        end
        shading flat;
        set(gca, 'layer', 'top');
        grid off;
        xlabel('Time [HH:MM UT]');
        ylabel('Altitude [km]');
        ylim([50, 700]);
        xlim(datenum([t1; t2]));
        datetick('x', 'HH:MM', 'keeplimits');
        
        % Inserting vertical and horizontal line marks.  Comment out if not
        % needed. SDB 11/29/18.
        ax = axis;
        if ~exist('scintstart') || ~exist('scintend')
            
        if floor(data_beam(1,1)) == datenum([2015 10 7])
%         ax = axis;
                scintstart = datenum([2015 10 7 6 17 0]); %[2015 10 7 6 2 0];
                scintend = datenum([2015 10 7 6 27 0]); %[2015 10 7 6 20 0];
        elseif floor(t1) == datenum([2015 3 17])
%            plot(repmat(datenum([2015 3 17 13 8 0]), 2, 1), ax(3:4), 'k', 'LineWidth', 2);
%            plot(repmat(datenum([2015 3 17 13 25 0]), 2, 1), ax(3:4), 'k', 'LineWidth',2);
        elseif floor(datenum(t1)) == datenum([2014 11 16])
            scintstart = datenum([2014 11 16 9 11 0]);
            scintend = datenum([2014 11 16 9 25 0]);
%            plot(repmat(datenum([2014 11 16 9 11 0]), 2, 1), ax(3:4), 'k', 'LineWidth', 2);
%            plot(repmat(datenum([2014 11 16 9 25 0]), 2, 1), ax(3:4), 'k', 'LineWidth', 2);
%             plot(repmat(datenum([2014 11 16 1 0 0]), 2, 1), ax(3:4), 'k', 'LineWidth', 2);
%             plot(repmat(datenum([2014 11 16 1 15 0]), 2, 1), ax(3:4), 'k', 'LineWidth', 2);
%             plot(repmat(datenum([2014 11 16 0 50 0]), 2, 1), ax(3:4), 'k--', 'LineWidth', 2);
%             plot(repmat(datenum([2014 11 16 1 28 0]), 2, 1), ax(3:4), 'k--', 'LineWidth', 2);
        elseif floor(t1) == datenum([2014 11 22])
            scintstart = datenum([2014 11 22 22 29 0]);
            scintend = datenum([2014 11 22 22 49 0]);
%            plot(repmat(datenum([2014 11 22 22 29 0]), 2, 1), ax(3:4), 'k', 'LineWidth', 2);
%            plot(repmat(datenum([2014 11 22 22 49 0]), 2, 1), ax(3:4), 'k', 'LineWidth', 2);
        elseif floor(t1) == datenum([2013 12 8])
            scintstart = datenum([2013 12 8 3 43 0]);
            scintend = datenum([2013 12 8 4 17 0]);
%            plot(repmat(datenum([2013 12 8 3 43 0]), 2, 1), ax(3:4), 'k', 'LineWidth', 2);
%            plot(repmat(datenum([2013 12 8 4 17 0]), 2, 1), ax(3:4), 'k', 'LineWidth', 2);
        elseif floor(t1) == datenum([2015 3 18])
            scintstart = datenum([2015 3 18 8 4 0]);
            scintend = datenum([2015 3 18 15 0]);
%            plot(repmat(datenum([2015 3 18 8 4 0]), 2, 1), ax(3:4), 'k', 'LineWidth', 2);
%            plot(repmat(datenum([2015 3 18 8 15 0]), 2, 1), ax(3:4), 'k', 'LineWidth', 2);
%            plot(repmat(datenum([2015 3 18 8 6 0]), 2, 1), ax(3:4), 'k--', 'LineWidth', 2);
%            plot(repmat(datenum([2015 3 18 8 21 0]), 2, 1), ax(3:4), 'k--', 'LineWidth', 2);
%            plot(repmat(datenum([2015 3 18 9 25 0]), 2, 1), ax(3:4), 'k', 'LineWidth', 2);
%            plot(repmat(datenum([2015 3 18 9 41 0]), 2, 1), ax(3:4), 'k', 'LineWidth', 2);
        end
        end
        
%        plot(ax(1:2), [alt_cutoff alt_cutoff], 'k', 'LineWidth', 2);
%        plot(ax(1:2), [150 150], 'k', 'LineWidth', 2);
%        plot(repmat(scintstart, 2, 1), ax(3:4), 'k', 'LineWidth',2);
%        plot(repmat(scintend, 2, 1), ax(3:4), 'k', 'LineWidth', 2);
        
%        if floor(t1) == datenum([2013 12 8])
%            plot_Lz('lzdata2013doy342.mat');
	    [year,~,~,~,~,~] = datevec(t1);
            plot_Lz(['~/mfiles/saga/' datestr(now, 'yymmdd') '_lzdata' ...
		datestr(t1, 'yyyy') 'doy' num2str(floor(t1) - datenum([year 0 0])) ...
		'.mat']);
%        end
        
        plotname = strjoin({'pfisr', datestr(t1,'yymmdd_HHMM'), num2str(beamid(ibeam)), flag, deltaflag, 'cutoff195&150'}, '_');
        plotpath = [op_path, plotname, '.png']
        saveas(gcf, plotpath, 'png');
        savefig(gcf, [plotpath, '.fig']);%, 'fig');
    end
else
    beamid = NaN;
    AZb = NaN;
    ELb = NaN;
end
% close all;
end
