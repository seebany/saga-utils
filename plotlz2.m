function [] = plotlz2(prnlist, tstt, vflag,n)
[prop, op_path0, MEGAVEST_path] = ver_chk;

close all;
init();
[~, op_path] = ver_chk();
dbstop if error;
warning off;
weimer_path = [MEGAVEST_path, filesep, 'Weimer', filesep, 'runs', filesep];
medium_blue = [0, 0.447, 0.741];
light_blue = [0.301, 0.745, 0.933];
format = 'png';
% format = 'epsc2';

if nargin == 0
    dbstop if error;
    %     vflag = 'vmag_vang';
    %     vflag = 've_vn';
    %     vflag = 'ar_psia_vctov';
    vflag = 'debug';
    prnlist = [23, 13, 10];
    tstt = [2013, 12, 8, 4, 3, 0];
    prnlist = [23];
    tstt = [2013, 12, 8, 3, 43, 0];
    %     prnlist = [8, 9, 26];
    %     tstt = [2013, 12, 8, 7, 0, 0];
    %     prnlist = [25, 29, 31];
    %     tstt = [2013, 12, 8, 16, 0, 0];
    %         prnlist = [29];
    %         tstt = [2014, 2, 20, 11, 20, 0];
end
year = num2str(tstt(1));
doy = num2str(floor(datenum(tstt)-datenum([tstt(1), zeros(1, 5)])), '%03i');
matfilestruct = dir([MEGAVEST_path, 'xcorr_*.mat']);
tau = 60;
col = [3, 4];
figbox = figure;
[sp, ~] = tight_subplot(length(col), 1, [0.05, 0], [0.18, 0.05], [0.11, 0.05]);
matfilesname = dir([MEGAVEST_path, ...
    'xcorr_', year, '_', doy, '_PRN', num2str(prnlist), ...
    datestr(tstt, '_'), num2str(tau, '*_%is.mat')]);

if ~isempty(matfilesname)
    for fi = 1:length(matfilesname)
        
        load([MEGAVEST_path, matfilesname(fi).name]);
        for kk = 1:length(sitenum_op)
            switch kk
                case 1
                    lcolor = [1, 0, 1];
                    lcolor2 = 0.66 * lcolor;
                    lcolor3 = 0.33 * lcolor;
                    marker = 'd';
                case 2
                    lcolor = [1, 0, 0];
                    lcolor2 = 0.66 * lcolor;
                    lcolor3 = 0.33 * lcolor;
                    marker = '^';
                case 3
                    lcolor = [0, 0, 0];
                    lcolor2 = 0.5 + lcolor;
                    lcolor3 = 0.75 + lcolor;
                    marker = 'o';
                case 4
                    lcolor = [0, 1, 0];
                    lcolor2 = 0.66 * lcolor;
                    lcolor3 = 0.33 * lcolor;
                    marker = 's';
                case 5
                    lcolor = [0, 0, 1];
                    lcolor2 = 0.66 * lcolor;
                    lcolor3 = 0.33 * lcolor;
                    marker = '>';
            end
            
            tlimprn(fi, :) = tlim / 24 / 3600 + init_time;
            
            
            filter = MEGA_LZ(:, kk, 5) >= 0;
            ngood = find(filter);
            nbad = find(~filter & ~isnan(MEGA_LZ(:, kk, 3)));
            
            for subi = 1:length(col)
                vest = MEGA_LZ(:, kk, [1:2, col(subi)]);
                tc{kk, fi} = mean(vest(:, :, 1:2), 3);
                t0{kk, fi} = vest(:, :, 1);
                tf{kk, fi} = vest(:, :, 2);
                plotconfig = {marker, 'LineWidth', 1, ...
                    'Markersize', 4,};
                
                x = MEGA_LZ(:, :, 1);
                y = MEGA_LZ(:, :, 2);
                for tt =1:size(x,1)
                if col(subi) == 3
                    %diy error bar with mean and std values of L,z over all
                    %rr and over all n
                    plot(sp(subi), linspace(x(tt,1),y(tt,1),100), ...
                        linspace(mean(MEGA_LZ(tt,:,col(subi))) + std(mean(MEGA_LZ(tt,:,col(subi)+2:col(subi)+2+n))),...
                        mean(MEGA_LZ(tt,:,col(subi))) + std(mean(MEGA_LZ(tt,:,col(subi)+2:col(subi)+2+n))),100), 'b-')
                    hold(sp(subi),'on')
                    plot(sp(subi), linspace(x(tt,1),y(tt,1),100), ...
                        linspace(mean(MEGA_LZ(tt,:,col(subi))) - std(mean(MEGA_LZ(tt,:,col(subi)+2:col(subi)+2+n))),...
                        mean(MEGA_LZ(tt,:,col(subi))) - std(mean(MEGA_LZ(tt,:,col(subi)+2:col(subi)+2+n))),100), 'b-')
                    hold(sp(subi),'on')                    
                    plot(sp(subi), linspace((y(tt,1)+x(tt,1))/2,(y(tt,1)+x(tt,1))./2,2), ...
                        linspace(mean(MEGA_LZ(tt,:,col(subi))) + std(mean(MEGA_LZ(tt,:,col(subi)+2:col(subi)+2+n))),...
                        mean(MEGA_LZ(tt,:,col(subi))) - std(mean(MEGA_LZ(tt,:,col(subi)+2:col(subi)+2+n))),2), 'b-')
                    hold(sp(subi),'on')
                    plot(sp(subi), (y(tt,1)+x(tt,1))/2, mean(MEGA_LZ(tt,:,col(subi))), 'b*')
                end
                
                if col(subi) == 4
                    % DIY version of boxplot
                    %plotBox(sp(subi), x(:, 1), y(:, 1), topheight); %Plots for cases with no error estimation
                    
                    plot(sp(subi), linspace(x(tt,1),y(tt,1),100), ...
                        linspace(mean(MEGA_LZ(tt,:,col(subi))) + std(mean(MEGA_LZ(tt,:,col(subi)+2:col(subi)+2+n))),...
                        mean(MEGA_LZ(tt,:,col(subi))) + std(mean(MEGA_LZ(tt,:,col(subi)+2:col(subi)+2+n))),100), 'b-')
                    hold(sp(subi),'on')
                    plot(sp(subi), linspace(x(tt,1),y(tt,1),100), ...
                        linspace(mean(MEGA_LZ(tt,:,col(subi))) - std(mean(MEGA_LZ(tt,:,col(subi)+2:col(subi)+2+n))),...
                        mean(MEGA_LZ(tt,:,col(subi))) - std(mean(MEGA_LZ(tt,:,col(subi)+2:col(subi)+2+n))),100), 'b-')
                    hold(sp(subi),'on')                    
                    plot(sp(subi), linspace((y(tt,1)+x(tt,1))/2,(y(tt,1)+x(tt,1))./2,2), ...
                        linspace(mean(MEGA_LZ(tt,:,col(subi))) + std(mean(MEGA_LZ(tt,:,col(subi)+2:col(subi)+2+n))),...
                        mean(MEGA_LZ(tt,:,col(subi))) - std(mean(MEGA_LZ(tt,:,col(subi)+2:col(subi)+2+n))),2), 'b-')
                    hold(sp(subi),'on')
                    plot(sp(subi), (y(tt,1)+x(tt,1))/2, mean(MEGA_LZ(tt,:,col(subi))), 'b*')
                end
                
                switch col(subi)
                    case 3
                        titlestr = ['SAGA Thickness'];
                        ylabel(sp(subi), '$\hat{L}$ [km]');
                    case 4
                        titlestr = ['SAGA Top Height'];
                        ylabel(sp(subi), '$\hat{z}$ [km]');
                end
                end
                title(sp(subi), [num2str(subi, '(%i)'), titlestr]);
            end
        end
    end
    for subi = 1:length(col)
        set(sp(subi), 'ytick', 100:200:1050, 'yticklabel', 100:200:1050);
        ylim(sp(subi), [50, 1050]);
    end
else
    return;
end

tmin = min(tlimprn(:, 1));
tmax = max(tlimprn(:, 2));
if (tmax - tmin) * 24 * 3600 <= 300
    ticklbl = 'HH:MM:SS';
    rotang = 0;
else
    ticklbl = 'HH:MM';
    rotang = 15;
end
set(sp, 'xlim', [tmin, tmax]);
for subi = 1:length(col)
    datetick(sp(subi), 'x', ticklbl, 'keeplimits');
    if subi ~= length(col)
        set(sp(subi), 'XTickLabel', []);
    else
        %         set(sp(subi), 'xticklabelrotation', rotang);
        xlabel(['Time [HH:MM UT] on: ', datestr(init_time, 'mm/dd/yyyy')]);
    end
end

%If Pfsir_NeTe data avaialble
flag_NeTe=0;
if flag_NeTe==1
    figcomp = figure;
    [~, ~, ~, op_path] = plotPFISR_NeTe(tmin, tmax, 'Ne');
    h_ = findobj(boxes_overlap, 'type', 'patch');
    copyobj(boxes_overlap, gca);
    set(gca, 'xtick', get(sp(end), 'xtick'));
    datetick('x', 'HH:MM', 'keepticks');
    set(gca, 'xlim', [tmin, tmax]);
    saveas(figcomp, [op_path, 'PFISR_SAGA_LZ_PRN', num2str(prnlist), '_', year, '_', doy], format);
end
saveas(figbox, [op_path, 'SAGA_LZ_PRN', num2str(prnlist), '_', year, '_', doy], format);
close all;
end