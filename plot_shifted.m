% Function plot_shifted.m plots CASES high-rate power and phase for as many
% receivers are given, time-shifting each receiver's time series so they
% are all aligned as much as possible.
%
% Created by Y. Su 2018
% Possibly modified by A. Lopez-Rubio 2019
% Commented by S. Datta-Barua 2 Apr 2020.  Adding output handles to figures
% so that calling function can save.
% 8 June 2020 SDB Just comment out plotting and keep only the ph_truncated and pwr_truncated output vars.

function [ph_truncated, pwr_truncated] = ...
    plot_shifted(xdata, tpeak, rcvr_op, sitenum_op, combos, Fs, fluct)
% load('lz.mat');
fig_each =[];
fig_all = [];

combos_auto = [1:size(rcvr_op, 1); 1:size(rcvr_op, 1)]';
[~, index] = sortrows([combos; fliplr(combos); combos_auto]);
tpeaks = [tpeak'; -tpeak'; zeros(size(rcvr_op, 1), 1)];
tpeaks_new = tpeaks(index);
tpeaks_array = reshape(tpeaks_new, [], size(rcvr_op, 1));

for rr = 1:size(rcvr_op, 1)
    tpeaks_rr = tpeaks_array(:, rr);
    leadnum(rr) = length(find(tpeaks_rr < 0));
    lagnum(rr) = size(rcvr_op, 1) - 1 - leadnum(rr);
end


rr_lead = find(leadnum == max(leadnum));
rr_lag = find(lagnum == max(lagnum));

leadsumsqr = sum(tpeaks_array(:, rr_lead).^2);
lagsumsqr = sum(tpeaks_array(:, rr_lag).^2);

rr_lead = rr_lead(leadsumsqr == min(leadsumsqr));
rr_lag = rr_lag(lagsumsqr == min(lagsumsqr));

% tmat(:,1) = tmat(:,1) - 2.79;
% tmat(:,2) = tmat(:,2) - 2.53;
% tmat(:,3) = tmat(:,3) - 1.62;
% tmat(:,4) = tmat(:,4);
% tmat(:,5) = tmat(:,5) - 2.13;

[ph_shifted, pwr_shifted, ph_truncated, pwr_truncated] = deal(cell(1, size(rcvr_op, 1)));
return

fig_all = figure;
for rr = 1:size(rcvr_op, 1)
    if (tpeaks_array(rr_lead, rr) * Fs + 1)<1
        return
    end
    ph_shifted{rr} = xdata{rr}(tpeaks_array(rr_lead, rr) * Fs + 1:end, 3);
    pwr_shifted{rr} = xdata{rr}(tpeaks_array(rr_lead, rr) * Fs + 1:end, 2);
end

[sp, ~] = tight_subplot(2, 1, [0, 0.03], [0.11, 0.05], [0.11, 0.05]);
for rr = 1:size(rcvr_op, 1)
if length(ph_shifted{rr_lag})>length(ph_shifted{rr})
    ph_truncated{rr} = [];
    pwr_truncated{rr} = [];
else
    ph_truncated{rr} = ph_shifted{rr}(1:length(ph_shifted{rr_lag}));
    pwr_truncated{rr} = pwr_shifted{rr}(1:length(pwr_shifted{rr_lag}));
    plot(sp(2), ph_truncated{rr}, 'color', rx_color(rcvr_op(rr, :)));
    %plot(sp(1), pwr_truncated{rr}, 'color', rx_color(rcvr_op(rr, :))); %NO DB
    plot(sp(1), 10*log10(pwr_truncated{rr}), 'color', rx_color(rcvr_op(rr, :))); %DB
    hold(sp(1), 'on');
    hold(sp(end), 'on');
end
end
set(sp(1:(end -1)), 'xticklabel', []);
legend(sp(1), sitenum_op, 'orientation', 'horizontal');
fig_all = tightfig;
%saveas(fig_all, '../../../all', 'png');
%saveas(fig_all, '/data1/home/alopez35/mfigures/all', 'png');
%close;

fig_each = figure;
if size(rcvr_op, 1) > 2
    set(gcf, 'papersize', [8, 2 * size(rcvr_op, 1)], ...
        'paperposition', [0, 0, 8, 2 * size(rcvr_op, 1)], ...
        'paperpositionmode', 'auto', ...
        'position', [0, 0, 8, 2 * size(rcvr_op, 1)]);
end
[sp, ~] = tight_subplot(size(rcvr_op, 1), 1, [0, 0.03], [0.11, 0.05], [0.11, 0.05]);
for rr = 1:size(rcvr_op, 1)
    %     plot(sp(rr), ph_shifted{rr});
    if fluct==0
    plot(sp(rr), ph_truncated{rr}, 'color', rx_color(rcvr_op(rr, :)));
    elseif fluct==1
    plot(sp(rr),10*log10(pwr_truncated{rr}), 'color', rx_color(rcvr_op(rr, :)));
    end
    legend(sp(rr), sitenum_op{rr});
end
set(sp(1:end-1), 'xticklabel', []);
set(sp(2:2:end), 'yticklabel', []);
fig_each = tightfig;
%saveas(fig_each, '../../../each', 'png');
%saveas(fig_each, '/data1/home/alopez35/mfigures/each', 'png');

%if ~exist('fig_each')
%	fig_each = [];
%end
%if ~exist('fig_all')
%	fig_all = [];
%end
%close;
