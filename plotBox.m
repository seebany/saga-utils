function [] = plotBox(ax, x, y, data)
w = y - x;
qmin = min(data, [], 2, 'omitnan');
qmax = max(data, [], 2, 'omitnan');
q1 = quantile(data, 0.25, 2);
q2 = quantile(data, 0.5, 2);
q3 = quantile(data, 0.75, 2);
dq = iqr(data, 2);
ql = q1 - 1.5 * dq;
qu = q3 + 1.5 * dq;

% errorbar(ax, (x + y)/2, q2, ...
%     min(q2-ql, q2-qmin), min(qmax-q2, qu-q2), ...
%     '.b', 'capsize', 12);
% errorbar(ax, (x + y)/2, q2, ...
%      min(q2-ql, q2-qmin), min(qmax-q2, qu-q2), ...
%      '.b');
errorbar_lon=0.2*10e-4; %For matlab versions without 'capsize'
for i = 1:size(x,1)
plot(ax, linspace((x(i) + y(i))/2-errorbar_lon,(x(i) + y(i))/2+errorbar_lon), ...
    linspace(q2(i) +  min(qmax(i)-q2(i), qu(i)-q2(i)),q2(i) +  min(qmax(i)-q2(i), qu(i)-q2(i))),'b'); 
plot(ax, linspace((x(i) + y(i))/2-errorbar_lon,(x(i) + y(i))/2+errorbar_lon), ...
    linspace(q2(i) - min(q2(i)-ql(i), q2(i)-qmin(i)),q2(i) - min(q2(i)-ql(i), q2(i)-qmin(i))),'b'); 
plot(ax,linspace((x(i) + y(i))/2, (x(i) + y(i))/2), ... 
    linspace(q2(i) - min(q2(i)-ql(i), q2(i)-qmin(i)), q2(i) +  min(qmax(i)-q2(i), qu(i)-q2(i))), 'b');
hold(ax, 'on');
for i = find(~isnan(data(:, 1)))'
%     rectangle(ax, 'position', [x(i), q1(i), w(i), dq(i)], ...
%         'facecolor', 'w', 'linewidth', 1.5);
%     rectangle(ax, 'position', [x(i), q2(i), w(i), 0], 'edgecolor', 'r');
    set(gcf,'CurrentAxes',ax)
    rectangle('position', [x(i), q1(i), w(i), dq(i)], ...
        'facecolor', 'w', 'linewidth', 1.5);
    rectangle('position', [x(i), q2(i), w(i), 0], 'edgecolor', 'r');
    %For matlab versions where cannot specify an axes handle 
end
end