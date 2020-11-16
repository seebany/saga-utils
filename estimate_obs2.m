%function [tauaarrn_, taucarrn_, ccvalarrn_, ccerrarrn_, ...
%    tauaarr_, taucarr_, ccvalarr_, ccerrarr_,A_nu,phi_nu] = ...
%    estimate_obs2(rcvr_op, xdata, combos, flag,scale,numsim)
% takes in high-rate data from all operational receivers, simulates
% noise on the signals, and estimates the cross-correlation for
% all pairs.
%
% Yang Su 2018
% Modified by Aurora Lopez 2019 to have the noise simulation at this stage.
% Commented by S. Datta-Barua 9 June 2020
% A_nu is the simulated noise amplitude matrix [npts x numrx x numsim].
% phi_nu is the simulated noise phase matrix [npts x numrx x numsim].

function [tauaarrn_, taucarrn_, ccvalarrn_, ccerrarrn_, ...
    tauaarr_, taucarr_, ccvalarr_, ccerrarr_,A_nu,phi_nu] = ...
    estimate_obs2(rcvr_op, xdata, combos, flag,scale,numsim)
% Construct observations for the linear system
dbstop if error;

% rng('default');
tic;
if nargin == 0
    disp('debug begins');
    close all;
    clear;
    load('testdata1.mat');
    %     load('../HRdataforDrBust/Bust_PRN23_2013_342');
    %     xdata = xdata_PRN;
else
    save('testdata.mat');
    %     exit;
end

% Specify which phase data are used
if strcmp(flag, 'phase')
    colnum = 3;
    dt = 0.01;
    fluct = 0;
elseif strcmp(flag, 'power')
    colnum = 2;  
    dt = 0.01;
    fluct = 1;
elseif strcmp(flag, 'JROph')
    colnum = 3;
    dt = 0.02;
    fluct = 0;
elseif strcmp(flag, 'JROpow')
    colnum = 2;
    dt = 0.02;
    fluct = 1;
end
% Generate normally distributed gaussian noise
%numsim = 100;

% SDB 7/6/20 Power seems to be detrended here with Matlab built-in linear trend removal.  Not sure why.
if colnum == 2
    for i = 1:size(combos,1)
        xdata{combos(i,1)}(:,colnum) = detrend(xdata{combos(i,1)}(:,colnum));
        xdata{combos(i,2)}(:,colnum) = detrend(xdata{combos(i,2)}(:,colnum));
    end
end
data = horzcat(xdata{:});
tmat = data(:, 1:size(xdata{1}, 2):end);
obsmat = data(:, colnum:size(xdata{1}, 2):end);

[rhoarr, lag] = xcorr(obsmat);
rhomax = max(rhoarr);

% Permutations of receiver pairs including auto-correlation
perms = sortrows([combos; fliplr(combos); repmat([1:size(rcvr_op, 1)]', 1, 2)]);
% [~, accols] = ismember(repmat([1:size(rcvr_op, 1)]', 1, 2), perms, 'rows');

% This scale factor is based on A. Lopez's experience of typical amplitude variations between receivers.  Hardcoded because if fluct is set to analyze phase scintillations, then scale will correspond to the sigma of the phase fluctuations rather than of the amplitude. SDB 6/9/20
% Also note that the Gaussian noise was originally coded to be on the power,
% but it's the amplitude that is Gaussian, so adjust noise to be on amplitude.
scale =0.1;%sqrt(0.3);

A_nu = scale* randn(size(xdata{1},1),size(rcvr_op,1),numsim); 
phi_nu = (6.7*pi/180)*randn(size(xdata{1},1),size(rcvr_op,1),numsim); 
% a =-pi/2;
% b = pi/2;
% phi_nu = a + (b-a).*rand(size(xdata{1},1),size(rcvr_op,1),numsim);
noise = A_nu.* exp(1i.*phi_nu);

% Some diagnostic plots.-----------------
fig =figure;
cont=1;
for i=1:size(noise,3) %Real and imag part of phasor noise of every simulation
    for ii=1:100:size(noise,1)
    y(cont,:) = imag(noise(ii,:,i));
    x(cont,:) = real(noise(ii,:,i));
    cont = 1 + cont;
    end
end

for rcvr_noiseplot=1:size(x,2)
    [color] = rx_color(rcvr_op(rcvr_noiseplot, :));
    hold on
    quiver(0*x(:,rcvr_noiseplot),0*y(:,rcvr_noiseplot),...
        x(:,rcvr_noiseplot),y(:,rcvr_noiseplot),'Color', color);
end
title('Complex noise added to signal (sample of 1/100)')
xlabel('Real noise')
ylabel('Imaginary noise')
legend(rcvr_op);
plotname = ['Complex_noise'];
[~, path] = ver_chk;
saveas(gcf, [path,plotname], 'png');
close

% SDB 7/7/20 Loops not needed for plotting.
%cont =1;
%for i=1:size(noise,3) %Real and imag part of phasor noise of every simulation
%    for ii=1:1:size(noise,1)
%    A_nu_plot(cont,:) = A_nu(ii,:,i);
%    phi_nu_plot(cont,:) = phi_nu(ii,:,i);
%    cont = 1 + cont;
%    end
%end
%clear cont

fighis=figure;
sp1 = subplot(2,1,1)
histogram(A_nu);%histogram(A_nu_plot);
sp2 = subplot(2,1,2)
histogram(phi_nu);%histogram(phi_nu_plot);
title(sp1,'Histogram');
ylabel(sp1, 'A_\nu') 
ylabel(sp2, 'phi_\nu')
saveas(gcf, [path,'Noise_histogram'], 'png');
close

%tmat2 = ones(size(tmat,1)*10,size(tmat,2));
%for i=1:numsim
%    tmat2((i-1)*size(tmat,1)+1:(i)*size(tmat,1),:) = tmat(:,:);
%end
%
%fighis2=figure;
%sp1 = subplot(2,1,1)
%
%hist3([A_nu_plot(:,1),tmat2(:,1)]);
%sp2 = subplot(2,1,2)
%hist3([phi_nu_plot(:,1),tmat2(:,1)]);
%title(sp1,'2D Histogram');
%xlabel(sp1, 'A_nu') 
%ylabel(sp1, 't')
%xlabel(sp2, 'phi_nu')
%ylabel(sp2, 't')
%saveas(gcf, [path,'Noise_histogram_2D'], 'png');
%close
%--------------------------------

% Initialize empty variables.
[taucarrn, tauaarrn, ccvalarrn, ccerrarrn] = deal(cell(length(combos), numsim));
[taucarr, tauaarr, ccvalarr, ccerrarr] = deal(cell(length(combos), 1));
%figobs = figure;

% Loop over k simulation number.
for k = 1:numsim
    
    [tobsn_k, tobs_k] = deal([]);
    %     if verLessThan('matlab', '9.1.0')
    %         noisearr = randn(size(obsmat));
    %     else
    %         noisearr = wgn(size(obsmat, 1), size(obsmat, 2), 0);
    %     end
    
    % SDB 7/2/20 Since noise is on the amplitude, the complex signals must
    % be added, and then the power is the signal times its complex conjugate.
    psi_tilde =sqrt(data(:, 2:size(xdata{1}, 2):end)).*exp(1i.*data(:, 3:size(xdata{1}, 2):end)) + ...
        noise(:, :, k); %noise is added to complex signal
    if fluct == 1
	% Look at correlation of power in decibels.
        obsmat2 = 10*log10((abs(psi_tilde)).^2);
    elseif fluct == 0
        obsmat2 = angle(psi_tilde);
        for ind1=1:size(rcvr_op,1)
            ind = find(((obsmat(:,ind1)>pi)&(obsmat2(:,ind1)<0))...
                |(obsmat(:,ind1)-obsmat2(:,ind1)>2));
            obsmat2(ind,ind1) = obsmat2(ind,ind1) + 2*pi;
            ind2 = find(((obsmat(:,ind1)<(-pi))&(obsmat2(:,ind1)>0))...
                |(obsmat(:,ind1)-obsmat2(:,ind1)<-2));
            obsmat2(ind2,ind1) = obsmat2(ind2,ind1) - 2*pi;
            clear ind ind2    
        end
    end
    
    rhoarrn = xcorr(obsmat2); 
    rhomaxn = max(rhoarrn);
    %loop through i pairs of receivers
    for i = 1:length(combos)
        [~, ijcol] = ismember(combos(i, :), perms, 'rows');
        %         [~, jicol] = ismember(combos(i,:), perms, 'rows');
        [~, iicol] = ismember([combos(i, 1), combos(i, 1)], perms, 'rows');
        [~, jjcol] = ismember([combos(i, 2), combos(i, 2)], perms, 'rows');
        
        %Normalization
        normscale = sqrt(rhomax(iicol)*rhomax(jjcol));
        %         abc = rhoarr(:,iicol) + ...
        %             xcorr(noisearr(:,combos(i, 1),k), noisearr(:,combos(i, 1),k)) + ...
        %             xcorr(noisearr(:,combos(i, 1),k), obsmat(:,combos(i, 1))) + ...
        %             xcorr(obsmat(:,combos(i, 1)), noisearr(:,combos(i, 1),k));
        normscale1 = rhomax(iicol);
        normscale1n = rhomaxn(iicol);
        normscale2 = rhomax(ijcol);
        normscale2n = rhomaxn(ijcol);
        normscalen = sqrt(rhomaxn(iicol)*rhomaxn(jjcol));
        cc = rhoarr(:, ijcol) / normscale;
        %     cc2(:, i) = rhomat(:,jicol) / normscale;
        ac = rhoarr(:, iicol) / normscale;
        %     ac2(:, i) = rhomat(:,jjcol) / normscale;
        
%keyboard
        ccn = rhoarrn(:, ijcol) / normscale;
        %     cc2n(:, i) = rhonmat(:,jicol) / normscale;
        acn = rhoarrn(:, iicol) / normscale;
        %         acn(lag==0,i) = NaN;
        %         ac(lag==0,i) = NaN;
        %     ac2n(:, i) = rhonmat(:,jjcol) / normscale;
        
        %     B = 1000e-4;
        %     rhoiin = rhomaxn(iicol) * sinc(2*pi*B*lag'*dt);
        %     plot(lag*dt, rhon(:, iicol), lag*dt, rhoiin)
        
        %     %average auto-correlation of pairwise receiver signals
        %     ac(:, i) = (ac(:, i) + ac2(:, i)) / 2;
        
        %     hold on;
        
        
        if nargin == 0 && i == 1 && k == 1
            figrawcorr = figure;
            hold on;
            h1 = plot(lag*dt, cc, 'k', lag*dt, ac, 'c');
            [ccl, ccr] = findmainlobe(cc, lag);
            [ccnl, ccnr] = findmainlobe(ccn, lag);
            title(['Receiver ', num2str(combos(i, 1)), ' \& ', num2str(combos(i, 2))]);
            h2 = plot(lag*dt, ccn, 'g.', lag*dt, acn, 'r.', 'markersize', 10);
            xlabel('Lag [s]');
            [acl, acr] = findmainlobe(ac, lag);
            [acnl, acnr] = findmainlobe(acn, lag);
            legend('$\rho_{ij}$', '$\rho_{ii}$', ...
                '$\tilde{\rho_{ij}}$', '$\tilde{\rho_{ii}}$');
           xlim([min([acl, ccl, acnl, ccnl]), max([acr, ccr, acnr, ccnr])]*dt);
           ylim([-0.08, 1.02]);
            tightfig;
            %saveas(figrawcorr, '../ccexp.pdf');
            %             keyboard;
            close(figrawcorr);
        end
        
        if nargin == 0 && k == 1 && i == 1
            debugflag = 0;
        else
            debugflag = 1;
        end

	% Simulation 1 is the noise-free version?  Or only need to do the first time. SDB 7/7/20
        if k == 1
            [rowscc{i}, tobs, tau_c, tau_a, ccval, ccerr] = ...
                findrows(cc, ac, lag, dt, 'original', [], debugflag);
        end
        % For all simulations find the tau values.
%keyboard        
        [~, tobsn, tau_cn, tau_an, ccvaln, ccerrn] = ...
            findrows(ccn, acn, lag, dt, 'noisy', rowscc{i}, debugflag);
        
        tobs_k = [tobs_k, tobs];
        tobsn_k = [tobsn_k, tobsn];
        taucarrn{i, k} = tau_cn;
        tauaarrn{i, k} = tau_an;
        ccvalarrn{i, k} = ccvaln;
        ccerrarrn{i, k} = ccerrn;
        
	% SDB 7/7/20 Only need to do the original signal the first time.
        if k == 1
            taucarr{i} = tau_c;
            tauaarr{i} = tau_a;
            ccvalarr{i} = ccval;
            ccerrarr{i} = ccerr;
        end
    end
    
    if nargin == 1
        figure(figobs);
        hold on;
        if k == 1
            tobsarr = zeros(numsim, length(tobs_k));
            errarr = zeros(numsim, length(tobs_k));
            plot(tobs_k, 'k.-', 'linewidth', 2);
        end
        errarr(k, :) = (tobsn_k - tobs_k)';
        plot(tobsarr(k, :), 'g');
        plot(errarr(k, :), 'r.');
        drawnow;
    end
end
% close(figobs);
[taucarrn_, tauaarrn_, ccvalarrn_, ccerrarrn_, ...
    taucarr_, tauaarr_, ccvalarr_, ccerrarr_] = deal(cell(length(combos), 1));
for i = 1:length(combos)
    taucarrn_{i, :} = vertcat(taucarrn{i, :})';
    tauaarrn_{i, :} = vertcat(tauaarrn{i, :})';
    ccvalarrn_{i, :} = vertcat(ccvalarrn{i, :})';
    ccerrarrn_{i, :} = vertcat(ccerrarrn{i, :})';
    taucarr_{i, :} = horzcat(taucarr{i, :})';
    tauaarr_{i, :} = horzcat(tauaarr{i, :})';
    ccvalarr_{i, :} = horzcat(ccvalarr{i, :})';
    ccerrarr_{i, :} = horzcat(ccerrarr{i, :})';
end
toc;
end

% function findrows----------------------------------------
function [tmprowscc, tobs, tau_c, tau_a, ccval, errs] = ...
    findrows(cc, ac, lag, dt, flag, rowscc, debugflag)
% figpks = figure;
% findpeaks(cc, lag, 'annotate', 'extents', 'widthreference', 'halfheight');
% [pks, loc, ~, ~] = ...
%     findpeaks(cc, lag, 'annotate', 'extents', 'widthreference', 'halfheight');

%above a cut-off, i.e. 0.6
rhocutoff = 0.65;
nearzero = 1e-2;
% nearzero = 0.1;

lagmax = lag(cc == max(cc));
[laglcc, lagrcc] = findmainlobe(cc, lag);  
[laglac, lagrac] = findmainlobe(ac, lag);
if isempty(lagrac)
    tobs = [];
    tau_c = [];
    tau_a = [];
    ccval = [];
    errs = [];
    tmprowscc = [];
    return
end
% Find the rows of the auto-correlation that are on the main lobe to the right of the peak.
tmprowsac = find(ac' > 0 & lag > 0 & lag < lagrac);


if ~strcmp(flag, 'noisy')
    %take the lag positive side of the main lobe
    %tmprowscc = find(cc' >= rhocutoff & lag >= lagmax & lag < lagrcc);
    tmprowscc = find(lag >= lagmax & lag < lagrcc);
else
    tmprowscc = rowscc;
end

if debugflag == 0
    figdebug = figure;
    %     ac(lag==0) = NaN;
    %     plot(lag*dt, cc, 'k', lag*dt, ac, 'k');
    hold on;
    %     plot(laglac*dt,ac(lag==laglac),'o');
    %     plot(lagrac*dt,ac(lag==lagrac),'o');
    hc = plot(lag*dt, cc, 'r');
    %     plot(lag(tmprowscc_)*dt, cc(tmprowscc_), 'rs-');
    ha = plot(lag(tmprowsac)*dt, ac(tmprowsac), 'c');
    if strcmp(flag, 'noisy')
        rhostr = '\tilde{\rho}';
    else
        rhostr = '\rho';
    end
    title([flag, ', $\rho_{cutoff}$ = ', num2str(rhocutoff), ...
        ', $\min |', rhostr, '_{ij}-', rhostr, '_{ii}| \le$', num2str(nearzero)]);
    legend([hc, ha], {['$', rhostr, '_{ij}$'], ['$', rhostr, '_{ii}$']});
    xlabel('Lag [s]');
    if ~isempty(tmprowscc)
        xlim(dt*[min(lag(tmprowscc(1)), 0), max(lagrcc, lagrac)]);
    else
        xlim(dt*[min(lagmax, 0), max(lagrcc, lagrac)]);
    end
    %     keyboard;
    %         close;
end

if ~isempty(tmprowscc)
    
    ccmat = repmat(cc(tmprowscc)', size(ac(tmprowsac)));
    acmat = repmat(ac(tmprowsac), size(cc(tmprowscc)'));
    [errs, rowsmin] = min(abs(ccmat-acmat), [], 1);
    
    %remove any cross-auto pair with large errors
    invalid = errs >= nearzero;
    %     rowsmin(errs >= nearzero) = NaN;
    %     tmprowscc(errs >= nearzero) = NaN;
    %     errs(errs >= nearzero) = NaN;
    
    
    %make sure one tau_a corresponds to one tau_c
    if debugflag == 0
        for indcc = find(~isnan(tmprowscc) & ~invalid)
            plot([lag(tmprowscc(indcc)) * dt, lag(tmprowsac(rowsmin(indcc))) * dt], ...
                [cc(tmprowscc(indcc)), ac(tmprowsac(rowsmin(indcc)))], 'k.-', ...
                'linewidth', 0.5);
        end
        %         for rowac = unique(rowsmin, 'stable')
        %             rowunique = find(rowsmin == rowac);
        %             h1 = plot(lag(tmprowscc(rowunique))*dt, cc(tmprowscc(rowunique)), 'k.');
        %             rowminerr = rowunique(errs(rowunique) == min(errs(rowunique)));
        %             h2 = plot(lag(tmprowscc(rowminerr))*dt, cc(tmprowscc(rowminerr)), 'k*');
        %             h3 = plot([lag(tmprowscc(rowminerr))*dt lag(tmprowsac(rowac))*dt], ...
        %                 [cc(tmprowscc(rowminerr)) ac(tmprowsac(rowac))], 'k.-');
        %             errs(rowunique(rowunique ~= rowminerr)) = [];
        %             rowsmin(rowunique(rowunique ~= rowminerr)) = [];
        %             tmprowscc(rowunique(rowunique ~= rowminerr)) = [];
        %         end
        tightfig;
        saveas(figdebug, ['../getobs', flag, '.pdf']);
        %         keyboard;
        close(figdebug);
    end
    
    if any(invalid)
        %         disp('There are NaNs in $tau_a$');
    end
    tau_a = lag(tmprowsac(rowsmin)) * dt;
% SDB 7/7/20 Commenting this change to Nans out because it seems to make
% amplitude estimation with ensemble simulation give nan drift estimation. 
%    tau_a(invalid) = NaN;
    if isempty(tau_a)
        tau_a = NaN(1, length(tmprowscc));
        errs = tau_a;
    end
    tau_c = lag(tmprowscc) * dt;
    ccval = cc(tmprowscc)';
    tobs = tau_a.^2 - tau_c.^2;
    if any(abs(tau_a-tau_c) > 10)
        [laglcc, lagrcc] = findmainlobe(cc, lag);
        [laglac, lagrac] = findmainlobe(ac, lag);
    end
    if debugflag == 0
        figcc = figure;
        hold on;
        plot(lag*dt, cc, lag*dt, ac);
        plot(lag(tmprowscc)*dt, cc(tmprowscc), 'r', 'LineWidth', 2);
        plot(lag(tmprowsac)*dt, ac(tmprowsac), 'c', 'LineWidth', 2);
        plot(tau_c, cc(tmprowscc), 'k.');
        plot(tau_a, ac(tmprowsac(rowsmin)), 'k.');
        close(figcc);
    end
else
    tobs = [];
    tau_c = [];
    tau_a = [];
    ccval = [];
    errs = [];
end
end

% function findmainlobe------------------------------
function [lagmainl, lagmainr] = findmainlobe(cc, lag)
% This function identifies the main lobe of the cross-correlation.
% For comparison to the autocorrelation.
% S. Datta-Barua commenting 7/7/20

% Find the time lag corresponding to the peak cross-correlation value. 
if all(abs(flip(cc)-cc) < 1e-13)
    lagmax = 0;
else
    lagmax = lag(cc == max(cc));
end

% Find times earlier than the peak lag time that cross from negative to positive valued, or vice versa.
lagzerosl = lag((circshift(cc, 1) .* cc < 0 | circshift(cc, -1) .* cc < 0) & lag' < lagmax);
% Find times later than the peak lag time that cross from negative to positive valued, or vice versa.
lagzerosr = lag((circshift(cc, 1) .* cc < 0 | circshift(cc, -1) .* cc < 0) & lag' > lagmax);

% Choose the minimum zero crossing to the right of the peak lag.
[~, indminr] = min(lagzerosr-lagmax);
lagmainr = lagzerosr(indminr);
% Choose the maximum zero crossing to the left of the peak lag.
[~, indminl] = max(lagzerosl-lagmax);
lagmainl = lagzerosl(indminl);
end
