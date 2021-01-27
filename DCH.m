% This script loops through the Detect, Classify, Hypothesize (DCH) method
% for estimating the scattering layer as described in Sreenivash et al., 2020.
%

clear all
clc
%frequency=  1 for L1; 2 for L2
%addpath(genpath('~/'));
outdir = '~/matfiles/saga_dch_2014-2019/';
signal_type = [0,2];
lronlyflag = 1;

tic
% Loop over frequency, year, doy.
for frequency=1:2
    sigstr = getSignal(signal_type(frequency));
    for year=2020:2020
	yearstr = num2str(year);
	% Make a directory for the files if the directory doesn't already exist
	inputdir = ['~/matfiles/saga_dch_2014-2019/', yearstr, '/'];
	if ~(exist(inputdir, 'dir'))
		system(['mkdir ' inputdir]);
	end
        % If some days' events have already been DCHed & stored as .mat, load.
        %if exist([outdir, 'SP_',yearstr,'_', sigstr, '_events.mat'],'file')
	%	load([outdir, 'SP_',yearstr,'_', sigstr, '_events.mat'],'phase_only');
        %	load([outdir, 'S4_',yearstr,'_', sigstr,'_events.mat'],'amplitude_only');
        %	load([outdir, 'SPS4_',yearstr,'_', sigstr,'_events.mat'],'phase_amplitude');
	%else
		phase_only=[];
        	amplitude_only=[];
        	phase_amplitude=[];
	%end % if exist SP_yearstr_sigstr_events.mat

	% If data were pre-loaded, check the last date that had any events
	if ~isempty(phase_only)
		lastdoy = phase_only(end,2);
	else
		lastdoy = 0;
	end

	% Start from the date the previous run left off at.
        endofyear = 365 + (mod(year,4) == 0);
	for doy=lastdoy+1:endofyear %2014 - 1 to 327; 2015 - 44 to 312
	    doystr = num2str(doy, '%03i');

	    % For 2014-15 data which were unpacked already for reproducing
	    % Vaishnavi's results, create links to save storage space.
	    %lncmd = ['ln ~/matfiles/saga_dch/' yearstr '/lrdata_' sigstr ...
	    %	'_' yearstr '_' doystr '.mat ~/matfiles/saga_dch_2014-2019/' ...
	    %	yearstr '/lrdata_' sigstr '_' yearstr '_' doystr '.mat'];
	    % For 2016-19 data which were unpacked already for
	    % Pau's results, create links to save storage space.
	    lncmd = ['ln ~/../pllado/matfiles/lrdata_' num2str(signal_type(frequency))  ...
		'_' yearstr '_' doystr '.mat ~/matfiles/saga_dch_2014-2019/' ...
		yearstr '/lrdata_' sigstr '_' yearstr '_' doystr '.mat'];
	    system(lncmd);


            disp(['Executing main to detect scintillation for ', ...
		yearstr,', ', doystr, ...
		', freq ', sigstr]);
            % If the day's data were previously stored, detection step is done. 
            if (exist([inputdir, sigstr,'_S4_',yearstr,'_', ...
		doystr,'.mat']) ==2 & ...
		exist([inputdir, sigstr,'_SP_',yearstr,'_', ...
		doystr,'.mat'])==2)
                
                display([sigstr,'_S4 and SP_',yearstr,'_', doystr,'.mat already existed']);
                %              display('Calling AC_LP_Vaishnavi_scint')
                
            else
                display([inputdir, sigstr,'_S4_',yearstr,'_', ...
			doystr,'.mat', ' not found and/or ']);
                display([inputdir, sigstr,'_SP_',yearstr,'_', ...
			doystr,'.mat', ' not found']);

		% Detect events using main.m
    		for fluct = 0:1
                	main(year,doy,lronlyflag,signal_type(frequency),fluct);
                	%main(year,doy,lronlyflag,0,signal_type);
    		end % for fluct
                
            end % if S4 and SP_*.mat files exist, detection step is done.
        %end % for doy

	%if redo_flag(2)
	%	lastdoyc = 0;
	%else
	%	lastdoyc = lastdoy;
	%end
        %for doy=lastdoy+1:365 %2014 - 1 to 327; 2015 - 44 to 312
	    % Next, classify layer.
            if (exist([inputdir, sigstr,'_S4_',yearstr,'_', ...
		doystr,'.mat']) ==2 & ...
		exist([inputdir, sigstr,'_SP_',yearstr,'_', ...
		doystr,'.mat'])==2)
                
                [s4sp,s4,sp,s4spregion,s4region,spregion] = ...
			classify_hypothesize(year,doy, ...
			signal_type(frequency), inputdir);
            spcat=horzcat(sp,spregion); %concatenating sp and sp region
            s4cat=horzcat(s4,s4region);
            s4spcat=horzcat(s4sp,s4spregion);
	    phase_only = [phase_only; spcat];
	    amplitude_only = [amplitude_only; s4cat];
	    phase_amplitude = [phase_amplitude; s4spcat];
            end % if detected events have been stored in L1CA_S*_yyyy_doy.mat.
        %end % for doy
	%end

	%if redo_flag(3)
	%	lastdoyh = 0;
	%	phase_only = phase_only(:,1:nspcols);
	%	amplitude_only = amplitude_only(:,1:nspcols);
	%	phase_amplitude = phase_amplitude(:,1:nspcols);
%
%	else
%		lastdoyh = lastdoy;
%	end
%	% Next, hypothesize layer.
%        for doy=lastdoyh+1:365 %2014 - 1 to 327; 2015 - 44 to 312
%	
%	    if exist(s4sp) & exist(s4) & exist(sp)
%		[s4spregion, s4region, spregion] = pfisr_hypothesis();
%	    end

%            phase_only=[phase_only;spregion];
%            amplitude_only=[amplitude_only;s4region];
%            phase_amplitude=[phase_amplitude;s4spregion];
                
        end % for doy
        save([outdir, 'SP_',yearstr,'_', sigstr, '_events.mat'],'phase_only');
        save([outdir, 'S4_',yearstr,'_', sigstr,'_events.mat'],'amplitude_only');
        save([outdir, 'SPS4_',yearstr,'_', sigstr,'_events.mat'],'phase_amplitude');

	if ~isempty(phase_only)
	xlswrite([outdir, 'SP_', yearstr, '_', sigstr, '.csv'], phase_only);
	end
	if ~isempty(amplitude_only)
	xlswrite([outdir, 'S4_', yearstr, '_', sigstr, '.csv'], amplitude_only);
	end
	if ~isempty(phase_amplitude)
	xlswrite([outdir, 'SPS4_', yearstr, '_', sigstr, '.csv'], phase_amplitude);
	end
%ACR are without the at least 4 receiver filter
% save(['ACR_LP_SP_',num2str(year),'_L',num2str(frequency),'_events_updated1.mat'],'phase_only');
%         save(['ACR_LP_S4_',num2str(year),'_L',num2str(frequency),'_events_updated1.mat'],'amplitude_only');
%         save(['ACR_LP_SPS4_',num2str(year),'_L',num2str(frequency),'_events_updated1.mat'],'phase_amplitude');

    end % for year
end % for frequency

toc
