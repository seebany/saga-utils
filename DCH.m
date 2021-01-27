% This script loops through the Detect, Classify, Hypothesize (DCH) method
% for estimating the scattering layer as described in Sreenivash et al., 2020.
%

clear all
clc
%frequency=  1 for L1; 2 for L2
%addpath(genpath('~/'));
outdir = '~/matfiles/saga_dch/';
signal_type = [0,2];
lronlyflag = 1;
tic
% Loop over frequency, year, doy.
for frequency=2:2
    sigstr = getSignal(signal_type(frequency));
    for year=2015:2015
	yearstr = num2str(year);
        
        if exist([outdir, 'SP_',yearstr,'_', sigstr, '_events.mat'],'file')
		load([outdir, 'SP_',yearstr,'_', sigstr, '_events.mat'],'phase_only');
        	load([outdir, 'S4_',yearstr,'_', sigstr,'_events.mat'],'amplitude_only');
        	load([outdir, 'SPS4_',yearstr,'_', sigstr,'_events.mat'],'phase_amplitude');
	else
		phase_only=[];
        	amplitude_only=[];
        	phase_amplitude=[];
	end

	inputdir = ['~/matfiles/saga_dch/', yearstr, '/'];
	if ~(exist(inputdir, 'dir'))
		system(['mkdir ' inputdir]);
	end

	if ~isempty(phase_only)
		lastdoy = phase_only(end,2);
	else
		lastdoy = 0;
	end
        for doy=lastdoy+1:365 %2014 - 1 to 327; 2015 - 44 to 312
            disp(['Executing main to detect scintillation for ',yearstr,', ',num2str(doy,'%03i'),', freq ',num2str(frequency)]);
            
            if (exist([inputdir, sigstr,'_S4_',yearstr,'_',num2str(doy,'%03i'),'.mat']) ==2 & exist([inputdir, sigstr,'_SP_',yearstr,'_',num2str(doy,'%03i'),'.mat'])==2)
                
                display([sigstr,'_S4 and SP_',yearstr,'_',num2str(doy,'%03i'),'.mat already existed']);
                %              display('Calling AC_LP_Vaishnavi_scint')
                
            else
                display([inputdir, sigstr,'_S4_',yearstr,'_',num2str(doy,'%03i'),'.mat', ' not found and/or ']);
                display([inputdir, sigstr,'_SP_',num2str(year,'%04i'),'_',num2str(doy,'%03i'),'.mat', ' not found']);
                %if frequency==1
                %    signal_type=0;
                %else signal_type=2;
                %end
    		for fluct = 0:1
                	main(year,doy,lronlyflag,signal_type(frequency),fluct);
                	%main(year,doy,lronlyflag,0,signal_type);
    		end % for fluct
                
            end % if S4 and SP_*.mat files exist
	    % Next, classify and hypothesize layer.
            if (exist([inputdir, sigstr,'_S4_',yearstr,'_',num2str(doy,'%03i'),'.mat']) ==2 & exist([inputdir, sigstr,'_SP_',yearstr,'_',num2str(doy,'%03i'),'.mat'])==2)
                
                [s4sp,s4,sp,s4spregion,s4region,spregion]=classify_hypothesize(year,doy,signal_type(frequency), inputdir);
%                 s4region=region2num(s4region);
%                 spregion=region2num(spregion);
%                 s4spregion=region2num(s4spregion);
                spcat=horzcat(sp,spregion); %concatenating sp and sp region
                s4cat=horzcat(s4,s4region);
                s4spcat=horzcat(s4sp,s4spregion);
                phase_only=[phase_only;spcat];
                amplitude_only=[amplitude_only;s4cat];
                phase_amplitude=[phase_amplitude;s4spcat];
            else
                
            end
        end % for doy
        save([outdir, 'SP_',yearstr,'_', sigstr, '_events.mat'],'phase_only');
        save([outdir, 'S4_',yearstr,'_', sigstr,'_events.mat'],'amplitude_only');
        save([outdir, 'SPS4_',yearstr,'_', sigstr,'_events.mat'],'phase_amplitude');

	xlswrite([outdir, 'SP_', yearstr, '_', sigstr, '.csv'], phase_only);
	xlswrite([outdir, 'S4_', yearstr, '_', sigstr, '.csv'], amplitude_only);
	xlswrite([outdir, 'SPS4_', yearstr, '_', sigstr, '.csv'], phase_amplitude);
%ACR are without the at least 4 receiver filter
% save(['ACR_LP_SP_',num2str(year),'_L',num2str(frequency),'_events_updated1.mat'],'phase_only');
%         save(['ACR_LP_S4_',num2str(year),'_L',num2str(frequency),'_events_updated1.mat'],'amplitude_only');
%         save(['ACR_LP_SPS4_',num2str(year),'_L',num2str(frequency),'_events_updated1.mat'],'phase_amplitude');

    end % for year
end % for frequency

toc
