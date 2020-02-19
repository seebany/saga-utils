year=2015; 
frequency=2; % 1 for L1; 2 for L2
phase_only=[];
amplitude_only=[];
phase_amplitude=[];
for doy=44:312 %2014 - 1 to 327; 2015 - 44 to 312
    if (exist(['L',num2str(frequency),'_S4_',num2str(year,'%04i'),'_',num2str(doy,'%03i'),'.mat']) ==2 & exist(['L',num2str(frequency),'_SP_',num2str(year,'%04i'),'_',num2str(doy,'%03i'),'.mat'])==2)
    [s4sp,s4,sp,s4spregion,s4region,spregion]=AC_LP_Vaishnavi_Scint(year,doy,frequency);
    spcat=horzcat(sp,spregion); %concatenating sp and sp region
    s4cat=horzcat(s4,s4region);
    s4spcat=horzcat(s4sp,s4spregion);
    phase_only=[phase_only;spcat];
    amplitude_only=[amplitude_only;s4cat];
    phase_amplitude=[phase_amplitude;s4spcat];
    else 
        continue;
    end
    
end
save('AC_LP_SP_2015_L2_events_updated1.mat','phase_only');
save('AC_LP_S4_2015_L2_events_updated1.mat','amplitude_only');
save('AC_LP_SPS4_2015_L2_events_updated1.mat','phase_amplitude');
% save('trial11.mat','phase_only');
% save('trial22.mat','amplitude_only');
% save('trial33','phase_amplitude');
