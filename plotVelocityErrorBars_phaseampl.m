function varargout = plotVelocityErrorBars_phaseampl(year,doy,prn,init_time_in,xtime_in)
close all;
dbstop if error;
init();
global sep;
sep = filesep;
%where cases data are stored
global cases_folder;
cases_folder = '/data1/public/Data/cases/pfrr/';

%where output data and plots are stored
global home_dir;
global mat_dir;
[~, dummy] = system('echo $HOME');
home_dir = [dummy(1:end-1), sep];

if strcmp(cases_folder(end-4:end-1), 'pfrr')
    %path for Poker Flat data
    mat_dir = ['PFRR_Data', sep];
    mat_dir = ['matfiles', sep];
else
    %folder_path for 2013 Calgary data
    mat_dir = ['Calgary_Data', sep];
end

init_time = datenum(init_time_in);
xtime = xtime_in;
tt = init_time + xtime / 24 / 3600;
dt = diff(tt) * 1.1;
tspan = tt;

for ts = tspan(1):dt:tspan(end) - mod(diff(tspan), dt)
    if ts == tspan(end) - mod(diff(tspan), dt)
        te = ts + mod(diff(tspan), dt);
    else
        te = ts + dt;
    end
    tt = [ts; te];
end

tspan_utc = datevec(tt);
tstt = tspan_utc(1, :);

close all
figure
plotSAGAvsPFISR3(prn, tstt, 'debug');
end
