function [tstt, tend] = zoomin_hrplot_iq_carrier(irow, tspan, dt, signal_type, home_dir, cases_folder, year, doy, prn, rcvr_op,fluc)
%zoomin_hrplot_iq_carrier plot and save zoom-in figures for high-rate IQ power &
%phase and carrier phase
%tspan: datenum, time interval for hr data processing
%dt: length of the zoom-in data segment

zcounter = irow;
set_plot = 'A'
for ts = tspan(1):dt:tspan(end) - mod(diff(tspan), dt)
    if ts == tspan(end) - mod(diff(tspan), dt)
        te = ts + mod(diff(tspan), dt);
    else
        te = ts + dt;
    end
    tt = [ts; te];
    datevec(tt)
    zcounter;
%          [tstt, tend] = hrplot_iq_carrier(signal_type,home_dir,cases_folder,year,doy,prn, ...
%         tt,rcvr_op,zcounter,'B'); 
     [tstt, tend] = hrplot_iq_carrier2(signal_type, home_dir, cases_folder, year, doy, prn, ...
          tt, rcvr_op, zcounter, set_plot,fluc);  %noise in lz estimation
    %[tstt, tend] = hrplot_iq_carrier(signal_type, home_dir, cases_folder, year, doy, prn, ...
        tt, rcvr_op, zcounter, set_plot,fluc);  %NO NOISE LZ SIMULATION
    %     zcounter = zcounter + 1;
end

end
