function [mega_t, TSP_hr, TSP_hrv] = find_general_times(MSP, rcvr_op, spth_hr)
%function [mega_t, TSP_hr, TSP_hrv] = find_general_times(MSP, rcvr_op, spth_hr)
% Finds a list of times for which sigma_phi exceeds the threshold.
% Input data matrix MSP contains all receivers.
% rcvr_op is cell array of receiver names
% spth_hr is the sigma_phi threshold desired.
% Outputs mega_t: a num_prns cell array containing the list of common times.
% TSP_hr: for data from 3 or more rxs falling within the common time of a duration of at least 10 elements (at 0.01 s cadence, this means 0.1 s), an array [prn, datenum tstart, datenum tend,duration (number of points), median time?, length(TSP), numrxs]. 
% TSP_hrv: for data from 3 or more rxs falling within the common time of a duration of at least 10 elements, an array [prn, starthr, startmin, endhr, endmin, duration, median, length, numrxs].

% load test.mat;
%[MSP] = [time, scintillaiton indices, prn, rx]
%daily average sigmaphi
% This appears to be mean over all rxs, all prns. SDB 11/16/20.
dailymean = mean(MSP(:, 2));
MSP = sortrows(MSP, 1);
spth_hr0 = spth_hr;
% mean(MSP(:,2))
% MSP(:,1) = (MSP(:,1)-datenum([2014 2 20 0 0 0]))*24*3600;

TSP_hr = [];
TSP_hrv = [];

prnlist = 1:32;
for prn = prnlist
    prnmean(prn) = mean(MSP(MSP(:, 3) == prn, 2));
    %disp(['The average sigmaphi for ', num2str(prn), ' is ', num2str(prnmean(prn))])
    n(prn) = length(find(MSP(MSP(:, 3) == prn, 2) > dailymean));
end

%identify scintillating satellites
prnscint = prnlist(prnmean >= dailymean);
prnscint = 1:32;
figure('visible', 'off');
for prn = prnscint
    %     %specify a threshold for each PRN data, for example mean
    spth_hr = max([prnmean(prn), spth_hr0]);
    %or a static threshold 0.6
    %     spth_hr = spth_hr0;
    for rr = 1:size(rcvr_op)
        MSP_rr = MSP(MSP(:, 4) == rr & MSP(:, 3) == prn, :);
        diff(MSP(:, 1));
        if isempty(MSP_rr)
            T = [];
        else
            plot(MSP_rr(:, 1), MSP_rr(:, 2), 'k.-');
            hold on;
            T = [];
            ii = 1;
            line([min(MSP_rr(:, 1)), max(MSP_rr(:, 1))], [spth_hr, spth_hr]);
            while ii < size(MSP_rr, 1)
                %     while MSP_rr(ii,2)<spth_hr && ii < size(MSP_rr,1)
                %         ii = ii + 1;
                %     end
                %     ts = ii + 1;
                %     while MSP_rr(ii,2)>=spth_hr && ii < size(MSP_rr,1)
                %         ii = ii + 1;
                %     end
                %     te = ii;
                ts = ii;
                while ii < size(MSP_rr, 1) && (MSP_rr(ii+1, 1) - MSP_rr(ii, 1)) * 24 * 3600 <= 600 && ((MSP_rr(ii, 2) < spth_hr && MSP_rr(ii+1, 2) >= spth_hr) || (MSP_rr(ii, 2) >= spth_hr && MSP_rr(ii+1, 2) < spth_hr) || (MSP_rr(ii, 2) >= spth_hr && MSP_rr(ii+1, 2) >= spth_hr))
                    %                 while ii < size(MSP_rr,1) && (MSP_rr(ii,2)>=spth_hr && MSP_rr(ii+1,2)>=spth_hr)
                    ii = ii + 1;
                end
                te = ii;
                if (MSP_rr(te, 1) - MSP_rr(ts, 1)) * 24 * 3600 > 0
                    T = [T; MSP_rr(ts, 1); MSP_rr(te, 1)];
                    plot(MSP_rr(ts:te, 1), MSP_rr(ts:te, 2), 'r.-');
                end
                ii = ii + 1;
            end
        end
        datevec(T);
        TRX{rr} = T';
        datetick('x', 'HH:MM');
    end
    % Find common times (across array?)
    t = find_common_times(TRX);
    mega_t{prn, :} = t;
    
    tslist = t(1:2:end);
    telist = t(2:2:end);
    
    for ti = 1:length(telist)
        TSP = MSP(MSP(:, 3) == prn & MSP(:, 1) <= telist(ti) & MSP(:, 1) >= tslist(ti), :);
        rxop = length(unique(TSP(:, 4)));
        tmin = (telist(ti) - tslist(ti)) * 24 * 60;
        TSP = sortrows(TSP, 1);
        if ~isempty(TSP) && tmin >= 10
            tv = datevec(TSP([1; end], 1));
            TSP_hr = [TSP_hr; [prn, TSP([1; end], 1)', tmin, median(TSP(:, 2)), size(TSP, 1), rxop]];
            TSP_hrv = [TSP_hrv; [prn, tv(1, 4:5), tv(end, 4:5), tmin, median(TSP(:, 2)), size(TSP, 1), rxop]];
        end
    end
end
clear TRX
try
    TSP_hr = sortrows(TSP_hr, -4-1);
    TSP_hrv = sortrows(TSP_hrv, -6-1);
catch
    return;
end

% Here is where the requirement that there be three or more rxs is enforced.
% Does not exist in the version used to run Vaishnavi's result.
if length(rcvr_op(:,1))<4
    TSP_hr=[];
    TSP_hrv=[];
else
	TSP_hr = [TSP_hr, repmat(length(rcvr_op(:,1)), size(TSP_hr,1),1)];
	TSP_hrv = [TSP_hrv, repmat(length(rcvr_op(:,1)), size(TSP_hrv,1),1)];
end
close;
end
