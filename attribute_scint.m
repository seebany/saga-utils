function layer= attribute_scint(year,month,day,starttime,endtime)

% This function is a modification of Updated_Test_Regiondetection.m as I (SDB) go 
% through, check that the algorithm is correct, and comment the code prior 
% to manuscript submission.
% Created by V. Sreenivash, 2018
%
% Inputs are scalars: year, month, day, starttime is UT decimal hour,
% endtime is UT decimal hour.
%
% Output is: 
% layer, a single char 'N' (no PFISR data), 'I' (inconclusive), 
% 'E' (E layer), 'F' (F layer).
%
% S. Datta-Barua
% 20 Feb 2019 
% 5 Apr 2019 Putting on the server and trying to merge Vaishnavi's edits
% for AC and LP analysis with my commented version. 

% Load the .mat data for the desired date. kindat is 5950 for lp, 5951 for
% ac.
Madrigal= load_pfisr_matfile(year, month, day, 5950);
Madrigalac = load_pfisr_matfile(year, month, day, 5951);
%load(['Madrigal','_',num2str(day),'_',num2str(month),'_',num2str(year),'.mat']);
%load('Madrigal19Mar2015.mat');

% Select times of PFISR data within scintillation interval.
timestamp=unique(Madrigal(:,18));
timestampac = unique(Madrigalac(:,18));

% Arrays of [datatimestamp, altitude, peakdensity] are saved for
% comparison.
maxNe_lp=[];
maxNe_ac = [];

scintrows=find(timestamp>=starttime & timestamp<endtime);
scinttime=timestamp(scintrows);
scintrowsac = find(timestampac >= starttime & timestampac < endtime);
scinttimeac = timestampac(scintrowsac);

E=0;
F=0;
I = 0;
T = 0;

% Desired beam's az, el in degrees.
az = -154.3;%14.04;
el = 77.5;%90;

% Find peak density's altitude at each time.
% LP
for i=1:length(scinttime)
    
    % Search for the beam by choice of azimuth and elevation
    
    rows=find(Madrigal(:,1)==az & Madrigal(:,14)==el & Madrigal(:,18)==scinttime(i) & Madrigal(:,4) >= 195);
    %rows=find(Madrigal(:,1)==-34.69&Madrigal(:,14)==66.09&Madrigal(:,18)==scinttime(i));% for 09Feb2014
    %rows=find(Madrigal(:,18)==scinttime(i));
    
    % Find peak density at that time.
    Maxlp = max(Madrigal(rows,8));
    % Find the altitude of the peak density.
    maxrowlp = rows(find(Madrigal(rows,8)==Maxlp));
    altitude=Madrigal(maxrowlp,4);
    % Check there are no nans.
    if isfinite(Maxlp)
    % Save the [time, peak altitude, peak density]
    maxNe_lp = [maxNe_lp; 
        [Madrigal(maxrowlp,18), altitude, Maxlp]];
    end
    
%    data=[Madrigal(rows,18),Madrigal(rows,4),Madrigal(rows,8)];
%    MaxNe=[MaxNe;data];
%    
%    for j=1:length(altitude)
%        if(altitude(j)<=200)
%            E=E+1;
%        elseif(altitude(j)>200)
%            F=F+1;
%            
%        end
%    end
end

% AC
for i = 1:length(scinttimeac)
    rowsac = find(Madrigalac(:,1) == az & Madrigalac(:,14) == el & Madrigalac(:,18) == scinttimeac(i) & Madrigalac(:,4)< 195);
    Maxac = max(Madrigalac(rowsac,8));
    maxrowac = rowsac(find(Madrigalac(rowsac,8) == Maxac));
    altitudeac = Madrigalac(maxrowac, 4);
    if isfinite(Maxac)
    % Save the [time, peak altitude, peak density]
    maxNe_ac = [maxNe_ac; 
        [Madrigalac(maxrowac, 18), altitudeac, Maxac]];
    end
end

% Compare LP and AC
if (size(maxNe_ac, 1) <= size(maxNe_lp, 1))
    % Loop through times of the shorter data set.
    for i = 1:size(maxNe_ac, 1)
        % Find the row of the nearest corresponding time in the other data set.
        nearesttrow = find(abs(maxNe_ac(i,1) - maxNe_lp(:,1)) ...
        == min(abs(maxNe_ac(i,1) - maxNe_lp(:,1))));
    
        % Compare the densities for those two nearest times.
        % If LP has the higher peak density then it's definitely F layer.
        if maxNe_ac(i,3) < maxNe_lp(nearesttrow,3)
            F = F+1;
            % If AC has the higher peak density then it's either E or T
            % layer.
        elseif maxNe_ac(i,3) >= maxNe_lp(nearesttrow,3) 
            if maxNe_ac(i,2) < 150
                E = E+1;
            else
                T = T+1;
            end
        end
    end
else
    for i = 1:size(maxNe_lp,1)
        % Find the row of the nearest corresponding time in the longer data
        % set.
        nearesttrow = find(abs(maxNe_lp(i,1) - maxNe_ac(:,1)) ...
            == min(abs(maxNe_lp(i,1)) - maxNe_ac(:,1)));
        
        % Compare the densities for those two nearest times.
        % If LP has the higher peak density then it's definitely F layer.
        if maxNe_lp(i,3) >= maxNe_ac(nearestrow,3)
            F = F+1;
        elseif maxNe_lp(i,3) < maxNe_ac(nearesttrow,3)
            % Check if it's E (<150 km) or T (150-195km) layer
            if maxNe_ac(nearesttrow,2) < 150
                E = E+1;
            else
                T = T+1;
            end
        end
    end
end


%Detecting region for a scintillation event on a day
if(E==0 && F==0 && T == 0)
    disp('No PFISR data available'); % N for "No PFISR data"
    layer='N';
    
elseif(E > F && E > T)
    disp('E region scintillation');  % E for "E region"
    layer='E';
    
elseif(F > E && F > T)
    disp('F region scintillation'); % F for "F region"
    layer='F';
    
elseif(T > E && T > F)
    disp('Transition region scintillation');
    layer = 'T';
else
    disp('region cannot be determined'); % I for "Inconclusive"
    layer='I';
end
end