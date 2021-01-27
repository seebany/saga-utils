% This function splits a multi-day Madrigal PFISR file downloaded 
% into one .mat file per day. Adapted from Vaishnavi's script
% Automated_Madrigal_File_for_a_day_for1year.m to use with her
% downloaded data.
%
% Created by Vaishnavi Sreenivash
% 27 March 2019
% Modified by S. Datta-Barua

function [] = split_madrigal_data_daily(inputfile, outputpath)

data=[];
tic
data=load(inputfile);
%data=load('Madrigal13-15.txt');
toc
disp(['Madrigal data loaded from ', inputfile])
% Check for code type as listed in filename.
if ~isempty(strfind(inputfile,'ac.txt'))
	kindstr = 'ac';
else
	kindstr = '';
end

% Loop through year, doy.
for year=2014:2015
	for doy=1:365
	    % Convert year, doy to a Matlab datenum.
	    daten = datenum([year 0 doy]);
	    eval(['Madrigal' kindstr '=[];'])
	    [yr, month, day] = datevec(daten);
%	    RINEX.year=year;
%	    RINEX.day=doy;
%[gps,abs,utc]=convert_rinex(RINEX,0); %leap =0 as 2014 and 2015 are not leap years
%month=utc.mon;
%day=utc.day;
%Oct 7 2015
%rows=find(data(:,12)==7&data(:,17)==2015&data(:,19)==10);
%rows=find(data(:,12)==19 &data(:,17)==2015&data(:,19)==3);
	    rows = find(data(:,12) == day & data(:,17) == year & data(:,19) == month);
	    if ~isempty(rows)
		eval(['Madrigal' kindstr '=data(rows,:);'])
	    end
	    save([outputpath, 'Madrigal', kindstr,'_', datestr(daten, 'yyyy_mm_dd'),'.mat'],['Madrigal',kindstr]);
	end
end
