function textfile = download_madrigal_pfisr(year, month, day, kindat, downloaddir)

% Created by V. Sreenivash 2018
% Commented by S. Datta-Barua
% 20 Feb 2019
% 20 Feb 2019 Converting to a function.
% 4 Apr 2019 Additional input kindat is numerical code for kind of data,
% where 5950 is the PFISR long pulse and 5951 is the PFISR alternating
% code.
%
% This script is based on the original Madrigal_13-15.m.  It downloads 
% PFISR data from the Madrigal website using
% globalIsprint.m provided by Madrigal maintainers (possibly 
% Bill Rideout).

% Instrument: Poker Flat IS Radar
% Kinds of Data: Long Pulse (480), Long Pulse (480)
% Experiment Name: All experiment names accepted
% StartDate = 7/12/2013
% EndDate = 1/1/2016
% Seasonal filter = 1/1 - 31/12 (no seasonal filter)
% Data filters:
% No filters entered
% 
% Parameters displayed: AZM, MIN, DAYNO, GDALT, GDLAT, HOUR, TR, NE, SEC, TI, TE, DAY, GLON, ELM, VO, RANGE, YEAR, UT, MONTH, KP

%lladodir='/data1/home/pllado/matfiles/madrigal_data/';
switch kindat
	case 5950
		kindstr = '';
	case 5951
		kindstr = 'ac';
end
textfile = [downloaddir,'Madrigal' kindstr datestr(datenum([year, month, day]), 'yymmdd') '.txt'];
%textfile = 'Madrigal13-15_ac_part2.txt';
%keyboard
globalIsprint('http://isr.sri.com/madrigal', ...
 'AZM,MIN,DAYNO,GDALT,GDLAT,HOUR,TR,NE,SEC,TI,TE,DAY,GLON,ELM,VO,RANGE,YEAR,UT,MONTH,KP', ...
 textfile, ...
 'Seebany Datta-Barua', ...
 'sdattaba@iit.edu', ...
 'IIT', ...
  datenum([year, month, day]),...% datenum('30-Apr-2014 00:00:00'), ... % datenum('07-Dec-2013 00:00:00'), ... % 
  datenum([year, month, day+1]),... %datenum('01-Jan-2016 23:59:59'), ... % 
 61, ...        % Instrument code is 61 for PFISR
 '', ...        % Optional Filters.  To select see comments in globalIsprint.m.
 kindat, ...%5951, ...%[5950], ...    % Optional Kinds of data.  ID 5950 is "long pulse (480)". 5951 is "Alternating Code (AC16-30)"
 '', ...        % Optional Experiment name.
 '');        % Optional File description.
