function Madrigal= load_pfisr_matfile(year, month, day, kindat, mypfisrdir)

% function Madrigal= load_pfisr_matfile(year, month, day, kindat, mypfisrdir)
% Created by V. Sreenivash 2018
% Commented by S. Datta-Barua
% 20 Feb 2019
% 4 Apr 2019 Once V's paper is done, change the file naming convention to
% be yy_mm_dd.mat instead of dd_mm_yy.mat.  Also use datestr to create, so
% you don't have to zero-pad the strings.
% Those changes have been done.
%
% This function returns a variable containing the PFISR data for the
% user-specified UT date.  If the .mat file from which it is loaded doesn't
% exist, then the function calls download_madrigal_pfisr.m which calls
% globalIsprint.m to download the data from the website.
%
% 17 Nov 2020 I accessed Vaishnavi's Madrigal downloads to reproduce her result.
% 8 Dec 2020 Commented to now access Pau's downloads to reproduce his result.

data=[];

switch kindat
    case 5950
        kindstr = '';
    case 5951
        kindstr = 'ac';
end
%keyboard
% Check if the matfile already exists in someone else's folder.
try
	% First try working only with Vaishnavi's original files. SDB 11/17/20
    	%studentdir = '/data1/home/vsreeni1';
	%infilename = [studentdir, 'Madrigal', kindstr, '_', ...
	%	datestr(datenum([year, month, day]), 'dd_mm_yyyy'), '.mat'];
    	%%studentdir='/data1/home/pllado/matfiles/madrigal_data/';
        %%infilename=[studentdir, 'Madrigal', kindstr, '_', ...
	%%	datestr(datenum([year, month, day]), 'yyyy_mm_dd'), '.mat'];
	%load(infilename);
    	%cmdstr = ['Madrigal = Madrigal' kindstr ';'];
    	%eval(cmdstr);

	% Next, try working only with Pau's original files. SDB 12/8/20
    	studentdir='/data1/home/pllado/matfiles/madrigal_data/';
        infilename=[studentdir, 'Madrigal', kindstr, '_', ...
		datestr(datenum([year, month, day]), 'yyyy_mm_dd'), '.mat'];
	load(infilename);
    	cmdstr = ['Madrigal = Madrigal' kindstr ';'];
    	eval(cmdstr);
catch
	disp('Not available in student folder.');

	% Otherwise extract from the already downloaded text files I split up.
	try %catch
		
        	infilename=[mypfisrdir, 'Madrigal', kindstr, '_', ...
		datestr(datenum([year, month, day]), 'yyyy_mm_dd'), '.mat'];
		load(infilename);
    		cmdstr = ['Madrigal = Madrigal' kindstr ';'];
		eval(cmdstr);
	% else download from madrigal website for the first time, then save.
	catch
	disp('Not available in specified folder.  Downloading from Madrigal.')
		pubpfisrdir = '/data2/public/Data/pfisr/';
        	filename = download_madrigal_pfisr(year, month, day, kindat, pubpfisrdir);
        	data = load(filename);
        	rows=find(data(:,12)==day&data(:,17)==year&data(:,19)==month);
            	eval(['Madrigal' kindstr '=data(rows,:);'])
    		cmdstr = ['Madrigal = Madrigal' kindstr ';'];
		eval(cmdstr);
        	outfilename=[mypfisrdir, 'Madrigal', kindstr, '_', ...
		datestr(datenum([year, month, day]), 'yyyy_mm_dd'), '.mat'];
    		save(outfilename, ['Madrigal', kindstr]);
    	%save(filesave, ['Madrigal', kindstr]);
	end % try loading from elsewhere, catch download for th first time.
end
