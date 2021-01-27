function [] = read_iono(doy, year, sep, folder_path, op_path, signal)

if strcmp(folder_path(end-7:end-1), 'ASTRArx') || strcmp(folder_path(end-16:end-10), 'ASTRArx')
    ionofilestruct = dir([folder_path, 'txt', sep, 'iono*.log']);
else
    ionofilestruct = dir([folder_path, 'txt', sep, 'ion*o*', year, '_', doy, '_*.log']);
end

if ~isempty(ionofilestruct)
    IONO = [];
    for jj = 1:size(ionofilestruct, 1)
        ionofile = ionofilestruct(jj, :);
        ionofilename = dlmread([folder_path, 'txt', sep, ionofile.name]);
	% SDB 12/9/20 Put in logic to check that timestamps in the data are correct.
	if ionofilename(1,1) == year
	        st = datestr(gps2utc(ionofilename(1, 1:2)), 'HHMM-');
	        se = datestr(gps2utc(ionofilename(end, 1:2)), 'HHMM UT');
	        %         disp([ionofile.name,' actually has data for ',st,se]);
	        IONO = [IONO; ionofilename];
	end
    end
    % SDB 12/9/20 Check that IONO is not empty because data had valid timestamps.
    if isempty(IONO)
	IONODATA = [];
    else    
    IONO = IONO(IONO(:, 1) <= 3640, :);
    IONO = sortrows(IONO, [1, 2, 3]);
    ORTW = IONO(:, 1);
    ORTS = IONO(:, 2) + IONO(:, 3);
    STEC = IONO(:, 4);
    DSTEC = IONO(:, 5);
    PRN = IONO(:, 6);
    
    IONODATA = [ORTW, ORTS, STEC, DSTEC, PRN];
    %     unique(PRN,'stable')
    end
else
    IONODATA = [];
end

outfilename = strcat('prn_files_', signal, sep, 'ionodata.mat');
save([op_path, outfilename], 'IONODATA');
end
