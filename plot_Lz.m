function [h] = plot_Lz(Lz_filename)
%plot_Lz plots rectangles enclosing the estimated thickness L about
%estimated height z for L_mean, z_mean, ts, te datenums provided in the
%matfile given by path string Lz_filename.  This is for use with the data
%Aurora Lopez provided for the 12/8/13 event in file lzdata2013doy342.mat.
%
%   S. Datta-Barua
% 17 Sept 2019

load(Lz_filename);
% Loop through each estimation time.
for i = 1:numel(ts)
    x = [ts(i) te(i) te(i) ts(i) ts(i)]';
    y = [z_mean(i) - L_mean(i);
        z_mean(i) - L_mean(i);
        z_mean(i);
        z_mean(i);
        z_mean(i) - L_mean(i)];
    plot(x, y, 'k')
    hold on
    clear x y

end

