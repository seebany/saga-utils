function [R_rytov, k_v, index_valid, R_rytov_Bust] = Lz(V, Psi_v, az, ze, thickness, height, f,signal_type)
%Lz rytov approximation
% function [R_rytov, k_v, index_valid, R_rytov_Bust] = Lz(V, Psi_v, az, ze, thickness, height, f)
% Given: az, azimuth angle; ze, zenith angle; thickness; height;
% f array of frequencies; (SDB not sure what this pre-existing text means: nyquist frequencies corresponding to given sampling rate)
% thickness is an array of possible thicknesses L in m.
% height is an array of possible heights z in m.
% signal_type is 0 or 2 to indicate whether L1 or L2C.
% Return: R_rytov, Rytov ratio vector; k_v, corresponding wavenumber array 
% along the velocity direction.
% See also spectral
%
% Yang Su 2018
% Commented by S. Datta-Barua 9 June 2020

%convert L,z to km
% thickness = thickness;
% height = height;

%speed of light in vacuum
c = 299792458;
% Wavelength of the incident signal
switch signal_type
	case 0
		lambda = c/1575.42e6;
	case 2
		lambda = c/1227.6e6;
end		

%frequency ranges from 0 to half the sampling rate 100Hz
%     ff = [linspace(0,Fs,NFFT/2+1)]';
%     ff = [0:Fs/NFFT:Fs]';
k_v = 2 * pi * f / V;

%limit the wavenumber within [10^-3, 0.5*10^-2]
%     k_par_ind = find(k_par<=10^-1&k_par>=10^-3);
% SDB 6/9/20 I can see why we need to limit based on Nyquist f, but
% don't understand why there is a limit on k_v here.
%index_valid = find(f > 0.2 & k_v <= 0.11); 
% index_valid = find(f > 0.2);

%incident wave number of the signal
%k = 2 * pi * 1575.42 * 10^6 / c;
k = 2 * pi / lambda;
%     az is measured from north in North East Down Coords
%     az= -az+90;
%     if az > 180
%     az = az-360;
%     end

% Dr. Bust's way
k_x = k_v * cos(Psi_v);
k_y = k_v * sin(Psi_v);
k_squared_Bust = repmat(k_v.^2, 1, length(az)) + ...
    (k_x.^2 * cos(az).^2 + k_y.^2 * sin(az).^2) .* repmat(tan(ze).^2, length(k_v), 1);

% Our way
k_squared = k_v.^2 * (sin(Psi_v+az).^2 .* sec(ze).^2 + cos(Psi_v+az).^2);

%effective thickness and height along signal ray path
thickness_effective = thickness * repmat(sec(ze), length(k_v), 1);
height_effective = height * repmat(sec(ze), length(k_v), 1);

% SDB try to clean up the representation for clarity. 6/29/20
alpha = k_squared;
L_e = thickness_effective;
z_e = height_effective;
x = L_e.*alpha./(2*k);
y = z_e.*alpha./k - x;
f = sinc(x/pi).*cos(y);
S_plus = 1+f;
S_minus = 1-f;
R_rytov = S_minus./S_plus;

ratio = k_squared ./ k;
%
ratio_Bust = k_squared_Bust ./ k;

%dummy = 2 ./ (ratio .* thickness_effective) ...
%    .* sin(.5*ratio.*thickness_effective) ...
%    .* cos(ratio.*(height_effective - .5 * thickness_effective));

dummy_Bust = 2 ./ (ratio_Bust .* thickness) ...
    .* sin(.5*ratio_Bust.*thickness) ...
    .* cos(ratio_Bust.*(height - thickness)/2);

%if thickness == 0
%    S_minus = sin(ratio.*height_effective).^2;
%    S_plus = cos(ratio.*height_effective).^2;
%else
%    S_minus = 1 - dummy;
%    S_plus = 1 + dummy;
%end

S_minus_Bust = 1 - dummy_Bust;
S_plus_Bust = 1 + dummy_Bust;

%R_rytov = S_minus ./ S_plus;
R_rytov_Bust = S_minus_Bust ./ S_plus_Bust;
%R_rytov_Bust = [];

% Trying a new way to get the valid fitting range of kappas.SDB 6/19/20
%kappa_min_ind = find(S_plus == max(S_plus));
%kappa_max_ind = find(S_minus == max(S_minus));
%keyboard
%kappa_min_ind = find(R_rytov == min(R_rytov));
%kappa_max_ind = find(R_rytov == max(R_rytov));
%index_valid = [kappa_min_ind:kappa_max_ind];
index_valid = [];

%[k_v(kappa_min_ind) k_v(kappa_max_ind)]
%[max(S_plus) max(S_minus)]
%kappa_min_ind = min(find(k_v > 2*1e-3));
%kappa_max_ind = min(find(k_v > 2*2e-2));
%index_valid = [kappa_min_ind:kappa_max_ind];
%[k_v(kappa_min_ind) k_v(kappa_max_ind)];

%     loglog(k_v, R_rytov, k_v, R_rytov_Bust);
% Difference between Bust and Ours
% diff_sqr = sumsqr(R_rytov-R_rytov_Bust);
%     fprintf('%f\n', diff_sqr);

end
