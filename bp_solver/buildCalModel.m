function contrast = buildCalModel(stp, prop)
% Requires SystemSetup + DP of Calibration Phantom (either 1 or 2).
center_x = 0;
center_y = 0;
X_radius = 199.2e-3 / 2;  % dimensions of cals.
Y_radius = 151.4e-3 / 2;  % dimensions of cals (not DOI).
eps_gama = ones(stp.Nx, stp.Ny) * stp.bg.eps_r;  % relative permit of CM everywhere.
sigma_gama = ones(stp.Nx, stp.Ny) * stp.bg.sig;  % sigma of CM everywhere.         
pho = (stp.axis_x - center_x) .^ 2 / X_radius ^ 2 + (stp.axis_y - center_y) .^ 2 / Y_radius ^ 2;  % Cal is ellipsoid.
eps_gama(pho <= 1) = prop.eps_r;  % replace cal area with relative permit of Cal1
sigma_gama(pho <= 1) = prop.sig;  % replace cal area with sigma of Cal1
eps_gama = eps_gama .* prop.eps_o;  % convert relative permitivities to absolute.
contrast = (eps_gama ./ stp.bg.eps_complex) - 1i * (sigma_gama ./ (prop.w * stp.bg.eps_complex)) - 1;   % generate contrast
contrast = contrast(:);  % vectorize the contrast.
end