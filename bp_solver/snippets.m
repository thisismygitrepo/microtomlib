
function [eps, sig] = getProperties(data_path, freq)
    % Takes in a path to a text or csv file (no headers!)
    % Returns dielctric properties at a given freq with units that
    % match the frequency units of the file passed.

    % Assumes that data_path has three columns
    % Freq, eps_p and eps_pp
    
    data = readmatrix(data_path);
    eps_func = fit(data(:, 1), data(:, 2), 'poly2');
    sig_func = fit(data(:, 1), data(:, 3), 'poly2');
    
    eps = eps_func(freq);
    sig = sig_func(freq);
end



pm = path_manager();
%%% Fit the EPs of coupling medium
[eps_r_b, sigma_b] = getProperties(pm.gdrive + "/emvision/Algorithm/toml_data/system/cm_props", freq);
eps_b = eps_r_b * eps_o - 1i * sigma_b / w;
kb = w * sqrt(uo * eps_b);

%%% Fit the EPs of Cal 1
[Cal1_eps, Cal1_sigma] = getProperties(pm.gdrive + "/emvision/Algorithm/toml_data/system/cal1_props", freq);

%%% Fit the EPs of Cal 2
[Cal2_eps, Cal2_sigma] = getProperties(pm.gdrive + "/emvision/Algorithm/toml_data/system/cal1_props", freq);
