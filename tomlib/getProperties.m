function [eps_p, eps_pp] = getProperties(data_path, freq) 
    freq = freq * 1e-9;
    % Takes in a path to a text or csv file (no headers!)
    % Returns dielctric properties at a given freq with units that
    % match the frequency units of the file passed.

    % Assumes that data_path has three columns
    % Freq, eps_p and eps_pp
    
    data = readmatrix(data_path);
    eps_p_func = fit(data(:, 1), data(:, 2), 'poly2');
    eps_pp_func = fit(data(:, 1), data(:, 3), 'poly2');
    
    eps_p = eps_p_func(freq);
    eps_pp = eps_pp_func(freq);
end


