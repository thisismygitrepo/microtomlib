function [eps_p, eps_pp] = getProperties(data_path, freq) 

    % Input:
    % 1- a path to a text or csv file (no headers!)
      % Assumes that data_path has three columns
      % Freq (GHz units), relative complex permitivity (eps', eps'').
    
    % 2- Frequency in Hz
      
    % Output:
    % Returns relative complex permitivity at a given freq.

    freq = freq * 1e-9;  % convert given Hz freq to GHz.

    data = readmatrix(data_path);
    eps_p_func = fit(data(:, 1), data(:, 2), 'poly2');
    eps_pp_func = fit(data(:, 1), data(:, 3), 'poly2');
    
    eps_p = eps_p_func(freq);
    eps_pp = eps_pp_func(freq);
end


