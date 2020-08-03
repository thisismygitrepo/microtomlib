function view_average_results(eps_data, sig_data, freq_low, freq_high)
% This function (previously called process_data) is supposed to do the following:
   % 1- Upscale the images (by calling a function internally)
   % 2- Average the results from different frequencies
   % 3- View the results nicely in one Figure using subplots.
   % 4- Adjsut the settings of colormap to show target clearly.    

[eps, sig] = view_per_freq_results(eps_data, sig_data, false);
 
x_dash = [-115 : 2 : 115 + 2] .* 1e-3;
y_dash = [-100 : 2 : 96 + 2] .* 1e-3;

[Axis_x, Axis_y] = ndgrid(x_dash, y_dash);

eps_bleeding = sum(eps(:, :, freq_low : freq_high), 3) ./ (freq_high - freq_low + 1);
sigma_bleeding = sum(sig(:, :, freq_low : freq_high), 3) ./ (freq_high - freq_low + 1);

eps_bleeding = imgaussfilt(eps_bleeding, 2);
sigma_bleeding = imgaussfilt(sigma_bleeding, 2);

figure
imagesc(Axis_y(:), Axis_x(:), eps_bleeding);axis image;view([0 -90]);axis off;

figure
imagesc(Axis_y(:), Axis_x(:), sigma_bleeding);axis image;view([0 -90]);axis off;caxis([0.15 0.24])


