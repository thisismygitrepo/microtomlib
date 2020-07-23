function process_data(eps_data, sigma_data, freq_low, freq_high, show_flag)

x_dash = [-115 : 2 : 115 + 2] .* 1e-3;
y_dash = [-100 : 2 : 96 + 2] .* 1e-3;

[Axis_x, Axis_y] = ndgrid(x_dash, y_dash);

eps_bleeding = sum(eps_data(:, :, freq_low : freq_high), 3) ./ (freq_high - freq_low + 1);
sigma_bleeding = sum(sigma_data(:, :, freq_low : freq_high), 3) ./ (freq_high - freq_low + 1);

eps_bleeding = imgaussfilt(eps_bleeding, 2);
sigma_bleeding = imgaussfilt(sigma_bleeding, 2);


if show_flag
	figure
	imagesc(Axis_y(:), Axis_x(:), eps_bleeding);axis image;view([0 -90]);axis off;
	
	figure
	imagesc(Axis_y(:), Axis_x(:), sigma_bleeding);axis image;view([0 -90]);axis off;caxis([0.15 0.24])
end

