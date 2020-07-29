function process_data(eps_data, sigma_data, freq_low, freq_high, show_flag)
% This function is supposed to do the following:
   % 1- Rescale
   % 2- 

%     aaa = imresize(XX_rec, 2);  It is not the job of this function to do
%     this. (Separation of concerns principal).
    

%     imagesc(eps_CSI(:, :, cnt - 1));axis image;axis off;colormap(jet);
% %     imagesc(sigma_CSI(:, :, cnt - 1));axis image;axis off;colormap(jet);caxis([0.1 0.3]);view([0 -90])
%     pause(0.001);
%     
%     toc

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

