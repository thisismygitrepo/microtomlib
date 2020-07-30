function [eps, sig] = view_per_freq_results(eps_data, sig_data, show_flag)

n, m, f = size(eps_data);
eps = zeros(n * 2, m * 2, f);
sig = zeros(n * 2, m * 2, f);
for i=1:f
    eps(:, :, i) = imresize(eps_data(:, :, i), 2);
    sig(:, :, i) = imresize(sig_data(:, :, i), 2);
end

if show_flag
    for i=1:f
        imagesc(eps(:, :, i));axis image;axis off;colormap(jet);
        imagesc(sig(:, :, i));axis image;axis off;colormap(jet);caxis([0.1 0.3]);view([0 -90])
        pause(0.001);
    end
end