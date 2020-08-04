function [eps, sig] = view_per_freq_results(eps_data, sig_data, show_flag)
% This function doubles the size of the image coming out from tomographic
% solvers for the purpose of better visualization.
% If show_flag is set to true, it can do an animation of results.

[n, m, f] = size(eps_data);
eps = zeros(n * 2, m * 2, f);
sig = zeros(n * 2, m * 2, f);
for i=1:f
    eps(:, :, i) = imresize(eps_data(:, :, i), 2);
    sig(:, :, i) = imresize(sig_data(:, :, i), 2);
end

if show_flag  % ==> animate the results (show them frequency by frequency)
    for i=1:f
        imagesc(eps(:, :, i));axis image;axis off;colormap(jet);
        imagesc(sig(:, :, i));axis image;axis off;colormap(jet);caxis([0.1 0.3]);view([0 -90])
        pause(0.001);
    end
end

end