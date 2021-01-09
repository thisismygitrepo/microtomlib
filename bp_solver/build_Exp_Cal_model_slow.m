function [eps_gama, sigma_gama] = build_Exp_Cal_model(x_dash, y_dash, eps_r, sigma, eps_r_b, sigma_b)

center_x = 0;
center_y = 0;

X_radius = 199.2e-3 / 2;
Y_radius = 151.4e-3 / 2;

Nx = length(x_dash);
Ny = length(y_dash);

[Axis_x, Axis_y] = ndgrid(x_dash, y_dash);

eps_gama = zeros(Nx, Ny);
sigma_gama = zeros(Nx, Ny);

for ii = 1 : Nx
    for jj = 1 : Ny
        
        xx = Axis_x(ii, jj);
        yy = Axis_y(ii, jj);
        
        pho = (xx - center_x) ^ 2 / X_radius ^ 2 + (yy - center_y) ^ 2 / Y_radius ^ 2;
        
        if pho <= 1
            eps_gama(ii, jj) = eps_r;
            sigma_gama(ii, jj) = sigma;
        else
            eps_gama(ii, jj) = eps_r_b;
            sigma_gama(ii, jj) = sigma_b;
        end
    end
end
        
        
