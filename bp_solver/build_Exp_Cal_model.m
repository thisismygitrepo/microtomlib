function [eps_gama, sigma_gama] = build_Exp_Cal_model_expedited(x_dash, y_dash, eps_r, sigma, eps_r_b, sigma_b)

center_x = 0;
center_y = 0;

X_radius = 199.2e-3 / 2;
Y_radius = 151.4e-3 / 2;

Nx = length(x_dash);
Ny = length(y_dash);

[Axis_x, Axis_y] = ndgrid(x_dash, y_dash);

eps_gama = ones(Nx, Ny) * eps_r_b;
sigma_gama = ones(Nx, Ny) * sigma_b;    

        
pho = (Axis_x - center_x) .^ 2 / X_radius ^ 2 + (Axis_y - center_y) .^ 2 / Y_radius ^ 2;


eps_gama(pho <= 1) = eps_r;
sigma_gama(pho <= 1) = sigma;

        
end
        
