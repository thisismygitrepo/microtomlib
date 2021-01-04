function [production] = gradient_inner_product(Ax, Ay, Bx, By, dx, dy)

production = dx * dy * (Ax.' * conj(Bx) + Ay.' * conj(By));