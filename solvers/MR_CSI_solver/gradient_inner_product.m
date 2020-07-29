function [production] = gradient_inner_product(Ax, Ay, Bx, By)

production = Ax.' * conj(Bx) + Ay.' * conj(By);