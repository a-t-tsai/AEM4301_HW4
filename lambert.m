r0_1 = [0.5; 0.6; 0.7];
r1_1 = [0.0; 1.0; 0.0];
delta_t_1 = 0.9667663;
[v0_1, v1_1] = lambert_solver(r0_1, r1_1, delta_t_1, -1);
fprintf('v0long: [%.3g, %.4g, %.4g], v1long: [%.3g, %.3g, %.3g]\n', v0_1(1), v0_1(2), v0_1(3), v1_1(1), v1_1(2), v1_1(3));
r0_2 = [1.0; 0.0; 0.0];
r1_2 = [1.0; 0.125; 0.125];
delta_t_2 = 0.125;
[v0_2, v1_2] = lambert_solver(r0_2, r1_2, delta_t_2, 1);
fprintf('v0short: [%.3g, %.4g, %.4g], v1short: [%.3g, %.3g, %.3g]\n', v0_2(1), v0_2(2), v0_2(3), v1_2(1), v1_2(2), v1_2(3));

function [v0, v1] = lambert_solver(r0_vec, r1_vec, delta_t, dir)
mu = 1.0; % Default gravitational parameter
% Calculate magnitudes of the position vectors
r0 = norm(r0_vec);
r1 = norm(r1_vec);
% Calculate the angle between the two position vectors
delta_theta = acos(dot(r0_vec, r1_vec) / (r0 * r1));
if dir == 1
%short way
sign = 1; %"direction of motion"
else
delta_theta = 2 * pi - delta_theta; 
%long way
sign = -1;
end
% Calculate parameter A using the Lambert's equation
A = sign * sqrt(r0 * r1 * (1 + cos(delta_theta)));
z0 = 0;
z = z0;
relErr = 1; 
tol = 1e-5; 
n = 0; 
nMax = 200;
while relErr > tol
if z == 0 
% No division by 0
S = 1/6; Sprime = -1/120;
C = 1/2; Cprime = -1/24;
else
C = (1 - cos(sqrt(z))) / z;
S = (sqrt(z) - sin(sqrt(z))) / (z^1.5);
Sprime = (1 / (2 * z)) * (C - 3 * S);
Cprime = (1 / (2 * z)) * (1 - z * S - 2 * C);
end
% Given Functions
% Calculate y(z) as per Lambert's problem
Y = r0 + r1 - A * (1 - z * S) / sqrt(C);
% Calculate X and y(z) after finding z
X = sqrt(Y / C);
% Calculate the u function
U = (1 / sqrt(mu)) * (X^3 * S + A * sqrt(Y)) - delta_t;
% Derivative of U and Calculates the v function
V = 1 / sqrt(mu) * (X^3 * (Sprime - (3 * S * Cprime) / (2 *C)) + (A / 8) * ((3 * S * sqrt(Y)) / C + (A / X)));
% Get new z value and relative error
z = z0 - U / V;
relErr = abs((z - z0) / z0);
z0 = z;
n = n + 1;
if n > nMax
break;
end
end
% Calculate the coefficients for velocity
f = 1 - Y / r0;
g = A * sqrt(Y / mu);
g_dot = 1 - Y / r1;
% Initial and final velocity vectors
v0 = (r1_vec - f * r0_vec) / g;
v1 = (g_dot * r1_vec - r0_vec) / g;
end
