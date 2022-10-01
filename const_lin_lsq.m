function x_sol = const_lin_lsq( mat, vec, rho, quadprog_options )
%Constrained linear least-squares, calls solver for:
%
%   min     || mat*x - vec ||_2^2 + rho*|| x ||_1
%   s.t.    x >= 0
%
%   x_0     is the initial point
%
% NOTE: Currently only option is "quadprog" solver 

num_var = size( mat, 2 );

H = 2*( transpose( mat )*mat);

f = -2*transpose( mat )*vec+ rho*ones( num_var, 1 );

[ x_sol, ~, exitflag ] = quadprog( H, f, [], [], [], [], zeros( num_var, 1 ), [], [], quadprog_options );

if exitflag < 0 &&  exitflag ~= -3
    disp(strcat('Warning: quadprog returned exitflag==', num2str( exitflag ), '\n' ));
end

end









