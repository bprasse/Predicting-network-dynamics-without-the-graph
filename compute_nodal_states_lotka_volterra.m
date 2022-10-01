function [ x, dx ] = compute_nodal_states_lotka_volterra( x_init, alpha, theta, B, t_observation, ode_options )

if nargin < 6
    [ ~, x ] = ode45(@(t,y) y.*( alpha - ( B + diag( theta ))*y ), t_observation, x_init);
else
    [ ~, x ] = ode45(@(t,y) y.*( alpha - ( B + diag( theta ))*y ), t_observation, x_init, ode_options);
end

x = transpose( x );

dx = x.*( alpha - ( B + diag( theta ))*x );

end











