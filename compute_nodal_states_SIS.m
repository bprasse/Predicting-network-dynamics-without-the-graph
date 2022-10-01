function [ x, dx ] = compute_nodal_states_SIS( x_init, delta, B, t_observation, ode_options )

if nargin < 5
    [ ~, x ] = ode45(@(t,y) -delta.*y + ( 1 - y ).*( B*y ), t_observation, x_init );
else
    [ ~, x ] = ode45(@(t,y) -delta.*y + ( 1 - y ).*( B*y ), t_observation, x_init, ode_options );
end

x = transpose( x );

dx = -delta.*x + ( 1 - x ).*( B*x );

end

