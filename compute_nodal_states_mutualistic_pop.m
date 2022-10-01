function [ x, dx ] = compute_nodal_states_mutualistic_pop( x_init, alpha, theta, B, t_observation, ode_options )

if nargin < 6
[ ~, x ] = ode45(@(t,y) y.*(( alpha - theta.*y ) + ( B*( y.^2./( 1 + y.^2 )))), t_observation, x_init);
else
[ ~, x ] = ode45(@(t,y) y.*(( alpha - theta.*y ) + ( B*( y.^2./( 1 + y.^2 )))), t_observation, x_init, ode_options );
end

x = transpose( x );

dx = x.*(( alpha - theta.*x ) + ( B*( x.^2./( 1 + x.^2 ))));

end











