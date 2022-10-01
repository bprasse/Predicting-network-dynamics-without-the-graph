function [ x, dx ] = compute_nodal_states_michaelis_menten( x_init, h, B, t_observation, ode_options )

if nargin < 5
    [ ~, x ] = ode45(@(t,y) -y + B*( y.^h./( 1 + y.^h )), t_observation, x_init);
else
    [ ~, x ] = ode45(@(t,y) -y + B*( y.^h./( 1 + y.^h )), t_observation, x_init, ode_options);
end

x = transpose( x );

dx = -x + B*( x.^h./( 1 + x.^h ));

end











