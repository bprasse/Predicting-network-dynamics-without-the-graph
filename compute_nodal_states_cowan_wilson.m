function [ x, dx ] = compute_nodal_states_cowan_wilson( x_init, tau, mu, B, t_observation, ode_options )

if nargin < 6
    [ ~, x ] = ode45(@(t,y) -y + B*(1./(1+exp(-tau*( y - mu )))), t_observation, x_init);
else
    [ ~, x ] = ode45(@(t,y) -y + B*(1./(1+exp(-tau*( y - mu )))), t_observation, x_init, ode_options);
end

x = transpose( x );

dx = -x + B*(1./(1+exp(-tau*( x - mu ))));

end











