function [ x, dx ] = compute_nodal_states_kuramoto( x_init, omega, B, t_observation, ode_options )

%using the identity sin( x - y ) = sin( x )*cos( y ) - cos( x )*sin( x ) to
%speed up computation:, ode_options )

if nargin < 5
[ ~, x ] = ode45(@(t,y) omega + sin( y ).*(B*cos( y )) - cos( y ).*(B*sin( y )), t_observation, x_init);
else
[ ~, x ] = ode45(@(t,y) omega + sin( y ).*(B*cos( y )) - cos( y ).*(B*sin( y )), t_observation, x_init, ode_options);
end

x = transpose( x );

dx = omega + sin( x ).*(B*cos( x )) - cos( x ).*(B*sin( x ));

end











