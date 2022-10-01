function [ V_i, F_i ] = linear_system_kuramoto( x, dx, omega_i, node )

n = size( x, 2 );

V_i = transpose( dx( node, : ) );

V_i = V_i - omega_i*ones( n, 1 );

F_i = transpose( sin( x( node, : ) - x ));

end

