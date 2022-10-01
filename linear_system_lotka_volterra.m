function [ V_i, F_i ] = linear_system_lotka_volterra( x, dx, alpha_i, theta_i, node )

N = size( x, 1 );

V_i = transpose( dx( node, : ) );
f_i = transpose( x( node, : ).*( alpha_i - theta_i*x( node, : ) ));
V_i = V_i - f_i;

mat_tmp = transpose( x( node, : ) )*ones( 1, N );

F_i = -transpose( x ).*mat_tmp;

end

