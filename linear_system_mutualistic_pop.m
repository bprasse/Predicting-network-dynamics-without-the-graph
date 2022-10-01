function [ V_i, F_i ] = linear_system_mutualistic_pop( x, dx, alpha_i, theta_i, node )

N = size( x, 1 );

V_i = transpose( dx( node, : ) );
f_i = transpose( x( node, : ).*( alpha_i - theta_i*x( node, : ) ));

V_i = V_i - f_i;

mat_tmp = transpose( x( node, : ) )*ones( 1, N );

F_i = (transpose( x.^2 )./transpose( 1 + x.^2 )).*mat_tmp;

end

