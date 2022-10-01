function [ V_i, F_i ] = linear_system_cowan_wilson( x, dx, tau, mu, node )

V_i = transpose( dx( node, : ) + x( node, : ) );

F_i = 1./transpose( 1 + exp( -tau*(x - mu ) ) );

end

