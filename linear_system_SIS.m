function [ V_i, F_i ] = linear_system_SIS( I, dI, delta_i, node )

N = size( I, 1 );

V_i = transpose( dI( node, : )  ) + delta_i*transpose( I( node, : ) );

F_i = ( transpose( 1 - I( node, : ) )*ones( 1, N ) ).*transpose( I );

end

