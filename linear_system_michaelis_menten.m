function [ V_i, F_i ] = linear_system_michaelis_menten( x, dx, h, node )

V_i = transpose( dx( node, : ) );
V_i = V_i + transpose( x( node, : ) );
F_i = transpose( x.^h )./transpose( 1 + x.^h );

end

