function [ V_i, F_i ] = linear_system_g_het( x_mat, dx, f_coeff, g_coeff, node )

N = size( x_mat, 1 );
n = size( x_mat, 2 );

V_i = transpose( dx( node, : ) );
F_i = nan( n, N );

for k = 1:n
    
    x = x_mat( :, k );
    x_i = x( node );
    f_i = [ 1, x_i, x_i.^2, x_i.^3 ]*f_coeff;
    
    xy_vec = [ ones( N, 1 ),...
        ones( N, 1 )*x_i, x,...
        ones( N, 1 )*x_i^2, x_i.*x, x.^2,...
        ones( N, 1 )*x_i^3, x_i.^2.*x, x_i.*x.^2, x.^3 ];
    
    
    V_i( k ) = V_i( k ) - f_i;
    
    F_i( k, : ) = transpose( xy_vec*g_coeff( :, node ) );
    
end
end

