function [ x, dx ] = compute_nodal_states_g_het( x_init, f_coeff, g_coeff, B, t_observation )

[ ~, x ] = ode45(@(t,y) f_dyn( y, f_coeff, g_coeff, B ), t_observation, x_init);

x = transpose( x );

dx = f_dyn( x, f_coeff, g_coeff, B );

end


function dx = f_dyn( x_mat, f_coeff, g_coeff, B )

N = size( x_mat, 1 );
n = size( x_mat, 2 );
dx = nan( N, n );

for k = 1:n
    
    x = x_mat( :, k );
    f_i = [ ones( N, 1 ), x, x.^2, x.^3 ]*f_coeff;
    
    for i = 1:N
        
        x_i = x( i );
        
        xy_vec = [ ones( N, 1 ),...
            ones( N, 1 )*x_i, x,...
            ones( N, 1 )*x_i^2, x_i.*x, x.^2,...
            ones( N, 1 )*x_i^3, x_i.^2.*x, x_i.*x.^2, x.^3 ];
        
        dx( i, k ) = f_i( i ) + B( i,: )*xy_vec*g_coeff( :, i );     
        
    end
end
end









