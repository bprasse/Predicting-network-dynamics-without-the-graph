function Fig_S1_S2( seed )

%set the random number generator for reproducibility
if nargin < 1
    seed = 1;
end

rng( seed );

%the model abbreviations
model_name = { 'LV', 'MP', 'MM', 'SIS', 'KUR', 'WC' };

%the number of nodes
N = 100;

%the number of observations
num_observations = 1e3;

%the random graph model
random_network_model = 'BA'; %'ER'

%the number of networks
num_networks = 100;

%the filename to save/read the results
filename = strcat( './results/Fig_S1_S2_seed_',  num2str( seed ), '_N_', num2str( N ) );

%only run simulations if filename does not exist
if exist( strcat( filename, '.mat' ), 'file' ) == 2
    load( filename )
else
    
    %pre-allocate the contribution Phi_cp of the p-th agitation mode
    norm_cp = nan( num_networks,N, length( model_name ) );
    
    %pre-allocate the error Phi_PODm of the POD with m agitation modes
    error_POD = nan( num_networks,N, length( model_name ) );
    
    %loop over all models
    parfor model_i = 1:length( model_name )
        
        %redefine some variables for parfor loop
        norm_cp_this_parfor = nan( num_networks, N );
        error_this_parfor = nan( num_networks, N );
        
        parameters = set_parameters;
        m0_BA = parameters.m0_BA;
        m_BA = parameters.m_BA;
        
        %loop over all num_networks networks
        for run_i = 1:num_networks
            
            %the unweighted network
            A = create_network( N, random_network_model, [], [], [], m0_BA, m_BA, true );
            
            %the link weights
            link_weights = parameters.min_link_weight + ( parameters.max_link_weight - parameters.min_link_weight )*rand( nnz( A ), 1);
            
            %the weighted network
            B = zeros( N );
            B( A > 0 ) = link_weights;
            
            %the observation time samples
            t_obs = linspace( 0, parameters.T_max( model_i ), num_observations );
            
            switch model_i
                case 1
                    %Lotka-Volterra
                    
                    %the initial nodal state
                    x_init = parameters.LV.x_init_max*rand( N, 1 );
                    
                    %the parameters of the model
                    alpha = 1 + parameters.LV.sigma_alpha*( 2*rand( N, 1 ) - 1);
                    theta = 1 + parameters.LV.sigma_theta*( 2*rand( N, 1 ) - 1);
                    
                    %generate the nodal state sequence
                    [ x, dx ] = compute_nodal_states_lotka_volterra( x_init, alpha, theta, B, t_obs );
                case 2
                    %Mutualistic population
                    
                    %the initial nodal state
                    x_init = parameters.MP.x_init_max*rand( N, 1 );
                    
                    %the parameters of the model
                    alpha = 1 + parameters.MP.sigma_alpha*( 2*rand( N, 1 ) - 1 );
                    theta = 1 + parameters.MP.sigma_theta*( 2*rand( N, 1 ) - 1 );
                    
                    %generate the nodal state sequence
                    [ x, dx ] = compute_nodal_states_mutualistic_pop( x_init, alpha, theta, B, t_obs );
                case 3
                    %Michaelis-Menten
                    
                    %the initial nodal state
                    x_init = parameters.MM.x_init_max*rand( N, 1 );
                    
                    %generate the nodal state sequence
                    [ x, dx ] = compute_nodal_states_michaelis_menten( x_init, parameters.MM.hill_coeff, B, t_obs  );
                case 4
                    %SIS
                    
                    %the initial nodal state
                    x_init = parameters.SIS.x_init_max*rand( N, 1 );
                    
                    %the parameters of the model
                    delta_init = 1 + parameters.SIS.sigma_delta*( 2*rand( N, 1 ) - 1);
                    R_0_init = eigs( diag( 1./sqrt( delta_init ))*B*diag( 1./sqrt( delta_init )), 1 );
                    delta = R_0_init./parameters.SIS.R_0_SIS*delta_init; %then it holds eigs(W, 1) ==parameters.SIS.R_0_SIS, where W = diag(1./results.SIS.delta  )*results.B
                    
                    %generate the nodal state sequence
                    [ x, dx] = compute_nodal_states_SIS( x_init, delta, B, t_obs );
                case 5
                    %Kuramoto
                    
                    %the initial nodal state
                    x_init = ( parameters.kuramoto.x_init_max - parameters.kuramoto.x_init_min )*rand( N, 1 ) - ( parameters.kuramoto.x_init_max - parameters.kuramoto.x_init_min )/2;
                    
                    %the parameters of the model
                    omega = parameters.kuramoto.sigma_omega*randn( N, 1 );
                    
                    %generate the nodal state sequence
                    [ x, dx] = compute_nodal_states_kuramoto( x_init, omega, B, t_obs);
                case 6
                    %Cowan-Wilson
                    
                    %the initial nodal state
                    x_init = parameters.cw.x_init_max*rand( N, 1 );
                    
                    %generate the nodal state sequence
                    [ x, dx] = compute_nodal_states_cowan_wilson( x_init, parameters.cw.tau, parameters.cw.mu, B, t_obs  );
            end
            
            %compute the POD/low-rank approximation of the nodal state
            [ x_apx, norm_cp_past, error_x_apx ] = obtain_low_rank_apx( x );
            
            %save the results
            norm_cp_this_parfor( run_i, : ) = norm_cp_past;
            error_this_parfor( run_i, : ) = error_x_apx;
        end
        
        %save the results
        norm_cp( :, :,  model_i ) =  norm_cp_this_parfor;
        error_POD( :, :,  model_i ) =  error_this_parfor;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%% save file  %%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    save( filename )
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% plot results, corresponding to Fig_S1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%open a new large figure
figure
fig = gcf;
fig.Units = 'normalized';
fig.OuterPosition = [ 0 0 1 1 ];

%loop over all models
for model_i = 1:6
    mean_norm_cp_i = mean( norm_cp(:, :, model_i), 1 );
    
    subplot( 3, 2, model_i )
    semilogy( mean_norm_cp_i, 'b')
    yline(eps, '--', 'color', [.5 .5 .5])
    title( strcat( 'Norm c_p; ', model_name{model_i} ) )
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% plot results, corresponding to Fig_S2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%open a new large figure
figure
fig = gcf;
fig.Units = 'normalized';
fig.OuterPosition = [ 0 0 1 1 ];

%loop over all models
for model_i = 1:6
    mean_error_i = mean( error_POD(:, :, model_i), 1 );
    
    subplot( 3, 2, model_i )
    semilogy( mean_error_i, 'b')
    yline(eps, '--', 'color', [.5 .5 .5])
    title( strcat( 'Error x_{apx}; ', model_name{model_i} ) )
end



end


function [ x_apx, norm_cp, error_x_apx ]= obtain_low_rank_apx( x )
%this function approximates the nodal state x by the POD,
%based on the observations of x

%network size
N = size( x, 1 );

%number of time points
n = size( x, 2 );

%obtain the agitation modes from the singular value decomposition of the
%nodal state sequence x( :, t_ind ) until the observation time
[ singular_vectors, ~, ~] = svds( x, N );

%the low-rank POD approximation x_apx of the nodal state x
x_apx = zeros( N, n );

%the contribution Phi_cp of the p-th agitation mode
norm_cp = nan( 1, N  );

%the error Phi_PODm of the POD with m agitation modes
error_x_apx = nan( 1, N  );

for p = 1:N
    %agitation mode y_p
    y_p = singular_vectors( :, p );
    
    %scalar function c_p
    c_p = transpose( y_p ) * x;
    
    %add the contribution of the p-th agitation mode to the POD x_apx
    x_apx = x_apx + y_p * c_p;
    
    error_x_apx( p ) = norm( x_apx - x )/norm( x );
    
    norm_cp( p ) = norm( c_p )/norm( x );
end

end


