function Fig_1( seed )

%set the random number generator for reproducibility
if nargin < 1
    seed = 1;
end

rng( seed );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% set the simulation parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

parameters = set_parameters;

%set the ratio t_obs/T_max:
T_obs_fraction = 1/5;

%the number of observations n
num_observations = 1e2;

%the indices of the sequence 1:num_observations that are past the
%observation time
t_future_ind = ceil( num_observations*T_obs_fraction )+1:num_observations;

%the number of agitation modes m
POD_m_max = 15;

%only for plotting: num_N_plot is the number of nodal traces that are plotted
num_N_plot = 6;

%only for plotting: num_k_plot is the time points that are plotted
num_k_plot = 50;

%save the results of the different models here:
results_LV = [];
results_SIS = [];
results_MP = [];
results_MM = [];
results_kuramoto = [];
results_cw = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Model 1: Lotka-Volterra
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%the zero-one (unweighted) network
A_LV = create_network( [], 'little_rock_lake',[], [],[], [], [], [] );

%the network size/number of nodes
N_LV = size( A_LV, 1 );

%the link weights
link_weights = parameters.min_link_weight + ( parameters.max_link_weight - parameters.min_link_weight )*rand( nnz( A_LV ), 1);

%the unweighted network
results_LV.B = zeros( N_LV );
results_LV.B( A_LV > 0 ) = link_weights;

%the parameters of the model
results_LV.alpha = 1 + parameters.LV.sigma_alpha*( 2*rand( N_LV, 1 ) - 1);
results_LV.theta = 1 + parameters.LV.sigma_theta*( 2*rand( N_LV, 1 ) - 1);

%the initial nodal state
results_LV.x_init = parameters.LV.x_init_max*rand( N_LV, 1 );

%the observation times samples
t_obs_LV = linspace( 0, parameters.T_max( 1 ), num_observations );

%generate the nodal state sequence
[ x_LV, ~ ] = compute_nodal_states_lotka_volterra( results_LV.x_init, results_LV.alpha, results_LV.theta, results_LV.B, t_obs_LV );

%obtain the POD of the nodal state from observation until the
%relative time T_obs_fraction
x_apx_LV = obtain_low_rank_apx( x_LV, T_obs_fraction, POD_m_max );

%the approximation error of the POD for future times (after T_obs_fraction)
error_LV =  mean( mean( abs( x_LV( :, t_future_ind ) - x_apx_LV( :, t_future_ind ) )));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Model 2: Mutualistic population
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%the zero-one (unweighted) network
A_MP = create_network( [], 'mutualistic',[], [],[], [], [], [] );

%the network size/number of nodes
N_MP = size( A_MP, 1 );

%the link weights
link_weights = parameters.min_link_weight + ( parameters.max_link_weight - parameters.min_link_weight )*rand( nnz( A_MP ), 1);

%the unweighted network
results_MP.B = zeros( N_MP );
results_MP.B ( A_MP > 0 ) = link_weights;

%the parameters of the model
results_MP.alpha = 1 + parameters.MP.sigma_alpha*( 2*rand( N_MP, 1 ) - 1 );
results_MP.theta = 1 + parameters.MP.sigma_theta*( 2*rand( N_MP, 1 ) - 1 );

%the initial nodal state
results_MP.x_init = parameters.MP.x_init_max*rand( N_MP, 1 );

%the observation times samples
t_obs_MP = linspace( 0, parameters.T_max( 2 ), num_observations );

%generate the nodal state sequence
[ x_MP, ~ ] = compute_nodal_states_mutualistic_pop( results_MP.x_init, results_MP.alpha, results_MP.theta, results_MP.B, t_obs_MP );

%obtain the POD of the nodal state from observation until the
%relative time T_obs_fraction
x_apx_MP = obtain_low_rank_apx( x_MP, T_obs_fraction, POD_m_max );

%the approximation error of the POD for future times (after T_obs_fraction)
error_MP =  mean( mean( abs( x_MP( :, t_future_ind ) - x_apx_MP( :, t_future_ind ) )));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Model 3: Michaelis-Menten
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%the zero-one (unweighted) network
A_MM = create_network( [], 'yeast_regulatory',[], [],[], [], [], [] );

%the network size/number of nodes
N_MM = size( A_MM, 1 );

%the link weights
link_weights = parameters.min_link_weight + ( parameters.max_link_weight - parameters.min_link_weight )*rand( nnz( A_MM ), 1);

%the unweighted network
results_MM.B = zeros( N_MM );
results_MM.B( A_MM > 0 ) = link_weights;

%the initial nodal state
results_MM.x_init = parameters.MM.x_init_max*rand( N_MM, 1 );

%the observation times samples
t_obs_MM = linspace( 0, parameters.T_max( 3 ), num_observations );

%generate the nodal state sequence
[ x_MM, ~ ] = compute_nodal_states_michaelis_menten( results_MM.x_init, parameters.MM.hill_coeff, results_MM.B, t_obs_MM  );

%obtain the POD of the nodal state from observation until the
%relative time T_obs_fraction
x_apx_MM = obtain_low_rank_apx( x_MM, T_obs_fraction, POD_m_max );

%the approximation error of the POD for future times (after T_obs_fraction)
error_MM =  mean( mean( abs( x_MM( :, t_future_ind ) - x_apx_MM( :, t_future_ind ) )));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Model 4: SIS epidemics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%the zero-one (unweighted) network
A_SIS = create_network( [], 'infectious',[], [],[], [], [], [] );

%the network size/number of nodes
N_SIS = size( A_SIS, 1 );

%the link weights
link_weights = parameters.min_link_weight + ( parameters.max_link_weight - parameters.min_link_weight )*rand( nnz( A_SIS ), 1);

%the unweighted network
results_SIS.B = zeros( N_SIS );
results_SIS.B( A_SIS > 0 ) = link_weights;

%the parameters of the model
delta_init = 1 + parameters.SIS.sigma_delta*( 2*rand( N_SIS, 1 ) - 1);
R_0_init = eigs( diag( 1./sqrt( delta_init ))*results_SIS.B*diag( 1./sqrt( delta_init )), 1 );
multiplicity_delta = R_0_init./parameters.SIS.R_0_SIS;
results_SIS.delta = multiplicity_delta*delta_init; %then it holds eigs(W, 1) ==parameters.SIS.R_0_SIS, where W = diag(1./results.SIS.delta  )*results.B

%the initial nodal state
results_SIS.x_init = parameters.SIS.x_init_max*rand( N_SIS, 1 );

%the observation times samples
t_obs_SIS = linspace( 0, parameters.T_max( 4 ), num_observations );

%generate the nodal state sequence
[ x_SIS, ~ ] = compute_nodal_states_SIS( results_SIS.x_init, results_SIS.delta, results_SIS.B, t_obs_SIS );

%obtain the POD of the nodal state from observation until the
%relative time T_obs_fraction
x_apx_SIS = obtain_low_rank_apx( x_SIS, T_obs_fraction, POD_m_max );

%the approximation error of the POD for future times (after T_obs_fraction)
error_SIS =  mean( mean( abs( x_SIS( :, t_future_ind ) - x_apx_SIS( :, t_future_ind ) )));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Model 5: Kuramoto oscillators
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%the zero-one (unweighted) network
A_kuramoto = create_network( [], 'human_connectome',[], [],[], [], [], [] );

%the network size/number of nodes
N_kuramoto = size( A_kuramoto, 1 );

%the link weights
link_weights = parameters.min_link_weight + ( parameters.max_link_weight - parameters.min_link_weight )*rand( nnz( A_kuramoto ), 1);

%the unweighted network
results_kuramoto.B = zeros( N_kuramoto );
results_kuramoto.B( A_kuramoto > 0 ) = link_weights;

%the parameters of the model
results_kuramoto.omega = parameters.kuramoto.sigma_omega*randn( N_kuramoto, 1 );

%the initial nodal state
results_kuramoto.x_init = ( parameters.kuramoto.x_init_max - parameters.kuramoto.x_init_min )*rand( N_kuramoto, 1 ) - ( parameters.kuramoto.x_init_max - parameters.kuramoto.x_init_min )/2;

%the observation times samples
t_obs_kuramoto = linspace( 0, parameters.T_max( 5 ), num_observations );

%generate the nodal state sequence
[ x_kuramoto, ~ ] = compute_nodal_states_kuramoto( results_kuramoto.x_init, results_kuramoto.omega, results_kuramoto.B, t_obs_kuramoto );

%obtain the POD of the nodal state from observation until the
%relative time T_obs_fraction
x_apx_kuramoto = obtain_low_rank_apx( x_kuramoto, T_obs_fraction, POD_m_max );

%the approximation error of the POD for future times (after T_obs_fraction)
error_kuramoto =  mean( mean( abs( x_kuramoto( :, t_future_ind ) - x_apx_kuramoto( :, t_future_ind ) )));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Model 6: Cowan-Wilson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%the zero-one (unweighted) network
A_cw = create_network( [], 'c_elegans_neuronal',[], [],[], [], [], [] );

%the network size/number of nodes
N_cw = size( A_cw, 1 );

%the link weights
link_weights = parameters.min_link_weight + ( parameters.max_link_weight - parameters.min_link_weight )*rand( nnz( A_cw ), 1);

%the unweighted network
results_cw.B = zeros( N_cw );
results_cw.B( A_cw > 0 ) = link_weights;

%the initial nodal state
results_cw.x_init = parameters.cw.x_init_max*rand( N_cw, 1 );

%the observation times samples
t_obs_cw = linspace( 0, parameters.T_max( 6 ), num_observations );

%generate the nodal state sequence
[ x_cw, ~ ] = compute_nodal_states_cowan_wilson( results_cw.x_init, parameters.cw.tau, parameters.cw.mu, results_cw.B, t_obs_cw  );

%obtain the POD of the nodal state from observation until the
%relative time T_obs_fraction
x_apx_cw = obtain_low_rank_apx( x_cw, T_obs_fraction, POD_m_max );

%the approximation error of the POD for future times (after T_obs_fraction)
error_cw =  mean( mean( abs( x_cw( :, t_future_ind ) - x_apx_cw( :, t_future_ind ) )));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% plot the results 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%the time steps to be plotted of the POD x_apx (choosing all time
%steps 1:num_observations would make the plot too busy)
k_plot = ceil( linspace( 1, num_observations, num_k_plot ));

%open a new large figure 
figure
fig = gcf;
fig.Units = 'normalized';
fig.OuterPosition = [ 0 0 1 1 ];

%loop over all models and plot the respective accuracy of the POD
for model_count = 1:6
    
    switch model_count
        case 1
            x = x_LV;
            x_apx = x_apx_LV;
            t_obs = t_obs_LV;
            error = error_LV;
            model_name = 'LV';
        case 2
            x = x_MP;
            x_apx = x_apx_MP;
            t_obs = t_obs_MP;
            error = error_MP;
            model_name = 'MP';
        case 3
            x = x_MM;
            x_apx = x_apx_MM;
            t_obs = t_obs_MM;
            error = error_MM;
            model_name = 'MM';
        case 4
            x = x_SIS;
            x_apx = x_apx_SIS;
            t_obs = t_obs_SIS;
            error = error_SIS;
            model_name = 'SIS';
        case 5
            x = x_kuramoto;
            x_apx = x_apx_kuramoto;
            t_obs = t_obs_kuramoto;
            error = error_kuramoto;
            model_name = 'KUR';
        case 6
            x = x_cw;
            x_apx = x_apx_cw;
            t_obs = t_obs_cw;
            error = error_cw;
            model_name = 'WC';
    end
    
    % obtain N_plot, the nodes whose state is plotted 
    N_plot = nan( 1, num_N_plot+1 );
    
    %x_end_quantised has N_plot+1 values, ranging from the minimum to the 
    %maximum nodal state at the last time point 
    x_end_quantised = linspace( min( x( :, end ) ), max( x( :, end ) ), num_N_plot+1 );
    
    %the nodes in N_plot are obtain by finding the nodes in x( :, end ) 
    %that are closes to the values in x_end_quantised
    for curr_N_plot = 1:(num_N_plot+1)
        [~, indd] = min( abs( x( :, end ) - x_end_quantised( curr_N_plot ) ));
        N_plot( curr_N_plot ) = indd;
    end
    
    %remove the node with the smallest nodal state (which is almost 
    %zero and has hardly visible dynamics on the plot)
    N_plot( 1 ) = [];
    
    %plot the results of the current model
    subplot( 3, 2, model_count )    
    plot( t_obs, x( N_plot, : )', 'b')
    hold
    plot( t_obs( k_plot ), x_apx( N_plot, k_plot )', 'or')
    xline( t_obs( t_future_ind( 1 ) ) )
    title( [ model_name, '; error = ' num2str( error )] )
end


end


function [ x_apx ]= obtain_low_rank_apx( x, T_obs, POD_m_max )
%this function approximates the nodal state x by the POD with POD_m_max
%agitation modes, based on the observations of x until T_obs

%network size
N = size( x, 1 );

%number of time points
n = size( x, 2 );

%the time points that are observed
t_ind = 1:ceil( n*T_obs );

%obtain the agitation modes from the singular value decomposition of the
%nodal state sequence x( :, t_ind ) until the observation time
[ singular_vectors, ~, ~] = svds( x( :, t_ind ), N );

%the low-rank POD approximation x_apx of the nodal state x
x_apx = zeros( N, n );

%loop over all agitation modes
for p = 1:POD_m_max
    
    %agitation mode y_p
    y_p = singular_vectors( :, p );
    
    %scalar function c_p
    c_p = transpose( y_p ) * x;
    
    %add the contribution of the p-th agitation mode to the POD x_apx
    x_apx = x_apx + y_p * c_p;
end

end





