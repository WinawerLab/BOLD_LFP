function NS = neural_sim_defaults()
% Define default parameters for neural simulation, and return these in the
% struct NS.params
%
% NS = neural_sim_defaults()

% see ns_set for description of parameters
NS = ns_set([], 'simulate_method', 'FAST'); 
NS = ns_set(NS, 'num_neurons', 600);      
NS = ns_set(NS, 'bb_fraction', .5);  
NS = ns_set(NS, 'num_conditions', 8);       
NS = ns_set(NS, 'num_averages', 30);       
NS = ns_set(NS, 'num_experiments', 1);
NS = ns_set(NS, 'dt', 0.001); 
NS = ns_set(NS, 'trial_length', 1);
NS = ns_set(NS, 'gamma_range', [40 60]); 
NS = ns_set(NS, 'alpha_range', [8 15]); 
NS = ns_set(NS, 'poisson_baseline', 1);
NS = ns_set(NS, 'poisson_bb_rg', [0 1]);
NS = ns_set(NS, 'poisson_g_rg', [0 1]);
NS = ns_set(NS, 'poisson_a_rg', [0 1]);
NS = ns_set(NS, 'gamma_coh', 1);
NS = ns_set(NS, 'alpha_coh', 1);