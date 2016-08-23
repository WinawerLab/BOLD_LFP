function NS = neural_sim_defaults()
% Define default parameters for neural simulation, and return these in the
% struct NS.params
%
% NS = neural_sim_defaults()

% see ns_set for description of parameters
NS = [];
NS = ns_set(NS, 'save_inputs', 0);
NS = ns_set(NS, 'trials_save_inputs', []);
NS = ns_set(NS, 'num_neurons', 600);      
NS = ns_set(NS, 'num_conditions', 8);       
NS = ns_set(NS, 'num_averages', 30);       
NS = ns_set(NS, 'num_experiments', 1);
NS = ns_set(NS, 'dt', 0.001); 
NS = ns_set(NS, 'trial_length', 1);
NS = ns_set(NS, 'gamma_range', [50 60]); 
NS = ns_set(NS, 'alpha_range', [9 12]); 
NS = ns_set(NS, 'poisson_baseline', 1);

% NS = ns_set(NS, 'poisson_bb', [0 1]);
% NS = ns_set(NS, 'coherence_bb', [0]);
% 
% NS = ns_set(NS, 'poisson_g', [0 1]);
% NS = ns_set(NS, 'coherence_g', [1]);
% 
% NS = ns_set(NS, 'poisson_a', [0 1]);
% NS = ns_set(NS, 'coherence_a', [0]);
