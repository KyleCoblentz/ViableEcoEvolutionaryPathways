clear; clc; clf;

%% NEED THE FOLLOWING FUNCTIONS
% V4_pick_individuals
% V4_medians_and_cis
% GEM_MacArthurRosenzweig
% rep
% MacArthurRosenzweig_bdLogistic

%% first define parameters and run the deterministic ODE and 
%% and the GEM version to get a feel for the process

% parameters for resource logistic growth

b = 4.5;                      params(1) = b; % prey birth rate

d = 1;                         params(2) = d; % prey death rate

qb = 0.00875;                      params(3) = qb; % birth density dependence

qd = 0.00875;                      params(4) = qd; % death density dependence

% parameters for b, d, qb, and qd give a K of 200

% parameters for the interaction and predator

a = 0.03;                       params(5) = a; % space clearance rate

h = 0.01;                        params(6) = h; % handling time

e = 0.3;                        params(7) = e; % conversion efficiency

m = 0.6;                        params(8) = m; % predator per capita mortality

% set up for GEM

cv_vector = [0 0 0 0 0 0 0 0];

h2_vector = [0 0 0 0 0 0 0 0];

y0 = [50 100];

num_replicates = 10;

t_max = 300;

tspan = [0 t_max];

% specify which parameters go with which states

state_parameter_match = zeros(length(y0), length(params)); % start matrix with all zeros

state_parameter_match(1, [1:4]) = 1; % parameters associated with R

state_parameter_match(2,[5:8]) = 1; % parameters associated with the predator

% ode to examine how dynamics play out at different parameter values

ode = @(t,y) MacArthurRosenzweig_bdLogistic(t,y,b,d,qb,qd,a,e,h,m);

[t1,y1] = ode45(ode, tspan, y0);
    figure(1);clf(1);
    hold on;
    plot(t1,y1, '-');
    legend({'R', 'C'});
    shg;

    averagestart = zeros(1,2);
    averagestart  = [round(mean(y1(round(size(y1,1)/2):end,1))), round(mean(y1(round(size(y1,1)/2):end,2)))];
    

% call and run GEM

[x_stand, pop_stand, stand_times, pop_data_out, x_data_out, x_var_data_out] = GEM_MR_bdLogistic(params, state_parameter_match, cv_vector, h2_vector, num_replicates, averagestart, t_max);

%% Now, define a loop to cover the space clearance rate-handling time plane
%% getting the proportion of persitent populations at each parameter value combination

%% Note that this can take a very long time to run i.e. days

h_vec = 0.01:0.02:0.65;

h_vec = rep(h_vec', 15);

a_vec = 0.01:0.01:0.15; 

a_vec = rep(a_vec, 33);

% set up data frame to store data

prop_ext_highint_lowmag = nan(length(a_vec), 5);

% store a and h values in the matrix

prop_ext_highint_lowmag(:,1) = h_vec;

prop_ext_highint_lowmag(:,2) = a_vec;

for i = 1:size(prop_ext_highint_lowmag,1)
    
    i % just print which step the thing is on. need to go back and get rid of all the updates from GEM algorithm
    
   h = prop_ext_highint_lowmag(i,1);      params(6) = h; % get handling time
   
   a = prop_ext_highint_lowmag(i,2);      params(5) = a; % get space clearance rate
   
   % set other parameter values and things needed for the GEM function
   
   % parameters for resource logistic growth

   b = 4.5;                      params(1) = b; % prey birth rate

   d = 1;                      params(2) = d; % prey death rate
   
   qb = 0.00875;                params(3) = qb; % density dependence in birth rate
   
   qd = 0.00875;                params(4) = qd; % density dependence in death rate
   
   
   % parameters for the interaction and predator
   
    e = 0.2;                        params(7) = e; % conversion efficiency

    m = 0.3;                        params(8) = m; % predator per capita mortality

    % set up for GEM

    cv_vector = [0 0 0 0 0 0 0 0];

    h2_vector = [0 0 0 0 0 0 0 0];

    num_replicates = 100;

    t_max = 300;
    
    t_max_ode = 2000;

    tspan = [0 t_max];
    
    tspan_ode = [0 t_max_ode];
    
    y0_ode = [100 50];
    
    % run ode use average densities after transient dynamics as starting values
    
    ode = @(t,y) MacArthurRosenzweig_bdLogistic(t,y,b,d,qb,qd,a,e,h,m);

    [t1,y1] = ode45(ode, tspan_ode, y0_ode);
    
    y0 = zeros(1,2);
    
    y0(1,1:2) = [round(mean(y1(round(size(y1,1)/2):end,1))), round(mean(y1(round(size(y1,1)/2):end,2)))];
    
    % store minimum deterministic resource and consumer densities after
    % transient dynamics
    
    prop_ext_highint_lowmag(i,4) = min(y1(round(size(y1,1)/2):end, 1));
    
    prop_ext_highint_lowmag(i,5) = min(y1(round(size(y1,1)/2):end, 2));
    
    
    % specify which parameters go with which states

    state_parameter_match = zeros(length(y0), length(params)); % start matrix with all zeros

    state_parameter_match(1, [1:4]) = 1; % parameters associated with R

    state_parameter_match(2,[5:8]) = 1; % parameters associated with the predator
    
    
    % call and run GEM

    
    [x_stand, pop_stand, stand_times, pop_data_out, x_data_out, x_var_data_out] = GEM_MR_bdLogistic(params, state_parameter_match, cv_vector, h2_vector, num_replicates, y0, t_max);


    % from GEM results calculate the proportion of populations that went
    % extinct and save into matrix
    
    prop_ext_highint_lowmag(i,3) = num_replicates - sum(isnan(pop_stand(2,end,:)));
    
end

% plot of heat map for probability to not go extinct ... 

table_prop_ext_highint = array2table(prop_ext_highint_lowmag, 'VariableNames', {'HandlingTime', 'AttackRate', 'NumberPersist', 'MinPrey', 'MinPred'});

h = heatmap(table_prop_ext_highint, 'HandlingTime', 'AttackRate', 'ColorVariable', 'NumberPersist');

h.YDisplayData = flipud(h.YDisplayData);

h.Colormap = parula;

