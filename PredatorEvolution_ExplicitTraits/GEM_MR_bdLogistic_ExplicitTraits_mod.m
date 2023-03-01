function [x_stand, pop_stand, stand_times, pop_data_out, x_data_out, x_var_data_out] = GEM_MR_bdLogistic(params, state_parameter_match, cv_vector, h2_vector, num_replicates, y0, t_max)
% This function is a GEM for the MacArthur-Rosenzweig predator prey model
%   with b,d-Logistic growth in the prey

%% standardized time steps for storing time series
% import the time span (t_max) but decide on step lengths
stand_times = 0:0.25:t_max;
num_time_steps = length(stand_times); % calculate the number of standardized time steps

%% preallocate space in which to log data at standardized times
% need one for each population, one for each evolving trait, and one for
% the variance in each evolving trait
no_species = length(y0); % the number of populations in the model = number of starting sizes
no_params = length(params);
no_columns = length(params)+1; % user must update this

pop_stand = zeros(no_species,num_time_steps,num_replicates); % population size
x_stand = nan(no_params,num_time_steps,num_replicates); % trait
x_var_stand = nan(no_params,num_time_steps,num_replicates); % trait variance

%% run for loop for each replicate simulation
parfor i = 1:num_replicates % for parallel computing
%for i = 1:num_replicates % not parallel
    t = 0; % initial time
    rng('shuffle'); % change random number seed    
    i % display replicate in the Command Window
    % pre-allocate for the sliced standardized variables
        blank_comm_matrix = nan(sum(y0),no_columns); % set up a NaN matrix for all individuals and traits
        pop_slice = zeros(no_species,num_time_steps);
        x_slice = nan(length(params),num_time_steps);
        x_var_slice = nan(length(params),num_time_steps);  
        
    % assign initial population sizes to sliced abundance variables    
        pop_slice(:,1) = y0';

    % pull initial distribution of trait (repeat for all populations)
        end_row = cumsum(y0);
        starting_row = [1 1+end_row(1:length(end_row)-1)];
    % determine how many of each state to pull
        ind_to_assign = state_parameter_match.*y0';
        
    for qq = 1:length(y0) % loop through states
        blank_comm_matrix(starting_row(qq):end_row(qq),1) = qq;
        for zz = 1:length(params) % loop through parameters   
            temp  = V4_pick_individuals(params(zz),cv_vector(zz)*params(zz),ind_to_assign(qq,zz));
            if isempty(temp) == 0
                blank_comm_matrix(starting_row(qq):end_row(qq),1+zz) = temp;               
            end
        end
    end
        
    x_dist_init = blank_comm_matrix;
    x_slice(1:no_columns-1,1) = nanmean(x_dist_init(:,2:no_columns),1)'; % initial mean trait vector
    x_var_slice(:,1) = var(x_dist_init(:,2:no_columns),1, 'omitnan')'; % initial variances in trait

    %% Initiate core GEM algorithm
    time_step_index = 2; % start a counter for standard times
    time_step = stand_times(time_step_index); % assign first standard time step

    % assign initial pop sizes to reduction variables (updated abundances)
        R = y0(1);
        C = y0(2);        
        x_dist = x_dist_init;
        
    while t < t_max
        % loop through species to find individuals for each state
        params_next = nan(no_species,no_params);
        whosnext = nan(size(y0));
        for zz = 1:no_species
            inds_in_state = find(x_dist(:,1)==zz);
            which_params = find(state_parameter_match(zz,:));
            if isempty(inds_in_state) == 0
                which_row = randi(length(inds_in_state));
                whosnext(zz) = inds_in_state(which_row);
                params_next(zz,which_params) = x_dist(whosnext(zz),1+which_params);
            else
                params_next(zz,which_params) = 0; % if pop is gone, set parameters to 0
            end
        end
       
        % pull out and re-assign parameters with names 
        
        % state 1 parameters (R1) [1:4]
        
        b = params_next(1,1);
        
        d = params_next(1,2);
        
        qb = params_next(1,3);
        
        qd = params_next(1,4);
        
        % state 2 parameters (C) [5:14]
        
        xa = params_next(2,5);
        
        xh = params_next(2,6);
        
        amax = params_next(2,7);
        
        hmin = params_next(2,8);
        
        theta_a = params_next(2,9);
        
        theta_h = params_next(2,10);
        
        tau_a = params_next(2,11);
        
        tau_h = params_next(2,12);
        
        e = params_next(2,13);
        
        m = params_next(2,14);
        
        % get a and h from the parameters that define them
        
        a = amax*exp(-((xa - theta_a)^2)/(2*tau_a^2));
        
        h = (1 + hmin) - exp(-((xh - theta_h)^2)/(2*tau_h^2));
        
        % set up rates of each possible event
        
            % 1: birth of prey
            
            b_R = R*max((b - qb*R), 0);
            
            % 2: natural death of prey
            
            d_Rd = R*(d + qd*R)
            
            % 3: prey dies from predation
            
            d_RC = a*R*C/(1 + a*h*R);
            
            % 4: predator is born
            
            b_C = e * d_RC;
            
            % 5: predator dies
            
            d_C =  m * C ;
            
        % sum the events to create the wheel of fortune
            
            CS_vector = cumsum([b_R d_Rd d_RC b_C d_C]);
            Slice_widths = CS_vector./CS_vector(end);
            LI = rand < Slice_widths;
            Event_index = find(LI, 1, 'first');
    
        % event outcomes 
        
        if Event_index == 1 % birth of Resource
            if h2_vector(1) == 0
                x_parent = (1-h2_vector).*nanmean(x_dist_init(find(x_dist_init(:,1)==1),2:no_columns)) + h2_vector.*x_dist(whosnext(1),2:no_columns);
            elseif h2_vector(1) > 0
                x_parent = (1-h2_vector).*nanmean(x_dist(find(x_dist(:,1)==1),2:no_columns)) + h2_vector.*x_dist(whosnext(1),2:no_columns);
            end
            off_std = sqrt(1-h2_vector.^2).*((1-h2_vector).*nanstd(x_dist_init(find(x_dist_init(:,1)==1),2:no_columns))+h2_vector.*nanstd(x_dist(find(x_dist(:,1)==1),2:no_columns))); % offspring trait distribution std
            x_dist(size(x_dist,1)+1,2:no_columns) = V4_pick_individuals(x_parent,off_std,1); % return trait
            x_dist(end,1) = 1;
            
        elseif Event_index == 2 % death of resource 
            x_dist(whosnext(1),:) = []; % reduce dist by lost individual
            
        elseif Event_index == 3 % death of resource by being eaten by the predator
            x_dist(whosnext(1),:) = [] % reduce dist by lost individual
            
        elseif Event_index == 4 % birth of predator
            if h2_vector(5) == 0
                 x_parent = (1-h2_vector).*nanmean(x_dist_init(find(x_dist_init(:,1)==2),2:no_columns)) + h2_vector.*x_dist(whosnext(2),2:no_columns);
            elseif h2_vector(5) > 0 
                 x_parent = (1-h2_vector).*nanmean(x_dist(find(x_dist(:,1)==2),2:no_columns)) + h2_vector.*x_dist(whosnext(2),2:no_columns);
            end
            off_std = sqrt(1-h2_vector.^2).*((1-h2_vector).*nanstd(x_dist_init(find(x_dist_init(:,1)==2),2:no_columns))+h2_vector.*nanstd(x_dist(find(x_dist(:,1)==2),2:no_columns))); % offspring trait distribution std
            x_dist(size(x_dist,1)+1,2:no_columns) = V4_pick_individuals(x_parent,off_std,1); % return trait
            x_dist(end,1) = 2;
            
        elseif Event_index == 5 % death of predator
            x_dist(whosnext(2),:) = []; % reduce dist by lost individual

        end
        
        R = sum(x_dist(:,1)==1); % count rows with State ID == 1
        C = sum(x_dist(:,1)==2); % count rows with State ID == 2
        
        if R == 0
            
            break
            
        end
        
        if C == 0
            
            break
            
        end
        
        
        if t > time_step
            time_step
            [R C]
            pop_slice(1:no_species,time_step_index) = [R C]; % assign current values to sliced standard times
            x_slice(1:no_params,time_step_index) = nanmean(x_dist(:,2:no_columns),1);
            x_var_slice(1:no_params,time_step_index) = var(x_dist(:,2:no_columns),1, 'omitnan');
            
            time_step_index = time_step_index + 1; % advance to next standardized time
            time_step = stand_times(time_step_index);
        end
            
        time_advance = exp(-1/CS_vector(end))/(CS_vector(end));
        if isnan(time_advance) == 0
            t = t + time_advance;
        else
            break
        end
         
    end  
    
    
    
    % pass in current values to the final standardized time
        pop_slice(1:no_species,time_step_index) = [R C];
        x_slice(1:no_params,time_step_index) = nanmean(x_dist(:,2:no_columns),1);
        x_var_slice(1:no_params,time_step_index) = var(x_dist(:,2:no_columns),1, 'omitnan');

    % pass sliced variable to standard time matrix
        pop_stand(:,:,i) = pop_slice;
        x_stand(:,:,i) = x_slice;
        x_var_stand(:,:,i) = x_var_slice;

        
end

    %% calculate ci's for time series

pop_stand2 = pop_stand; pop_stand(pop_stand ==0) = NaN;
x_stand2 = x_stand; x_stand(x_stand ==0) = NaN;
x_var_stand2 = x_var_stand; x_var_stand(x_var_stand ==0) = NaN;

    upper_ci_level = 75; % choose ci levels
    lower_ci_level = 25; % choose ci levels    
    pop_data_out = V4_medians_and_cis(upper_ci_level,lower_ci_level,pop_stand); % abundance      
    x_data_out = V4_medians_and_cis(upper_ci_level,lower_ci_level,x_stand2); % trait
    x_var_data_out = V4_medians_and_cis(upper_ci_level,lower_ci_level,x_var_stand2); % variance in trait 
    

