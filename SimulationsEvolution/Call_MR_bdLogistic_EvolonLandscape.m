%%% Rosenzweig-MacArthur BD-Logistic Evolution on Landacape GEM

%%% needs the following 
    % GEM_MacArthurRosenzweig
    % V4_pick_individuals
    % V4_medians_and_cis
    % MacArthurRosenzweig_bdLogistic
   
% parameters for resource logistic growth

b = 4.5;                      params(1) = b; % prey birth rate

d = 1;                         params(2) = d; % prey death rate

qb = 0.00875;                      params(3) = qb; % birth density dependence

qd = 0.00875;                      params(4) = qd; % death density dependence

% parameters for b, d, qb, and qd give a K of 200

% parameters for the interaction and predator

a = 0.04;                       params(5) = a; % space clearance rate

h = 0.25;                        params(6) = h; % handling time

e = 0.3;                        params(7) = e; % conversion efficiency

m = 0.6;                        params(8) = m; % predator per capita mortality

% set up for GEM

cv_vector = [0 0 0 0 0.75 1.2 0 0];

h2_vector = [0 0 0 0 0.35 0.85 0 0];

y0 = [50 100];

num_replicates = 100;

t_max = 200;

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
    averagestart(1,1:2) = [round(mean(y1(120:145,1))), round(mean(y1(120:145,2)))];


% call and run GEM

[x_stand, pop_stand, stand_times, pop_data_out, x_data_out, x_var_data_out] = GEM_MR_bdLogistic(params, state_parameter_match, cv_vector, h2_vector, num_replicates, averagestart, t_max);

num_replicates - sum(isnan(pop_stand(2,end,:)))


% export data to .csv files 

% export median population size data

popmediandata = horzcat(pop_data_out(:,:,1).', pop_data_out(:,:,2).', stand_times.');

popmediandata = mat2dataset(popmediandata, 'VarNames', {'LowerPrey', 'MedPrey', 'UpperPrey','LowerPred', 'MedPred', 'UpperPred', 'Time'});

export(popmediandata, 'File','popmediandata.csv','Delimiter',',')

% export trait median and variance data

mediantraitdata = horzcat(x_data_out(:,:,5).', x_var_data_out(:,:,5).', x_data_out(:,:,6).', x_var_data_out(:,:,6).', stand_times.');

mediantraitdata = mat2dataset(mediantraitdata, 'VarNames', {'lowera', 'meda', 'uppera', 'loweravar', 'medavar', 'upperavar', 'lowerh', 'medh', 'upperh', 'lowerhvar', 'medhvar', 'upperhvar', 'time'});

export(mediantraitdata, 'File', 'mediantraitdata.csv', 'Delimiter', ',');

% individual run data

% individual run population sizes

indpops = reshape(pop_stand, 2, 80100).';  

% individual run space clearance rates and handling times

indtraits = reshape(x_stand(5:6, :, :), 2, 80100).';

% population id

populationnumber = rep(1:100, 801);

% whether each system goes extinct or not

endpop = squeeze(pop_stand(1, end,:));

endpop = isnan(endpop);

endpop = repelem(endpop, 801);

% times for the individual runs

indtimes = repmat(stand_times, 1, 100).';

% put all of the individual run data together and save to a .csv file

indrundata = horzcat(indpops, indtraits, populationnumber, endpop, indtimes);

% convert to dataset 

indrundata = mat2dataset(indrundata, 'VarNames', {'ResDens', 'PredDens', 'SCRate', 'HandlingTime', 'Pop', 'Extinct', 'Time'});

% save as .csv file

export(indrundata, 'File', 'indrundata.csv', 'Delimiter', ',');


