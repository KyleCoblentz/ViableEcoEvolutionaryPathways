%%% Rosenzweig-MacArthur BD-Logistic Evolution on Landscape GEM

%%% needs the following 
    % GEM_MacArthurRosenzweigh
    % V4_pick_individuals
    % V4_medians_and_cis
    % MacArthurRosenzweig_bdLogistic

colors(1:6,1:3) = [[1 0.47 0]; [0.4 0 1]; [0.07 1 0]; [0.9 0 1]; [0.96 1 0]; [1 0 0.2]];
fill_colors = colors.*0.8;
    
% parameters for resource logistic growth

b = 4.5;                      params(1) = b; % prey birth rate

d = 1;                         params(2) = d; % prey death rate

qb = 0.000875;                      params(3) = qb; % birth density dependence

qd = 0.000875;                      params(4) = qd; % death density dependence

% parameters for b, d, qb, and qd give a K of 2000

% parameters for the interaction and predator

xa = 24.8945;                              params(5) = xa;

xh = 11.3636;                              params(6) = xh;

amax = 0.05;                            params(7) = amax;

hmin = 0.005;                            params(8) = hmin;

theta_a = 8;                         params(9) = theta_a;

theta_h = 8;                         params(10) = theta_h;

tau_a = 10;                           params(11) = tau_a;

tau_h = 10;                           params(12) = tau_h;

e = 0.1;                        params(13) = e; % conversion efficiency

m = 0.6;                        params(14) = m; % predator per capita mortality

% set up for GEM

cv_vector = [0 0 0 0 0.12 0.125 0 0 0 0 0 0 0 0];

h2_vector = [0 0 0 0 0.2 0.2 0 0 0 0 0 0 0 0];

y0 = [50 100];

num_replicates = 100;

t_max = 50;

tspan = [0 t_max];

% specify which parameters go with which states

state_parameter_match = zeros(length(y0), length(params)); % start matrix with all zeros

state_parameter_match(1, [1:4]) = 1; % parameters associated with R

state_parameter_match(2,[5:14]) = 1; % parameters associated with the predator

% ode to examine how dynamics play out at different parameter values

ode = @(t,y) MacArthurRosenzweig_bdLogistic_ExplicitTraits(t,y,b,d,qb,qd,xa, amax, theta_a, tau_a, e, xh, hmin, theta_h, tau_h,m);

[t1,y1] = ode45(ode, tspan, y0);
    figure(1);clf(1);
    hold on;
    plot(t1,y1, '-');
    legend({'R', 'C'});
    shg;

    averagestart = zeros(1,2);
    averagestart(1,1:2) = [round(mean(y1(90:109,1))), round(mean(y1(90:109,2)))];

% call and run GEM

[x_stand, pop_stand, stand_times, pop_data_out, x_data_out, x_var_data_out] = GEM_MR_bdLogistic_ExplicitTraits_mod(params, state_parameter_match, cv_vector, h2_vector, num_replicates, averagestart, t_max);

% print the number of extant populations at the end of the sampling

num_replicates - sum(isnan(pop_stand(2,end,:)))

% mean population densities

subplot(1,3,1);
hold on;
    jbfill(stand_times,pop_data_out(1,:,2),pop_data_out(3,:,2),fill_colors(2,:),'w',1,0.2); hold on;
    plot(stand_times,pop_data_out(2,:,2),'-','Color',colors(2,:),'LineWidth',2); hold on;        
    jbfill(stand_times, pop_data_out(1,:,1), pop_data_out(3,:,1), fill_colors(1,:), 'w',1,0.2); hold on;
    plot(stand_times,pop_data_out(2,:,1),'-','Color',colors(1,:), 'LineWidth',2); hold on;
    xlabel('Time')
    ylabel('Density')
   
% mean space clearance rate changes
  
subplot(1,3,2);
    hold on;
    jbfill(stand_times,amax*exp(-(x_data_out(1,:,5) - theta_a).^2/(2*tau_a^2)),amax*exp(-(x_data_out(3,:,5)- theta_a).^2/(2*tau_a^2)),fill_colors(2,:),'w',1,0.2); hold on;
    plot(stand_times,amax*exp(-(x_data_out(2,:,5)-theta_a).^2/(2*tau_a^2)),'-','Color',colors(2,:),'LineWidth',2); hold on;
    xlabel('Time')
    ylabel('Mean Space Clearance Rate')

% mean handling time changes

subplot(1,3,3);
    hold on;
    jbfill(stand_times,1 + hmin - exp(-(x_data_out(1,:,6) - theta_h).^2/(2*tau_h^2)),1 + hmin - exp(-(x_data_out(3,:,6)-theta_h).^2/(2*tau_h^2)),fill_colors(2,:),'w',1,0.2); hold on;
    plot(stand_times,1+hmin-exp(-(x_data_out(2,:,6)-theta_h).^2/(2*tau_h^2)),'-','Color',colors(2,:),'LineWidth',2); hold on;
    xlabel('Time')
    ylabel('Mean Handling Time')


% want to export data

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

indpops = reshape(pop_stand, 2, 20100).';  

% individual run space clearance rates and handling times

indtraits = reshape(x_stand(5:6, :, :), 2, 20100).';

% population id

populationnumber = rep(1:100, 201);

% whether each system goes extinct or not

endpop = squeeze(pop_stand(1, end,:));

endpop = isnan(endpop);

endpop = repelem(endpop, 201);

% times for the individual runs

indtimes = repmat(stand_times, 1, 100).';

% put all of the individual run data together and save to a .csv file

indrundata = horzcat(indpops, indtraits, populationnumber, endpop, indtimes);

% convert to dataset 

indrundata = mat2dataset(indrundata, 'VarNames', {'ResDens', 'PredDens', 'SCRate', 'HandlingTime', 'Pop', 'Extinct', 'Time'});

% save as .csv file

export(indrundata, 'File', 'indrundata.csv', 'Delimiter', ',');





