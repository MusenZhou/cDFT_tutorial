% This program is used to illustrate the very basic idea of classic density
% functional theory (cDFT) in the ideal limit.


% Description: rho here is density.
clear;
clc;


%% constant define
kb = 1.38064852e-23;                                                        % unit: J/Klevin
T = 300;                                                                    % unit: Klevin
beta = 1./(kb.*T);                                                          % unit: 1./J
plank_constant = 6.626070040e-34;                                           % unit: J*m
NA = 6.022140857e23;                                                        % unit: mol-1
mass_neon = 20.18e-3/NA;                                                    % unit: kg
lambda = sqrt((beta.*plank_constant.^2)./(2*pi*mass_neon));                 % unit: m
epsilon = 1.8e3/NA;                                                         % unit: J
sigma = 2.6;                                                                % unit: Angstrom
dr = 1e-1;
r = 0.1:dr:10;
delta_coeff = 0.1;
rho_bulk = 0.033;

mu = log(rho_bulk.*(lambda.^3)) ./ beta;

%% analytic solution
v_ext_ana = v_ext_cal(epsilon, sigma, r);
rho_eq_ana = (exp(beta.*mu)./(lambda.^3)) .* exp(-beta.*v_ext_ana);


%% numerical solution
rho_initial = rho_bulk .* 0.6 .* ones(1,length(r));
mixing_parameter = 0.5;
torr = 1e-5;
rho_update(1,1) = rho_initial(1);

figure(1)
hold on;
axis([-1 10 0 1.5]);
p1 = plot(r, rho_eq_ana./rho_bulk, 'blue');


step = 0;
while (1)
    
    if (step == 0)
        rho_prev = rho_initial;
    end
    
    
    p2 = plot(r, rho_prev./rho_bulk, 'red');
    
    v_ext_ana = v_ext_cal(epsilon, sigma, r);
    
    rho_eq = (exp(beta.*mu)./(lambda.^3)) .* exp(-beta.*v_ext_ana);
    
    rho_new = (1-mixing_parameter) .* rho_prev + mixing_parameter.* rho_eq;
    
    rho_update((step+1),1) = rho_new(1);
    
    if (sqrt(sum((rho_new - rho_prev).^2)) < torr)
        break;
    else
        rho_prev = rho_new;
        step = step + 1;
    end
    
    pause(0.2);
    delete(p2);
    
end



% figure(1)
% plot(r, rho_eq_ana./rho_bulk);
% 
% legend('ideal');

% figure(1)
% plot(r, rho_eq_ana);


% %% numerical solution
% rho_bulk = 0.033;
% chemical_potential = ((lambda*1e10).^3) .* rho_bulk ./ beta;
% rho_iter0 = ones(1, length(r))*rho_bulk;








%% additional fxn
function v = v_ext_cal(epsilon_cal, sigma_cal, r)
    %% In ideal case, v_external equals to zero
    v = r.*0;
end