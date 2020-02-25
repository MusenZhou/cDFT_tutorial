% This program is used to illustrate the very basic idea of classic density
% functional theory (cDFT) in the ideal limit , which is without the
% consideration of the term of excess Helmholtz energy
% For more details, please refer to paper (DOI: 10.1021/ed500049m)
clear;
clc;

% Description: variable rho here represents density.
% dr is the step size of system and r is the length of the system and
% subject to appropriate changes from the user
dr = 0.5e-1;
r = 0.1:dr:10;
% Picard iteration parameters, which are subject to appropriate changes
% from the user (torr is the tolearance)
mixing_parameter = 0.5;
torr = 1e-5;



%% constant define
% Boltzmann constant, unit: J/Klevin
kb = 1.38064852e-23;
% Temperature, unit: Klevin
T = 300;
% 1/kBT, unit: 1./J
beta = 1./(kb.*T);
% plank constant, unit: J*m
plank_constant = 6.626070040e-34;
% Avogadro constant, unit: mol-1
NA = 6.022140857e23;
% mass of neon, unit: kg
mass_neon = 20.18e-3/NA;
% thermal wavelength, unit: m
lambda = sqrt((beta.*plank_constant.^2)./(2*pi*mass_neon));
% LJ energy parameter, unit: J
epsilon = 1.8e3/NA;
% LJ size parameter, unit: Angstrom
sigma = 2.6;
% bulk density
rho_bulk = 0.033;
% chemical potential
mu = log(rho_bulk.*(lambda.^3)) ./ beta;

%% analytic solution
v_ext_ana = v_ext_cal(epsilon, sigma, r);
rho_eq_ana = (exp(beta.*mu)./(lambda.^3)) .* exp(-beta.*v_ext_ana);


%% numerical solution
rho_initial = rho_bulk .* ones(1,length(r));

figure(1)
hold on;
axis([-1 10 0 3]);
p1 = plot(r, rho_eq_ana./rho_bulk, 'blue');
xlabel('r', 'FontSize', 15);
ylabel('\rho / \rho_{bulk}', 'FontSize', 15);
legend('analytical solution');

step = 0;
while (1)
    
    if (step == 0)
        rho_prev = rho_initial;
    end
    p2 = plot(r, rho_prev./rho_bulk, 'red');
    legend('analytical solution', 'numerical solution');
    
    rho_eq = (exp(beta.*mu)./(lambda.^3)) .* exp(-beta.*v_ext_ana);
    rho_new = (1-mixing_parameter) .* rho_prev + mixing_parameter.* rho_eq;
    
    if (sqrt(sum((rho_new - rho_prev).^2)) < torr)
        break;
    else
        rho_prev = rho_new;
        step = step + 1;
    end
    
    pause(0.3);
    delete(p2);
    
end










%% external potential calculated by LJ potential form
function v = v_ext_cal(epsilon_cal, sigma_cal, r)
    v = 4*epsilon_cal.*( (sigma_cal./r).^12 - (sigma_cal./r).^6 );
    for i = 1:length(r)
        if (v(i) >= 100)
            v(i) = 100;
        end
    end
end