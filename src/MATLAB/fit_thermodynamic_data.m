% Fit the density of solid, fluid and mass fraction of the non-volatile
% component as a function of P along an isothermal path.
%
% Import the data for the fitting as a matlab table with fields: T [C], P [kbar], rho of solid [kg/m3],
% rho of water [kg/m3], total weigth fraction of the non-volatile components [].
%
% Export the parameters of the 3 fittings as a text file.
% Params(1:8): fitting of density of solid
% Params(9:11): fitting of density of fluid
% Params(12:15); fitting of mass fraction of solid component.
% Params(16): Pressure of the transition [Pa] on which the pressure is
% rescaled.
% Important: the fitting functions give the solid, fluid and the mass fraction of the non-volatile
% component, respectively, as a function of the SCALED pressure.

clc, clear
Pfitting            = [16 30]; %set the pressure range for fitting [kbar]


%Load data
data                = load("data_constantT_Perplex.mat");
constantTdata       = data.constantTdata;

% Get the required properties from the imported data
data_for_fitting    = data.constantTdata(constantTdata.P>=Pfitting(1) & constantTdata.P<Pfitting(2),:); %restrict the pressure range for fitting
data_for_fitting    = data_for_fitting(1:1:end,:);
Rho_s_LU            = data_for_fitting.rhos;                     % Precalculated solid density as function of pressure
Rho_f_LU            = data_for_fitting.rhow;                     % Precalculated water density as function of water pressure
X_LU                = data_for_fitting.nonVolatile_wt./100;      % Precalculated mass fraction of MgO as function of fluid pressure
P_LU                = data_for_fitting.P;                        % Corresponding fluid pressure array [kbar]; 

%Convert P units from kbar to Pa
P_LU               = P_LU*1e8;
Pfitting           = Pfitting*1e8;
%Find the pressure of reaction as the discontinuity in solid density
[value,idx]         = max(abs(diff(Rho_s_LU)));
Preaction           = P_LU(idx);                               % Pressure of reaction [Pa]

% Scale pressure units
P_ini               = 1;                                       % Initial ambient pressure [Pa]
Pini_Pappl          = P_ini/Preaction;
P_LU_scaled         = P_LU*Pini_Pappl;                         % Rescale the look-up table P to the P of transition
Preaction_scaled    = Preaction*Pini_Pappl;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solid density

% Parameters for solid density
p_min_ana           = min(P_LU_scaled);
p_max_ana           = max(P_LU_scaled);                 
rho_s_max           = Rho_s_LU(1);
rho_s_min           = min(Rho_s_LU);
rho_s_dif           = rho_s_max-rho_s_min;

% Fitting of density of solid
Rho_s_func   = @(Rparams,P) -tanh(Rparams(1)*(P-Preaction_scaled))*(rho_s_dif/2+Rparams(2)) + (rho_s_dif/2-Rparams(2)) + rho_s_min + ((P-p_min_ana)./p_max_ana.*Rparams(3));
Rparams      = nlinfit(P_LU_scaled,Rho_s_LU,Rho_s_func,[1,1,1]);
Rho_s_param  = [Rparams(1) Preaction_scaled rho_s_dif Rparams(2) rho_s_min p_min_ana p_max_ana Rparams(3)];
Rho_s_ana    = Rho_s_func(Rparams,P_LU_scaled); %Rho_s_ana [kg/m3] as a function of the SCALED pressure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fluid density

% Parameters for fluid density
rho_f_max_ana       = max(Rho_f_LU);

% Fitting of density of fluids
Rho_f_func   = @(Rparams,P) rho_f_max_ana*log(P+Rparams(1)).^Rparams(2);
Rparams      = nlinfit(P_LU_scaled,Rho_f_LU,Rho_f_func,[1,1]);
Rho_f_param  = [rho_f_max_ana Rparams];
Rho_f_ana    = Rho_f_func(Rparams,P_LU_scaled); %Rho_f_ana [kg/m3] as a function of the SCALED pressure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mass fraction

% Parameters for mass fraction
x_max               = max(X_LU);
x_min_ana           = min(X_LU);
x_dif_ana           = x_max-x_min_ana;

% Fitting of mass fraction
X_func   = @(Rparams,P) -tanh(Rparams(1)*(P-Preaction_scaled))*x_dif_ana/2 + x_dif_ana/2 + x_min_ana;
Rparams  = nlinfit(P_LU_scaled,X_LU,X_func,[1]);
X_param  = [Rparams(1) Preaction_scaled x_dif_ana x_min_ana];
X_ana    = X_func(Rparams,P_LU_scaled);  %X_ana [] as a function of the SCALED pressure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rescale back to physical P
P_ana      = P_LU_scaled/Pini_Pappl; %scale P back to physical value in Pa

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save fitting parameters and the P of transition
Tab = table;
Tab.rhos_params     = Rho_s_param;
Tab.rhof_params     = Rho_f_param;
Tab.X_params        = X_param;
Tab.Preaction_Pa    = Preaction; 
writetable(Tab, 'fitting_parameters.txt')


% Plotting
figure
subplot(311)
hold on
plot(P_ana, Rho_s_LU, 'o','DisplayName','Thermodynamic data')
plot(P_ana, Rho_s_ana,'-r','DisplayName','Analytical fit')
%plot(P_ana, Rho_s_ana(end)+100*erfc(5000*(P_ana-Preaction)),'-g','DisplayName','Analytical fit')
ylabel('Density [kg/m^3]')
xlabel('Pressure [Pa]')
title('A) Solid density vs fluid pressure')
set(gca,'ytick',[2000:250:3500])
axis([Pfitting(1) Pfitting(2) 2400 3500])
grid on
plot([Preaction Preaction],[0 4000],'--k')
text(17,2500,'fo + en + H2O','fontangle','italic','fontsize',10,'Rotation',90)
text(22,3000,'atg','fontangle','italic','fontsize',10)

subplot(312)
hold on
plot(P_ana, Rho_f_LU, 'o','DisplayName','Thermodynamic data')
plot(P_ana, Rho_f_ana, '-r','DisplayName','Analytical fit')
ylabel('Density [kg/m^3]')
xlabel('Pressure [Pa]')
title('B) Fluid density vs fluid pressure')
set(gca,'ytick',[250:250:1500])
axis([Pfitting(1) Pfitting(2)  750 1500])
grid on
plot([Preaction Preaction],[0 4000],'--k')



subplot(313)
hold on
plot(P_ana, X_LU, 'o', 'DisplayName','Thermodynamic data')
plot(P_ana,X_ana, '-r','DisplayName','Analytical fit (manual)')
%plot(P_ana,X_ana_fit, '-c','DisplayName','Analytical fit')
ylabel('X_s [ ]')
xlabel('Pressure [Pa]')
title('C) Mass fraction of non-volatile component vs fluid pressure')
axis([Pfitting(1) Pfitting(2) 0.65 1.1])
set(gca,'ytick',[0.5:0.1:1])
grid on
set(gcf,'renderer','opengl')
plot([Preaction Preaction],[0.4 1.1],'--k')

text(17,0.70,'fo + en + H2O','fontangle','italic','fontsize',10,'Rotation',90)
text(22,0.75,'atg','fontangle','italic','fontsize',10)


set(gcf,'position',[419.4000 66.6000 629.6000 696.4000])

