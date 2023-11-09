% Import data from Perplex, plot them and save the isothermal profile to a
% Matlab file for further fitting.
%
% The data imported from Perplex are:
% 1) a P-T table with: T, P, rho solid, rho fluid, id of assemblage
% 2) an isothermal profile at a given T with increasing P with variables:
% T, P, rho solid, rho fluid, weight percent of each non-volatile component
% 3) an isothermal profile for water at a given T with increasing P with
% variables: T, P, rho of water
%
% A table is created and saved to a .mat file for further fitting, with
% fields: T, P, rhos, rhow, sum of weight fractions of the non-volatile components

clc, clear
Tprof           = 650;          %temperture of isothermal profile [C]
Trange_plot     = [400, 800];   %T range for the plots [C]
Prange_plot     = [1 30];       %P range for the plots [kbar]
%% Import data generated from Perplex
PTfile              = readtable('atg_PT.tab', 'FileType', 'text', 'VariableNamingRule', 'preserve');
constantTdata       = readtable('antigorite_650C.tab', 'FileType', 'text', 'VariableNamingRule', 'preserve');
constantTwater      = readtable('water_650C.tab', 'FileType', 'text', 'VariableNamingRule', 'preserve');       %Thermo properties of H2O
constantTdata.Properties.VariableNames  = {'T', 'P', 'rhos', 'rhof', 'SiO2_wt', 'MgO_wt', 'fo_wt'};            %This needs to be set by the user
constantTwater.Properties.VariableNames = {'Tf', 'Pf', 'rhow'};
% Rearrange data
constantTdata.nonVolatile_wt = constantTdata.MgO_wt + constantTdata.SiO2_wt;                                    %This needs to be set by the user
constantTdata = [constantTdata,constantTwater];

%% Reshape data for 2D plot
PTdata      = PTfile.Variables; %convert to array
TK          = PTdata(:,1);
Pbar        = PTdata(:,2);
rhos        = PTdata(:,3);
rhow        = PTdata(:,4);
assemblage  = PTdata(:,5);
nx          = sqrt(length(PTdata(:,1))); %number of colums (the grid is square)

%Reshape columns into matrices
T           = reshape(TK, nx, nx);      
P           = reshape(Pbar, nx, nx);    
rhos   = reshape(rhos, nx, nx);
rhow     = reshape(rhow, nx, nx);
assemblage  = reshape(assemblage, nx, nx);

% Convert units to C and kbar
T                       = T-273.15; %T in C
P                       = P*1e-3; %P in kbar
constantTdata.T         = constantTdata.T-273.15;
constantTdata.P         = constantTdata.P*1e-3;
constantTdata.Pf     = constantTwater.Pf*1e-3;
constantTdata.Tf     = constantTwater.Tf-273.15;

%Save data to mat files
save("data_constantT_Perplex.mat","constantTdata")
%% Find the pressure of reaction as the discontinuity in solid density
[value,idx]         = max(abs(diff(constantTdata.rhos)));
Preaction           = constantTdata.P(idx);                               % Pressure of reaction [kbar]
%% Plot 2D data
figure
subplot(2,2,1)
contourf(T,P,rhos), colorbar, shading interp %plot
caxis([2.6e3 3.2e3]);
hold on
contour(T,P,assemblage, 'r')
hold on
plot([Tprof,Tprof], [0, 520], '--r')
hold on
xlim(Trange_plot), ylim(Prange_plot)
title('Density solid (kg/m3)')
xlabel('T(C)'), ylabel('P(kbar)')
grid on

subplot(2,2,3)
contourf(T,P,rhow), colorbar, shading flat %plot
hold on
contour(T,P,assemblage, 'r')
xlim(Trange_plot), ylim(Prange_plot)
title('Density H2O (kg/m3)')
xlabel('T(C)'), ylabel('P(kbar)')

subplot(2,2,2)
plot(constantTdata.P, constantTdata.rhos, 'DisplayName', 'rho solid')
hold on
plot(constantTdata.Pf, constantTdata.rhow, 'g', 'DisplayName', 'rho water')
hold on
plot(constantTdata.P, constantTdata.rhof, '--r','DisplayName', 'rho fluid')
hold on
hold on
plot([Preaction ,Preaction], [0, 5000], '-.k')
xlim([10, 30]), ylim([0, 4000])
text(12,1250,'fo + ta + H2O','fontangle','italic','fontsize',14,'Rotation',90)
text(17,1250,'fo + en + H2O','fontangle','italic','fontsize',14,'Rotation',90)
text(22,1500,'atg','fontangle','italic','fontsize',14)
title(sprintf('Solid and fluid density at %.1f °C', Tprof))
ylabel('Density (kg/m3)'), xlabel('P(kbar)')
legend

subplot(2,2,4)
% plot(constantTdata.P, constantTdata.MgO_wt./100, 'DisplayName', 'MgO wt')
% hold on
% plot(constantTdata.P, constantTdata.SiO2_wt./100, 'DisplayName', 'SiO2 wt')
hold on
plot(constantTdata.P, constantTdata.nonVolatile_wt./100, 'DisplayName', 'Non-volatile component wt')
hold on
plot(constantTdata.P, constantTdata.fo_wt./100, 'DisplayName', 'fo wt')
hold on
hold on
plot([Preaction ,Preaction], [0, 5000], '-.k')
xlim([10, 30]), ylim([0.4, 1.05])
text(12,0.65,'fo + ta + H2O','fontangle','italic','fontsize',14,'Rotation',90)
text(17,0.65,'fo + en + H2O','fontangle','italic','fontsize',14,'Rotation',90)
text(22,0.65,'atg','fontangle','italic','fontsize',14)
title(sprintf('Mass fraction at %.1f °C', Tprof))
ylabel('X'), xlabel('P(kbar)')
grid on
legend

