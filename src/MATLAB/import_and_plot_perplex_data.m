clc, clear
Tprof           = 650;          %temperture of isothermal profile [C]
Trange_plot     = [400, 800];
Prange_plot     = [1 30];
%% Import data generated from Perplex
PTfile              = readtable('atg_PT.tab', 'FileType', 'text', 'VariableNamingRule', 'preserve');
constantTdata       = readtable('antigorite_650C.tab', 'FileType', 'text', 'VariableNamingRule', 'preserve');
constantTwater      = readtable('water_650C.tab', 'FileType', 'text', 'VariableNamingRule', 'preserve');       %Thermo properties of H2O
constantTdata.Properties.VariableNames  = {'T', 'P', 'rho_solid', 'rho_fluid', 'SiO2_wt', 'MgO_wt', 'fo_wt'};
constantTwater.Properties.VariableNames = {'T_h2o', 'P_h2o', 'rho_h2o'};
% Rearrange data
constantTdata.nonVolatile_wt = constantTdata.MgO_wt + constantTdata.SiO2_wt;   %This needs to be set by the user
constantTdata = [constantTdata,constantTwater];

%% Reshape data for 2D plot
PTdata      = PTfile.Variables; %convert to array
TK          = PTdata(:,1);
Pbar        = PTdata(:,2);
rho_solid   = PTdata(:,3);
rho_h2o     = PTdata(:,4);
assemblage  = PTdata(:,5);
nx          = sqrt(length(PTdata(:,1))); %number of colums (the grid is square)

%Reshape columns into matrices
T           = reshape(TK, nx, nx);      
P           = reshape(Pbar, nx, nx);    
rho_solid   = reshape(rho_solid, nx, nx);
rho_h2o     = reshape(rho_h2o, nx, nx);
assemblage  = reshape(assemblage, nx, nx);

% Convert units to C and kbar
T                       = T-273.15; %T in C
P                       = P*1e-3; %P in kbar
constantTdata.T         = constantTdata.T-273.15;
constantTdata.P         = constantTdata.P*1e-3;
constantTdata.P_h2o     = constantTwater.P_h2o*1e-3;
constantTdata.T_h2o     = constantTwater.T_h2o-273.15;

%Save data to mat files
save("data_constantT.mat","constantTdata")
%% Find the pressure of reaction as the discontinuity in solid density
[value,idx]         = max(abs(diff(constantTdata.rho_solid)));
Preaction           = constantTdata.P(idx);                               % Pressure of reaction [kbar]
%% Plot 2D data
figure
subplot(2,2,1)
contourf(T,P,rho_solid), colorbar, shading interp %plot
caxis([2.6e3 3.2e3]);
hold on
contour(T,P,assemblage, 'r')
hold on
plot([Tprof,Tprof], [0, 520], '--r')
hold on
%plot([0, 1000], [Ptrans2 , Ptrans2 ], '--k')
xlim(Trange_plot), ylim(Prange_plot)
title('Density solid (kg/m3)')
xlabel('T(C)'), ylabel('P(kbar)')
grid on

subplot(2,2,3)
contourf(T,P,rho_h2o), colorbar, shading flat %plot
hold on
contour(T,P,assemblage, 'r')
xlim(Trange_plot), ylim(Prange_plot)
title('Density H2O (kg/m3)')
xlabel('T(C)'), ylabel('P(kbar)')

subplot(2,2,2)
plot(constantTdata.P, constantTdata.rho_solid, 'DisplayName', 'rho solid')
hold on
plot(constantTdata.P_h2o, constantTdata.rho_h2o, 'g', 'DisplayName', 'rho water')
hold on
plot(constantTdata.P, constantTdata.rho_fluid, '--r','DisplayName', 'rho fluid')
hold on
%plot([Ptrans2 ,Ptrans2 ], [0, 5000], '--k')
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
%plot([Ptrans2 ,Ptrans2], [0, 5000], '--k', 'DisplayName', 'atg -> fo + en +H2O')
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


save("data_constantT.mat","constantTdata")
clear constantTdata