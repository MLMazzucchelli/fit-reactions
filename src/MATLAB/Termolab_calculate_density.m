% Calculate the desity of solid and fluid along an isothermal path with
% Thermolab.
% Save the data to a .mat file for further fitting.
clc, clear
runname = 'antigorite';
%% Calculate Gibbs energy
Tprofile  = 650+273.15;
Pprofile  = linspace(16,30,1000) .*1e8; %P [bar]
% Numerics
delP      = 1e5;    % for numerical differentiation
delc      = 1e-5;  


X               = {'Si','Mg','H','O'};
load tl_dataset
fluid           = 'H2O,tc-ds55';
phases          = {'atg,tc-ds55','br,tc-ds55','fo,tc-ds55','ta,tc-ds55','en,tc-ds55','anth,tc-ds55', fluid};
[T2d,P2d]       = ndgrid(Tprofile,Pprofile);
[G,Nphs,pc_id]  = tl_gibbs_energy(T2d(:),P2d(:),phases);
Nsys = Nphs(:,1); 
LB  = zeros(1,size(G,1)); % do not look for alph below zero, because stable phase amount cannot be negative.
for iPT = 1:length(T2d(:))
    alph(iPT,:) = linprog(G(:,iPT),[],[],Nphs,Nsys,LB); % The Gibbs energy minimization       
end
[G_dP,Nphs,~]  = tl_gibbs_energy(T2d(:),P2d(:)+delP,phases);% get Gibbs energy at P+dP for derivatives

%% Save Thermolab calculation
save(sprintf("%s_isothermal_%.2fK_Thermolab_linprog.mat", runname, Tprofile));
%% Load Thermolab calculation
clc, clear
load("antigorite_isothermal_923.15K_Thermolab_linprog.mat");
%% Postprocessing
molm     = molmass_fun(X);                                                          % get molar mass of the components
Nphases  = Nphs([14,12,1,8],:);
fluid_id =  strcmp(phases(pc_id),fluid);                                            % Find index of fluid
solid_id = ~fluid_id;                                                               % Find index of solids

Mmol     = Nphases'*molm;                                                           % Equation 56, molar mass of each phase
Cwt      = Nphases.*repmat(molm,1,size(Nphs,2))./repmat(Mmol',size(Nphases,1),1);   % Equation 58

% Calculate density variation with P
for i = 1:length(Pprofile)
    phim     = alph(i,:)/sum(alph(i,:));                            % Equation 52, molar fraction of each stable phase
    phiw     = phim.*Mmol'./(Mmol'*phim');                          % Equation 55, weigth fraction of each phase
    Vmol     = (G_dP(:,i)-G(:,i))/delP;                             % Equation 54, molar volume of each phase    
    phi      = phim.*Vmol'./(Vmol'*phim');                          % Equation 53, volume fraction of each phase
    rho      = Mmol./Vmol;                                          % Equation 57, density of each phase [kg/m3]
    rhos(i)  = rho(solid_id)'*phi(solid_id)'/sum(phi(solid_id));    % Equation 59, density of solid [kg/m3]
    rhof(i)  = rho(fluid_id)'*phi(fluid_id)'/sum(phi(fluid_id));    % Equation 59, fluid density [kg/m3]
end
rhow         = rho_H2O(Tprofile,Pprofile,'ZD05');                   % Calculate density of water in the same P,T range
%% Save data for fitting
constantTdata       = table();
constantTdata.T     = T2d'-273.15;
constantTdata.P     = P2d'.*1e-8;
constantTdata.rhos  = rhos';
constantTdata.rhow  = rhow';
save("data_constantT_Thermolab.mat","constantTdata")
