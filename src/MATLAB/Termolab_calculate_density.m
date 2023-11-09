%% 
clc, clear
% runname = 'antigorite'
% load([runname '_linprog']);
%%
Tprofile         = 650+273.15;
Pprofile         = linspace(16,30,1000) .*1e8; %P [bar]
% Numerics
delP     = 1e5;                                                          % for numerical differentiation
delc     = 1e-5;  


X    = {'Si','Mg','H','O'};
load tl_dataset
fluid  = 'H2O,tc-ds55';
phases = {'atg,tc-ds55','br,tc-ds55','fo,tc-ds55','ta,tc-ds55','en,tc-ds55','anth,tc-ds55','H2O,tc-ds55'};
[T2d,P2d] = ndgrid(Tprofile,Pprofile);
[G,Nphs,pc_id] = tl_gibbs_energy(T2d(:),P2d(:),phases);
Nsys = Nphs(:,1); 
LB  = zeros(1,size(G,1)); % do not look for alph below zero, because stable phase amount cannot be negative.
for iPT = 1:length(T2d(:))
    alph(iPT,:) = linprog(G(:,iPT),[],[],Nphs,Nsys,LB); % The Gibbs energy minimization       
end

[G_dP,Nphs,~]  = tl_gibbs_energy(T2d(:),P2d(:)+delP,phases);% get Gibbs energy at P+dP
%%
molm = molmass_fun(X);                                 % get molar mass of the components

for i = 1:length(Pprofile)
Nphases = Nphs([14,12,1,8],:)
% Postprocessing
Vmol     = (G_dP(:,i)-G(:,i))/delP;                                           % Equation 54
Mmol     = Nphases'*molm;                                         % Equation 56
rho      = Mmol./Vmol;                                                        % Equation 57
phim     = alph(i,:)/sum(alph(i,:));                                          % Equation 52
phi      = phim.*Vmol'./(Vmol'*phim');                                        % Equation 53
phiw     = phim.*Mmol'./(Mmol'*phim');                                        % Equation 55
Cwt      = Nphases.*repmat(molm,1,size(Nphs,2))./repmat(Mmol',size(Nphases,1),1);  % Equation 58
fluid_id =  strcmp(phases(pc_id),fluid);                                     % Find index of fluid
solid_id = ~fluid_id;                                                        % Find index of solids
rhos(i)     = rho(solid_id)'*phi(solid_id)'/sum(phi(solid_id));                  % Equation 59
rhof(i)     = rho(fluid_id)'*phi(fluid_id)'/sum(phi(fluid_id));                  % Fluid density

end


figure
subplot(311)
hold on
plot(Pprofile, rhos, 'o','DisplayName','Thermodynamic data')
ylabel('Density [kg/m^3]')
xlabel('Pressure [kbar]')
title('A) Solid density vs fluid pressure')
set(gca,'ytick',[2000:250:3500])
axis([Pfitting(1) Pfitting(2) 2400 3500])
grid on
plot([Preaction Preaction],[0 4000],'--k')
text(17,2500,'fo + en + H2O','fontangle','italic','fontsize',10,'Rotation',90)
text(22,3000,'atg','fontangle','italic','fontsize',10)
