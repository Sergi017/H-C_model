clear variables, close all, clc
%% Hydro-Chemical code for hydration of zero-porosity rock
% Modified by Sérgio Medeiros
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Look up related stuff; approxmation of LU variables with tanh, and Pf as
% "reverse" ln-function of total density
% Look up table for atg + bru = fo + H2O reaction
load LOOK_UP_atg
Pf_lu       = P_vec*1e8/1e9;    % Pf in GPa
Pf_re       = 1.265;            % Reaction Pf in GPa
rhos_lp     = Rho_s(1);
rhos_hp     = Rho_s(end);
rhof        = mean(Rho_f);
xs_lp       = X_s_vec(1);
xs_hp       = X_s_vec(end);

% Fitting rhos and xs from LU with tanh; rhof is assumed constant
rw              = 0.06; % Width of step in tanh
rhos            = rhos_lp + (rhos_hp-rhos_lp) *(0.5*(tanh((Pf_lu-Pf_re)*2*pi/rw)+1));
xs              = xs_lp   + (xs_hp-xs_lp)     *(0.5*(tanh((Pf_lu-Pf_re)*2*pi/rw)+1));
rhostot         = min(rhos)*(1-0)*min(xs);
Phi_lu          = 1-(rhostot./(rhos.*xs));
rhoTOT_lu       = (1-Phi_lu).*rhos + Phi_lu.*rhof;
rT_lp           = rhoTOT_lu(1);
rT_hp           = rhoTOT_lu(end);
rhoTOT          = rT_lp   + (rT_hp-rT_lp)     *(0.5*(tanh((Pf_lu-Pf_re)*2*pi/rw)+1));
Pf_reverse      = real( (rw.*log((rhoTOT-rT_lp)./(-rhoTOT+rT_hp))+4*pi.*Pf_re)./4/pi );

% Start of HC model, configuration and intialization
% Physics
Lx              = 1;
betaf           = 1e11/1e9;
k0_etaf         = 1e-19/(1e-3/1e9);
n_poro          = 3;
% Numerics
nx              = 100;
dx              = Lx/(nx-1);
X               = 0:dx:Lx;
nt              = 1e10;
CFL             = 2e-2/betaf/1;
rel             = 0.25;
% Initial values, in detail. First define ambient zero-porosity values and
% then add pertubation at boundary assuming total solid mass is constant.
np                      = 1;
Phi_amb                 = zeros(size(X));
Delta_Pf                = 0.1;
Pf_amb                  = Pf_re*ones(size(X)) + Delta_Pf/2;
rhos_amb                = rhos_lp + (rhos_hp-rhos_lp)*(0.5*(tanh( (Pf_amb-Pf_re)*2*pi/rw ) + 1));
xs_amb                  = xs_lp + (xs_hp-xs_lp)*(0.5*(tanh( (Pf_amb-Pf_re)*2*pi/rw ) + 1));
rhostot                 = max( rhos_amb.*(1-Phi_amb).*xs_amb );
Pf_amb(1:np)            = Pf_amb(1) - Delta_Pf;
Pf_ini                  = Pf_amb;
rhos_ini                = rhos_lp + (rhos_hp-rhos_lp)*(0.5*(tanh( (Pf_ini-Pf_re)*2*pi/rw ) + 1));
xs_ini                  = xs_lp + (xs_hp-xs_lp)*(0.5*(tanh( (Pf_ini-Pf_re)*2*pi/rw ) + 1));
Phi_ini                 = 1-(rhostot./(rhos_ini.*xs_ini));
phi0                    = max(Phi_ini);
rhoTOT                  = (1-Phi_ini).*rhos_ini + Phi_ini.*rhof;
rhoTOT_ini              = rhoTOT;
Phi                     = Phi_ini;
rhos                    = rhos_ini;
Pf_reverse              = real( (rw.*log((rhoTOT-rT_lp)./(-rhoTOT+rT_hp))+4*pi.*Pf_re)./4/pi );
Pf                      = Pf_reverse;
Pf_ini                  = Pf_reverse;
% figure(2)
% subplot(211)
% plot(X,Pf_ini,'-k'), hold on
% plot(X,Pf_reverse,'o')
% plot([X(1) X(end)],[1 1]*Pf_re,'--xr')
% subplot(212)
% plot(X,rhoTOT,'-k')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time loop and solver.
Time                    = [0];
Dehy_front              = [ ];
Max_rhoTOT              = [ ];
for it = 1:nt
    Dehy_front(it)      = X(max(find(Pf<Pf_re)));
    Max_rhoTOT(it)      = max(real(rhoTOT));
    rhoTOT_old          = rhoTOT;
    Diffmax             = max(k0_etaf .* (Phi/phi0).^n_poro) ./ betaf;
    dt                  = CFL*dx^2*min(rhof)/max(real(rhoTOT))/Diffmax;
    if it==1;dt_ini = dt; end
    error_it            = 1;
    iter                = 0;
    while error_it > 1e-6
        iter            = iter+1;
        rhoTOT_it       = rhoTOT;
        Perm            = rhof*k0_etaf .* (Phi/phi0).^n_poro;
        Flux            = (Perm(1:end-1) + Perm(2:end))/2 .* diff(Pf)/dx;
        rhoTOT(2:end-1) = real( rhoTOT_old(2:end-1) + dt*diff(Flux)/dx );
        error_rho       = max(abs(rhoTOT_it-rhoTOT));
        Pfr             = real( (rw.*log((rhoTOT-rT_lp)./(-rhoTOT+rT_hp))+4*pi.*Pf_re)./4/pi );
        Pf              = (1-rel)*Pf + rel*Pfr;
        Pf(1:np)        = Pf_ini(1:np);
        error_Pf        = max(abs(Pf-Pfr));
        error_it        = max([error_rho error_Pf]);
        rhos            = rhos_lp + (rhos_hp-rhos_lp)*(0.5*(tanh( (Pf-Pf_re)*2*pi/rw ) + 1));
        xs              = xs_lp   + (xs_hp-xs_lp)    *(0.5*(tanh( (Pf-Pf_re)*2*pi/rw ) + 1));
        Phi             = 1-rhostot./(rhos.*xs);
    end
    if it>1, Time = [Time Time(end)+dt];  end
    t_vein=sqrt(Time);
    if it==1, Flux_ini = Flux; end
    if mod(it,10000)== 0 || it==1
        figure(7)
        subplot(311), hold off
        plot(X,Pf,'-b',X,Pf_ini,'--b'),hold on
        plot([X(1) X(end)],[1 1]*Pf_re,'--r')
        legend('P_f','P_{f0}','P_{f reaction}'), set(gca,'FontSize',12)
        title('Fluid pressure evolution')
        ylabel('P_f'),xlabel('X'),axis([-0.01 Lx min(Pf_ini)*0.98 max(Pf_ini)*1.02])
        subplot(312)
        plot(X,Phi,'-k',X,Phi_ini,'--k')
        legend('\phi','\phi_0'), set(gca,'FontSize',12)
        ylabel('\phi'),xlabel('X'),axis([-0.01 Lx -0.01 max(Phi_ini)*1.2]), grid on
        title('Porosity evolution')
%         subplot(312)
%         plot(X,rhoTOT,'-ok',X,rhoTOT_ini,'--k',X,rhos,'-ob',X,rhos_ini,'--b')
%         legend('\rho_T','\rho_{T 0}','\rho_s','\rho_{s 0}'), set(gca,'FontSize',12)
%         ylabel('\rho'),xlabel('X'),axis([-0.01 Lx 2200 3300]), grid on
%         title(['Density evolution. Maximal \rho_T: ',num2str(max(rhoTOT))])
%         if Time(end)<5e3;
%         subplot(413)
%         Phi_plot                = zeros(length(X),3);
%         Phi_plot(:,3)           = Phi;
%         Phi_plot(Pf>Pf_re,2)    = 1-Phi(Pf>Pf_re);
%         Phi_plot(Pf<=Pf_re,1)   = 1-Phi(Pf<=Pf_re);
%         ar = area(X,Phi_plot);ar(1).FaceColor=[10 141 10]/255;ar(2).FaceColor=[186 252 228]/255;ar(3).FaceColor=[14 176 246]/255;
%         xlabel('X'),ylabel('Proportion'),legend('Olivine', 'Serpentinite', 'H2O'), set(gca,'FontSize',12)
%         end
        %title('Volumetric proportions of phases')
        %         plot((X(1:end-1)+X(2:end))/2,-Flux,'-k',(X(1:end-1)+X(2:end))/2,-Flux_ini,'--k')
        %         legend('q','q_{0}'), set(gca,'FontSize',12)
        %         ylabel('Total mass flux'),xlabel('X'),axis([-0.01 Lx -50 500]), grid on
        %         title('Total mass flux evolution')
        subplot(313)
        Phi_plot                = zeros(length(X),3);
        Phi_plot(:,3)           = Phi;
        Phi_plot(Pf>Pf_re,2)    = 1-Phi(Pf>Pf_re);
        Phi_plot(Pf<=Pf_re,1)   = 1-Phi(Pf<=Pf_re);
        ar = area(X,Phi_plot);ar(1).FaceColor=[10 141 10]/255;ar(2).FaceColor=[186 252 228]/255;ar(3).FaceColor=[14 176 246]/255;
        xlabel('X'),ylabel('Proportion'),legend('Olivine', 'Serpentinite', 'H2O'), set(gca,'FontSize',12)
        set(gcf,'Position',[553.8000 50.6000 968.0000 722.4000])
        drawnow

%                        % search of slope and constant  
%         line_equation1     = polyfit(t_vein1, Dehy_front1,1);
%         line_equation2    = polyfit(t_vein2, Dehy_front2,1);
%         line_equation3    = polyfit(t_vein3, Dehy_front3,1);
%         line_equation4    = polyfit(t_vein4, Dehy_front4,1);
%         line_equation5    = polyfit(t_vein5, Dehy_front5,1);
%         line_equation6    = polyfit(t_vein6, Dehy_front6,1);
%         line_equation7    = polyfit(t_vein7, Dehy_front7,1);
%         line_equation8    = polyfit(t_vein8, Dehy_front8,1);
%         line_equation9    = polyfit(t_vein9, Dehy_front9,1);
%         line_equation10    = polyfit(t_vein10, Dehy_front10,1);
% 
%         % creation of an equation for the fitted line with the slope and
%         % constant found before
%         y_line_eq1   =  line_equation1(1)*t_vein1 + line_equation1(2);    
%         y_line_eq2  =  line_equation2(1)*t_vein2 + line_equation2(2);    
%         y_line_eq3  =  line_equation3(1)*t_vein3 + line_equation3(2);    
%         y_line_eq4  =  line_equation4(1)*t_vein4 + line_equation4(2);    
%         y_line_eq5  =  line_equation5(1)*t_vein5 + line_equation5(2);    
%         y_line_eq6  =  line_equation6(1)*t_vein6 + line_equation6(2);    
%         y_line_eq7  =  line_equation7(1)*t_vein7 + line_equation7(2);    
%         y_line_eq8  =  line_equation8(1)*t_vein8 + line_equation8(2);    
%         y_line_eq9  =  line_equation9(1)*t_vein9 + line_equation9(2);    
%         y_line_eq10  =  line_equation10(1)*t_vein10 + line_equation10(2);    
% 
%         % calculation of effective diffusivity (m^2/s) 
%         K_eff1       =   (line_equation1(1))^2; 
%         K_eff2      =   (line_equation2(1))^2; 
%         K_eff3      =   (line_equation3(1))^2; 
%         K_eff4      =   (line_equation4(1))^2; 
%         K_eff5      =   (line_equation5(1))^2; 
%         K_eff6      =   (line_equation6(1))^2; 
%         K_eff7      =   (line_equation7(1))^2; 
%         K_eff8      =   (line_equation8(1))^2; 
%         K_eff9      =   (line_equation9(1))^2; 
%         K_eff10      =   (line_equation10(1))^2;
%         K_eff_tot   = [K_eff1,K_eff2,K_eff3,K_eff4,K_eff5,K_eff6,K_eff7,K_eff8,K_eff9,K_eff10];
%         n_por       = (1:10);

        if  Time(end)>1e5; break
        end
%           save('error_Dehy_front_nx800', 'Dehy_front')
        
        
    end
end

% % resolution test
% error_nx100= abs(Dehy_front100(end)-Dehy_front800(end))/Dehy_front800(end);
% error_nx400=abs(Dehy_front400(end)-Dehy_front800(end))/Dehy_front800(end);
% error_nx600=abs(Dehy_front600(end)-Dehy_front800(end))/Dehy_front800(end);
% error_tot=[error_nx100,error_nx400,error_nx600];
% nodes=[100 400 600];

% figure(5)
% subplot(521),plot(t_vein1,Dehy_front1,'k',t_vein1, y_line_eq1,"r"),xlabel('Time(s)^{1/2}'), ylabel('vein width [m]'),title('Front reaction width with porosity exponent = 1'),equ = sprintf('y = %.2ex %.4f', line_equation1(1), line_equation1(2));text(100,0.17,equ,'FontSize',11,'Color','r')
% subplot(522),plot(t_vein2,Dehy_front2,'k',t_vein2, y_line_eq2,"r"),xlabel('Time(s)^{1/2}'), ylabel('vein width [m]'),title('With porosity exponent = 2'),equ2 = sprintf('y = %.2ex %.4f', line_equation2(1), line_equation2(2));text(100,0.17,equ2,'FontSize',11,'Color','r')
% subplot(523),plot(t_vein3,Dehy_front3,'k',t_vein3, y_line_eq3,"r"),xlabel('Time(s)^{1/2}'), ylabel('vein width [m]'),title('With porosity exponent = 3'),equ3 = sprintf('y = %.2ex %.4f', line_equation3(1), line_equation3(2));text(100,0.17,equ3,'FontSize',11,'Color','r')
% subplot(524),plot(t_vein4,Dehy_front4,'k',t_vein4, y_line_eq4,"r"),xlabel('Time(s)^{1/2}'), ylabel('vein width [m]'),title('With porosity exponent = 4'),equ4 = sprintf('y = %.2ex %.4f', line_equation4(1), line_equation4(2));text(100,0.17,equ4,'FontSize',11,'Color','r')
% subplot(525),plot(t_vein5,Dehy_front5,'k',t_vein5, y_line_eq5,"r"),xlabel('Time(s)^{1/2}'), ylabel('vein width [m]'),title('With porosity exponent = 5'),equ5 = sprintf('y = %.2ex %.4f', line_equation5(1), line_equation5(2));text(100,0.17,equ5,'FontSize',11,'Color','r')
% subplot(526),plot(t_vein6,Dehy_front6,'k',t_vein6, y_line_eq6,"r"),xlabel('Time(s)^{1/2}'), ylabel('vein width [m]'),title('With porosity exponent = 6'),equ6 = sprintf('y = %.2ex %.4f', line_equation6(1), line_equation6(2));text(100,0.17,equ6,'FontSize',11,'Color','r')
% subplot(527),plot(t_vein7,Dehy_front7,'k',t_vein7, y_line_eq7,"r"),xlabel('Time(s)^{1/2}'), ylabel('vein width [m]'),title('With porosity exponent = 7'),equ7 = sprintf('y = %.2ex %.4f', line_equation7(1), line_equation7(2));text(100,0.17,equ7,'FontSize',11,'Color','r')
% subplot(528),plot(t_vein8,Dehy_front8,'k',t_vein8, y_line_eq8,"r"),xlabel('Time(s)^{1/2}'), ylabel('vein width [m]'),title('With porosity exponent = 8'),equ8 = sprintf('y = %.2ex %.4f', line_equation8(1), line_equation8(2));text(100,0.17,equ8,'FontSize',11,'Color','r')
% subplot(529),plot(t_vein9,Dehy_front9,'k',t_vein9, y_line_eq9,"r"),xlabel('Time(s)^{1/2}'), ylabel('vein width [m]'),title('With porosity exponent = 9'),equ9 = sprintf('y = %.2ex %.4f', line_equation9(1), line_equation9(2));text(100,0.17,equ9,'FontSize',11,'Color','r')
% subplot(5,2,10),plot(t_vein10,Dehy_front10,'k',t_vein10, y_line_eq10,"r"),xlabel('Time(s)^{1/2}'), ylabel('vein width [m]'),title('With porosity exponent = 10'),equ10 = sprintf('y = %.2ex %.4f', line_equation10(1), line_equation10(2));text(100,0.17,equ10,'FontSize',11,'Color','r')
% drawnow
% %porosity exponent vs diffusivity
%         figure(6)
%         plot(n_por,log(K_eff_tot)),title("Diffusivities in function of porositiy exponents"),xlabel("Porosity exponents"),ylabel("Log(diffusivity [m^{2}/s])")
% 
% 
% 
% figure(9)
% plot(log10(nodes),log10(error_tot))
% xlabel('Logarithm of number of nodes (nx)')
% ylabel('Logarithm of error [(nx-nx800)/nx800]')







