clc;
close all;
clear all;
constant = m_constant;

[tmp,str.main]=system('pwd');

iseed = 1;
ikslice = 1;


SED = load('./SEDavg.mat');
NMD = load('./NMD.mat');

% PT_PERC = 0.01;
% c01 = linspace(1,15,1);%(10,30,1)
% c02 = linspace(0.001,0.08,1);%(0.08,0.10,1)
%% Loop_SED_fig
if exist(['./Loop_SED_fig'], 'file')~=0
    system(['rm -r ' ' Loop_SED_fig']);
end

str.cmd = [ 'mkdir Loop_SED_fig'];
system(str.cmd);

fig_seed = 'Loop_SED_fig/'; % remove after use

%% All peaks for a particular k-point -- expected to be in increasing order

VOLUME =(NMD.Nx*NMD.Ny*NMD.Nz)*det(NMD.latvec)/((NMD.LJ.sigma*(1e+10))^3);% Ang^3 --> LJ (NMD.Nx*NMD.Ny*NMD.Nz)*det(NMD.latvec)/(NMD.LJ.sigma*(1e+10)^3)
Temp = 300; %in K
x_fac = NMD.constant.hbar/(NMD.constant.kb*Temp*NMD.LJ.tau); %hbar/(Kb*T)
factor = NMD.constant.kb/(NMD.LJ.sigma*NMD.LJ.tau*VOLUME);

str_write=strcat(str.main, '/', int2str(iseed),'/kappa_mode_lifetime_omega.txt');
if exist(['./' int2str(iseed) '/kappa_mode_lifetime_omega.txt'], 'file')~=0
    system(['rm -f ./' int2str(iseed) '/kappa_mode_lifetime_omega.txt']);
end


cnt = 1;
for ikpt = 1:size(NMD.kptlist(:,1:3,ikslice),1)
    
    
    tic
    %     PT_PERC_LEFT = 0.1; PT_PERC_RIGHT = 0.1;
    INV_PERC = 1.0;
    groupvel_i = (ikpt-1)*NMD.NUM_ATOMS_UCELL*3;
    
    for i = 1:NMD.NUM_MODES
        
           % ZA
        if (i==1)
            PT_PERC = 0.01;
            c01 = linspace(1,3,1);
            c02 = linspace(0.0278,0.035,1);
            % TA
        elseif (i==2)
            PT_PERC = 0.01;
            c01 = linspace(1,8,1);
            c02 = linspace(0.0278,0.025,1);
            % LA
        elseif (i ==3)
            PT_PERC = 0.01;
            c01 = linspace(1,35,1);
            c02 = linspace(0.0278,0.038,1);
            % Optical
        else % Other than acoustic modes will have diffrent initial parameters
            PT_PERC = 0.01;
            c01 = linspace(1,10,1);
            c02 = linspace(0.0278,0.035,1);
        end
        
        
        start =1;
        [Imax,Jmax] = max(SED.SED(ikpt,:,i));
        
        %--
        if Jmax < (2) || Jmax > (NMD.NUM_TSTEPS/2 - 1)
            i = i + 1;
            continue
        end
        
        %--
        %#####################################################################
        
%         shift = 200;
        
        %#####################################################################
        
        
        %Find wleft
        [I,J] = find(SED.SED(ikpt,start:start+Jmax,i) <...
            PT_PERC*SED.SED(ikpt,start+Jmax,i) );% find --> I willl be either 1 or 0. J will give cell no. of those values.
        
        
        
        %################################################################################
        %         wleft = Jmax - shift;
        wleft = (length(I));
        %                 wleft = start+(length(I));
        
        
        %################################################################################
        if wleft < 1
            i = i + 1;
            continue
            
        end
        %------NEW----------
        %Find wright
        [I,J] = find(SED.SED(ikpt,start+Jmax:end,i) <...
            PT_PERC*SED.SED(ikpt,start+Jmax,i) );
        
        %################################################################################
        %         wright = Jmax + shift;
        wright = length(SED.SEDavg(ikpt,:,i))-(length(I));%+50;
        
        
        %                 wright = start+Jmax + (length(I));%+50;
        
        %################################################################################
        
        if wright > NMD.NUM_TSTEPS/2
            i = i + 1;
            continue
        end
        
        %% Creating a 2D array for initial guess mode_guess_c01 & mode_guess_c02
        
        %         c01 = linspace(10,30,5);
        %         c02 = linspace(0.08,0.10,5);
        count = 1;
        
        for x = 1:size(c01,2)
            for y = 1:size(c02,2)
                mode_guess_c01 = c01(x);
                mode_guess_c02 = c02(y);
                
                c0 = [ mode_guess_c01*Imax, mode_guess_c02, SED.omega(Jmax)];%start+
                
                lb(1:length(c0)) = 0.0; ub(1:3:length(c0)) = 100000*Imax;
                ub(2:3:length(c0)) = 1000*1e15;
                ub(3:3:length(c0)) = 1000*SED.omega(length(SED.omega));
                
                % To begin, define the parameters in terms of one variable c
                % Then define the curve as a function of the parameters c and the data w:
                lor_func = @(c,w)(c(1))./(1 + ( (w - c(3))./ c(2) ).^2 );
                
                options = optimset('MaxIter',5000,'MaxFunEvals',5000,'TolFun',1e-5,'TolX',1e-5);
                
                %         [c_fit] = lsqcurvefit(lor_func,c0,SED.omega(wleft:wright),...
                %             SED.SEDpoi_kpt_avg(ikpt,wleft:wright,i),lb,ub,options);
                
                %                 [c_fit,resnorm,residual,exitflag,output]  = lsqcurvefit(lor_func,c0,SED.omega(wleft:wright),...
                %                     SED.SEDpoi_kpt_avg(ikpt,wleft:wright,i),lb,ub,options);
                %
                [c_fit,resnorm,residual,exitflag,output]  = lsqcurvefit(lor_func,c0,SED.omega(wleft:wright),...
                    SED.SEDavg(ikpt,wleft:wright,i),lb,ub,options);
                
                
                mode_RESNORM(count,1) = mode_guess_c01;
                mode_RESNORM(count,2) = mode_guess_c02;
                mode_RESNORM(count,3) = c_fit(1);
                mode_RESNORM(count,4) = c_fit(2);
                mode_RESNORM(count,5) = c_fit(3);
                mode_RESNORM(count,6) = resnorm;
                
                count = count+1;
            end
        end
        %       %Taking the minimum row based upon the min. resnorm
        %                 [value,ind] = min(mode_RESNORM(:,6));
        %       %Taking the maximum row based upon the min. resnorm
        [value,ind] = max(mode_RESNORM(:,6));
        
        
        
        c_fit(1) = mode_RESNORM(ind,3);
        c_fit(2) = mode_RESNORM(ind,4);
        c_fit(3) = mode_RESNORM(ind,5);
        
        ikpt_imode_c01_c02(cnt,:) = [ikpt,i,mode_RESNORM(ind,1),mode_RESNORM(ind,2)];
        cnt = cnt+1;
        %% SED figures
        
%         figure(1);        
%        
%         semilogy( SED.omega(1,wleft:wright),SED.SEDavg(ikpt,wleft:wright,i),'ob');
%         hold on;
%         semilogy(SED.omega(wleft:wright),lor_func([c_fit],SED.omega(wleft:wright)),'*-k','Markersize',1);
%         hold off;
% %         plot( SED.omega(1,wleft:wright),SED.SEDavg(ikpt,wleft:wright,i),'ob');
% %         hold on;
% %         plot(SED.omega(wleft:wright),lor_func([c_fit], SED.omega(wleft:wright)),'*-k','Markersize',1);
% %         hold off;
%         
% %         xt = [0 10 20 30 40 50 60 70]./NMD.constant.w_LJ2THz;
% %         set(gca, 'XTick', xt, 'XTickLabel', xt * NMD.constant.w_LJ2THz);
%         xticks([]);
%         fq_left = SED.omega(1,wleft);
%         fq_right = SED.omega(1,wright);
%         fq_mid = (fq_left+fq_right)/2;
%         xticks([fq_left, fq_mid, fq_right]); % Set custom x-axis ticks
%         xticklabels({num2str(fq_left*(NMD.constant.w_LJ2THz)), ...
%                                 num2str(fq_mid*(NMD.constant.w_LJ2THz)), ...
%                                     num2str(fq_right*(NMD.constant.w_LJ2THz))});
%         
%         xlabel('$\omega$(THz)','Interpreter','latex','fontsize',12);
%         ylabel('SED(LJ)','Interpreter','latex','fontsize',12);
%         w_0 = num2str(c_fit(3)*NMD.constant.w_LJ2THz);
%         xline(c_fit(3),'--k',{w_0});                
%         legend off;
        %% Store data
        
        %Store separate liftimes and frequencies for single and MULTIPLE FITS
        center=c_fit(3); lifetime=1/(2*c_fit(2));
        
        SED_FIT.freq(ikpt,i) = center;
        SED_FIT.life(ikpt,i) = lifetime;
        
        % At Gamma : Frequency of LA,TA  ZA = 0 
        if ((NMD.kpt.cart_uvw(ikpt,1) == 0 && NMD.kpt.cart_uvw(ikpt,2) == 0 && NMD.kpt.cart_uvw(ikpt,3) == 0) && ...
                 (i == 1 || i == 2 || i == 3))
            flag=1
            SED_FIT.freq(ikpt,i) = 0.0;
            SED_FIT.life(ikpt,i) = 0.0;
        end        
        
        x_cv = x_fac*NMD.freq(ikpt,i);
        cv_ph = (x_cv^2)*(exp(x_cv)/((exp(x_cv)-1)^2)); %specific heat per mode divided by Kb
        groupvel = NMD.vel(groupvel_i+i,1:3); %LJ
        kappa_ph = lifetime*(groupvel.^2); %LJ
        %         SED.kappa(ikpt,i,1:3) = factor*kappa_ph; % W/mK --> checked
        SED.kappa(ikpt,i,1:3) = cv_ph*factor*kappa_ph; % W/mK --> checked
        
        %first 5 elements in array below are in LJ units, kappa components are in W/mK
        output_life = [SED_FIT.freq(ikpt,i),SED_FIT.life(ikpt,i),...
            groupvel(1),groupvel(2),groupvel(3),SED.kappa(ikpt,i,1),...
            SED.kappa(ikpt,i,2),SED.kappa(ikpt,i,3),ikpt, i,cv_ph];
        dlmwrite(str_write,output_life,'-append','delimiter',' ');
               
        disp(i);
                
    end %end of i
%     mean_cv_ph = mean(output_life(:,11));
%     mean_gv = sqrt(output_life(:,11));
%     mean_tau = mean(output_life(:,2));
%     [kpt,mean_cv_ph,mean_gv,mean_tau]=[ikpt,mean(i),tau];
    
    ikpt
    toc
end %end of kpt

SED.omega_tau_gv_kappa = dlmread(str_write);
save(strcat(str.main,'/SEDomega_tau_gv_kappa.mat'), '-struct', 'SED');

%% tau vs omega

str_write=strcat(str.main, '/', int2str(iseed),'/kappa_mode_lifetime_omega.txt');
SED.omega_tau_gv_kappa = dlmread(str_write);

omega = SED.omega_tau_gv_kappa(:,1)*NMD.constant.w_LJ2THz;
tau = SED.omega_tau_gv_kappa(:,2)*NMD.constant.tau_LJ2ps;
figure;
loglog(omega,tau,"r*");%,'MarkerSize',10,'MarkerFaceColor','g','MarkerEdgeColor','g');

%% LD vs SED frquency

f_LD = NMD.freq(1,:)'*NMD.constant.w_LJ2THz;
f_SED = SED.omega_tau_gv_kappa(1:NMD.NUM_ATOMS_UCELL*3,1)*NMD.constant.w_LJ2THz;

figure;
plot(1:size(f_LD),f_LD,'.k','MarkerSize',8);
hold on;
plot(1:size(f_SED),f_SED,'or','MarkerSize',10);
legend('LD','SED');
xlabel('mode number');
ylabel('$\omega$ (THz)','Interpreter','latex');

figure;
f_diff = f_SED - f_LD;% Because of anharmonicity the frequency of SED is decreased
avg_diff = mean(f_diff);
plot(1:size(f_diff),f_diff,'om','MarkerSize',8);
yline(avg_diff,'-k');
xlabel('mode number');
ylabel('$\Delta \omega$ (THz)','Interpreter','latex');