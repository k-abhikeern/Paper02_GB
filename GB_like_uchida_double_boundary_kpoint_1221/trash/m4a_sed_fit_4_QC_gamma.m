clc;
close all;
clear all;
constant = m_constant;

[tmp,str.main]=system('pwd');

iseed = 1;
ikslice = 1;


SED = load('./SEDavg.mat');
NMD = load('./NMDavg.mat');
r = NMD.NUM_MODES;
c = 3;
NMD.vel = zeros(r,c);
% NMD = load('./NMD.mat');

PT_PERC = 0.01;
c01 = linspace(1,15,1);%(10,30,1)
c02 = linspace(0.001,0.08,1);%(0.08,0.10,1)
% PT_PERC = 0.01;
% c01 = linspace(1,30,1);%(10,30,1)
% c02 = linspace(0.001,0.05,1);%(0.08,0.10,1)
%% Loop_SED_fig
if exist(['./Loop_SED_fig'], 'file')~=0
    system(['rm -r ' ' Loop_SED_fig']);
end

str.cmd = [ 'mkdir Loop_SED_fig'];
system(str.cmd);

fig_seed = 'Loop_SED_fig/'; % remove after use

%% All peaks for a particular k-point -- expected to be in increasing order
% for ikpt = 574:574%1:size(NMD.kptlist(:,1:3,ikslice),1)
%     for i = 1:NMD.NUM_MODES
%             plot(SED.omega(1,:),SED.SEDpoi_kpt_avg(ikpt,:,i));
%             hold on;
%             name = strcat(NMD.str.main,fig_seed,'SED','_',int2str(ikpt),'_',int2str(i)); % will create FIG1, FIG2,...
%             saveas(gcf,name,'jpg');
%     end
% end
% %% -------------------------------
VOLUME =(NMD.Nx*NMD.Ny*NMD.Nz)*det(NMD.latvec)/((NMD.LJ.sigma*(1e+10))^3);% Ang^3 --> LJ (NMD.Nx*NMD.Ny*NMD.Nz)*det(NMD.latvec)/(NMD.LJ.sigma*(1e+10)^3)
Temp = 300; %in K
x_fac = NMD.constant.hbar/(NMD.constant.kb*Temp*NMD.LJ.tau); %hbar/(Kb*T)
factor = NMD.constant.kb/(NMD.LJ.sigma*NMD.LJ.tau*VOLUME);

str_write=strcat(str.main, '/', int2str(iseed),'/kappa_mode_lifetime_omega.txt');
if exist(['./' int2str(iseed) '/kappa_mode_lifetime_omega.txt'], 'file')~=0
    system(['rm -f ./' int2str(iseed) '/kappa_mode_lifetime_omega.txt']);
end

str_write_omega_tau=strcat(str.main, '/', int2str(iseed),'/omega_tau.txt');
if exist(['./' int2str(iseed) '/omega_tau.txt'], 'file')~=0
    system(['rm -f ./' int2str(iseed) '/omega_tau.txt']);
end

%########################################################
w_LJ2THz = (10^(-12))/(2*pi*NMD.LJ.tau);
% str_read = strcat('./1/irrkpt_symm_point.txt');
% kpoints = dlmread(str_read);
% NMD.kpt.kpt_poi = kpoints(:,1);
%#########################################################
cnt = 1;
for ikpt = 1:size(NMD.kptlist(:,1:3,ikslice),1)%size(NMD.kpt.kpt_poi,1)%size(NMD.kptlist(:,1:3,ikslice),1) size(NMD.kpt.symm_points_index,1)
    
    %     ikpt = NMD.kpt.kpt_poi(ikpt)% NMD.kpt.symm_points_index(ikpt)
    tic
    
    %     PT_PERC_LEFT = 0.1; PT_PERC_RIGHT = 0.1;
    %     PT_PERC = 0.01;
    INV_PERC = 1.0;
    %     groupvel_i = (ikpt-1)*NMD.NUM_ATOMS_UCELL*3;
    %
    
    groupvel_i = (ikpt-1)*NMD.gamma.NUM_ATOMS_UCELL*3;
    
    for i = 1:NMD.NUM_MODES
        
        
        
        
        %         % For not counting firt 3 modes of kpt = [0,0,0]
        % %         ikpt
        % %         i
        %         if (ikpt == 200) && (i ==1 || i==2 || i ==3)
        %                 continue;
        %         end
        
        
        
        
        
        
        start =1;
        %         [Imax,Jmax] = max(SED.SEDpoi_kpt_avg(ikpt,:,i));
        [Imax,Jmax] = max(SED.SED(ikpt,:,i));
        %--
        if Jmax < (2) || Jmax > (NMD.NUM_TSTEPS/2 - 1)
            i = i + 1;
            continue
        end
        
        
        
        %         %Find wleft
        %         [I,J] = find(SED.SEDpoi_kpt_avg(ikpt,start:start+Jmax,i) <...
        %             PT_PERC_LEFT*SED.SEDpoi_kpt_avg(ikpt,start+Jmax,i) );% find --> I willl be either 1 or 0. J will give cell no. of those values.
        %Find wleft
        %         [I,J] = find(SED.SEDpoi_kpt_avg(ikpt,start:Jmax,i) <...
        %             PT_PERC*SED.SEDpoi_kpt_avg(ikpt,Jmax,i) );% find --> I willl be either 1 or 0. J will give cell no. of those values.
        %
        
        %         [I,J] = find(SED.SED(ikpt,start:start+Jmax,i) <...
        %             PT_PERC_LEFT*SED.SED(ikpt,start+Jmax,i) );% find --> I willl be either 1 or 0. J will give cell no. of those values.
        %
        [I,J] = find(SED.SED(ikpt,start:start+Jmax,i) <...
            PT_PERC*SED.SED(ikpt,start+Jmax,i) );% find --> I willl be either 1 or 0. J will give cell no. of those values.
        
        %################################################################################
        %         wleft = Jmax - shift;
        
        %         w(:,1)=(1:length(SED.irrkpt.sedavg(:,kpt_list(kpt_cnt))));
        %         wleft = w(I(length(I)));
        %                         wleft = start+(length(I));
        wleft = (length(I));
        
        %################################################################################
        if wleft < 1
            i = i + 1;
            continue
            
        end
        %------NEW----------
        %Find wright
        %         [I,J] = find(SED.SEDpoi_kpt_avg(ikpt,start+Jmax:end,i) <...
        %             PT_PERC_RIGHT*SED.SEDpoi_kpt_avg(ikpt,start+Jmax,i) );
        %         [I,J] = find(SED.SEDpoi_kpt_avg(ikpt,start+Jmax:end,i) <...
        %             PT_PERC*SED.SEDpoi_kpt_avg(ikpt,Jmax,i) );
        
        %         [I,J] = find(SED.SED(ikpt,start+Jmax:end,i) <...
        %             PT_PERC_RIGHT*SED.SED(ikpt,start+Jmax,i) );
        [I,J] = find(SED.SED(ikpt,start+Jmax:end,i) <...
            PT_PERC*SED.SED(ikpt,start+Jmax,i) );
        %################################################################################
        %         wright = Jmax + shift;
        %         wright = w_guess(imode) + I(length(I));
        
        %                         wright = start+Jmax + (length(I));%+50;
        %         wright = length(SED.SEDpoi_kpt_avg(ikpt,:,i))-(length(I));%+50;
        wright = length(SED.SEDavg(ikpt,:,i))-(length(I));%+50;
        %Find wleft
        %                     [I,J] = find(SED.irrkpt.sedavg(1:w_guess(imode),kpt_list(kpt_cnt),imode)<PT_PERC*SED.irrkpt.sedavg(Jpeak,kpt_list(kpt_cnt),imode));
        %                     wleft = w(I(length(I)));
        %             %Find wright
        %                     [I,J] = find(SED.irrkpt.sedavg(w_guess(imode):length(w),kpt_list(kpt_cnt),imode)>PT_PERC*SED.irrkpt.sedavg(Jpeak,kpt_list(kpt_cnt),imode));
        %                     wright = w_guess(imode) + I(length(I));
        
        
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
        %% --
                %################################################################################
        %                                 plot(xdata,ydata,'ko',xdata,fun(x,xdata),'b-')
        shift = 50;
        kpt = strcat(int2str(ikpt),'/',...
            int2str(i),'/',int2str(shift)); % will create FIG1, FIG2,...
        %         kpt = num2str(NMD.kpt.cart_uvw(ikpt,:)); % giving u,v,w points as the 'legend' of kpoints in the plot %NMD.kpt.cart(ikpt,:)
        figure(1);
        
        
%         semilogy( SED.omega(1,:),SED.SEDavg(ikpt,:,i),'ob');
%         hold on;
%         semilogy(SED.omega(wleft:wright),lor_func([c_fit], SED.omega(wleft:wright)),'*-k','Markersize',1);
%         hold on;

%         semilogy( SED.omega(wleft:wright),SED.SEDavg(ikpt,wleft:wright,i),'ob');
%         hold on;
%         semilogy(SED.omega(wleft:wright),lor_func([c_fit], SED.omega(wleft:wright)),'*-k','Markersize',1);
%         hold off;
        
        plot( SED.omega(1,:),SED.SEDavg(ikpt,:,i),'ob');
        hold on;
        plot(SED.omega(wleft:wright),lor_func([c_fit], SED.omega(wleft:wright)),'*-k','Markersize',1);
        hold off;
              
        
%         xt = get(gca, 'XTick');
%         set(gca, 'XTick',xt, 'XTickLabel',xt*w_LJ2THz)
        
        xt = [0 10 20 30 40 50 60 70]./w_LJ2THz;
        set(gca, 'XTick', xt, 'XTickLabel', xt * w_LJ2THz);
%         xticklabels({'0','10','20','30','40','50','60','70'});
        
        %         xlabel('omega (LJ)','fontsize',12);
        xlabel('omega','fontsize',12);
        ylabel('SED(LJ)','fontsize',12);
        grid on;
        w_0= num2str(c_fit(3)*w_LJ2THz);
        xline(c_fit(3),'--k',{w_0});
        
        w_NMD= num2str(NMD.freq(i)*w_LJ2THz);
        xline(NMD.freq(i),'--r',{w_NMD});
        
        name = strcat(NMD.str.main,fig_seed,'SED','_',int2str(ikpt),'_',...
            int2str(i),'_',int2str(shift)); % will create FIG1, FIG2,...
        legend off;
% %         saveas(gcf,name,'jpg');
%         
        %################################################################################
        %% --
        
        %Store separate liftimes and frequencies for single and MULTIPLE FITS
        center=c_fit(3); lifetime=1/(2*c_fit(2));
        
        SED_FIT.freq(ikpt,i) = center;
        SED_FIT.life(ikpt,i) = lifetime;
        
        
        % For not counting firt 3 modes of kpt = [0,0,0]
        if (ikpt == 200) && (i ==1 || i==2 || i ==3)
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
        output_life = [SED_FIT.freq(ikpt,i),  SED_FIT.life(ikpt,i), ...
            groupvel(1), groupvel(2),  groupvel(3), SED.kappa(ikpt,i,1), ...
            SED.kappa(ikpt,i,2), SED.kappa(ikpt,i,3),ikpt, i];
        dlmwrite(str_write,output_life,'-append','delimiter',' ');
        
        
        %         output_omega_tau = [ikpt, i, SED_FIT.freq(ikpt,i),  SED_FIT.life(ikpt,i) ];
        %         dlmwrite(str_write_omega_tau,output_omega_tau,'-append','delimiter',' ');
    end
    
    ikpt
    
    toc
    
end %end of kpt

%%
SED.omega_tau_gv_kappa = dlmread(str_write);
save(strcat(str.main,'/SEDomega_tau_gv_kappa.mat'), '-struct', 'SED');

% Comparing the frequencies between GULP(NMD.freq) & SED() 
% f_NMD(:,1)=NMD.freq.*w_LJ2THz;
% f_SED(:,1)=SED.omega_tau_gv_kappa(:,1).*w_LJ2THz;
% f_SED(744,1)=0; % As NMD.freq is 744 and SED.omega_tau_gv_kappa is 743
% f_NMD_SED(:,1)=f_NMD;
% f_NMD_SED(:,2)=f_SED;
% f_NMD_SED=f_NMD_SED.*w_LJ2THz;
%
% SED.omega_tau = dlmread(str_write_omega_tau);
% save(strcat(str.main,'/SEDomega_tau.mat'), '-struct', 'SED');
%%
%
% % close all;
% w_LJ2THz = (10^(-12))/(2*pi*NMD.LJ.tau);
% tau_LJ2ps = NMD.LJ.tau*(10^(12));
%
% % SED.omega_tau_gv_kappa(1:3,2) = 0.0;
%
% figure;
% mark = ["ro","g+","b*","cs","md","k^"];
% for i = 1:NMD.NUM_MODES
%     index = find(SED.omega_tau_gv_kappa(:,10) == i );
%     omega = SED.omega_tau_gv_kappa(index,1)*w_LJ2THz;
%     tau = SED.omega_tau_gv_kappa(index,2)*tau_LJ2ps;
%     loglog(omega,tau,mark(i));%,'MarkerSize',10,'MarkerFaceColor','g','MarkerEdgeColor','g');
% % plot(omega,tau,mark(i));%,'MarkerSize',10,'MarkerFaceColor','g','MarkerEdgeColor','g');
%     hold on;
% end
% hold off;
% legend("ZA","TA","LA","ZO","TO","LO");
% xlabel('omega(THz)','fontsize',12);
% ylabel('tau(ps)','fontsize',12);

%%
% close all;

str_write=strcat(str.main, '/', int2str(iseed),'/kappa_mode_lifetime_omega.txt');
SED.omega_tau_gv_kappa = dlmread(str_write);
w_LJ2THz = (10^(-12))/(2*pi*NMD.LJ.tau);
tau_LJ2ps = NMD.LJ.tau*(10^(12));

% index = find(SED.omega_tau_gv_kappa(:,10) == i );
omega = SED.omega_tau_gv_kappa(:,1)*w_LJ2THz;
tau = SED.omega_tau_gv_kappa(:,2)*tau_LJ2ps;
loglog(omega,tau,"r*");%,'MarkerSize',10,'MarkerFaceColor','g','MarkerEdgeColor','g');

%% LD vs SED frquency

f_LD = NMD.freq(4:end)'*w_LJ2THz;
f_SED = SED.omega_tau_gv_kappa(:,1)*w_LJ2THz;

figure;
plot(1:size(f_LD),f_LD,'.k','MarkerSize',8);
hold on;
plot(1:size(f_SED),f_SED,'or','MarkerSize',10);

figure;
f_diff = f_SED - f_LD;% Because of anharmonicity the frequency of SED is decreased
avg_diff = mean(f_diff);
plot(1:size(f_diff),f_diff,'om','MarkerSize',8);
yline(avg_diff,'-k');