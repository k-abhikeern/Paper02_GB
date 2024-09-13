 clc;
close all;
clear all;
constant = m_constant;

[tmp,str.main]=system('pwd');

iseed = 1;
ikslice = 1;

NMD = load('./NMD.mat');
SED = load('./SEDavg.mat');


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

%########################################################
str_read = strcat('./1/irrkpt_symm_point.txt');
kpoints = dlmread(str_read);
NMD.kpt.kpt_poi = kpoints(:,1);
%########################################################
cnt_kpt = 1;
cnt_all = 1;
for ikpt = 2%1:size(NMD.kpt.kpt_poi,1)
    
    ikpt = NMD.kpt.kpt_poi(ikpt)
    tic
    %     PT_PERC_LEFT = 0.1; PT_PERC_RIGHT = 0.1;
    INV_PERC = 1.0;
    groupvel_i = (ikpt-1)*NMD.NUM_ATOMS_UCELL*3;
    
    for i = 300:320%1:NMD.NUM_MODES
        
          
        %--
        PT_PERC = 0.01;
            c01 = 31;
            c02 = 0.06;
        % %--
        % PT_PERC = 0.01;
        %     c01 = 31;
        %     c02 = 0.06;
        % %--
            
            
        start =1;
        [Imax,Jmax] = max(SED.SEDpoi_kpt_avg(ikpt,:,i));
        
        %--
        if Jmax < (2) || Jmax > (NMD.NUM_TSTEPS/2 - 1)
            i = i + 1;
            continue
        end
        
        
        %Find wleft
        [I,J] = find(SED.SEDpoi_kpt_avg(ikpt,start:start+Jmax,i) <...
            PT_PERC*SED.SEDpoi_kpt_avg(ikpt,start+Jmax,i) );% find --> I willl be either 1 or 0. J will give cell no. of those values.
        
        
        
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
        [I,J] = find(SED.SEDpoi_kpt_avg(ikpt,start+Jmax:end,i) <...
            PT_PERC*SED.SEDpoi_kpt_avg(ikpt,start+Jmax,i) );
        
        %################################################################################
        %         wright = Jmax + shift;
        wright = length(SED.SEDpoi_kpt_avg(ikpt,:,i))-(length(I));%+50;
        
        
        %                 wright = start+Jmax + (length(I));%+50;
        
        %################################################################################
        
        if wright > NMD.NUM_TSTEPS/2
            i = i + 1;
            continue
        end
        
        %% Initial guess mode_guess_c01 & mode_guess_c02
        
        mode_guess_c01 = c01;
        mode_guess_c02 = c02;
        
        c0 = [ mode_guess_c01*Imax, mode_guess_c02, SED.omega(Jmax)];%start+
        
        lb(1:length(c0)) = 0.0; ub(1:3:length(c0)) = 100000*Imax;
        ub(2:3:length(c0)) = 1000*1e15;
        ub(3:3:length(c0)) = 1000*SED.omega(length(SED.omega));
        
        lor_func = @(c,w)(c(1))./(1 + ( (w - c(3))./ c(2) ).^2 );
        
        options = optimset('MaxIter',5000,'MaxFunEvals',5000,'TolFun',1e-5,'TolX',1e-5);
        
        [c_fit,resnorm,residual,exitflag,output]  = lsqcurvefit(lor_func,c0,SED.omega(wleft:wright),...
            SED.SEDpoi_kpt_avg(ikpt,wleft:wright,i),lb,ub,options);
        
        
        %% SED figures
        
                figure;

                semilogy( SED.omega(1,:),SED.SEDpoi_kpt_avg(ikpt,:,i),'ob');
                hold on;
                semilogy(SED.omega(wleft:wright),lor_func([c_fit],SED.omega(wleft:wright)),'*-k','Markersize',1);
                hold off;
                % plot( SED.omega(1,wleft:wright),SED.SEDpoi_kpt_avg(ikpt,wleft:wright,i),'ob');
                % hold on;
                % plot(SED.omega(wleft:wright),lor_func([c_fit], SED.omega(wleft:wright)),'*-k','Markersize',1);
                % hold off;

                xticksPT_PERC = 0.01;
            c01 = 35;
            c02 = 0.065;([]);
                fq_left = SED.omega(1,wleft);
                fq_right = SED.omega(1,wright);
                fq_mid = (fq_left+fq_right)/2;
                xticks([fq_left, fq_mid, fq_right]); % Set custom x-axis ticks
                xticklabels({num2str(fq_left*(NMD.constant.w_LJ2THz)), ...
                                        num2str(fq_mid*(NMD.constant.w_LJ2THz)), ...
                                            num2str(fq_right*(NMD.constant.w_LJ2THz))});

                xlabel('$\omega$(THz)','Interpreter','latex','fontsize',12);
                ylabel('SED(LJ)','Interpreter','latex','fontsize',12);
                w_0 = num2str(c_fit(3)*NMD.constant.w_LJ2THz);
                xline(c_fit(3),'--k',{w_0});
                legend off;
        %% Store data
        
        %Store separate liftimes and frequencies for single and MULTIPLE FITS
        center=c_fit(3); lifetime=1/(2*c_fit(2));
        
        SED_FIT.freq(ikpt,i) = center;
        SED_FIT.life(ikpt,i) = lifetime;
        
        % At Gamma : Frequency of LA,TA & ZA = 0
        if ((NMD.kpt.cart_uvw(ikpt,1) == 0 && ...
                NMD.kpt.cart_uvw(ikpt,2) == 0 && ...
                NMD.kpt.cart_uvw(ikpt,3) == 0) && ...
                (i == 1 || i == 2 || i == 3))
            flag = 1
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
        
        SED.omega_tau_gv_kappa(cnt_all,:) = output_life;
        cnt_all = cnt_all+1;
        disp(i);
        
    end %end of i
    disp(cnt_kpt);
    
    
    NMD.time_perkpt(cnt_kpt) = toc;
    cnt_kpt = cnt_kpt+1;
end %end of kpt

% %save(strcat(NMD.str.main,'NMD.mat'), '-struct', 'NMD');

%% Plots
% tau vs omega
clc;
close all;
NMD = load('./NMD.mat');
str_write=strcat('./1/kappa_mode_lifetime_omega.txt');
SED.omega_tau_gv_kappa = dlmread(str_write);

omega = SED.omega_tau_gv_kappa(:,1)*NMD.constant.w_LJ2THz;
tau = SED.omega_tau_gv_kappa(:,2)*NMD.constant.tau_LJ2ps;
figure;
loglog(omega,tau,"r*");
hold on;
NMD.tau_mean = mean(SED.omega_tau_gv_kappa(:,2))*NMD.constant.tau_LJ2ps;
yline(NMD.tau_mean,'--k');
xlabel('$\omega$ (THz)','Interpreter','latex');
ylabel('$\tau$ (ps)','Interpreter','latex');
hold off;

% Plot LD vs SED frequencies for the last kpt
last_kpt = (SED.omega_tau_gv_kappa(end,9));
f_LD = NMD.freq(last_kpt,:)'*NMD.constant.w_LJ2THz;

index = find(SED.omega_tau_gv_kappa(:,9)==last_kpt);
f_SED = SED.omega_tau_gv_kappa(index,1)*NMD.constant.w_LJ2THz;

figure;
plot(1:size(f_LD),f_LD,'.k','MarkerSize',8);
hold on;
plot(1:size(f_SED),f_SED,'.r','MarkerSize',10);
legend('LD','SED');
xlabel('mode number','Interpreter','latex');
ylabel('$\omega$ (THz)','Interpreter','latex');

% % Plot frequency difference of LD vs SED
% figure;
% f_diff = abs(f_SED - f_LD);% Because of anharmonicity the frequency of SED is decreased
% NMD.w_avg_diff = mean(f_diff);
% plot(1:size(f_diff),f_diff,'om','MarkerSize',8);
% yline(NMD.w_avg_diff,'--k');
% xlabel('mode number');
% ylabel('$\Delta \omega$ (THz)','Interpreter','latex');

% Plot cv
cv = SED.omega_tau_gv_kappa(index,11);
figure;
plot(1:size(cv),cv,'.r','MarkerSize',10);
xlabel('mode number','Interpreter','latex');
ylabel('cv','Interpreter','latex');

%% -----------------New patch from previous codes---------------------

% Bar chart
clc;
SED.omega_tau_gv_kappa = dlmread(str_write);
% Frequency from LJ to THz
row = size(SED.omega_tau_gv_kappa,1);
column = 10;
data = zeros(row,column);

% Converting all units from LJ --> required units
data(:,1) = SED.omega_tau_gv_kappa(:,1)*NMD.constant.w_LJ2THz; %THz
data(:,2) = SED.omega_tau_gv_kappa(:,2)*NMD.constant.tau_LJ2ps; %ps
data(:,3) = SED.omega_tau_gv_kappa(:,3)*NMD.constant.LJ2km_s_inv; %km/s
data(:,4) = SED.omega_tau_gv_kappa(:,4)*NMD.constant.LJ2km_s_inv; %km/s
data(:,5) = SED.omega_tau_gv_kappa(:,5)*NMD.constant.LJ2km_s_inv; %km/s
data(:,6) = SED.omega_tau_gv_kappa(:,6); % W/mK
data(:,7) = SED.omega_tau_gv_kappa(:,7); % W/mK
data(:,8) = SED.omega_tau_gv_kappa(:,8); % W/mK
data(:,9) = SED.omega_tau_gv_kappa(:,9); % kpoint
data(:,10) = SED.omega_tau_gv_kappa(:,10); % mode

data(:,11) = sqrt(data(:,3).^2).*data(:,2)*10; % Angs % mfpx
data(:,12) = sqrt(data(:,4).^2).*data(:,2)*10; % Angs % mfpy

% Sort according to the frequncy
sorted_data = sortrows(data, 1);


%% 1st Quadrant Group velocity together
figure(21);
NMD.grp_vel(:,1) = sorted_data(:,1);
NMD.grp_vel(:,2) = sqrt(sorted_data(:,3).^2+sorted_data(:,4).^2+sorted_data(:,5).^2);
plot(NMD.grp_vel(:,1),NMD.grp_vel(:,2),'ok','MarkerSize',5);
xlim([0,55]);

NMD.grp_vel_vgx_vgy_vgz_Q1 = [mean(sqrt(sorted_data(:,3).^2)),mean(sqrt(sorted_data(:,4).^2)),mean(sqrt(sorted_data(:,5).^2))];

NMD.grp_vel_avg = mean(NMD.grp_vel(:,2),1); % averageRow = mean(matrix, 1);
yl = yline(NMD.grp_vel_avg,'--r','LineWidth',3);

xticks([]);
new_xticks = [0, 25, 50]; % Define your new x-tick positions
xticks(new_xticks);

yticks([]);
new_xticks = [0,10,20,30,40,50]; % Define your new x-tick positions
yticks(new_xticks);

xlabel('$\omega$ (THz)','Interpreter','latex');
ylabel('$v_g$ (km/s)','Interpreter','latex');

setFigureProperties();
%save_fig('fig9g_SLG_wo_GB');

%% All Group velocity together
figure(31);

% Assuming your original 64x84 matrix is called 'originalMatrix'
NMD.grp_vel_all(:,1) = reshape((NMD.freq).', [], 1);
NMD.grp_vel_all(:,1) = NMD.grp_vel_all(:,1).*NMD.constant.w_LJ2THz;% LJ--> THz

NMD.grp_vel_all(:,2) = sqrt(NMD.vel(:,1).^2+NMD.vel(:,2).^2+NMD.vel(:,3).^2);
NMD.grp_vel_all(:,2) = NMD.grp_vel_all(:,2).*NMD.constant.LJ2km_s_inv;

plot(NMD.grp_vel_all(:,1),NMD.grp_vel_all(:,2),'ok','MarkerSize',5);
xlim([0,55]);

NMD.grp_vel_all_avg = mean(NMD.grp_vel_all(:,2),1); % averageRow = mean(matrix, 1);
yl = yline(NMD.grp_vel_all_avg,'--r','LineWidth',3);


xlabel('$\omega$ (THz)','Interpreter','latex');
ylabel('$v_g$ (km/s)','Interpreter','latex');

setFigureProperties();
%save_fig('fig9h_SLG_wo_GB');

%% MFP

% mfpx
count_mfpx = 0;
figure;
semilogy(sorted_data(:,1),sorted_data(:,11),'sk');
xlabel('\omega (THz)','Interpreter','latex');ylabel('MFP_x (A)','Interpreter','latex');
yline(30, '--r', 'Y = 30', 'LineWidth', 2);
hold on;
% mfpx mean
NMD.avg_mfpx = mean(sorted_data(:,11));
avg_mfpx_str = num2str(NMD.avg_mfpx);
yline(NMD.avg_mfpx,'--b',['<mfp-x> = ', avg_mfpx_str], 'LineWidth', 2);

for i = 1:size(sorted_data(:,11),1)
    if sorted_data(i,11) <= 30
        count_mfpx = count_mfpx +1;
    end
end
NMD.percent_mfpx = 100*count_mfpx/size(sorted_data(:,11),1);


%mfpy
count_mfpy = 0;
figure;
semilogy(sorted_data(:,1),sorted_data(:,12),'sk');
xlabel('\omega(THz)','Interpreter','latex');ylabel('MFP_y (A)','Interpreter','latex');
yline(30, '--r', 'Y = 30', 'LineWidth', 2);
hold on;
% mfpx mean
NMD.avg_mfpy = mean(sorted_data(:,12));
avg_mfpy_str = num2str(NMD.avg_mfpy);
yline(NMD.avg_mfpy,'--b',['<mfp-y> = ', avg_mfpy_str], 'LineWidth', 2);

for i = 1:size(sorted_data(:,12),1)
    if sorted_data(i,12) <= 30
        count_mfpy = count_mfpy +1;
    end
end
NMD.percent_mfpy = 100*count_mfpy/size(sorted_data(:,12),1);


% Write in a file
if exist(['./' int2str(iseed) '/omega_lx_ly.txt'], 'file')~=0
    system(['rm -f ./' int2str(iseed) '/omega_lx_ly.txt']);
end

str_write = strcat('./1/omega_lx_ly.txt');
output_mfp = [sorted_data(:,1),sorted_data(:,11),sorted_data(:,12)];
dlmwrite(str_write,output_mfp,'-append','delimiter',' ');
%%  TC of Each kpt of 1st Quadrant contributes
count = 1;
clear NMD.kappa.kpt_kx_ky_kz_vg_tau_1st;
fid = fopen('./1/irrkpt_symm_point.txt','r');
while ~feof(fid)    
    tline = fgetl(fid);
    index = str2num(tline);
    disp(index);
    batch_Q1m_kpoints = 1; %only 1st quadrant --> batch_Q1m_kpoints will be only 1
    
% Multiply k_x,k_y,k_z by no. of similar points
    [imode,ind] = find(SED.omega_tau_gv_kappa(:,9)==index(1));
    kpt = SED.omega_tau_gv_kappa(imode,9);
    kappa_x_estimated = sorted_data(imode,6)*batch_Q1m_kpoints;
    kappa_y_estimated = sorted_data(imode,7)*batch_Q1m_kpoints;
    kappa_z_estimated = sorted_data(imode,8)*batch_Q1m_kpoints;
    
    vgx_sq = sorted_data(imode,3).^2;
    vgy_sq = sorted_data(imode,4).^2;
    vgz_sq = sorted_data(imode,5).^2;    
    vgx = sqrt(vgx_sq)*batch_Q1m_kpoints;
    vgy = sqrt(vgy_sq)*batch_Q1m_kpoints;
    vgz = sqrt(vgz_sq)*batch_Q1m_kpoints;
    
    tau = sorted_data(imode,2)*batch_Q1m_kpoints;
% Each kpt of 1st Quadrant contributes 
    NMD.kappa.kpt_kx_ky_kz_vg_tau_1st(count,:) = [index(1),sum(kappa_x_estimated),sum(kappa_y_estimated),sum(kappa_z_estimated),...
                                                    sum(vgx),sum(vgy),sum(vgz),sum(tau)];
    count = count+1;
end
fclose(fid);

%
NMD.kappa.kappax_estm = sum(NMD.kappa.kpt_kx_ky_kz_vg_tau_1st(:,2));
NMD.kappa.kappay_estm = sum(NMD.kappa.kpt_kx_ky_kz_vg_tau_1st(:,3));
NMD.kappa.kappaz_estm = sum(NMD.kappa.kpt_kx_ky_kz_vg_tau_1st(:,4));
[NMD.kappa.kappax_estm,NMD.kappa.kappay_estm,NMD.kappa.kappaz_estm] % in W/mK

%%  TC of whole
%  All corresponding kpt similar to 1st Quadrant kpt contributes
count = 1;
clear NMD.kappa.kpt_kx_ky_kz_vg_tau_all;
fid = fopen('./1/irrkpt_symm_point.txt','r');
while ~feof(fid)    
    tline = fgetl(fid);
    index = str2num(tline);
    disp(index);
    batch_Q1m_kpoints = size(index,2);
    
% Multiply k_x,k_y,k_z by no. of similar points
    [imode,ind] = find(sorted_data(:,9)==index(1));
    kpt = sorted_data(imode,9);
    kappa_x_estimated = sorted_data(imode,6)*batch_Q1m_kpoints;
    kappa_y_estimated = sorted_data(imode,7)*batch_Q1m_kpoints;
    kappa_z_estimated = sorted_data(imode,8)*batch_Q1m_kpoints;

vgx_sq = sorted_data(imode,3).^2;
    vgy_sq = sorted_data(imode,4).^2;
    vgz_sq = sorted_data(imode,5).^2;    
    vgx = sqrt(vgx_sq)*batch_Q1m_kpoints;
    vgy = sqrt(vgy_sq)*batch_Q1m_kpoints;
    vgz = sqrt(vgz_sq)*batch_Q1m_kpoints;
    
    tau = sorted_data(imode,2)*batch_Q1m_kpoints;
% Each kpt of 1st Quadrant contributes 
    NMD.kappa.kpt_kx_ky_kz_vg_tau_all(count,:) = [index(1),sum(kappa_x_estimated),sum(kappa_y_estimated),sum(kappa_z_estimated),...
                                                    sum(vgx),sum(vgy),sum(vgz),sum(tau)];
    count = count+1;
end
fclose(fid);

%
NMD.kappa.kappax_estm = sum(NMD.kappa.kpt_kx_ky_kz_vg_tau_all(:,2));
NMD.kappa.kappay_estm = sum(NMD.kappa.kpt_kx_ky_kz_vg_tau_all(:,3));
NMD.kappa.kappaz_estm = sum(NMD.kappa.kpt_kx_ky_kz_vg_tau_all(:,4));
[NMD.kappa.kappax_estm,NMD.kappa.kappay_estm,NMD.kappa.kappaz_estm] % in W/mK
%% %save the structure
save(strcat('NMD.mat'), '-struct', 'NMD'); %saving NMD structure variables as a file with name "NMDavg.mat"
save(strcat('SEDomega_tau_gv_kappa.mat'), '-struct', 'SED');
