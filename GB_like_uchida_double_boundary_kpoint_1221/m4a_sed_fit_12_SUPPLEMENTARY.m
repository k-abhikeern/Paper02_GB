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
for ikpt = 4%1:size(NMD.kpt.kpt_poi,1)

    ikpt = NMD.kpt.kpt_poi(ikpt)
    tic
    %     PT_PERC_LEFT = 0.1; PT_PERC_RIGHT = 0.1;
    INV_PERC = 1.0;
    groupvel_i = (ikpt-1)*NMD.NUM_ATOMS_UCELL*3;

    for i = 219%1:NMD.NUM_MODES


        % %--
        % PT_PERC = 0.01;
        %     c01 = 24;
        %     c02 = 0.035;
        % % %--
        PT_PERC = 0.01;
        c01 = 11;
        c02 = 0.03;
        %--


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

        semilogy( SED.omega(1,wleft:wright),SED.SEDpoi_kpt_avg(ikpt,wleft:wright,i),'pb','MarkerSize',12);
        hold on;
        semilogy(SED.omega(wleft:wright),lor_func([c_fit],SED.omega(wleft:wright)),'-k','LineWidth',3);
        hold off;
        % plot( SED.omega(1,wleft:wright),SED.SEDpoi_kpt_avg(ikpt,wleft:wright,i),'ob');
        % hold on;
        % plot(SED.omega(wleft:wright),lor_func([c_fit], SED.omega(wleft:wright)),'*-k','Markersize',1);
        % hold off;



        fq_left = floor(SED.omega(1,wleft));
        fq_right = ceil(SED.omega(1,wright));
        fq_mid = (fq_left+fq_right)/2;

        xticks([fq_mid]); % Set custom x-axis ticks

        formatted_text = [num2str(fq_mid*(NMD.constant.w_LJ2THz), '%.1f')];
        xticklabels(formatted_text);

        xlabel('$\omega$ (THz)','Interpreter','latex','fontsize',16);
        ylabel('SED (LJ)','Interpreter','latex','fontsize',16);
        w_0 = num2str(c_fit(3)*NMD.constant.w_LJ2THz);
        xline(c_fit(3),'--k');%,{w_0})

        % Define the position where you want to place the text
        x_pos = c_fit(3) + 0.2; % Adjust 'offset' to move the text slightly to the right
        % y_pos = ylim-20800; % Get current y-axis limits
        y_pos(2) = 500;
        text(x_pos, y_pos(2), '$\omega_0$', 'FontSize', 20, 'Interpreter', 'latex');
        str = {'SED Data','Lorentzian Fit'};
        lgd = legend(str);
        set(lgd, 'Interpreter', 'latex');
        ylim([50,25000]);
        setFigureProperties1(lgd);
        save_fig('figS02a');
        

    end %end of i
    disp(cnt_kpt);


    NMD.time_perkpt(cnt_kpt) = toc;
    cnt_kpt = cnt_kpt+1;
end %end of kpt

