clc;
close all;
clear all;
constant = m_constant;

[tmp,str.main]=system('pwd');

iseed = 1;
ikslice = 1;

SED = load('./SEDomega_tau_gv_kappa.mat');
NMD = load('./NMD.mat');

y1 = isnan(SED.omega_tau_gv_kappa(:,6));
z1 = find(y1 == 0);

kappax = sum(SED.omega_tau_gv_kappa(z1,6));
kappay = sum(SED.omega_tau_gv_kappa(z1,7));
kappaz = sum(SED.omega_tau_gv_kappa(z1,8));
factor = (NMD.constant.kb/(NMD.LJ.sigma*NMD.LJ.tau));%factor units will come out to be W/mK --> checked
[kappax kappay kappaz] % in W/mK

f = fopen('kappa.txt', 'wt');      % 'wt' - write file, text mode
formspec = '%f  %f  %f %s \n';      % three values in a row
fprintf(f, formspec, kappax,kappay,kappaz,"W/mK");
fclose(f);

%[kappax kappay kappaz]*factor %i W/mK
% for i = 1:size(NMD.kptlist(:,1:3,ikslice),1)
%     inti = (ikpt-1)*NMD.NUM_MODES;
%     intf = inti + NMD.NUM_MODES;
%     omega = SED.omega_tau_gv_kappa(inti:intf,1);
%     mfp = SED.omega_tau_gv_kappa(inti:intf,2).*SED.omega_tau_gv_kappa(inti:intf,3);
% end

%###########################################################
w_LJ2cm_inv = NMD.constant.w_LJ2THz*33.356;
%###########################################################

for i = 1:size(SED.omega_tau_gv_kappa,1)%size(NMD.kptlist(:,1:3,ikslice),1)*NMD.NUM_MODES % size = (no. of kpoints * no. of modes)
    %###########################################################
    omega(i) = SED.omega_tau_gv_kappa(i,1); %--> LJ
    %###########################################################
    velg = sqrt(sum(SED.omega_tau_gv_kappa(i,3:5).^2)); %--> LJ
    %###########################################################
    mfp(i) = velg*SED.omega_tau_gv_kappa(i,2); %--> LJ
    %###########################################################
    tau(i) = SED.omega_tau_gv_kappa(i,2); %--> LJ
    %###########################################################
end

%####################################################################################
% omega = omega*w_LJ2cm_inv;
% tau = tau*NMD.constant.tau_LJ2ps;
% %%
% xt = get(gca, 'XTick');
% yl = get(gca, 'YLim');
% % str = cellstr( num2str(xt(:),'2^{%d}') );      %# format x-ticks as 2^{xx}
% str.xticks = cellstr( num2str(2.^xt(:)) );      %# format x-ticks as 2^{xx}
% hTxt = text(xt, yl(ones(size(xt))), str.xticks, ...   %# create text at same locations
%     'Interpreter','tex', ...                   %# specify tex interpreter
%     'VerticalAlignment','top', ...             %# v-align to be underneath
%     'HorizontalAlignment','center');           %# h-aligh to be centered

%%
% xlabel('omega'); 
% ylabel('tau');
%% ####################################################################################

[mfp_sort I] = sort(mfp(z1));
mfp_sort = mfp_sort*NMD.LJ.sigma*1E9; %in nm
%kappax_sort = SED.omega_tau_gv_kappa(I,6); % in LJ units
kappax_sort = SED.omega_tau_gv_kappa(I,6)*factor; % in W/mK
cum_kappax_sort = cumsum(kappax_sort);
percent_cum_kappax_sort =100*cum_kappax_sort/max(cum_kappax_sort);
%#####################################################################################
figure;
semilogx(mfp_sort,percent_cum_kappax_sort,'-*g');
xlabel('phonon mfp(nm)'); ylabel('accumulated thermal conductivity(%)')'
saveas(gcf,'percent_kappa','jpg');
%#####################################################################################

str_write=strcat(str.main, '/', int2str(iseed),'/accumulated_kappa.txt');
if exist(['./' int2str(iseed) '/accumulated_kappa.txt'], 'file')~=0
    system(['rm -f ./' int2str(iseed) '/accumulated_kappa.txt']);
end

dlmwrite(str_write,[mfp_sort'  percent_cum_kappax_sort],'-append','delimiter',' ');

