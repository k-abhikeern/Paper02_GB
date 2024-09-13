close all;
constant = m_constant;

[tmp,str.main]=system('pwd');

iseed = 1;
ikslice = 1;

SED = load('./SEDomega_tau_gv_kappa.mat');
NMD = load('./NMD.mat');
size(SED.omega_tau_gv_kappa,1)
y1 = isnan(SED.omega_tau_gv_kappa(:,6));
z1 = find(y1 == 0);

kappax = sum(SED.omega_tau_gv_kappa(z1,6));
kappay = sum(SED.omega_tau_gv_kappa(z1,7));
kappaz = sum(SED.omega_tau_gv_kappa(z1,8));
factor = (NMD.constant.kb/(NMD.LJ.sigma*NMD.LJ.tau));%factor units will come out to be W/mK
[kappax kappay kappaz] % in W/mK


%% Open file in 1

fid = fopen('./1/irrkpt_symm_point.txt','r');

while ~feof(fid)    
    tline = fgetl(fid);
    index = str2num(tline);
    disp(index);
    batch_symm_kpoints = size(index,2);
    
% Multiply k_x,k_y,k_z by no. of similar points
    [imode,ind] = find(SED.omega_tau_gv_kappa(:,9)==index(1));
    kappa_x_estimated = SED.omega_tau_gv_kappa(imode,6)*batch_symm_kpoints;
    kappa_y_estimated = SED.omega_tau_gv_kappa(imode,7)*batch_symm_kpoints;
    kappa_z_estimated = SED.omega_tau_gv_kappa(imode,8)*batch_symm_kpoints;

    imode_kx_ky_kz(imode,:) = [imode,kappa_x_estimated,kappa_y_estimated,kappa_z_estimated];
end
fclose(fid);

y1 = isnan(imode_kx_ky_kz(:,2));
z1 = find(y1 == 0);

kappax_estm = sum(imode_kx_ky_kz(z1,2));
kappay_estm = sum(imode_kx_ky_kz(z1,3));
kappaz_estm = sum(imode_kx_ky_kz(z1,4));
[kappax_estm,kappay_estm,kappaz_estm]


f = fopen('kappa.txt', 'wt');      % 'wt' - write file, text mode
%formspec = '%f  %f\n';            % two values in a row (\n - line break)
formspec = '%f  %f  %f %s \n';      % three values in a row
fprintf(f, formspec, kappax_estm,kappay_estm,kappaz_estm,"W/mK");
fclose(f);

