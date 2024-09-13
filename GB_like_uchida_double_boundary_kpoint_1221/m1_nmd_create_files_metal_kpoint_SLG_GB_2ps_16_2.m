%--------GRAPHENE+BRENNER---------------------------
clc;
clear all;
close all;
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%INPUT
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
NMD.LJ.eps = 0.00296*1.6022e-19; %   *4.13041118e-21; %4.55018144e-22 ; % J       %(TABLE II)0.02578 eV --> 0.00284*1.6022e-19 J
NMD.LJ.sigma = 3.34e-10;      % m       % (TABLE II) 3.40 Angs --> 3.40e-10 m
NMD.LJ.mass =  12.01*1.66054e-27; %  1.9944235e-26;  % kg      % 12.01 amu * 6.6326e-26 kg
NMD.LJ.tau = sqrt((NMD.LJ.mass*(NMD.LJ.sigma^2))/NMD.LJ.eps); % s
NMD.constant.kb = 1.3806e-23;  % J/K
NMD.constant.hbar = 1.054e-34; % Js
NMD.constant.i = sqrt(-1);
NMD.constant.c = 29979245800.00019; % cm/s
NMD.constant.ps2s = 1e-12;
NMD.constant.ang2m = 1e-10;
NMD.constant.eV2J = 1.60217646e-19;

NMD.constant.w_LJ2THz = (10^(-12))/(2*pi*NMD.LJ.tau);
NMD.constant.LJ2km_s_inv = (NMD.LJ.sigma/NMD.LJ.tau)*1e-03;
NMD.constant.tau_LJ2ps = NMD.LJ.tau*(10^(12));
%% Files 
NMD.str.main =  '/home/abhikeern/research2/graphene_with_GB/GB_like_uchida_double_boundary_kpoint_1221/';%pwd + "/" ;
NMD.gulp.path = '/usr/local/bin/gulp';
NMD.gulp.name = 'gulp_disp_lj_conv_num.tmp';
NMD.gulp.freq_name = 'gulp_lj_conv_num.tmp';

% Prepare the .tmp file from the following "del" files
NMD.gulp.name_del = 'gulp_disp_lj_conv_num_del.tmp';
NMD.gulp.freq_name_del = 'gulp_lj_conv_num_del.tmp';
NMD.lmp.name_del = 'lmp_SLG_read_replace.tmp';

%% First run the LAMMPS for getting the tmp_dump_NEMD.0.lammpstrj file (recangular shaped) and then extract the coordinates
NMD.Nx = 3;NMD.Ny = 5;NMD.Nz = 1; % Considering only 1 unit-cell 
NMD.cvt.N = 1;%  (n,m) = (0,2) , (n,m) = (1,2) , (n,m) = (2,3) 
NMD.cvt.M = 2; 
NMD.N = [NMD.cvt.N NMD.cvt.M
         NMD.cvt.M NMD.cvt.N]; % N1=[1, 2], N2=[2, 1],

bond_length = 1.3968418;%NMD.lmp.a/sqrt(3);
NMD.lmp.a = bond_length*sqrt(3);%2.46;% NMD.lmp.a*sqrt() lattice vector length = 2.46 for REBO potential


% Inputs to CVT code
x_periodic = 60;

% GB length on left and right --> which should be almost equal
NMD.l1 = bond_length*sqrt(3*(NMD.N(1,1)^2+NMD.N(1,1)*NMD.N(1,2)+NMD.N(1,2)^2)); %left GB-length
NMD.l2 = bond_length*sqrt(3*(NMD.N(2,1)^2+NMD.N(2,1)*NMD.N(2,2)+NMD.N(2,2)^2)); %right GB-length
y_periodic = min(NMD.l1,NMD.l2);
z_periodic = 3.34;

NMD.m(1) = 1.0; NMD.m(2) = 1.0; NMD.m(3) = 1.0; NMD.m(4) = 1.0; NMD.ATOM_TYPE = [1];
NMD.mass_basis = 12.01;
NMD.mass = [12.01];
NMD.mass_gulp = 12.01;
NMD.seed.lj = 1;
NMD.seed.initial = 1:1;

%% Angle relation: theta_ophus = (60 + theta_uc/2)

% Angle of misorientation for left side of GB theta_1
% NMD.theta_1_op = atan((2*NMD.N(1,1) + NMD.N(1,2))/(sqrt(3)*NMD.N(1,2))); % theta of ophus paper
% NMD.theta_1 = -acos((NMD.cvt.N^2+4*NMD.cvt.N*NMD.cvt.M+NMD.cvt.M^2)/(2*(NMD.cvt.N^2+NMD.cvt.N*NMD.cvt.M+NMD.cvt.M^2)));
% NMD.theta_1_uc = -2*(NMD.theta_1_op - pi*(60/180) );
% Change of axis 
% e1 = [     1                   0                   0
%         cos(pi/3)           sin(pi/3)               0
%           0                   0                        6.500000000000000/(NMD.lmp.a)];     
% R1 = [cos(NMD.theta_1_uc/2) -sin(NMD.theta_1_uc/2) 0
%       sin(NMD.theta_1_uc/2) cos(NMD.theta_1_uc/2)  0
%       0                 0                1];
% e1p = R1*e1;

% NMD.vec = (NMD.lmp.a)*e1p; 


% NMD.vec = [NMD.l1 0 0
%            NMD.l1*cos(pi/3) NMD.l1*sin(pi/3) 0
%            0 0 6.5];
% % Angle of misorientation for right side of GB i.e. theta_2
% NMD.theta_2 = atan((2*NMD.N(2,1) + NMD.N(2,2))/(sqrt(3)*NMD.N(2,2)));
% 
% % Change of axis 
% e2 = [     1                   0                   0
%         cos(pi/3)           sin(pi/3)               0
%           0                   0                        6.500000000000000/(NMD.lmp.a)];
% % NMD.vec = (NMD.lmp.a)*e;      
% R2 = [cos(NMD.theta_2) -sin(NMD.theta_2) 0
%       sin(NMD.theta_2) cos(NMD.theta_2)  0
%             0               0               1];
% e2p = R2*e2;



% %NMD.vec = v1,v2,v3
% NMD.vec = (NMD.lmp.a)*[     1                   0                   0
%                         cos(pi/3)           sin(pi/3)               0
%                             0                   0                        6.500000000000000/(NMD.lmp.a)];

%% Unit cell lattice vectors
% NMD.vec = [NMD.Nx*x_periodic 0                  0
%            0                NMD.Ny*y_periodic   0
%            0                0                   6.5];

NMD.ucell.vec = [x_periodic 0                  0
           0                y_periodic   0
           0                0                   z_periodic];

%% Run CVT - Unit cell atom coordinates       
NMD.x0_cvt = f_dump_coordinates_cvt(x_periodic,y_periodic,bond_length,NMD.ATOM_TYPE,NMD.cvt.M,NMD.cvt.N);
NMD.NUM_ATOMS = size(NMD.x0_cvt(:,1),1);

figure;
plot(NMD.x0_cvt(:,3),NMD.x0_cvt(:,4),".b",'MarkerSize',8);% xq,yq
title('Rectangular Unit Cell from CVT');
hold on;
axis equal;
hold off;
%% Unit cell vectors 
% Required by GULP
% NMD.gamma.vec(1,:) = NMD.Nx * NMD.vec(1,:);
% NMD.gamma.vec(2,:) = NMD.Ny * NMD.vec(2,:);
% NMD.gamma.vec(3,:) = NMD.Nz * NMD.vec(3,:);
%% Minimize with LAMMPS - Stablize unit cell coordinates 
% Make the NMD.x0 more stable by mnimizing it by airebo potential by LAMMPS

NMD.x0_GULP = f_min_x0(NMD.ucell.vec,NMD.NUM_ATOMS,NMD.ATOM_TYPE(1),NMD.mass(1),NMD.x0_cvt);
%% Prepare the two GULP files
% GULP requires only: Unit cell (a)Vectors (b) Coordinates (c) KPOINT(which we calculate later)
% f_gamma_get_cartesian() function: It solves 3 purpose by getting (a)  vector (b) cartesian (c) Making GULP type coorinate structure

f_gamma_get_cartesian(NMD.gulp.name_del,NMD.gulp.name,NMD.ucell.vec,NMD.x0_GULP,NMD.ATOM_TYPE);
f_gamma_get_cartesian(NMD.gulp.freq_name_del,NMD.gulp.freq_name,NMD.ucell.vec,NMD.x0_GULP,NMD.ATOM_TYPE);

%% --------------------------------------------------------------------------

%SED PARAMETERS------------------------------------------------------------

%ISEED---------------------------------------------------------------------
NMD.NUM_SEEDS = size(NMD.seed.initial,2);


%--Creating the required no. of folders
for   iseed = 1:NMD.NUM_SEEDS
    
    %     if exist(['./' int2str(iseed)], 'file')~=0
    %         system(['rm -r ./' int2str(iseed)]);
    %     end
    
    system([ 'mkdir ' int2str(iseed)]);
    
end

%--------------------------------------------------------------------------

%---IKSLICE----------------------------------------------------------------
NMD.NUM_KSLICES = 1;
%---------------------------FURTHER CHANGES FROM HERE-----------------------------------------------

%TIMES---------------------------------------------------------------------
%############################################################################
NMD.t_total = 2^16; NMD.t_fft = 2^16; NMD.t_step = 2^2; NMD.dt = (2*1e-15)/NMD.LJ.tau; % lj_final_big_step LJ (0.5 fs = (0.5*1e-15)/NMD.LJ.tau)
% NMD.t_total = 2^20; NMD.t_fft = 2^20; NMD.t_step = 2^5; NMD.dt = (1e-15)/NMD.LJ.tau; % lj_final LJ (0.001 fs = (1.0*1e-15)/NMD.LJ.tau)

%###############################################################################
NMD.NUM_TSTEPS = NMD.t_fft/NMD.t_step;
%--------------------------------------------------------------------------

%IFFT----------------------------------------------------------------------
NMD.NUM_FFTS = NMD.t_total/NMD.t_fft;
%--------------------------------------------------------------------------

%FREQS---------------------------------------------------------------------
NMD.w_step = 2*pi/(NMD.t_fft*NMD.dt); % LJ everything
NMD.w_max = 2*pi/(NMD.t_step*NMD.dt*2); % its not in THz,because (2*pi) is involved here. THEORY : "(2s)*(dt)< 1/(max. freq)" --> (max. freq) = 1/(s*dt*2). But here 2*pi is extra multiplied
NMD.NUM_OMEGAS = NMD.t_fft/(2*NMD.t_step);
%--------------------------------------------------------------------------

% NMD.alat = It is the lattice constant.
% Physical Meaning: Which remains outside the bracket and inside the
% bracket is the fractional coordinate.
%######################################################################################
NMD.alat = NMD.lmp.a;       %300K
%######################################################################################
% alat = 5.341/3.4;       %30K
% alat = 5.370/3.4;       %40K
% alat = 5.401/3.4;       %50K
% alat = 5.436/3.4;       %60K
% alat = 5.476/3.4;       %70K
% alat = 5.527/3.4;       %80K
% NMD.mass_basis = 12.01;
% NMD.mass = [12.01];
% NMD.mass_gulp = 12.01;

% Physical meaning : 1st cell extreme positions in fractional form
%######################################################################################
NMD.cartesian = NMD.x0_cvt;
usr_latvec = (1/NMD.alat)*NMD.ucell.vec;
usr_basis = (1/NMD.alat)*NMD.cartesian(:,3:5);

NMD.NUM_ATOMS_UCELL_BLG = size(usr_basis,1);
%######################################################################################
dummy = [usr_latvec;usr_basis]; % LJ

%--
%Define box size and conventional cell lattice parameters
NMD.latvec(1,1) = dummy(1,1); NMD.latvec(1,2) = dummy(1,2);
NMD.latvec(1,3) = dummy(1,3);
NMD.latvec(2,1) = dummy(2,1); NMD.latvec(2,2) = dummy(2,2);
NMD.latvec(2,3) = dummy(2,3);
NMD.latvec(3,1) = dummy(3,1); NMD.latvec(3,2) = dummy(3,2);
NMD.latvec(3,3) = dummy(3,3);

NMD.latvec = NMD.alat*NMD.latvec; % metal
%% User Reciprocal latt vec has to be computed from the original unit cell lattice vector
b1 = (cross(usr_latvec(2,:),usr_latvec(3,:)))/(dot(usr_latvec(1,:),cross(usr_latvec(2,:),usr_latvec(3,:))));
b2 = (cross(usr_latvec(3,:),usr_latvec(1,:)))/(dot(usr_latvec(2,:),cross(usr_latvec(3,:),usr_latvec(1,:))));
b3 = (cross(usr_latvec(1,:),usr_latvec(2,:)))/(dot(usr_latvec(3,:),cross(usr_latvec(1,:),usr_latvec(2,:))));
NMD.latvec_rec = [b1;b2;b3];

%first 3 rows are the lattice vectors
NMD.x.direct = dummy(4:length(dummy),:);

NMD.x.cart = NMD.alat * NMD.x.direct ; % orthogonal system (x-y)

%--------------------------------------------------------------------------

% build supercell
N_cnt = 1;
for iNx = 0:NMD.Nx-1
    for iNy = 0:NMD.Ny-1
        for iNz = 0:NMD.Nz-1  % covering z-direction first then y and then x
            inti = (N_cnt-1)*size(NMD.x.direct,1)+1;
            intf = (N_cnt)*size(NMD.x.direct,1);
            NMD.x0( inti:intf ,1) = inti:intf;
            NMD.x0( inti:intf ,2) = NMD.ATOM_TYPE(1);

            
            NMD.x0( inti:intf ,3) = NMD.x.cart(:,1) + iNx * NMD.latvec(1,1) +...
                iNy*NMD.latvec(2,1) + iNz*NMD.latvec(3,1); %--> all T_x componenets
            NMD.x0( inti:intf ,4) = NMD.x.cart(:,2) + iNx * NMD.latvec(1,2) +...
                iNy*NMD.latvec(2,2) + iNz*NMD.latvec(3,2); %--> all T_y componenets
            NMD.x0( inti:intf ,5) = NMD.x.cart(:,3) + iNx * NMD.latvec(1,3) +...
                iNy*NMD.latvec(2,3) + iNz*NMD.latvec(3,3); %--> all T_z componenets
            N_cnt =N_cnt+1;
        end
    end
end

%% --
figure
plot(NMD.x0(:,3),NMD.x0(:,4),'ok','MarkerSize',5); % change this in an appropriate way for zig-zag and arm chair
xlabel('x');
ylabel('y');
title('Atom positions initially');
axis equal;
grid on;

X = NMD.alat*usr_latvec; % the componenets should be edge vectors and not

%% 1. a_rot is the the rotation matrix [R] in --> x = R.X
NMD.a_rot(1,1) = norm( X(1,:) ); % metal everything
NMD.a_rot(1,2) = 0;
NMD.a_rot(1,3) = 0;
NMD.a_rot(2,1) = dot( X(2,:),X(1,:)/norm(X(1,:)) );
NMD.a_rot(2,2) = norm( cross(X(1,:)./norm(X(1,:)),X(2,:)) );
NMD.a_rot(2,3) = 0;
NMD.a_rot(3,1) = dot( X(3,:),X(1,:)./norm(X(1,:)) );
NMD.a_rot(3,2) = (dot( X(2,:),X(3,:)) - NMD.a_rot(2,1)*NMD.a_rot(3,1))/NMD.a_rot(2,2);
NMD.a_rot(3,3) = sqrt( norm(X(3,:))^2 - NMD.a_rot(3,1)^2 - NMD.a_rot(3,2)^2);

%% 2. Volume of supercell box
%VOLUME = (NMD.Nx*NMD.Ny*NMD.Nz)*det(X)
Vol = dot(X(1,:), cross(X(2,:),X(3,:)) ); % metal
%% 3. Matrix = MAT

%usr_latvec = NMD.alat*usr_latvec; % making the fractional into distance units
MAT(1,:) =  cross(X(2,:),X(3,:));
MAT(2,:) =  cross(X(3,:),X(1,:));
MAT(3,:) =  cross(X(1,:),X(2,:));

%% 4. original fractional coordinate = usr_latvec

%% Result

NMD.Q = NMD.a_rot'*(1/Vol)*MAT; % metal
NMD.x0(:,5) = zeros(size(NMD.x0,1),1);
NMD.x0(:,3:5) = (NMD.Q*NMD.x0(:,3:5)')';  % metal , non-orthogonal system (parallelopiped directions or a1,a2 & a3 directions)

% --
figure;
plot(NMD.x0(:,3),NMD.x0(:,4),'or','MarkerSize',5); % change this in an appropriate way for zig-zag and arm chair
xlabel('x');
ylabel('y');
title('Atom positions finally');
axis equal;
grid on;
hold on;

% --
if exist(['./' int2str(NMD.seed.lj) '/x0.data'], 'file')~=0
    system(['rm -f ./' int2str(NMD.seed.lj) '/x0.data']);
end
str.write=strcat(NMD.str.main,int2str(NMD.seed.lj),'/x0.data');
dlmwrite(str.write,NMD.x0,'-append','delimiter',' ');
NMD.NUM_ATOMS = size(NMD.x0,1);

%% Making the lmp_metal.in.x0.1 file
%NMD.alat has already been considered in NMD.latvec_new

% Simulation box
% xlo = 0.0; ylo = 0.0; zlo = 0.0; % non-orthogonal system (parallelopiped directions or a1,a2 & a3 directions)
xhi = (NMD.Nx)*NMD.a_rot(1,1); % changes NMD.Nx-1
yhi = (NMD.Ny)*NMD.a_rot(2,2);
% zhi = (NMD.Nz)*NMD.a_rot(3,3);
xy  = (NMD.Ny)*NMD.a_rot(2,1);
xz  = (NMD.Nz)*NMD.a_rot(3,1);
yz  = (NMD.Nz)*NMD.a_rot(3,2);


%[xlo_bound, ylo_bound, zlo_bound ] = min(NMD.x0(:,3:5));
lo_bound = min(NMD.x0(:,3:5));
xlo = lo_bound(1) - min( [0.0 xy xz xy+xz] );
ylo = lo_bound(2) - min( [0.0 yz] );
% zlo = lo_bound(3);

zlo = -20.0;zhi = 20.0;% just for making the z value high
%-------------------------------------------------------
if exist(['./lmp_metal.in.x0.1'], 'file')~=0
    system(['rm -f ./lmp_metal.in.x0.1']);
end
fileID = fopen("lmp_metal.in.x0.1","w");
fprintf(fileID, "#\n%d\t atoms\n %d\t atom types\n\n\n",NMD.NUM_ATOMS, size(NMD.ATOM_TYPE,2));

fprintf(fileID, " %f %f\t xlo xhi\n %f %f\t ylo yhi\n %f %f\t zlo zhi\n %f %f %f\t xy xz yz\n\n\n", .....
    xlo,xhi,ylo,yhi,zlo,zhi,xy,xz,yz);

fprintf(fileID, "Masses\n\n 1 %f\n\n",NMD.mass(1));
fprintf(fileID, "Atoms\n\n");
fclose(fileID);
str.write=strcat(NMD.str.main,'lmp_metal.in.x0.1');
dlmwrite(str.write,NMD.x0,'-append','delimiter',' ');

%% kptlist limits

% Eq 9 (PAPER = "alternate_SED" in LD folder)
% Approach : if even no. of cells in x-direction then take even limits ,
%            else the odd limits. Silimarly for other directions.

%                           x               y               z
% NMD.kpt.cell_lim = 2*[    -(NMD.Nx-2)/2   -(NMD.Ny-2)/2   -(NMD.Nz-2)/2       % even lower limit
%     NMD.Nx/2        NMD.Ny/2        NMD.Nz/2            % even upper limit
%     -(NMD.Nx-1)/2   -(NMD.Ny-1)/2   -(NMD.Nz-1)/2       % odd lower limit
%     (NMD.Nx-1)/2   (NMD.Ny-1)/2   (NMD.Nz-1)/2      ];  % odd upper limit
% 
% if  mod( NMD.Nx , 2 ) == 0 % checking if it is even no. ?
%     x_low = NMD.kpt.cell_lim(1,1);
%     x_high = NMD.kpt.cell_lim(2,1);
% else
%     x_low = NMD.kpt.cell_lim(3,1);
%     x_high = NMD.kpt.cell_lim(4,1);
% end
% 
% if  mod( NMD.Ny , 2 ) == 0
%     y_low = NMD.kpt.cell_lim(1,2);
%     y_high = NMD.kpt.cell_lim(2,2);
% else
%     y_low = NMD.kpt.cell_lim(3,2);
%     y_high = NMD.kpt.cell_lim(4,2);
% end
% 
% if  mod( NMD.Nz , 2 ) == 0
%     z_low = NMD.kpt.cell_lim(1,3);
%     z_high = NMD.kpt.cell_lim(2,3);
% else
%     z_low = NMD.kpt.cell_lim(3,3);
%     z_high = NMD.kpt.cell_lim(4,3);
% end

NMD.kpt.cell_lim = [    -(NMD.Nx-2)/2   -(NMD.Ny-2)/2   -(NMD.Nz-2)/2       % even lower limit
    NMD.Nx/2        NMD.Ny/2        NMD.Nz/2            % even upper limit
    -(NMD.Nx-1)/2   -(NMD.Ny-1)/2   -(NMD.Nz-1)/2       % odd lower limit
    (NMD.Nx-1)/2   (NMD.Ny-1)/2   (NMD.Nz-1)/2      ];  % odd upper limit

if  mod( NMD.Nx , 2 ) == 0 % checking if it is even no. ?
    x_low = NMD.kpt.cell_lim(1,1);
    x_high = NMD.kpt.cell_lim(2,1);
else
    x_low = NMD.kpt.cell_lim(3,1);
    x_high = NMD.kpt.cell_lim(4,1);
end

if  mod( NMD.Ny , 2 ) == 0
    y_low = NMD.kpt.cell_lim(1,2);
    y_high = NMD.kpt.cell_lim(2,2);
else
    y_low = NMD.kpt.cell_lim(3,2);
    y_high = NMD.kpt.cell_lim(4,2);
end

if  mod( NMD.Nz , 2 ) == 0
    z_low = NMD.kpt.cell_lim(1,3);
    z_high = NMD.kpt.cell_lim(2,3);
else
    z_low = NMD.kpt.cell_lim(3,3);
    z_high = NMD.kpt.cell_lim(4,3);
end
%% User Reciprocal latt vec has to be computed from the original lattice vector

% NMD.gamma.usr_latvec(1,:) = NMD.Nx * usr_latvec(1,:);
% NMD.gamma.usr_latvec(2,:) = NMD.Ny * usr_latvec(2,:);
% NMD.gamma.usr_latvec(3,:) = NMD.Nz * usr_latvec(3,:);
% 
% b1 = (cross(NMD.gamma.usr_latvec(2,:),NMD.gamma.usr_latvec(3,:)))/(dot(NMD.gamma.usr_latvec(1,:),cross(NMD.gamma.usr_latvec(2,:),NMD.gamma.usr_latvec(3,:))));
% b2 = (cross(NMD.gamma.usr_latvec(3,:),NMD.gamma.usr_latvec(1,:)))/(dot(NMD.gamma.usr_latvec(2,:),cross(NMD.gamma.usr_latvec(3,:),NMD.gamma.usr_latvec(1,:))));
% b3 = (cross(NMD.gamma.usr_latvec(1,:),NMD.gamma.usr_latvec(2,:)))/(dot(NMD.gamma.usr_latvec(3,:),cross(NMD.gamma.usr_latvec(1,:),NMD.gamma.usr_latvec(2,:))));
% NMD.latvec_rec = [b1;b2;b3];


% the k-points are to be given in fractional form to GULP i.e. 2*pi/a(0.5 0.5 0.0) OR
% 2*pi/a(0.5 0.0 0.0) types --> so there is no need to worry about the
% factor pi/a . It is not required.
cnt=1;
for ix=x_low:1:x_high
    for iy=y_low:1:y_high
        for iz= z_low:1:z_high %(-(NMD.Nz/2)+1):1:(NMD.Nz/2)
            NMD.kpt.integer(cnt,:) = [ix iy iz];
            NMD.kpt.cart(cnt,1:3) = ix/NMD.Nx*NMD.latvec_rec(1,:) +...
                iy/NMD.Ny*NMD.latvec_rec(2,:) +...
                iz/NMD.Nz*NMD.latvec_rec(3,:) ; % double checked -- correct
            
            NMD.kpt.NUM_KPTS = cnt;
            cnt=cnt+1;
        end
    end
end

%% B-Z points covered initially

figure;
plot(NMD.kpt.cart(:,1),NMD.kpt.cart(:,2),'o','MarkerSize',10,'MarkerEdgeColor','b') % change this in an appropriate way for zig-zag and arm chair
xlabel('$k_x/(2*\pi/a$)','Interpreter','latex');
ylabel('$k_y/(2*\pi/a$)','Interpreter','latex');
title('BZ before cut');
grid on;
hold on;

%% User input required. This input is from the website where in we put the "vec" values into the website.
% Define corner points of BZ

% tol = 1e-3;
% corn= 1/(2*pi/NMD.alat)*[0.851371133786863, 1.4746491990188664, 0.0
%     1.7027692349853263, 0, 00.0
%     0.851371133786863, -1.4746491990188664, 00.0
%     -0.851371133786863, -1.4746491990188664, 00.0
%     -1.7027692349853263, 0, 00.0
%     -0.851371133786863, 1.4746491990188664, 00.0 ];  %K6
% corn(:,1) = corn(:,1)+tol*ones(6,1);
% 
% % Plot the corner points
% plot(corn(:,1),corn(:,2),'-k');
% hold on;
% plot([corn(6,1),corn(1,1)],[corn(6,2),corn(1,2)],'-k');
% hold on;
%% Equation of line
% cnt_new = 1;
% 
% for i = 1:NMD.kpt.NUM_KPTS
%     [in,on] = inpolygon(NMD.kpt.cart(i,1),NMD.kpt.cart(i,2),corn(:,1),corn(:,2));
%     if in == 1
%         NMD.kpt.cart_new(cnt_new,:) = NMD.kpt.cart(i,:);
%         cnt_new = cnt_new+1;
%     end
% end
% 
% % Plot the in & on points
% % plot(NMD.kpt.cart_new(:,1),NMD.kpt.cart_new(:,2),'d','MarkerSize',3,'MarkerFaceColor','r','MarkerEdgeColor','r') % change this in an appropriate way for zig-zag and arm chair
% plot(NMD.kpt.cart_new(:,1),NMD.kpt.cart_new(:,2),'o','MarkerSize',5,'MarkerEdgeColor','k') % change this in an appropriate way for zig-zag and arm chair
% xlabel('k_x/(2*pi/a)');
% ylabel('k_y/(2*pi/a)');
% title('Red BZ after cut');
% grid on;

%%

% NMD.kpt.NUM_KPTS = cnt_new-1;
% max_ky_coord = max(NMD.kpt.cart_new(:,2));
% min_ky_coord = min(NMD.kpt.cart_new(:,2));
% cnt_new = 1;
% for i = 1:NMD.kpt.NUM_KPTS
%     if(((NMD.kpt.cart_new(i,2) == max_ky_coord)||(NMD.kpt.cart_new(i,2) == min_ky_coord))&& NMD.kpt.cart_new(i,1) < 0)
%         continue;
%     end
%     NMD.kpt.cart_new2(cnt_new,:) = NMD.kpt.cart_new(i,:);
%     cnt_new = cnt_new+1;
% end


%#################################################
% NMD.kpt.cart=NMD.kpt.cart_new;
%#################################################

%% B-Z points covered finally

figure;
plot(NMD.kpt.cart(:,1),NMD.kpt.cart(:,2),'p','MarkerSize',3,'MarkerFaceColor','k','MarkerEdgeColor','k') % change this in an appropriate way for zig-zag and arm chair
xlabel('$k_x/(2*\pi/a$)','Interpreter','latex');
ylabel('$k_y/(2*\pi/a$)','Interpreter','latex');
title('BZ after cut');
grid on;
hold off;
name = strcat(NMD.str.main,'BZ_finally','_',int2str(NMD.Nx),'_',...
    int2str(NMD.Ny),'_',int2str(NMD.Nz)); % will create FIG1, FIG2,...
saveas(gcf,name,'jpg');


%% --------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%INPUT
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%GULP
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

NMD.kptmaster = bsxfun(@times,NMD.kpt.cart,[NMD.Nx NMD.Ny NMD.Nz]);
NMD.NUM_KPTS = size(NMD.kptmaster(:,1:3),1);
NMD.kptmaster_index = 1:NMD.NUM_KPTS;

%--------------------------------------------------------------------------
%This is required to get the kpt + dk to properly input to GULP
format long
%--------------------------------------------------------------------------

%Define number of atoms
% NMD.NUM_ATOMS =size( NMD.x0,1);
% NMD.NUM_ATOMS_UCELL = size(NMD.x.cart,1) ;
% NMD.NUM_UCELL_COPIES=NMD.NUM_ATOMS/NMD.NUM_ATOMS_UCELL;
% NMD.NUM_MODES = 3*NMD.NUM_ATOMS_UCELL;

NMD.NUM_ATOMS = size(NMD.x0,1);
NMD.NUM_ATOMS_UCELL = size(NMD.x.cart,1);
NMD.NUM_UCELL_COPIES = NMD.NUM_ATOMS/NMD.NUM_ATOMS_UCELL;
NMD.NUM_MODES = 3*NMD.NUM_ATOMS_UCELL;



%Define finite difference increment
dk = 10e-5;
%--------------------------------------------------------------------------
%-------------------------KPTS---------------------------------------------
%--------------------------------------------------------------------------

if exist(['./' int2str(NMD.seed.lj) '/eigvec.dat'], 'file')~=0
    system(['rm -f ./' int2str(NMD.seed.lj) '/eigvec.dat']);
    system(['rm -f ./' int2str(NMD.seed.lj) '/freq.dat']);
    system(['rm -f ./' int2str(NMD.seed.lj) '/vel.dat']);
end


% NMD.NUM_ATOMS_UCELL = NMD.Nx*NMD.Ny*NMD.Nz*NMD.NUM_ATOMS_UCELL;

for ikpt=1:size(NMD.kptmaster,1)
    tic
    kpt = NMD.kptmaster(ikpt,:)./[NMD.Nx NMD.Ny NMD.Nz];
    
    
    
    %--
    
    % we want [u v w] in the reduced system.[k_x k_y k_z] = [u v w]*[b_1x b_1y b_1z
    %                                                                   b_2x b_2y b_2z
    %                                                                   b_3x b_3y b_3z]
    % So [K] = [R]*[B]-->and we can get [K][B]^-1 = [R]
    
    %#####################################################################################
    
    kpt =  kpt*inv(NMD.latvec_rec); % so that it displays the value on command window
    disp(kpt);
    NMD.kpt.cart_uvw(ikpt,:) = kpt;%(always correct --double checked)
    
    %#####################################################################################
    
    eigvec =...
        m_gulp_lj_eig(kpt,NMD.NUM_ATOMS_UCELL,NMD.mass_gulp,...
        NMD.alat,NMD.str.main,NMD.gulp.path,NMD.gulp.name);
    freq =...
        m_gulp_lj_freq(kpt,NMD.NUM_ATOMS_UCELL,NMD.mass_gulp,...
        NMD.alat,NMD.LJ,NMD.str.main,NMD.gulp.path,NMD.gulp.freq_name); % LJ
    %vel =...
    %   gulp_lj_vel(kpt,NMD.NUM_ATOMS_UCELL,NMD.mass_gulp,...
    %   NMD.alat*NMD.LJ.sigma/NMD.constant.ang2m,NMD.LJ,NMD.str.main,NMD.gulp.path,NMD.gulp.freq_name);
    vel =...
        m_gulp_lj_groupvel(kpt,NMD.NUM_ATOMS_UCELL,NMD.mass_gulp,...
        NMD.alat,NMD.LJ,NMD.str.main,NMD.gulp.path,NMD.gulp.freq_name); % LJ
%     
    %Output formatted eigvec
    str.write=strcat(NMD.str.main,int2str(NMD.seed.lj),'/eigvec.dat');
    dlmwrite(str.write,eigvec,'-append','delimiter',' ');
    %Output formatted freqs
    str.write=strcat(NMD.str.main,int2str(NMD.seed.lj),'/freq.dat');
    dlmwrite(str.write,freq,'-append','delimiter',' ');
    %Output velocities
    str.write=strcat(NMD.str.main,int2str(NMD.seed.lj),'/vel.dat');
    dlmwrite(str.write,vel,'-append','delimiter',' ');
    toc
end

%#####################################################################################
figure
plot(NMD.kpt.cart_uvw(:,1),NMD.kpt.cart_uvw(:,2),'o','MarkerSize',8,'MarkerFaceColor','k','MarkerEdgeColor','k'); % change this in an appropriate way for zig-zag and arm chair
xlabel('u');
ylabel('v');
title('uvw of BZ after cut');
grid on;
name = strcat(NMD.str.main,'BZ_finally','_','u','_',...
    'v','_','w'); % will create FIG1, FIG2,...
saveas(gcf,name,'jpg');
%
%#####################################################################################

NMD.eigvec =...
    dlmread(strcat(NMD.str.main,int2str(NMD.seed.lj),'/eigvec.dat'));
NMD.freq =...
    dlmread(strcat(NMD.str.main,int2str(NMD.seed.lj),'/freq.dat'));
NMD.vel =...
    dlmread(strcat(NMD.str.main,int2str(NMD.seed.lj),'/vel.dat'));

%% Make the frequencu 0 if its negative
% NMD.freq(NMD.freq < 0) = 0;
%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%GULP
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------


%KPT LISTS-----------------------------------------------------------------

slice_length = size(NMD.kptmaster,1)/NMD.NUM_KSLICES;
% remainder_length = size(NMD.kptlist,1) - slice_length*(NMD.NUM_KSLICE-1);
for ikslice = 1:NMD.NUM_KSLICES
    NMD.kptlist(:,1:3,ikslice) =...
        NMD.kptmaster( (ikslice-1)*slice_length+1:(ikslice)*slice_length,1:3);
    NMD.kpt_index(:,ikslice) =...
        NMD.kptmaster_index( (ikslice-1)*slice_length+1:(ikslice)*slice_length);
end




%% All Group velocity together
if exist(['./figures/vg'], 'dir')~=0
    system(['rmdir ./figures/vg']);
end

str.cmd = 'mkdir -p ./figures/vg';
system(str.cmd);

figure;
NMD.grp_vel_all(:,1) = reshape((NMD.freq).', [], 1);
NMD.grp_vel_all(:,1) = NMD.grp_vel_all(:,1).*NMD.constant.w_LJ2THz;% LJ--> THz

NMD.grp_vel_all(:,2) = sqrt(NMD.vel(:,1).^2+NMD.vel(:,2).^2+NMD.vel(:,3).^2);
NMD.grp_vel_all(:,2) = NMD.grp_vel_all(:,2).*NMD.constant.LJ2km_s_inv;

plot(NMD.grp_vel_all(:,1),NMD.grp_vel_all(:,2),'ok','MarkerSize',5);
hold on;    
NMD.grp_vel_all_avg = mean(NMD.grp_vel_all(:,2),1); % averageRow = mean(matrix, 1);
yl = yline(NMD.grp_vel_all_avg,'--r','LineWidth',3);


xlabel('$\omega$ (THz)','Interpreter','latex');
ylabel('$v_g$ (km/s)','Interpreter','latex');

name = strcat(NMD.str.main,'figures/vg/','vg_all.fig'); % will create FIG1, FIG2,...
saveas(gcf,name);
hold off;

num_modes = NMD.NUM_ATOMS_UCELL*3;

for i=0:(NMD.Nx*NMD.Ny-1)       
    
    val = num_modes*i+1:num_modes*(i+1);
    
    figure;
    plot(NMD.grp_vel_all(val,1),NMD.grp_vel_all(val,2),'ok','MarkerSize',5);
    kpt_name = strcat(num2str(NMD.kpt.cart_uvw(i+1,1)),',',num2str(NMD.kpt.cart_uvw(i+1,2)),',',num2str(NMD.kpt.cart_uvw(i+1,3)));
    legend(kpt_name);
    
    name = strcat(NMD.str.main,'figures/vg/',kpt_name,'.fig'); % will create FIG1, FIG2,...
    saveas(gcf,name);
    close(gcf);
end

% setFigureProperties();
% save_fig('fig9h_all_tau_omega_1221');

%SAVE NMD structure--------------------------------------------------------
save(strcat(NMD.str.main,'NMD.mat'), '-struct', 'NMD');