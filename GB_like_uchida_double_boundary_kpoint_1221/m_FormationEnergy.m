clc;
clear all;
close all;

%% Constants
E_coh = -7.3949;%-7.96;
NMD.ATOM_TYPE = 1;
NMD.mass(1) = 12.01;

NMD.lmp.a = 2.46;% NMD.lmp.a*sqrt() lattice vector length = 2.46 for REBO potential
bond_length = NMD.lmp.a/sqrt(3);

%% Inputs to CVT code
NMD.cvt.N = 1;% (n,m) = (1,2) , (n,m) = (2,3)
NMD.cvt.M = 2;
NMD.N = [NMD.cvt.N NMD.cvt.M
    NMD.cvt.M NMD.cvt.N]; % N1=[1, 2], N2=[2, 1],

x_periodic_DFT = 24.8; % initial guess of most stable length from DFT
y_periodic = NMD.lmp.a*sqrt((NMD.N(1,1)^2+NMD.N(1,1)*NMD.N(1,2)+NMD.N(1,2)^2)); % y-length = GB length

n_points = 20; % no. of points in graph around "x_periodic_DFT"
diff = 0.3;% units = Angs --> difference of d between two points


%% Store the values
filename = 'FE.txt';
if exist(['./' filename], 'file')~=0
    system(['rm -f ./' filename]);
end

fidout = fopen(filename,'a');
fprintf(fidout,"#N \t X_per \t Y_per \t E_min \t E_pristine \t FE \n");

count = 1;
for i = -n_points/2:1:n_points/2
    x_periodic = x_periodic_DFT + i*diff;
    
    % Run CVT
    NMD.x0_cvt = f_dump_coordinates_cvt(x_periodic,y_periodic,bond_length,NMD.ATOM_TYPE,NMD.cvt.M,NMD.cvt.N);
    NMD.NUM_ATOMS = size(NMD.x0_cvt(:,1),1);
    
%     figure;
%     plot(NMD.x0_cvt(:,3),NMD.x0_cvt(:,4),".b",'MarkerSize',8);% xq,yq
%     axis equal;
    
    % Run LAMMPS
    % Unit cell vectors
    NMD.vec = [x_periodic       0                       0
        0          y_periodic                    0
        0                0                       6.5];
    
    NMD.ucell.vec(1,:) = NMD.vec(1,:);
    NMD.ucell.vec(2,:) = NMD.vec(2,:);
    NMD.ucell.vec(3,:) = NMD.vec(3,:);
    
    % Make the NMD.x0 more stable by minimizing it by airebo potential
    NMD.x0 = f_min_x0(NMD.ucell.vec,NMD.NUM_ATOMS,NMD.ATOM_TYPE(1),NMD.mass(1),NMD.x0_cvt);
    
    % Extract the start, penultimate and end energy
    str.cmd = strcat("grep -A 1 'Energy initial' log.lammps | awk 'NR==2' > energy_val.txt");
    system(str.cmd);
    
    % Store these values in cloumn format    
    E_val(count,:) = dlmread('energy_val.txt');    
    
    E_pristine = E_coh*NMD.NUM_ATOMS;    
    FE = E_val(count,3)-E_pristine;% Formation Energy
    FE = FE/(2*y_periodic);
    
     output_val(count,:) = [size(NMD.x0,1), x_periodic, ...
            y_periodic, E_val(count,3), E_pristine, FE];
        
    fprintf(fidout,"%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n",output_val(count,:));
    count = count+1;
end

fclose(fidout);
%% 
figure;
plot(output_val(:,2),output_val(:,6),'-ks','MarkerSize',8);
xlabel('two GB distance (\AA)','Interpreter','latex');
ylabel('Formation Energy (eV)','Interpreter','latex');
