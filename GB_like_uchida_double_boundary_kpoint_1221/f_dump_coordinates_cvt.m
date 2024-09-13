function x0 = f_dump_coordinates_cvt(width,y_length,bond_length,atom_type,cvt_M,cvt_N)


% %% Running the NEMD LAMMPS SLG file.
% % The rectagular shaped structure is of 100X100 size which hopefully will cover our requirements of extracting required atoms.
% str.cmd = strcat("mpirun -np 64 lmp_mpi -in in.SLG_NEMD_big_cell");
% system(str.cmd);
%
% % Delete starting 9 lines and saving into new file
% newfile = "newfile.txt";
% str.cmd = strcat("tail -n +10 tmp_dump_NEMD.0.lammpstrj >"," ",newfile);
% system(str.cmd);
%
% % Read the new file and store its value in a variable
% fidin = fopen(newfile,"r");% tmp_dump_NEMD.*.lammpstrj
% dump_coord = textscan(fidin,"%d\t%d\t%f\t%f\t%f");
% fclose(fidin);
% dump_x = dump_coord{3}(:);
% dump_y = dump_coord{4}(:);
% dump_z = dump_coord{5}(:);
%
% figure;
% plot(dump_x,dump_y,"sb");% xq,yq
% axis equal;
% hold on;

%% CVT code
% Running the cvt code first and then extracting the 5 columns.The cell
% side has to be big so taht a GB can be extracted in the trapezoid

% Run CVT code and generate initial configuration file
% M = 2; % (m,n) = (1,1) is not AA-BLG
% N = 1;

loc_cvt = ('/home/abhikeern/Paper02_GB/cvt');
% % Inputs to CVT code
% width = 30;
% y_length = 20;
% a = 1.3968418;
L_cell = bond_length*sqrt(3*(cvt_N^2+cvt_M^2+cvt_M*cvt_N));
fac_Ly = y_length/L_cell;

str.cvt = strcat(loc_cvt,'/','periodic_2.py');
if exist(str.cvt, 'file')~=0
    system(['rm -rf ' str.cvt ]);
end

fid = fopen(str.cvt,'w');
fprintf(fid,"import grainBdr as gb \n");
fprintf(fid,"import polyCrystal as pc \n");
fprintf(fid,"import ase.io \n");
fprintf(fid,"import numpy \n");

% fprintf(fid,"cr = gb.onePeriodicGB(N1=[%s, %s], N2=[%s, %s],fac_ly=%s, cell_width=%s, verbose=False) \n",num2str(cvt_M),num2str(cvt_N),num2str(cvt_N),num2str(cvt_M),num2str(fac_Ly),num2str(width));
fprintf(fid,"cr = gb.twoPeriodicGB(%s,N1=[%s, %s], N2=[%s, %s],cell_width=%s, verbose=False) \n",num2str(bond_length),num2str(cvt_N),num2str(cvt_M),num2str(cvt_M),num2str(cvt_N),num2str(width));
fprintf(fid,"ase.io.write('GB_1Periodic.pdb', cr) \n");
fprintf(fid,"ase.io.write('GB_1Periodic.cfg', cr) \n");
fprintf(fid,"gb.writeLammpsData(cr, 'GB_1Periodic.lammps') \n");

fclose(fid);

% Remove the existing GB_1Periodic.lammps
if exist(['./GB_1Periodic.lammps'], 'file')~=0 
    system(['rm -rf ./GB_1Periodic.lammps' ]);
end

% Run the CVT code now which generates GB_1Periodic.lammps in this folder
str_run = strcat("python"," ",str.cvt);
system(str_run);


% str_cp = strcat("cp"," ",loc_cvt,'/','GB_1Periodic.lammps'," ",".");
% system(str_cp);
% str_mv = strcat("cp"," ",'GB_1Periodic.lammps'," ","lmp_metal.in.x0.1");
% system(str_mv);




% Delete starting 11 lines and saving into new file
newfile = "newfile.txt";
str.cmd = strcat("tail -n +12 GB_1Periodic.lammps >"," ",newfile);
system(str.cmd);

% Read the new file and store its value in a variable
fidin = fopen(newfile,"r");% tmp_dump_NEMD.*.lammpstrj
dump_coord = textscan(fidin,"%d\t%d\t%f\t%f\t%f");
fclose(fidin);
dump_x = dump_coord{3}(:);
dump_y = dump_coord{4}(:);
% dump_z = dump_coord{5}(:);
dump_z = zeros(size(dump_y,1),1);


%% vec from previous matlab

% vec_supercell(1,:) = NMD_Nx * vec(1,:);
% vec_supercell(2,:) = NMD_Ny * vec(2,:);
% vec_supercell(3,:) = NMD_Nz * vec(3,:);

% corn = [0.0                                     0.0
%     vec_supercell(1,1)                      vec_supercell(1,2)
%     vec_supercell(1,1)+vec_supercell(2,1)   vec_supercell(1,2)+vec_supercell(2,2)
%     vec_supercell(2,1)                      vec_supercell(2,2)];%go ant-clockwise

% corn = [0.0                                     0.0
%     vec_supercell(1,1)                      vec_supercell(1,2)
%     vec_supercell(1,1)+vec_supercell(2,1)   vec_supercell(1,2)+vec_supercell(2,2)
%     vec_supercell(2,1)                      vec_supercell(2,2)];%go ant-clockwise

% corn = [0.0                                     0.0
%     vec_supercell(1,1)                          0.0
%     vec_supercell(1,1)+vec_supercell(2,1)   vec_supercell(2,2)
%     vec_supercell(2,1)                      vec_supercell(2,2)];%go ant-clockwise

%% Shifting the trapezoid towards left and then down in order to accomodate the required atoms inside the the trapezoid
% tol = 0;%1e-3;
% corn(:,1) = corn(:,1)-tol*ones(4,1);
% corn(:,2) = corn(:,2)-tol*ones(4,1);
% 
% plot(corn(:,1),corn(:,2),"-k");
% hold on;
% plot([corn(1,1),corn(4,1)],[corn(1,2),corn(4,2)],"-k");

%% Make a polyhedra and extract the required atoms inside the shape

% cnt_new = 1;
% 
% for i = 1:size(dump_x,1)
%     [in,on] = inpolygon(dump_x(i,1),dump_y(i,1),corn(:,1),corn(:,2));
%     if in == 1
%         cart_new(cnt_new,1) = dump_x(i,1);
%         cart_new(cnt_new,2) = dump_y(i,1);
%         cnt_new = cnt_new+1;
%     end
% end
% 
% plot(cart_new(:,1),cart_new(:,2),"og","MarkerEdgeColor","g","MarkerSize",8);
% hold off;
%% Plot only the extarcted part now
% figure;
% plot(corn(:,1),corn(:,2),"-k");
% hold on;
% plot([corn(1,1),corn(4,1)],[corn(1,2),corn(4,2)],"-k");
% plot(cart_new(:,1),cart_new(:,2),"og","MarkerFaceColor","g","MarkerSize",8);
% axis equal;
% hold off;


%% Make the x0 file where the z-coordinate will be 0.0 as it is SLG
% Sorting the cart_new in increasing order of x-coordinates

% cart_new_2 = sortrows(cart_new, 1);
% x0(:,1) = [1:size(cart_new_2)];
% x0(:,2) = atom_type;
% x0(:,3) = cart_new_2(:,1);
% x0(:,4) = cart_new_2(:,2);
% x0(:,5) = 0.0;

dump_xyz(:,1) = dump_x;
dump_xyz(:,2) = dump_y;
dump_xyz(:,3) = dump_z;
dump_x_sorted = sortrows(dump_xyz, 1); % sorted in increasing oredr to x-coordinates

x0(:,1) = [1:size(dump_xyz,1)]; % 1st column to give number 
x0(:,2) = atom_type;
x0(:,3:5) = dump_x_sorted;

end
