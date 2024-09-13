function x0_AIREBO = f_min_x0(gamma_vec,num_atoms,atom_types,mass,x0_cvt)

% Making the lmp_metal_min.in.x0.1 file
xhi = gamma_vec(1,1); % changes NMD.Nx-1
yhi = gamma_vec(2,2);
zhi = gamma_vec(3,3);
xy  = gamma_vec(2,1);
xz  = gamma_vec(3,1);
yz  = gamma_vec(3,2);

xlo = 0.0;%lo_bound(1) - min( [0.0 xy xz xy+xz] );
ylo = 0.0;%lo_bound(2) - min( [0.0 yz] );
zlo = 0.0;%lo_bound(3);

% zlo = -50.0;zhi = 50.0;% just for making the z value high
%-------------------------------------------------------
if exist(['./lmp_metal_min.in.x0.1'], 'file')~=0
    system(['rm -f ./lmp_metal_min.in.x0.1']);
end
fileID = fopen("lmp_metal_min.in.x0.1","w");
fprintf(fileID, "#\n%d\t atoms\n %d\t atom types\n\n\n",num_atoms,atom_types);

fprintf(fileID, " %f %f\t xlo xhi\n %f %f\t ylo yhi\n %f %f\t zlo zhi\n %f %f %f\t xy xz yz\n\n\n", .....
    xlo,xhi,ylo,yhi,zlo,zhi,xy,xz,yz);

fprintf(fileID, "Masses\n\n 1 %f\n\n",mass);
fprintf(fileID, "Atoms\n\n");
fclose(fileID);
str.write=strcat('./lmp_metal_min.in.x0.1');
dlmwrite(str.write,x0_cvt,'-append','delimiter',' ');

%% Now call the lammps with AIREBO potential and minimize the above coordinates

% cmd = strcat("lmp_serial <"," ","in.min_AIREBO");
cmd = strcat("/share/apps/lammps/lmp_mpi <"," ","in.min_AIREBO");
system(cmd);

newfile = "x0_AIREBO.txt";
cmd = strcat("tail -n"," ",num2str(num_atoms)," ","dump_min_AIREBO.lammpstrj > ",newfile);
system(cmd);


% Read the new file and store its value in a variable
fidin = fopen(newfile,"r");
data = textscan(fidin,"%f\t%f\t%f\t%f\t%f");
fclose(fidin);
x0_AIREBO = [data{1},data{2},data{3},data{4},data{5}];

end