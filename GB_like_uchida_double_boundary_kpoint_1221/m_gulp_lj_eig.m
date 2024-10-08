function eigvec = m_gulp_lj_eig(kpt,NUM_ATOMS_UCELL,MASS,ALAT,str_main,gulp_path,name)
%--------------------------------------------------------------------------
%eigvec = gulp_lj_eig(kpt,NUM_ATOMS_UCELL)
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%This is required to get the kpt + dk to properly input to GULP
    format long
%--------------------------------------------------------------------------

orig(1).str = 'KPT';
change(1).str = strcat( num2str(kpt(1)),'\t',num2str(kpt(2)),'\t',num2str(kpt(3)) );
orig(2).str = 'MASS';
change(2).str = num2str(MASS);
orig(3).str = 'ALAT';
change(3).str = num2str(ALAT);

m_change_file_strings(...
    [str_main name],...
    orig,...
    [str_main 'disp.gin'],...
    change);
%  str.cmd = [gulp_path '<disp.gin']; system(str.cmd);
 str.cmd = [gulp_path ' < disp.gin > disp.gout']; 
 system(str.cmd);%str.cmd = [gulp_path ' disp']; system(str.cmd);

%grep out eigenvectors
    eigvec = zeros(3*NUM_ATOMS_UCELL,3*NUM_ATOMS_UCELL);
    str1 = 'grep -A ';
    str2 = [int2str(3*NUM_ATOMS_UCELL),...
        ' " 1 x" ' str_main 'disp.gout > ' str_main 'eigvec_grep.dat'];
    str.cmd = [str1,str2]; system(str.cmd);
    str.cmd =...
        ['sed ''s/x//g'' ' str_main 'eigvec_grep.dat > '...
        str_main 'eigvec2.dat'];
    system(str.cmd);
    system(['rm ' str_main 'eigvec_grep.dat']);
    str.cmd = ['sed ''s/y//g'' ' str_main 'eigvec2.dat > '...
        str_main 'eigvec3.dat']; 
    system(str.cmd); system(['rm ' str_main 'eigvec2.dat']);  
    str.cmd = ['sed ''s/z//g'' ' str_main 'eigvec3.dat > '...
        str_main 'eigvec4.dat']; 
    system(str.cmd); system(['rm ' str_main 'eigvec3.dat']);

%read in eigvec to sort properly		
    str.read= [str_main 'eigvec4.dat'];
    fid=fopen(str.read);
    dummy = textscan(fid,'%f%f%f%f%f%f%f','Delimiter','\t',...
        'commentStyle', '--'); 
    fclose(fid);
    system(['rm ' str_main 'eigvec4.dat']); 
    
if kpt(1) == 0 & kpt(2) == 0 & kpt(3) == 0 
%Gamma has only real components	
    eigvec = zeros(3*NUM_ATOMS_UCELL,3*NUM_ATOMS_UCELL); 
    for imode = 1:(3*NUM_ATOMS_UCELL/3/2)
    eigvec(:,(imode-1)*6+1) =...
        dummy{2}((imode-1)*3*NUM_ATOMS_UCELL+1:(imode)*3*NUM_ATOMS_UCELL);
    eigvec(:,(imode-1)*6+2) =...
        dummy{3}((imode-1)*3*NUM_ATOMS_UCELL+1:(imode)*3*NUM_ATOMS_UCELL);
    eigvec(:,(imode-1)*6+3) =...
        dummy{4}((imode-1)*3*NUM_ATOMS_UCELL+1:(imode)*3*NUM_ATOMS_UCELL);
    eigvec(:,(imode-1)*6+4) =...
        dummy{5}((imode-1)*3*NUM_ATOMS_UCELL+1:(imode)*3*NUM_ATOMS_UCELL);
    eigvec(:,(imode-1)*6+5) =...
        dummy{6}((imode-1)*3*NUM_ATOMS_UCELL+1:(imode)*3*NUM_ATOMS_UCELL);
    eigvec(:,(imode-1)*6+6) =...
        dummy{7}((imode-1)*3*NUM_ATOMS_UCELL+1:(imode)*3*NUM_ATOMS_UCELL);
    end
else
    
    
%Put Real and Imag in right place		
eigvec = zeros(3*NUM_ATOMS_UCELL,3*NUM_ATOMS_UCELL); 
    for imode = 1:(3*NUM_ATOMS_UCELL/3)
        eigvec(:,(imode-1)*3+1) =...
            dummy{2}((imode-1)*3*NUM_ATOMS_UCELL+1:(imode)*3*NUM_ATOMS_UCELL)...
            + i*dummy{3}((imode-1)*3*NUM_ATOMS_UCELL+1:(imode)*3*NUM_ATOMS_UCELL);
        eigvec(:,(imode-1)*3+2) =...
            dummy{4}((imode-1)*3*NUM_ATOMS_UCELL+1:(imode)*3*NUM_ATOMS_UCELL)...
            + i*dummy{5}((imode-1)*3*NUM_ATOMS_UCELL+1:(imode)*3*NUM_ATOMS_UCELL);
        eigvec(:,(imode-1)*3+3) =...
            dummy{6}((imode-1)*3*NUM_ATOMS_UCELL+1:(imode)*3*NUM_ATOMS_UCELL)...
            + i*dummy{7}((imode-1)*3*NUM_ATOMS_UCELL+1:(imode)*3*NUM_ATOMS_UCELL);
      end
   end

 end
