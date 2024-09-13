function vel = m_gulp_lj_groupvel(kpt,NUM_ATOMS_UCELL,MASS,ALAT,m_lj,str_main,gulp_path,name)
%--------------------------------------------------------------------------
%eigvec = gulp_lj_eig(kpt,NUM_ATOMS_UCELL)
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%This is required to get the kpt + dk to properly input to GULP
    format long
%--------------------------------------------------------------------------

constant = m_constant;
lj = m_lj;

gulp.dk = 1E-4;

if kpt(1) == 0 & kpt(2) == 0 & kpt(3) == 0 
    for idim=1:3
        kpt(idim) = kpt(idim) + gulp.dk;
    end
end



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

%str.cmd = [gulp_path 'disp']; system(str.cmd);
str.cmd = [gulp_path '< disp.gin > disp.gout']; system(str.cmd);%str.cmd = [gulp_path ' disp']; system(str.cmd);

%grep out eigenvectors
	vel = zeros(3*NUM_ATOMS_UCELL,3);		
    str1 = 'grep -A';
    str2 = [int2str(3*NUM_ATOMS_UCELL+3),...
        ' -m1 -e  "Group velocities" ' str_main 'disp.gout | sed "1,4d" > ' str_main 'vel_grep.dat'];
    str.cmd = [str1,str2]; system(str.cmd);
    
%read in eigvec to sort properly		
    str.read= [str_main 'vel_grep.dat'];
    fid=fopen(str.read);
    dummy = textscan(fid,'%d%f%f%f%s','Delimiter','\t',...
        'commentStyle', '--'); 
    fclose(fid); system(['rm ' str_main 'vel_grep.dat']); 
    
    for idim=1:3
        vel(:,idim) = dummy{idim+1}(:)*constant.c*2*pi*lj.tau*constant.ang2m/lj.sigma; % cm-1.Ang --> LJ (correct conversion--> checked)
    end
   

%--

for idim = 1:3
    if kpt(idim)==0.5
        freq = m_gulp_lj_freq(kpt,NUM_ATOMS_UCELL,MASS,ALAT,...
            lj,str_main,gulp_path,name);
        kpt(idim) = kpt(idim) - gulp.dk;
        freq_mdk = m_gulp_lj_freq(kpt,NUM_ATOMS_UCELL,MASS,ALAT,...
            lj,str_main,gulp_path,name);
        vel(:,idim) = (( freq - freq_mdk )/ gulp.dk / 4 );
    %Put kpt back to orig
        kpt(idim) = kpt(idim) + gulp.dk;

    elseif kpt(idim)==-0.5
        freq = m_gulp_lj_freq(kpt,NUM_ATOMS_UCELL,MASS,ALAT,...
           lj,str_main,gulp_path,name);
        kpt(idim) = kpt(idim) + gulp.dk;
        freq_pdk = m_gulp_lj_freq(kpt,NUM_ATOMS_UCELL,MASS,ALAT,...
           lj,str_main,gulp_path,name);
        vel(:,idim) = (( freq_pdk - freq )/ gulp.dk / 4 );
    %Put kpt back to orig
        kpt(idim) = kpt(idim) - gulp.dk;
    end
end

%--
end