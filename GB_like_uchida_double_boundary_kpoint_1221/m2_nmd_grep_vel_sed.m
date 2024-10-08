NMD = load('./NMD.mat');

%GREP LAMMPS DUMP----------------------------------------------------------

for iseed=1:NMD.NUM_SEEDS
    for ifft=1:NMD.NUM_FFTS
        cmd = ['grep -A ' int2str(NMD.NUM_ATOMS) ...
            ' "ITEM: ATOMS vx vy vz" dump_' int2str(iseed) '_' int2str(ifft) '.vel > vel.vel'];
        system(cmd);
        cmd = ['grep -v "ITEM: ATOMS vx vy vz" vel.vel > vel2.vel'];
        system(cmd);
        cmd = ['rm dump_' int2str(iseed) '_' int2str(ifft) '.vel'];
        system(cmd);
        cmd = ['mv vel2.vel dump_' int2str(iseed) '_' int2str(ifft) '.vel'];
        system(cmd);
    end
end

str.cmd = 'rm -rf vel.vel';
system(str.cmd);