clc;
clear all;
close all;
NMD = load('./NMD.mat');

[tmp,str.main]=system('pwd');

%--Good that it was earlier inside the loop. And so everytime for a new
%loop the inv was becoming back to original quantity. Now its all fine and good.
Qinv = inv(NMD.Q); % metal % required to convert velxprime back to velx


%% reverting back to orthogonal system (x-y)
% x_new back to x_old as x_old = [Q_inv]*x_new

NMD.x0(:,3:5) = (Qinv*NMD.x0(:,3:5)')'; % metal

%%
NMD.x0(:,3:5) = NMD.x0(:,3:5)/(NMD.LJ.sigma*(1e+10)); % metal --> LJ
NMD.alat = NMD.alat/(NMD.LJ.sigma*(1e+10)); % metal --> LJ
%--
% NMD.NUM_MODES = 3*NMD.NUM_ATOMS_UCELL;
%% --------------------------------------------------------------------------
SED.SEDavg(1:NMD.NUM_KPTS,1:(NMD.NUM_TSTEPS/2),1:NMD.NUM_MODES) = 0.0;
for   iseed = 1:NMD.NUM_SEEDS
    %--------------------------------------------------------------------------
    %--removing previously built SED_Phi_* files
    
    if exist(['./' int2str(iseed)], 'file')~=0
        system(['rm ./' int2str(iseed) '/SED_Phi_*']);
    end
    
    
    %--
    %--------------------------------------------------------------------------
    ikslice = 1;
    %--------------------------------------------------------------------------
    
    
    SED.SED(1:NMD.NUM_KPTS,1:(NMD.NUM_TSTEPS/2),1:NMD.NUM_MODES) = 0.0;
    
    %--------------------------------------------------------------------------
    %tic
    %--------------------------------------------------------------------------
    for ifft = 1:NMD.NUM_FFTS
        %VElOCITIES
        str_read=...
            strcat(...
            str.main ,'/dump_',int2str(iseed),'_',int2str(ifft),'.vel')
        fid=fopen(str_read,'r') % str_read
        dummy = textscan(fid,'%f%f%f','delimiter',' ','commentstyle', '--');
        
        %     %--
        %     % Converting velocity from metal units to LJ units as dump_1_1.vel is
        %     % currently in metal units (Ang/ps)
        %     factor = (NMD.LJ.tau/NMD.LJ.sigma) * (NMD.constant.ang2m /1e-12);
        %         for idim=1:3
        %             dummy{idim}(:) = dummy{idim}(:)*factor; % Ang/ps --> LJ
        %         end
        %
        %     %--
        
        fclose(fid);
        %Store velocity data of all atoms: subtract off the last time step
        velx = zeros(NMD.NUM_ATOMS,NMD.NUM_TSTEPS);
        vely = zeros(NMD.NUM_ATOMS,NMD.NUM_TSTEPS);
        velz = zeros(NMD.NUM_ATOMS,NMD.NUM_TSTEPS);
        velxprime = zeros(NMD.NUM_ATOMS,NMD.NUM_TSTEPS);
        velyprime = zeros(NMD.NUM_ATOMS,NMD.NUM_TSTEPS);
        velzprime = zeros(NMD.NUM_ATOMS,NMD.NUM_TSTEPS);
        %--------------------------------------------------------------------------
        %toc
        %--------------------------------------------------------------------------
        
        %--------------------------------------------------------------------------
        %tic
        %--------------------------------------------------------------------------
        
        %--
        %        Qinv = inv(NMD.Q); % metal
        
        %%
        
        % v_new back to v_old as v_old = [Q_inv]*v_new
        for iatom = 1:NMD.NUM_ATOMS
            velxprime(iatom,1:NMD.NUM_TSTEPS) =...
                dummy{1}...
                (iatom:NMD.NUM_ATOMS:(length(dummy{1}(:))-NMD.NUM_ATOMS));% metal
            velyprime(iatom,1:NMD.NUM_TSTEPS) =...
                dummy{2}...
                (iatom:NMD.NUM_ATOMS:(length(dummy{1}(:))-NMD.NUM_ATOMS));% metal
            velzprime(iatom,1:NMD.NUM_TSTEPS) =...
                dummy{3}...
                (iatom:NMD.NUM_ATOMS:(length(dummy{1}(:))-NMD.NUM_ATOMS));% metal
        end
        velx = Qinv(1,1)*velxprime+Qinv(1,2)*velyprime+Qinv(1,3)*velzprime;
        vely = Qinv(2,1)*velxprime+Qinv(2,2)*velyprime+Qinv(2,3)*velzprime;
        velz = Qinv(3,1)*velxprime+Qinv(3,2)*velyprime+Qinv(3,3)*velzprime;
        
        %--
        % Converting velocity from metal units to LJ units as dump_1_1.vel is
        % currently in metal units (Ang/ps)
        %     factor = (NMD.LJ.tau/NMD.LJ.sigma) * (NMD.constant.ang2m /1e-12);
        %         for idim=1:3
        %             dummy{idim}(:) = dummy{idim}(:)*factor; % Ang/ps --> LJ
        %         end
        factor = (NMD.LJ.tau/NMD.LJ.sigma) * (NMD.constant.ang2m /NMD.constant.ps2s);
        velx = velx*factor; % metal --> SI --> LJ
        vely = vely*factor;
        velz = velz*factor;
        %--
        
        clear velxprime;
        clear velyprime;
        clear velzprime;
        %--------------------------------------------------------------------------
        %toc
        %--------------------------------------------------------------------------
        %Remove dummy
        clear dummy;
        %Set mass array
        %     m = repmat(NMD.mass(:,1),1,NMD.NUM_TSTEPS);
        m = NMD.mass(:,1)/NMD.mass_basis; % amu --> LJ
        %EIGENVECTORS
        eigenvec = NMD.eigvec;
        %FREQUENCIES
        %    freq = NMD.freq;
        %Zero main SED FP: this gets averaged as you loop over the NUM_FFTS
        Q = zeros(1,NMD.NUM_TSTEPS);
        QDOT = zeros(1,NMD.NUM_TSTEPS);
        
        SED.SED(...
            size(NMD.kptlist(:,1:3,ikslice),1),...
            1:(NMD.NUM_TSTEPS/2),1:NMD.NUM_MODES) = 0.0;
        %--------------------------------------------------------------------------
%         tic
        %--------------------------------------------------------------------------
        
        %--
        %         NMD.x0(:,3:5) = NMD.x0(:,3:5)/(NMD.LJ.sigma*(1e+10)); % metal --> LJ
        %         NMD.alat = NMD.alat/(NMD.LJ.sigma*(1e+10)); % metal --> LJ
        %--
        
        for ikpt = 1:size(NMD.kptlist(:,1:3,ikslice),1)
            
            spatial = 2*pi*1i*(...
                NMD.x0(:,3)*( (NMD.kptlist(ikpt,1,ikslice))/(NMD.alat*NMD.Nx) ) +...
                NMD.x0(:,4)*( (NMD.kptlist(ikpt,2,ikslice))/(NMD.alat*NMD.Ny) ) +...
                NMD.x0(:,5)*( (NMD.kptlist(ikpt,3,ikslice))/(NMD.alat*NMD.Nz) ) );
            
            kindex = NMD.kpt_index(ikpt,ikslice);
            
            for imode = 1:NMD.NUM_MODES
                
                tic
                eigx = repmat(...
                    conj(...
                    eigenvec(...
                    ((NMD.NUM_ATOMS_UCELL*3)*(kindex-1)+1)...
                    :3:...
                    ((NMD.NUM_ATOMS_UCELL*3)*kindex),imode...
                    )...
                    ),NMD.NUM_UCELL_COPIES,1);
                
                eigy = repmat(...
                    conj(...
                    eigenvec(...
                    ((NMD.NUM_ATOMS_UCELL*3)*(kindex-1)+2)...
                    :3:...
                    ((NMD.NUM_ATOMS_UCELL*3)*kindex),imode...
                    )...
                    ),NMD.NUM_UCELL_COPIES,1);
                
                eigz = repmat(...
                    conj(...
                    eigenvec(...
                    ((NMD.NUM_ATOMS_UCELL*3)*(kindex-1)+3)...
                    :3:...
                    ((NMD.NUM_ATOMS_UCELL*3)*kindex),imode...
                    )...
                    ),NMD.NUM_UCELL_COPIES,1);
                
                QDOT = sum(...
                    bsxfun(@times,...
                    bsxfun(@times, velx, eigx) + ...
                    bsxfun(@times, vely, eigy) + ...
                    bsxfun(@times, velz, eigz) ...
                    , exp(spatial).*(sqrt(m/NMD.NUM_UCELL_COPIES)) )...
                    , 1 );
                
                KEFFT = real(fft(QDOT)).^2 + imag(fft(QDOT)).^2;%LJ
                
                SED.SED(ikpt,:,imode) =...
                    SED.SED(ikpt,:,imode)+KEFFT(1:(NMD.NUM_TSTEPS/2)) ;%LJ
                
                toc
                
                disp(imode);
            end %END imode
            
            
        end %END ikpt
        %--------------------------------------------------------------------------
        %toc
        %--------------------------------------------------------------------------
    end %END ifft
    
    %Average over FFTS
    factor = ((NMD.t_step^2))/(2*pi*NMD.NUM_TSTEPS);
    SED.SED = factor*SED.SED/NMD.NUM_FFTS;
    %Define frequencies
    SED.omega = (1:NMD.NUM_OMEGAS)*(NMD.w_max/NMD.NUM_OMEGAS);
    %Output SED
%     for ikpt = 1:size(NMD.kptlist(:,1:3,ikslice),1)
%         str_write_single=...
%             strcat(str.main, '/', int2str(iseed),'/SED_Phi_',...
%             num2str(NMD.kptlist(ikpt,1,ikslice)),...
%             num2str(NMD.kptlist(ikpt,2,ikslice)),...
%             num2str(NMD.kptlist(ikpt,3,ikslice)),...
%             '_',int2str(iseed),'.txt');
%         output(1:length(SED.omega),1) = SED.omega;
%         output(1:length(SED.omega),2:(NMD.NUM_MODES+1)) = SED.SED(ikpt,:,:);
%         dlmwrite(str_write_single,output,'delimiter',' ');
%         clear output
%     end %END ikpt
    SED.SEDavg = SED.SEDavg + SED.SED;
    
    iseed
end %END iseed
SED.SEDavg = SED.SEDavg/NMD.NUM_SEEDS;
SED.SED = SED.SEDavg; % SED.SED has now been replaced with SED.SEDavg value
save(strcat(NMD.str.main,'NMD.mat'), '-struct', 'NMD'); %saving NMD structure variables as a file with name "NMDavg.mat"
save(strcat(NMD.str.main,'SEDavg.mat'), '-struct', 'SED'); %saving SED structure variables as a file with name "SEDavg.mat"