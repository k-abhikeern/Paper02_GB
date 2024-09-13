clc;
clear all;
close all;

NMD = load('./NMD.mat');
SED = load('./SEDavg.mat');

iseedf = 1;
%% Make the "irrkpt_symm_point.txt" file

%% First quardrant kpoints
NMD.kpt.cart_1st_quard_index = find(NMD.kpt.cart(:,1)> -0.001 & NMD.kpt.cart(:,2)> -0.001);
kpt = NMD.kpt.cart(NMD.kpt.cart_1st_quard_index,:);

% Plot the corner points
% plot(corn(:,1),corn(:,2),'*-k');
% hold on;
% plot([corn(6,1),corn(1,1)],[corn(6,2),corn(1,2)],'*-k');
% hold on;

% Plot the quardrant points
% plot(NMD.kpt.cart_new(:,1),NMD.kpt.cart_new(:,2),'s','MarkerSize',3,'MarkerFaceColor','r','MarkerEdgeColor','r') % change this in an appropriate way for zig-zag and arm chair
% hold on;
% plot(kpt(:,1),kpt(:,2),'bs');
% hold off;
% xlabel('k_x/(2*pi/a)');
% ylabel('k_y/(2*pi/a)');
% title('BZ');
% grid on;

% name = strcat(NMD.str.main,'BZ_finally','_',int2str(NMD.Nx),'_',...
%     int2str(NMD.Ny),'_',int2str(NMD.Nz)); % will create FIG1, FIG2,...
% saveas(gcf,name,'jpg');

kpt =  kpt*inv(NMD.latvec_rec); % so that it displays the value on command window
NMD.kpt.cart_1st_quard_uvw = kpt;%(always correct --double checked)

figure;
plot(NMD.kpt.cart_uvw(:,1),NMD.kpt.cart_uvw(:,2),'^k','MarkerSize',8);
hold on;
plot(kpt(:,1),kpt(:,2),'or','MarkerFaceColor','r');
hold off;

% Saving the 1st quardrant symmetric kpoint names 
if exist(['./' int2str(iseedf) '/irrkpt_symm_point' '.txt'], 'file')~=0
    system(['rm -f ./' int2str(iseedf) '/irrkpt_symm_point' '.txt']);
end
str_write=strcat(NMD.str.main, '/', int2str(iseedf),'/irrkpt_symm_point','.txt');
fileID = fopen(str_write,'w');

index = NMD.kpt.cart_1st_quard_index;
for p = 1:size(index,1)
    irrkpt_index = index(p);
    kpt1=NMD.kpt.cart(irrkpt_index,:);
    
    fprintf(fileID,'%d \t',irrkpt_index);
    
    for q = 1:size(NMD.kpt.cart,1)
        
        kpt2=NMD.kpt.cart(q,:);
        
        istrue=issym(kpt1,kpt2); % issym function calling               
        
        if istrue == 1 && (irrkpt_index ~= q)            
            
            fprintf(fileID,'%d \t',q);           
            
        end        
    end
    fprintf(fileID,'\n');    
end


 
%% Coount the no. of lines in the irrkpt_symm_point.txt file
str.cmd = ("wc -l < ./1/irrkpt_symm_point.txt");
[tmp,lines] = system(str.cmd);
lines = str2num(lines);

%%
% SED.SEDpoi_kpt_avg(1:NMD.NUM_KPTS,1:(NMD.NUM_TSTEPS/2),1:NMD.NUM_MODES) = 0.0;

% for m = 1:lines % m = batch no.
%     % Now go line by line into the "irrkpt_symm_point.txt" file
%     
%     
% end

fid = fopen('./1/irrkpt_symm_point.txt','r');
% tline = fgetl(fid);
% index = str2num(tline);
% disp(index);

SED.SEDpoi_kpt_avg(1:NMD.NUM_KPTS,1:(NMD.NUM_TSTEPS/2),1:NMD.NUM_MODES) = 0.0;
while ~feof(fid)    
    tline = fgetl(fid);
    index = str2num(tline);
    disp(index);
    batch_symm_kpoints = size(index,2);
    
%     SED.SEDpoi_kpt_avg(1:NMD.NUM_KPTS,1:(NMD.NUM_TSTEPS/2),1:NMD.NUM_MODES) = 0.0;
    for i = 1:batch_symm_kpoints
%         disp(index(1,1));
        SED.SEDpoi_kpt_avg(index(1,1),:,:) = SED.SEDpoi_kpt_avg(index(1,1),:,:) + SED.SEDavg(index(i),:,:);
        
%         SED.SEDpoi_kpt_avg(index(1,1),:,imode) = SED.SEDpoi_kpt_avg(ikpt,:,imode)+SED.SEDpoi_kpt_avg(ikpt,:,imode)
    end
%     index(1,1)
    SED.SEDpoi_kpt_avg(index(1,1),:,:) = SED.SEDpoi_kpt_avg(index(1,1),:,:)/batch_symm_kpoints;
%     disp(index);
end
fclose(fid);


save(strcat(NMD.str.main,'NMD.mat'), '-struct', 'NMD'); %saving NMD structure variables as a file with name "NMDavg.mat"
save(strcat(NMD.str.main,'SEDavg.mat'), '-struct', 'SED'); %saving SED structure variables as a file with name "SEDavg.mat"