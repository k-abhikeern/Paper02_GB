function f_gamma_get_cartesian(gamma_gulp_name_del,gamma_gulp_name,vec,x0,atom_type)

%% Vectors
if exist(['gamma_gulp_vec.txt'], 'file')~=0
    system(['rm -f gamma_gulp_vec.txt']);
end

fid = fopen("gamma_gulp_vec.txt","w");
fprintf(fid,"%d \t %d \t %d \n %d \t %d \t %d \n %d \t %d \t %d",vec(1,:),vec(2,:),vec(3,:));

fclose(fid);

if exist([gamma_gulp_name], 'file')~=0
    str.cmd = strcat("rm -f"," ",gamma_gulp_name);
    system(str.cmd);
end

fidin = fopen(gamma_gulp_name_del,'r');
fidout = fopen(gamma_gulp_name,'w');
fidin_gulp_vec = fopen("gamma_gulp_vec.txt","r");

expression_1 = '@vectors';

while ~feof(fidin)
    tline = fgetl(fidin);
    fprintf(fidout,"%s \n",tline);
    
    if strcmp(tline,expression_1)==1
        
        while ~feof(fidin_gulp_vec)
            tline = fgetl(fidin_gulp_vec);
            fprintf(fidout,"%s \n",tline);
        end
    end
end

fclose(fidin);
fclose(fidin_gulp_vec);
fclose(fidout);


%% Cartesian
str.cmd = strcat("cp"," ",gamma_gulp_name," ","gulp_disp_lj_conv_num_del2.tmp");% Copy previous created file to this new tmp file
system(str.cmd);

gamma_gulp_name_del = "gulp_disp_lj_conv_num_del2.tmp";

% Lower layer
cnt = 1;
for i = 1:size(x0,1)
    if x0(i,2) == atom_type(1)
        info(cnt,:) = x0(i,:);
        cnt = cnt+1;
    end
end

if exist(['./gamma_gulp_lower.txt'], 'file')~=0
    system(['rm ./gamma_gulp_lower.txt']);
end

fileID = fopen('gamma_gulp_lower.txt','w');
for i = 1:size(info,1)
    fprintf(fileID,'C    core  %6.5f %6.5f %6.5f   0.0     1.0 \n',info(i,3),info(i,4),info(i,5));
end
fclose(fileID);

% Upper layer
% cnt = 1;
% for i = 1:size(x0,1)
%     if x0(i,2) == atom_type(2)
%         info(cnt,:) = x0(i,:);
%         cnt = cnt+1;
%     end
% end

if exist(['./gamma_gulp_upper.txt'], 'file')~=0
    system(['rm ./gamma_gulp_upper.txt']);
end

% fileID = fopen('gamma_gulp_upper.txt','w');
% for i = 1:size(info,1)
%     fprintf(fileID,'C2    core  %6.5f %6.5f %6.5f   0.0     1.0 \n',info(i,3),info(i,4),info(i,5));
% end
% fclose(fileID);

%
fidin = fopen(gamma_gulp_name_del,'r');
fidin_gulp_lower = fopen("gamma_gulp_lower.txt","r");
% fidin_gulp_upper = fopen("gamma_gulp_upper.txt","r");
fidout = fopen(gamma_gulp_name,'w');
expression_1 = strcat('@cartesian_r1'," ");
expression_2 = strcat('@cartesian_r2'," ");

while ~feof(fidin)
    tline = fgetl(fidin);
    fprintf(fidout,"%s \n",tline);
    
    if strcmp(tline,expression_1)==1
        
        while ~feof(fidin_gulp_lower)
            tline = fgetl(fidin_gulp_lower);
            fprintf(fidout,"%s \n",tline);
        end
        
%     elseif strcmp(tline,expression_2)==1
%         while ~feof(fidin_gulp_upper)
%             tline = fgetl(fidin_gulp_upper);
%             fprintf(fidout,"%s \n",tline);
%         end
        
    end
end

fclose(fidin);
fclose(fidin_gulp_lower);
% fclose(fidin_gulp_upper);
fclose(fidout);

%% Delete occurences of '@vectors' , '@cartesian_r1' & '@cartesian_r2' and other unwanted files
str.cmd = strcat("sed -i '/@vectors/d'"," ", gamma_gulp_name);
system(str.cmd);

str.cmd = strcat("sed -i '/@cartesian_r*/d'"," ", gamma_gulp_name);
system(str.cmd);

str.cmd = ("rm -f gulp_disp_lj_conv_num_del2.tmp ");
system(str.cmd);

end