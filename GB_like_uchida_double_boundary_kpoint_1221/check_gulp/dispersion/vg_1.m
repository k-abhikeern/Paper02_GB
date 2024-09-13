%load PbTe.disp
close all;
clear all;
nBranches = 152*3;
nk = 20;
ndata = nBranches*nk;
%length = sqrt(3)*0.25; %Length of the symmetry line
length = sqrt(0.5^2); %Length of the symmetry line
unitlength = length/nk;

%making matrices which will contain values of w,k and vg
k=zeros(nBranches,nk); %x=k
omega=zeros(nBranches,nk); %y=w
vg=zeros(nBranches,nk); %z=vg

%% GULP run

% gulp_path = '/home/kunwar/softwareskunwar/gulp-5.1/Src/gulp';
% str.cmd = [gulp_path '< disp.gin > disp.gout']; system(str.cmd);

finput = fopen('SLG.disp','r');
tline = fgets(finput);
tline = fgets(finput);
tline = fgets(finput);

A = fscanf(finput, '%25f %25f' ,[2 ndata]);
count = 0;
for i = 1:nk
    for j = 1:nBranches
        count = count + 1;
        k(j,i) = A(1,count);
        omega(j,i) = A(2,count);
    end
end

k = k*unitlength;
omega = omega*0.03; %Omega in THz

%The outer loop runs for each branch and the inner loop is for each k point
figure;

%% Dispersion curve
for i=1:nBranches
    %setting up the spline interpolation
    yy=spline(k(i,:),omega(i,:));% arguments- the data points for which spline is required
    y_val=ppval(yy,k(i,:));% calculates the value of spline polynomial at points xx
    dq=fnder(yy);%calculating derivative for spline
    vg(i,:)=ppval(dq,k(i,:));
   % figure(i)
    p1 = plot(k(i,:),omega(i,:) ,'-k','LineWidth',2);
    set(p1,'Color','black');
    hold on;
    p2 = plot(k(i,:),y_val);
    set(p2,'Color','black');
end

%%%%%%%%%%%%%%%
name = strcat('dispersion curve'); % will create FIG1, FIG2,...
    xlabel('scaled units (uvw)');
    ylabel('\omega');
    grid on;
    hold on;
saveas(gcf,name,'jpg');
