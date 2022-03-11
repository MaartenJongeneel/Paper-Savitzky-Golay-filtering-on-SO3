clearvars; clc; close all; addpath('functions');
set(groot,'defaulttextinterpreter','latex'); set(groot,'defaultAxesTickLabelInterpreter','latex'); set(groot,'defaultLegendInterpreter','latex');
%% ---------------- Savitzky-Golay Filtering on SO(3) ----------------- %%
%% Constants and settings
%User inputs
doSave = false;    %Boolean: set true if you want to save figures
Fc = 1;            %Signal frequency                  [Hz]
a  = 2;            %Signal amplitude                  [deg]
te = 2;            %Signal length                     [s]
Fs = 1000;         %Sampling frequency fine grid      [Hz]
m  = 5;            %Down-sampling rate                [-]
sigma = 0.06;      %Standard deviation of added noise [rad]
n  = 15;           %Window size SG-filter             [-]
p  = 2;            %Savitzky Golay filter order       [-]

%Computed values
dt1 = 1/Fs;        %Time step                         [s]
dt2 = m/Fs;        %Time step lower sampled           [s]
t1 = (0:dt1:te);   %Signal time vector                [s]
t2 = (0:dt2:te);   %Signal time vector lower sampled  [s]
N1 = length(t1);   %Number of samples                 [-]
N2 = length(t2);   %Number of samples lower sampled   [-]

%% Preallocate memory
omg = NaN(3,N1);    omg_FD = NaN(3,N2);     
domg = NaN(3,N1);   domg_FD = NaN(3,N2);    
R = NaN(3,3,N1);    R_noise = NaN(3,3,N2);  

phi = NaN(3,N1); dphi = NaN(3,N1); ddphi = NaN(3,N1); g_noise = NaN(3,N2);

%% Creating data on SO(3)
%Create a random sine wave in R3 with first and second order derivative
% lambda0 = randn(3,1);
% lambda1 = randn(3,1);

%Vectors below are created by randn(3,1) but placed here s.t. we can give
%the values in the paper and show the corresponding plots
lambda0 = [-0.4831; 0.6064; -2.6360];
lambda1 = [ 0.9792; 1.4699; -0.4283];

for ii = 1:N1
    freq= 2*pi*Fc;
    phi(:,ii) = lambda0 + lambda1*a*sin(freq*t1(ii)); 
    dphi(:,ii) = lambda1*a*(freq)*cos(freq*t1(ii)); 
    ddphi(:,ii) = -lambda1*a*(freq)^2*sin(freq*t1(ii)); 
      
    %Compute analytically the rotation matrices, ang. vel., and ang. acc.
    R(:,:,ii) = expSO3(phi(:,ii));
    omg(:,ii) = dexpSO3(phi(:,ii))*dphi(:,ii);
    domg(:,ii) = DdexpSO3(phi(:,ii),dphi(:,ii))*dphi(:,ii) +  dexpSO3(phi(:,ii))*ddphi(:,ii);
end

%Noisy, lower sampled signal ("measurement")
cnt = 1;
for ii = 1:m:N1
R_noise(:,:,cnt) = expSO3(phi(:,ii)+sigma*randn(3,1));
cnt=cnt+1;
end

%Finite differencing from noisy lower sampled signal ("measurement"):
for ii = 2:N2-1
    omg_FD(:,ii) = vee(1/(2*dt2)*(logm((R_noise(:,:,ii+1))/R_noise(:,:,ii))-logm((R_noise(:,:,ii-1))/R_noise(:,:,ii))));
end
for ii = 2:N2-1
    domg_FD(:,ii) = 1/(2*dt2)*(omg_FD(:,ii+1)-omg_FD(:,ii-1));
end

%% ---------------- Applying the Savitzky-Golay filter ----------------- %%
%Now, from the noisy lower sampled data, we want to get back the estimated
%rotation matrix, angular velocity and angular acceleration
[R_est,omg_est,domg_est,t3] = sgolayfiltSO3(R_noise,p,n,1/dt2);

%% ---------------- Computing errors, plotting results ----------------- %%
%Time indices of R for which we have a measurement:
close all;
tR1 = find(ismember(t1,t2)==1);
tR2 = find(ismember(single(t1),single(t3))==1);

for ii = 1:length(tR1)
eR_meas(:,:,ii) = logm(R(:,:,tR1(ii))\R_noise(:,:,ii));
NeR_meas(ii) = norm(eR_meas(:,:,ii));
eomg_FD(:,ii) = omg_FD(:,ii)-omg(:,tR1(ii));
edomg_FD(:,ii) = domg_FD(:,ii)-domg(:,tR1(ii)); 
end

for ii = 1:length(tR2)
eR_est(:,:,ii) = logm(R(:,:,tR2(ii))\R_est(:,:,ii));
NeR_est(ii) = norm(eR_est(:,:,ii));
eomg_est(:,ii) = omg_est(:,ii)-omg(:,tR2(ii));
edomg_est(:,ii) = domg_est(:,ii)-domg(:,tR2(ii));
end

Eomg_FD = vecnorm(eomg_FD);
Eomg_est= vecnorm(eomg_est);
Edomg_FD = vecnorm(edomg_FD);
Edomg_est = vecnorm(edomg_est);

%Mean error in rotation
mean_ER_est = mean(NeR_est);
mean_ER_meas = mean(NeR_meas);

%Mean errors in velocity
mean_Eomg_FD = mean(Eomg_FD,'omitnan');
mean_Eomg_est = mean(Eomg_est,'omitnan');

%Mean errors in acceleration
mean_Edomg_FD = mean(Edomg_FD,'omitnan');
mean_Edomg_est = mean(Edomg_est,'omitnan');

%% Figures
%Check if figures directory exists, if not, it will create one.
if ~isfolder('figures')
    mkdir('figures');
end

%Create a plot grid
sizex = 380;
sizey = 250;
px = (0:7)*(sizex+10)+10;
py = (0:4)*(sizey+40)+45;
for  ii = 1:length(px)
    for jj = 1:length(py)
        pp{jj,ii} = [px(ii) py(jj)];
    end
end

%Plot the orientation error
figure('rend','painters','pos',[pp{1,1} sizex 0.62*sizey]);
    ha = tight_subplot(1,1,[.08 .07],[.2 .08],[0.12 0.03]);  %[gap_h gap_w] [lower upper] [left right]
    axes(ha(1));
    plot(t2,NeR_meas); hold on; 
    plot(t3,NeR_est);
    xlim([0 2]);
    xlabel('Time [s]');
    ylabel('Orientation error [rad]');
    L1 = legend('$e_{\widetilde{\mathbf{R}}}$','$e_{\widehat{\mathbf{R}}}$','NumColumns',2);
    L1.FontSize = 9;
    grid on;
    
    if doSave
        fig = gcf;
        fig.PaperPositionMode = 'auto';
        fig_pos = fig.PaperPosition;
        fig.PaperSize = [fig_pos(3) fig_pos(4)];
        print(fig,'figures/norm_eR.pdf','-dpdf','-painters')
    end
    
% Try plotting the rotation
figure('rend','painters','pos',[pp{1,2} 450 200]); 
    ha = tight_subplot(1,3,[.08 -0.2],[-0.05 0],[-0.1 -0.1]);  %[gap_h gap_w] [lower upper] [left right]
    ax = gca;
    axes(ha(1));    
    [Xding,Yding,Zding] = sphere(20);  
    hSurface = surf(Xding,Yding,Zding);hold on;
    set(hSurface,'FaceColor',[0.9 0.9 0.9],'FaceAlpha',0.3,'FaceLighting','gouraud','EdgeColor','none','LineStyle','none','BackFaceLighting','unlit',...
                 'AmbientStrengt',0.3,'DiffuseStrength',0.7,'SpecularStrength',0.5,'SpecularExponent',5,'SpecularColorReflectance',1);             
    plot3(squeeze(R(1,1,1)),squeeze(R(2,1,1)),squeeze(R(3,1,1)),'*','MarkerSize',10,'color','k','LineWidth', 2);
    plot3(squeeze(R(1,1,2:1001)),squeeze(R(2,1,2:1001)),squeeze(R(3,1,2:1001)),'color',[0.9290 0.6940 0.1250],'linewidth',1.2);
    plot3(squeeze(R_est(1,1,1:201)),squeeze(R_est(2,1,1:201)),squeeze(R_est(3,1,1:201)),'color',[0.8500 0.3250 0.0980],'linewidth',1.2);
    plot3(squeeze(R_noise(1,1,1:201)),squeeze(R_noise(2,1,1:201)),squeeze(R_noise(3,1,1:201)),'color',[0 86 140]/255,'linewidth',1.2);   
    view(128,31) %Set the initial viewpoint
    light('Position',[-1 1 1])    
    axis vis3d %Allow to rotate without changing size
    view(-122,31)
    axis off
    text(0,0,-1.5,'x')

    axes(ha(2));    
    [Xding,Yding,Zding] = sphere(20);   
    hSurface = surf(Xding,Yding,Zding);hold on;
    set(hSurface,'FaceColor',[0.9 0.9 0.9],'FaceAlpha',0.3,'FaceLighting','gouraud','EdgeColor','none','LineStyle','none','BackFaceLighting','unlit',...
                 'AmbientStrengt',0.3,'DiffuseStrength',0.7,'SpecularStrength',0.5,'SpecularExponent',5,'SpecularColorReflectance',1);
    plot3(squeeze(R(1,2,1)),squeeze(R(2,2,1)),squeeze(R(3,2,1)),'*','MarkerSize',10,'color','k','LineWidth', 2);
    plot3(squeeze(R(1,2,2:1001)),squeeze(R(2,2,2:1001)),squeeze(R(3,2,2:1001)),'color',[0.9290 0.6940 0.1250],'linewidth',1.2);
    plot3(squeeze(R_est(1,2,1:201)),squeeze(R_est(2,2,1:201)),squeeze(R_est(3,2,1:201)),'color',[0.8500 0.3250 0.0980],'linewidth',1.2);
    plot3(squeeze(R_noise(1,2,1:201)),squeeze(R_noise(2,2,1:201)),squeeze(R_noise(3,2,1:201)),'color',[0 86 140]/255,'linewidth',1.2); 
    view(128,32) %Set the initial viewpoint
    light('Position',[-1 1 1])    
    axis vis3d %Allow to rotate without chaning size
    view(90,32)
    axis off
    text(0,0,-1.5,'y')

    axes(ha(3));    
    [Xding,Yding,Zding] = sphere(20);  
    hSurface = surf(Xding,Yding,Zding);hold on;
    set(hSurface,'FaceColor',[0.9 0.9 0.9],'FaceAlpha',0.3,'FaceLighting','gouraud','EdgeColor','none','LineStyle','none','BackFaceLighting','unlit',...
                 'AmbientStrengt',0.3,'DiffuseStrength',0.7,'SpecularStrength',1,'SpecularExponent',5,'SpecularColorReflectance',0);
    plot3(squeeze(R(1,3,1)),squeeze(R(2,3,1)),squeeze(R(3,3,1)),'*','MarkerSize',10,'color','k','LineWidth', 2);
    g1 = plot3(squeeze(R(1,3,2:1001)),squeeze(R(2,3,2:1001)),squeeze(R(3,3,2:1001)),'color',[0.9290 0.6940 0.1250],'linewidth',1.2);
    g2 = plot3(squeeze(R_est(1,3,1:201)),squeeze(R_est(2,3,1:201)),squeeze(R_est(3,3,1:201)),'color',[0.8500 0.3250 0.0980],'linewidth',1.2);
    g3 = plot3(squeeze(R_noise(1,3,1:201)),squeeze(R_noise(2,3,1:201)),squeeze(R_noise(3,3,1:201)),'color',[0 86 140]/255,'linewidth',1.2); 
    view(128,32) %Set the initial viewpoint
    axis vis3d %Allow to rotate without chaning size
    light('Position',[-1 1 1]) 
    axis off
    text(0,0,-1.5,'z'); 
    L1 = legend([g1 g3 g2],{'$\mathbf{R}$','$\widetilde{\mathbf{R}}$','$\widehat{\mathbf{R}}$'},'NumColumns',3,'location','northeast');
    L1.Position(2) = 0.88;
    L1.Position(1) = 0.5-(L1.Position(3)/2);
    L1.FontSize = 9;   
    
    if doSave
        fig = gcf;
        fig.PaperPositionMode = 'auto';
        fig_pos = fig.PaperPosition;
        fig.PaperSize = [fig_pos(3) fig_pos(4)];
        print(fig,'figures/Rotation.pdf','-dpdf','-painters')
    end
    
    
% Plot velocity in 1 plot
figure('rend','painters','pos',[pp{2,1} 2*sizex 0.8*sizey]);
    ha = tight_subplot(1,3,[.05 .04],[.18 .26],[0.06 0.03]);  %[gap_h gap_w] [lower upper] [left right] 
    axes(ha(1));
    g1=plot(t2,omg_FD(1,:),'color',[0 0.4470 0.7410 0.6]); hold on; grid on
    g2=plot(t3,omg_est(1,:),'linewidth',1.5);
    g3=plot(t1,omg(1,:),'linewidth',1.5); 
    xlim([0,te]);
    xlabel('Time [s]');
    ylim([-20 20]);
    ylabel('Angular velocity [rad/s]');
    t=text(0.5,0.5,'x-component','parent',ha(1),'Fontsize',9); 
    t.Position = [ha(1).XLim(1)+0.5*(abs(ha(1).XLim(1))+abs(ha(1).XLim(2)))-0.5*t.Extent(3) ha(1).YLim(1)+1.1*(abs(ha(1).YLim(1))+abs(ha(1).YLim(2)))];
    
    axes(ha(2));
    plot(t2,omg_FD(2,:),'color',[0 0.4470 0.7410 0.6]); hold on; grid on
    plot(t3,omg_est(2,:),'linewidth',1.5);
    plot(t1,omg(2,:),'linewidth',1.5); 
    xlim([0,te]);
    ylim([-20 20]);
    xlabel('Time [s]');
    t=text(0.5,0.5,'y-component','parent',ha(2),'Fontsize',9); 
    t.Position = [ha(2).XLim(1)+0.5*(abs(ha(2).XLim(1))+abs(ha(2).XLim(2)))-0.5*t.Extent(3) ha(2).YLim(1)+1.1*(abs(ha(2).YLim(1))+abs(ha(2).YLim(2)))];
    
    axes(ha(3));
    plot(t2,omg_FD(3,:),'color',[0 0.4470 0.7410 0.6]); hold on; grid on
    plot(t3,omg_est(3,:),'linewidth',1.5);
    plot(t1,omg(3,:),'linewidth',1.5); 
    xlim([0,te]);
    ylim([-20 20]);
    xlabel('Time [s]');
    t=text(0.5,0.5,'z-component','parent',ha(3),'Fontsize',9); 
    t.Position = [ha(3).XLim(1)+0.5*(abs(ha(3).XLim(1))+abs(ha(3).XLim(2)))-0.5*t.Extent(3) ha(3).YLim(1)+1.1*(abs(ha(3).YLim(1))+abs(ha(3).YLim(2)))];
    L1 = legend([g3 g1 g2],{'Analytical solution \boldmath${{\omega}}$','Finite differencing \boldmath$\breve{{\omega}}$',...
        'Savitzky-Golay \boldmath$\widehat{{\omega}}$'},'NumColumns',3,'location','northeast');
    L1.Position(2) = 0.88;
    L1.Position(1) = 0.5-(L1.Position(3)/2);
    L1.FontSize = 9;    
    
    if doSave
        fig = gcf;
        fig.PaperPositionMode = 'auto';
        fig_pos = fig.PaperPosition;
        fig.PaperSize = [fig_pos(3) fig_pos(4)];
        print(fig,'figures/omg.pdf','-dpdf','-painters')
    end
    
% Plot acceleration in 1 plot
figure('rend','painters','pos',[pp{3,1} 2*sizex 0.8*sizey]);
    ha = tight_subplot(1,3,[.05 .04],[.18 .26],[0.06 0.03]);  %[gap_h gap_w] [lower upper] [left right] 
    axes(ha(1));
    g1=plot(t2,domg_FD(1,:),'color',[0 0.4470 0.7410 0.6]); hold on; grid on
    g2=plot(t3,domg_est(1,:),'linewidth',1.5);
    g3=plot(t1,domg(1,:),'linewidth',1.5); 
    xlim([0,te]);
    ylim([-150,150]);
    yticks([-150 -100 -50 0 50 100 150])
    yticklabels({'-150','-100','-50','0','50','100','150'})
    xlabel('Time [s]');
    ylabel('Angular acceleration [rad/s$^2$]');
    t=text(0.5,0.5,'x-component','parent',ha(1),'Fontsize',9); 
    t.Position = [ha(1).XLim(1)+0.5*(abs(ha(1).XLim(1))+abs(ha(1).XLim(2)))-0.5*t.Extent(3) ha(1).YLim(1)+1.1*(abs(ha(1).YLim(1))+abs(ha(1).YLim(2)))];
    
    axes(ha(2));
    plot(t2,domg_FD(2,:),'color',[0 0.4470 0.7410 0.6]); hold on; grid on
    plot(t3,domg_est(2,:),'linewidth',1.5);
    plot(t1,domg(2,:),'linewidth',1.5); 
    xlim([0,te]);
    ylim([-300,300]);
    yticks([-300 -200 -100 0 100 200 300])
    yticklabels({'-300','-200','-100','0','100','200','300'})
    xlabel('Time [s]');
    t=text(0.5,0.5,'y-component','parent',ha(2),'Fontsize',9); 
    t.Position = [ha(2).XLim(1)+0.5*(abs(ha(2).XLim(1))+abs(ha(2).XLim(2)))-0.5*t.Extent(3) ha(2).YLim(1)+1.1*(abs(ha(2).YLim(1))+abs(ha(2).YLim(2)))];
    
    axes(ha(3));
    plot(t2,domg_FD(3,:),'color',[0 0.4470 0.7410 0.6]); hold on; grid on
    plot(t3,domg_est(3,:),'linewidth',1.5);
    plot(t1,domg(3,:),'linewidth',1.5); 
    xlim([0,te]);
    ylim([-300,300]);
    yticks([-300 -200 -100 0 100 200 300])
    yticklabels({'-300','-200','-100','0','100','200','300'})
    xlabel('Time [s]');
    t=text(0.5,0.5,'z-component','parent',ha(3),'Fontsize',9); 
    t.Position = [ha(3).XLim(1)+0.5*(abs(ha(3).XLim(1))+abs(ha(3).XLim(2)))-0.5*t.Extent(3) ha(3).YLim(1)+1.1*(abs(ha(3).YLim(1))+abs(ha(3).YLim(2)))];
    L1 = legend([g3 g1 g2],{'Analytical solution \boldmath$\dot{{\omega}}$','Finite differencing \boldmath$\breve{\dot{{\omega}}}$',...
        'Savitzky-Golay \boldmath$\widehat{\dot{{\omega}}}$'},'NumColumns',3,'location','northeast');
    L1.Position(2) = 0.88;
    L1.Position(1) = 0.5-(L1.Position(3)/2);
    L1.FontSize = 9;    
    
    if doSave
        fig = gcf;
        fig.PaperPositionMode = 'auto';
        fig_pos = fig.PaperPosition;
        fig.PaperSize = [fig_pos(3) fig_pos(4)];
        print(fig,'figures/domg.pdf','-dpdf','-painters')
    end
