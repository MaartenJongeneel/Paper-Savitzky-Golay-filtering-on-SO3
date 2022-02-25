clearvars; clc; close all; addpath('functions');
set(groot,'defaulttextinterpreter','latex'); set(groot,'defaultAxesTickLabelInterpreter','latex'); set(groot,'defaultLegendInterpreter','latex');
%% ---------------- Savitzky-Golay Filtering on SO(3) ----------------- %%
%% Constants and settings
%User inputs
doSave = true;     %Boolean: set true if you want to save figures
doPlot = false;    %Determine if you want to plot the trajectory in 3D plot
Fs = 1000;         %Sampling frequency fine grid      [Hz]
Fc = 1;            %Signal frequency                  [Hz]
a  = 2;            %Signal amplitude                  [deg]
n  = 15;           %Window size SG-filter             [-]
te = 2;            %Signal length                     [s]
order = 2;         %Savitzky Golay filter order       [-]
alpha = 5;         %Down-sampling rate                [-]
sigma = 0.06;      %Standard deviation of added noise [rad]
Fo_G = [1;0;1];    %position of G w.r.t. F            [m]
Go_E = [0;0;-0.4]; %position of E w.r.t. G            [m]
GR_E = eye(3);     %rotation of E w.r.t. G            [rad]

%Computed values
dt1 = 1/Fs;        %Time step                         [s]
dt2 = alpha/Fs;    %Time step lower sampled           [s]
t1 = (0:dt1:te);   %Signal time vector                [s]
t2 = (0:dt2:te);   %Signal time vector lower sampled  [s]
N1 = length(t1);   %Number of samples                 [-]
N2 = length(t2);   %Number of samples lower sampled   [-]
w = -n:n;          %Window for Golay
I = eye(3);        %Short hand notation
t3 = t2((n+1):(N2-(n+1)));       %Time vector filtered signal 
GH_E = [GR_E,Go_E;zeros(1,3),1]; %Transformtation from E to G

%% Preallocate memory
Fo_G = repmat(Fo_G,1,N1);   FH_E = NaN(4,4,N1);
FR_G = NaN(3,3,N1);         Fo_E = NaN(3,N1);
FH_G = NaN(4,4,N1);         FR_E = NaN(3,3,N1);

omg = NaN(3,N1);    omg_FD = NaN(3,N2);     est_omg = NaN(3,N2-length(w));
domg = NaN(3,N1);   domg_FD = NaN(3,N2);    est_domg = NaN(3,N2-length(w));
R = NaN(3,3,N1);    R_noise = NaN(3,3,N2);  est_R = NaN(3,3,N2-length(w));

g = NaN(3,N1); dg = NaN(3,N1); ddg = NaN(3,N1); g_noise = NaN(3,N2);

%% Creating data on SO(3)
%Create a random sine wave in R3 with first and second order derivative
% x0 = randn(3,1);
% x1 = randn(3,1);

%Vectors below are created by randn(3,1) but placed here s.t. we can give
%the values in the paper and show the corresponding plots
x0 = [-0.4831; 0.6064; -2.6360];
x1 = [ 0.9792; 1.4699; -0.4283];

for ii = 1:N1
    freq= 2*pi*Fc;
    g(:,ii) = x0 + x1*a*sin(freq*t1(ii)); 
    dg(:,ii) = x1*a*(freq)*cos(freq*t1(ii)); 
    ddg(:,ii) = -x1*a*(freq)^2*sin(freq*t1(ii)); 
      
    %Compute analytically the rotation matrices, ang. vel., and ang. acc.
    R(:,:,ii) = expSO3(g(:,ii));
    omg(:,ii) = dexpSO3(g(:,ii))*dg(:,ii);
    domg(:,ii) = DdexpSO3(g(:,ii),dg(:,ii))*dg(:,ii) +  dexpSO3(g(:,ii))*ddg(:,ii);
    
    %Create the homogeneous transformation matrices FH_E 
    Fo_E(:,ii)   = Fo_G(:,ii)+R(:,:,ii)*Go_E;
    FR_E(:,:,ii) = R(:,:,ii)*GR_E;
end

%Noisy, lower sampled signal ("measurement")
tel = 1;
for ii = 1:alpha:N1
g_noise(:,tel) = g(:,ii)+sigma*randn(3,1);
R_noise(:,:,tel) = expSO3(g_noise(:,tel));
tel=tel+1;
end

%Finite differencing from noisy lower sampled signal ("measurement"):
for ii = 2:N2-1
    omg_FD(:,ii) = vee(1/(2*dt2)*(logm((R_noise(:,:,ii+1))/R_noise(:,:,ii))-logm((R_noise(:,:,ii-1))/R_noise(:,:,ii))));
end
for ii = 2:N2-1
    domg_FD(:,ii) = 1/(2*dt2)*(omg_FD(:,ii+1)-omg_FD(:,ii-1));
end

%% Savitzky-Golay
%Now, from the noisy lower sampled data we want to get back the correct
%rotation matrix, angular vel. and angular acc. 

%For each time step (where we can apply the window)
cnt = 1;
for ii = (n+1):(N2-(n+1))
    %Build matrix A and vector b based on window size w
    row = 1;
    for jj = 1:length(w)
        %Time difference between 0^th element and w(jj)^th element
        Dt = (t2(ii+w(jj))-t2(ii)); 
        %Determine row of A matrix
        Ajj = I;
        for kk = 1:order
            Ajj = cat(2,Ajj,(1/kk)*Dt^kk*I); %catenation based on SG filter order
        end
        A(row:row+length(I)-1,:) = Ajj;
        b(row:row+length(I)-1,:) = vee(logm(R_noise(:,:,ii+w(jj))/R_noise(:,:,ii)));
        row = row+length(I); %Update to next row
    end
    eta_hat = (A'*A)\A'*b; %Solve the LS problem
    eta0 = eta_hat(1:3);
    eta1 = eta_hat(4:6);
    eta2 = eta_hat(7:9);
    
    %Compute analytically the rotation matrices, ang. vel., and ang. acc.
    est_R(:,:,cnt) = expSO3(eta0)*R_noise(:,:,ii);
    est_omg(:,cnt) = dexpSO3(eta0)*eta1;
    est_domg(:,cnt) = DdexpSO3(eta0,eta1)*eta1 +  dexpSO3(eta0)*eta2; 

    cnt = cnt+1;
end

%% Figures
%Create a plot grid
sizex = 380;
sizey = 250;
px = (0:7)*(sizex+10)+10;
py = (0:4)*(sizey+90)+45;
for  ii = 1:length(px)
    for jj = 1:length(py)
        pp{jj,ii} = [px(ii) py(jj)];
    end
end

tend = 2; %Define the end time to which we plot (tend should be smaller than te)
   

%% wee if we can do something with quaternions
close all
est_quat = compact(quaternion(est_R,'rotmat','frame'));
meas_quat = compact(quaternion(R_noise,'rotmat','frame'));
quat = compact(quaternion(R,'rotmat','frame'));

figure('rend','painters','pos',[pp{1,5} 1.1*sizex 2*sizey]);
    ha = tight_subplot(4,1,[.06 .06],[.08 .08],[0.08 0.03]);  %[gap_h gap_w] [lower upper] [left right]
    axes(ha(1));hold on;grid on;     
    g2=plot(t2,meas_quat(:,1),'color',[0 0.4470 0.7410]);
    g3=plot(t3,est_quat(:,1),'color',[0.8500 0.3250 0.0980],'linewidth',1.3);
    g1=plot(t1,quat(:,1),'color',[0.9290 0.6940 0.1250],'linewidth',1.3);    
    xlim([0 2]);ylim([0 1]);    
    xticks([0 0.25 0.5 0.75 1 1.25 1.5 1.75 2])
    xticklabels({'0','','0.5','','1','','1.5','','2'})
    yticks([0 0.25 0.5 0.75 1])
    yticklabels({'0','','0.5','','1'})
    ylabel('$q_0$');
   

    axes(ha(2)); hold on;grid on;
    plot(t2,meas_quat(:,2),'color',[0 0.4470 0.7410]);
    plot(t3,est_quat(:,2),'color',[0.8500 0.3250 0.0980],'linewidth',1.3);
    plot(t1,quat(:,2),'color',[0.9290 0.6940 0.1250],'linewidth',1.3);
    xlim([0 2]);ylim([-1 1]);
    xticks([0 0.25 0.5 0.75 1 1.25 1.5 1.75 2])
    xticklabels({'0','','0.5','','1','','1.5','','2'})
    yticks([-1 -0.5 0 0.5 1])
    yticklabels({'-1','','0','','1'})
    ylabel('$q_1$');

    axes(ha(3));hold on;grid on;
    plot(t2,meas_quat(:,3),'color',[0 0.4470 0.7410]);
    plot(t3,est_quat(:,3),'color',[0.8500 0.3250 0.0980],'linewidth',1.3);
    plot(t1,quat(:,3),'color',[0.9290 0.6940 0.1250],'linewidth',1.3);
    xlim([0 2]);ylim([-1 1]);
    xticks([0 0.25 0.5 0.75 1 1.25 1.5 1.75 2])
    xticklabels({'0','','0.5','','1','','1.5','','2'})
    yticks([-1 -0.5 0 0.5 1])
    yticklabels({'-1','','0','','1'})
    ylabel('$q_2$');

    axes(ha(4));hold on;grid on;
    plot(t2,meas_quat(:,4),'color',[0 0.4470 0.7410]);
    plot(t3,est_quat(:,4),'color',[0.8500 0.3250 0.0980],'linewidth',1.3);
    plot(t1,quat(:,4),'color',[0.9290 0.6940 0.1250],'linewidth',1.3);
    xlim([0 2]);ylim([-1 1]);
    xlabel('Time [s]')
    xticks([0 0.25 0.5 0.75 1 1.25 1.5 1.75 2])
    xticklabels({'0','','0.5','','1','','1.5','','2'})
    yticks([-1 -0.5 0 0.5 1])
    yticklabels({'-1','','0','','1'})
    ylabel('$q_3$');
    
    L1 = legend([g1 g2 g3],{'Analytical solution','Measurement',...
        'Savitzky-Golay'},'NumColumns',3,'location','northeast');
    L1.Position(2) = 0.95;
    L1.Position(1) = 0.5-(L1.Position(3)/2);
    L1.FontSize = 9;  

    if doSave
        fig = gcf;
        fig.PaperPositionMode = 'auto';
        fig_pos = fig.PaperPosition;
        fig.PaperSize = [fig_pos(3) fig_pos(4)];
        print(fig,'figures/quaternions.pdf','-dpdf','-painters')
    end

%%
% %plot angular velocity (\omega_x)
% figure('rend','painters','pos',[pp{1,1} sizex sizey]);
%     ha = tight_subplot(1,1,[.08 .07],[.18 .1],[0.12 0.03]);  %[gap_h gap_w] [lower upper] [left right]
%     axes(ha(1));
%     plot(t2,omg_FD(1,:)); hold on; grid on
%     plot(t3,est_omg(1,:));
%     plot(t1,omg(1,:)); 
%     xlim([0,2]);
% %     ylim([-45,30]);
%     xlabel('Time [s]');
%     ylabel('Angular velocity [rad/s]');
%     legend('$(^F$\boldmath$\widetilde{{\omega}}_{F,G})_x$','$(^F$\boldmath$\widehat{{\omega}}_{F,G})_x$','$(^F$\boldmath${{\omega}}_{F,G})_x$','location','southeast')
%     if doSave
%         fig = gcf;
%         fig.PaperPositionMode = 'auto';
%         fig_pos = fig.PaperPosition;
%         fig.PaperSize = [fig_pos(3) fig_pos(4)];
%         print(fig,'figures/omg.pdf','-dpdf','-painters')
%     end
%     
% %plot angular velocity (\omega_x)
% figure('rend','painters','pos',[pp{1,1} sizex sizey]);
%     ha = tight_subplot(1,1,[.08 .07],[.18 .1],[0.12 0.03]);  %[gap_h gap_w] [lower upper] [left right]
%     axes(ha(1));
%     plot(t2,omg_FD(1,:)); hold on; grid on
%     plot(t3,est_omg(1,:));
%     plot(t1,omg(1,:)); 
%     xlim([0,2]);
% %     ylim([-45,30]);
%     xlabel('Time [s]');
%     ylabel('Angular velocity [rad/s]');
%     legend('$(^F$\boldmath$\widetilde{{\omega}}_{F,G})_x$','$(^F$\boldmath$\widehat{{\omega}}_{F,G})_x$','$(^F$\boldmath${{\omega}}_{F,G})_x$','location','southeast')
%     if doSave
%         fig = gcf;
%         fig.PaperPositionMode = 'auto';
%         fig_pos = fig.PaperPosition;
%         fig.PaperSize = [fig_pos(3) fig_pos(4)];
%         print(fig,'figures/omg.pdf','-dpdf','-painters')
%     end    
%     
% %plot angular acceleration (\dot\omega_x)
% figure('rend','painters','pos',[pp{2,1} sizex sizey]);
%     ha = tight_subplot(1,1,[.08 .07],[.18 .1],[0.15 0.03]);  %[gap_h gap_w] [lower upper] [left right]
%     axes(ha(1));
%     plot(t2,domg_FD(1,:)); hold on; grid on
%     plot(t3,est_domg(1,:));
%     plot(t1,domg(1,:));
% %     ylim([-2000,4000]);
%     xlim([0,2]);
%     xlabel('Time [s]');
%     ylabel('Angular acceleration [rad/$s^2$]');
%     legend('$(^F$\boldmath$\widetilde{\dot{\omega}}_{F,G})_x$',...
%         '$(^F$\boldmath$\widehat{\dot{\omega}}_{F,G})_x$',...
%         '$(^F$\boldmath${\dot{\omega}}_{F,G})_x$',...
%         'location','northwest')
%     if doSave
%         fig = gcf;
%         fig.PaperPositionMode = 'auto';
%         fig_pos = fig.PaperPosition;
%         fig.PaperSize = [fig_pos(3) fig_pos(4)];
%         print(fig,'figures/domg.pdf','-dpdf','-painters')
%     end
%     
% %plot angular acceleration (\dot\omega_x) close up
% figure('rend','painters','pos',[pp{2,2} sizex sizey]);
%     ha = tight_subplot(1,1,[.08 .07],[.18 .1],[0.12 0.03]);  %[gap_h gap_w] [lower upper] [left right]
%     axes(ha(1));
%     plot(t2,domg_FD(1,:)); hold on; grid on
%     plot(t3,est_domg(1,:));
%     plot(t1,domg(1,:));
%     xlim([0.7,1.4]);
% %     ylim([-300,400]);
%     xlabel('Time [s]');
%     ylabel('Angular acceleration [rad/$s^2$]');
%     legend('$(^F$\boldmath$\widetilde{\dot{\omega}}_{F,G})_x$','$(^F$\boldmath$\widehat{\dot{\omega}}_{F,G})_x$','$(^F$\boldmath${\dot{\omega}}_{F,G})_x$','location','northeast')
%     if doSave
%         fig = gcf;
%         fig.PaperPositionMode = 'auto';
%         fig_pos = fig.PaperPosition;
%         fig.PaperSize = [fig_pos(3) fig_pos(4)];
%         print(fig,'figures/domg_zoom.pdf','-dpdf','-painters')
%     end
% 
%     
% %Plot trajectory in 3D plot
% if doPlot
%     figure('rend','painters','pos',[pp{1,3} sizex sizey]);
%     ha = tight_subplot(1,1,[.08 .07],[.18 .1],[0.12 0.03]);  %[gap_h gap_w] [lower upper] [left right]
%     axes(ha(1));
%     tel1 = 1;
%     plotBool = true;
%     for ii = (n*alpha+1):3*alpha:N1
%         if plotBool
%         %Plot the origin of the world coordinate frame
%         tip = [0.3*[1;0;0] 0.3*[0;1;0] 0.3*[0;0;1]];
%         plot3([0 tip(1,1)],[0 tip(2,1)],[0 tip(3,1)],'r'); hold on
%         plot3([0 tip(1,2)],[0 tip(2,2)],[0 tip(3,2)],'g');
%         plot3([0 tip(1,3)],[0 tip(2,3)],[0 tip(3,3)],'b');
%         
%         %Plot the origin of the robot joint with its unit vectors
%         tip = [Fo_G(:,ii)+0.3*R(:,1,ii) Fo_G(:,ii)+0.3*R(:,2,ii) Fo_G(:,ii)+0.3*R(:,3,ii)];
%         plot3([Fo_G(1,ii) tip(1,1)],[Fo_G(2,ii) tip(2,1)],[Fo_G(3,ii) tip(3,1)],'r'); hold on
%         plot3([Fo_G(1,ii) tip(1,2)],[Fo_G(2,ii) tip(2,2)],[Fo_G(3,ii) tip(3,2)],'g');
%         plot3([Fo_G(1,ii) tip(1,3)],[Fo_G(2,ii) tip(2,3)],[Fo_G(3,ii) tip(3,3)],'b');
%         
%         %Plot the origin of the end effector with its unit vectors
%         tip = [Fo_E(:,ii)+0.3*FR_E(:,1,ii) Fo_E(:,ii)+0.3*FR_E(:,2,ii) Fo_E(:,ii)+0.3*FR_E(:,3,ii)];
%         plot3([Fo_E(1,ii) tip(1,1)],[Fo_E(2,ii) tip(2,1)],[Fo_E(3,ii) tip(3,1)],'r'); hold on
%         plot3([Fo_E(1,ii) tip(1,2)],[Fo_E(2,ii) tip(2,2)],[Fo_E(3,ii) tip(3,2)],'g');
%         plot3([Fo_E(1,ii) tip(1,3)],[Fo_E(2,ii) tip(2,3)],[Fo_E(3,ii) tip(3,3)],'b');
%         
%         %Plot the noisy "measured" coordinate frame
%         tip = [Fo_G(:,ii)+0.3*est_R(:,1,tel1) Fo_G(:,ii)+0.3*est_R(:,2,tel1) Fo_G(:,ii)+0.3*est_R(:,3,tel1)];
%         plot3([Fo_G(1,ii) tip(1,1)],[Fo_G(2,ii) tip(2,1)],[Fo_G(3,ii) tip(3,1)],'r'); hold on
%         plot3([Fo_G(1,ii) tip(1,2)],[Fo_G(2,ii) tip(2,2)],[Fo_G(3,ii) tip(3,2)],'g');
%         plot3([Fo_G(1,ii) tip(1,3)],[Fo_G(2,ii) tip(2,3)],[Fo_G(3,ii) tip(3,3)],'b');
%         tel1=tel1+3*1;
%         
%         if tel1 >= length(t3)
%             plotBool = false;
%         end
%         
%         text(0,0-0.1,0-0.1,'F');
%         text(Fo_E(1,ii),Fo_E(2,ii)-0.1,Fo_E(3,ii)-0.1,'E');
%         text(Fo_G(1,1),Fo_G(2,1)-0.1,Fo_G(3,1)-0.1,'G');
%         grid on;axis equal;
%         axis([-0.1 1.3 -0.5 0.5 -0.1 1.5]);
%         xlabel('x [m]');
%         ylabel('y [m]');
%         zlabel('z [m]');
%         drawnow
%         hold off;
%         end
%     end
% end    
    
    
    
%% Try plotting the rotation

figure('rend','painters','pos',[pp{3,1} 450 200]); 
    ha = tight_subplot(1,3,[.08 -0.2],[-0.05 0],[-0.1 -0.1]);  %[gap_h gap_w] [lower upper] [left right]
    ax = gca;
    axes(ha(1));    
    [Xding,Yding,Zding] = sphere(20);   % draw the pi-ball
    hSurface = surf(Xding,Yding,Zding);hold on;
    set(hSurface,'FaceColor',[0.9 0.9 0.9],'FaceAlpha',0.3,'FaceLighting','gouraud','EdgeColor','none','LineStyle','none','BackFaceLighting','unlit',...
                 'AmbientStrengt',0.3,'DiffuseStrength',0.7,'SpecularStrength',0.5,'SpecularExponent',5,'SpecularColorReflectance',1);
             
    plot3(squeeze(R(1,1,1)),squeeze(R(2,1,1)),squeeze(R(3,1,1)),'*','MarkerSize',10,'color','k','LineWidth', 2);
    plot3(squeeze(R(1,1,2:1001)),squeeze(R(2,1,2:1001)),squeeze(R(3,1,2:1001)),'color',[0.9290 0.6940 0.1250],'linewidth',1.2);
    plot3(squeeze(est_R(1,1,1:201)),squeeze(est_R(2,1,1:201)),squeeze(est_R(3,1,1:201)),'color',[0.8500 0.3250 0.0980],'linewidth',1.2);
    plot3(squeeze(R_noise(1,1,1:201)),squeeze(R_noise(2,1,1:201)),squeeze(R_noise(3,1,1:201)),'color',[0 86 140]/255,'linewidth',1.2); 
    
    view(128,31) %Set the initial viewpoint
    light('Position',[-1 1 1])    
    axis vis3d %Allow to rotate without chaning size
    view(-122,31)
    axis off
    text(0,0,-1.5,'x')


    axes(ha(2));    
    [Xding,Yding,Zding] = sphere(20);   % draw the pi-ball
    hSurface = surf(Xding,Yding,Zding);hold on;
    set(hSurface,'FaceColor',[0.9 0.9 0.9],'FaceAlpha',0.3,'FaceLighting','gouraud','EdgeColor','none','LineStyle','none','BackFaceLighting','unlit',...
                 'AmbientStrengt',0.3,'DiffuseStrength',0.7,'SpecularStrength',0.5,'SpecularExponent',5,'SpecularColorReflectance',1);
    plot3(squeeze(R(1,2,1)),squeeze(R(2,2,1)),squeeze(R(3,2,1)),'*','MarkerSize',10,'color','k','LineWidth', 2);
    plot3(squeeze(R(1,2,2:1001)),squeeze(R(2,2,2:1001)),squeeze(R(3,2,2:1001)),'color',[0.9290 0.6940 0.1250],'linewidth',1.2);
    plot3(squeeze(est_R(1,2,1:201)),squeeze(est_R(2,2,1:201)),squeeze(est_R(3,2,1:201)),'color',[0.8500 0.3250 0.0980],'linewidth',1.2);
    plot3(squeeze(R_noise(1,2,1:201)),squeeze(R_noise(2,2,1:201)),squeeze(R_noise(3,2,1:201)),'color',[0 86 140]/255,'linewidth',1.2); 
    view(128,32) %Set the initial viewpoint
    light('Position',[-1 1 1])    
    axis vis3d %Allow to rotate without chaning size
    view(90,32)
    axis off
    text(0,0,-1.5,'y')

    axes(ha(3));    
    [Xding,Yding,Zding] = sphere(20);   % draw the pi-ball
    hSurface = surf(Xding,Yding,Zding);hold on;
    set(hSurface,'FaceColor',[0.9 0.9 0.9],'FaceAlpha',0.3,'FaceLighting','gouraud','EdgeColor','none','LineStyle','none','BackFaceLighting','unlit',...
                 'AmbientStrengt',0.3,'DiffuseStrength',0.7,'SpecularStrength',1,'SpecularExponent',5,'SpecularColorReflectance',0);
    plot3(squeeze(R(1,3,1)),squeeze(R(2,3,1)),squeeze(R(3,3,1)),'*','MarkerSize',10,'color','k','LineWidth', 2);
    g1 = plot3(squeeze(R(1,3,2:1001)),squeeze(R(2,3,2:1001)),squeeze(R(3,3,2:1001)),'color',[0.9290 0.6940 0.1250],'linewidth',1.2);
    g2 = plot3(squeeze(est_R(1,3,1:201)),squeeze(est_R(2,3,1:201)),squeeze(est_R(3,3,1:201)),'color',[0.8500 0.3250 0.0980],'linewidth',1.2);
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
% close all;
%plot angular velocity (\omega_x)
figure('rend','painters','pos',[pp{2,3} 2*sizex 0.8*sizey]);
    ha = tight_subplot(1,3,[.05 .04],[.18 .26],[0.06 0.03]);  %[gap_h gap_w] [lower upper] [left right] 
    axes(ha(1));
    g1=plot(t2,omg_FD(1,:),'color',[0 0.4470 0.7410 0.6]); hold on; grid on
    g2=plot(t3,est_omg(1,:),'linewidth',1.5);
    g3=plot(t1,omg(1,:),'linewidth',1.5); 
    xlim([0,tend]);
    xlabel('Time [s]');
    ylim([-20 20]);
    ylabel('Angular velocity [rad/s]');
    t=text(0.5,0.5,'x-component','parent',ha(1),'Fontsize',9); 
    t.Position = [ha(1).XLim(1)+0.5*(abs(ha(1).XLim(1))+abs(ha(1).XLim(2)))-0.5*t.Extent(3) ha(1).YLim(1)+1.1*(abs(ha(1).YLim(1))+abs(ha(1).YLim(2)))];
    
    axes(ha(2));
    plot(t2,omg_FD(2,:),'color',[0 0.4470 0.7410 0.6]); hold on; grid on
    plot(t3,est_omg(2,:),'linewidth',1.5);
    plot(t1,omg(2,:),'linewidth',1.5); 
    xlim([0,tend]);
    ylim([-20 20]);
    xlabel('Time [s]');
    t=text(0.5,0.5,'y-component','parent',ha(2),'Fontsize',9); 
    t.Position = [ha(2).XLim(1)+0.5*(abs(ha(2).XLim(1))+abs(ha(2).XLim(2)))-0.5*t.Extent(3) ha(2).YLim(1)+1.1*(abs(ha(2).YLim(1))+abs(ha(2).YLim(2)))];
    
    axes(ha(3));
    plot(t2,omg_FD(3,:),'color',[0 0.4470 0.7410 0.6]); hold on; grid on
    plot(t3,est_omg(3,:),'linewidth',1.5);
    plot(t1,omg(3,:),'linewidth',1.5); 
    xlim([0,tend]);
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
% close all;
%plot angular acceleration (\omega_x)
figure('rend','painters','pos',[pp{3,3} 2*sizex 0.8*sizey]);
    ha = tight_subplot(1,3,[.05 .04],[.18 .26],[0.06 0.03]);  %[gap_h gap_w] [lower upper] [left right] 
    axes(ha(1));
    g1=plot(t2,domg_FD(1,:),'color',[0 0.4470 0.7410 0.6]); hold on; grid on
    g2=plot(t3,est_domg(1,:),'linewidth',1.5);
    g3=plot(t1,domg(1,:),'linewidth',1.5); 
    xlim([0,tend]);
    ylim([-150,150]);
    yticks([-150 -100 -50 0 50 100 150])
    yticklabels({'-150','-100','-50','0','50','100','150'})
    xlabel('Time [s]');
    ylabel('Angular acceleration [rad/s$^2$]');
    t=text(0.5,0.5,'x-component','parent',ha(1),'Fontsize',9); 
    t.Position = [ha(1).XLim(1)+0.5*(abs(ha(1).XLim(1))+abs(ha(1).XLim(2)))-0.5*t.Extent(3) ha(1).YLim(1)+1.1*(abs(ha(1).YLim(1))+abs(ha(1).YLim(2)))];
    
    axes(ha(2));
    plot(t2,domg_FD(2,:),'color',[0 0.4470 0.7410 0.6]); hold on; grid on
    plot(t3,est_domg(2,:),'linewidth',1.5);
    plot(t1,domg(2,:),'linewidth',1.5); 
    xlim([0,tend]);
    ylim([-300,300]);
    yticks([-300 -200 -100 0 100 200 300])
    yticklabels({'-300','-200','-100','0','100','200','300'})
    xlabel('Time [s]');
    t=text(0.5,0.5,'y-component','parent',ha(2),'Fontsize',9); 
    t.Position = [ha(2).XLim(1)+0.5*(abs(ha(2).XLim(1))+abs(ha(2).XLim(2)))-0.5*t.Extent(3) ha(2).YLim(1)+1.1*(abs(ha(2).YLim(1))+abs(ha(2).YLim(2)))];
    
    axes(ha(3));
    plot(t2,domg_FD(3,:),'color',[0 0.4470 0.7410 0.6]); hold on; grid on
    plot(t3,est_domg(3,:),'linewidth',1.5);
    plot(t1,domg(3,:),'linewidth',1.5); 
    xlim([0,tend]);
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
    
%% Compute the error and plot the results
%Time indices of R for which we have a measurement:
% close all;
tR1 = find(ismember(t1,t2)==1);
tR2 = find(ismember(t1,t3)==1);

for ii = 1:length(tR1)
eR_meas(:,:,ii) = logm(R(:,:,tR1(ii))\R_noise(:,:,ii));
NeR_meas(ii) = norm(eR_meas(:,:,ii));
eomg_FD(:,ii) = omg_FD(:,ii)-omg(:,tR1(ii));
edomg_FD(:,ii) = domg_FD(:,ii)-domg(:,tR1(ii)); 
end


for ii = 1:length(tR2)
eR_est(:,:,ii) = logm(R(:,:,tR2(ii))\est_R(:,:,ii));
NeR_est(ii) = norm(eR_est(:,:,ii));
eomg_est(:,ii) = est_omg(:,ii)-omg(:,tR2(ii));
edomg_est(:,ii) = est_domg(:,ii)-domg(:,tR2(ii));
end

Eomg_FD = vecnorm(eomg_FD);
Eomg_est= vecnorm(eomg_est);
Edomg_FD = vecnorm(edomg_FD);
Edomg_est = vecnorm(edomg_est);

%Mean error in rotation
mean_ER_est = mean(NeR_est)
mean_ER_meas = mean(NeR_meas);

%Mean errors in velocity
mean_Eomg_FD = mean(Eomg_FD,'omitnan');
mean_Eomg_est = mean(Eomg_est,'omitnan')

%Mean errors in acceleration
mean_Edomg_FD = mean(Edomg_FD,'omitnan');
mean_Edomg_est = mean(Edomg_est,'omitnan')



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
    
%%
% close all;
figure('rend','painters','pos',[pp{1,4} sizex 2*sizey]);
    ha = tight_subplot(3,1,[.08 .07],[.08 .07],[0.1 0.03]);  %[gap_h gap_w] [lower upper] [left right]
    axes(ha(1));
    plot(t2,eomg_FD(1,:));hold on; grid on;
    plot(t3,eomg_est(1,:));    
    xlim([0 2]);
    ylim([-4 4]);
    xlabel('Time [s]');
    ylabel('Velocity error [rad/s]');

    axes(ha(2));
    plot(t2,eomg_FD(2,:));hold on; grid on;
    plot(t3,eomg_est(2,:));
    xlim([0 2]);
    ylim([-4 4]);
    xlabel('Time [s]');
    ylabel('Velocity error [rad/s]');
    
    axes(ha(3));
    plot(t2,eomg_FD(3,:));hold on; grid on;
    plot(t3,eomg_est(3,:));
    xlim([0 2]);
    ylim([-4 4]);
    xlabel('Time [s]');
    ylabel('Velocity error [rad/s]');    
    
    L1 = legend('$e_{\widetilde{\mathbf{R}}}$','$e_{\widehat{\mathbf{R}}}$','NumColumns',2);
    L1.Position(2) = 0.95;
    L1.Position(1) = 0.5-(L1.Position(3)/2);
    L1.FontSize = 9;


