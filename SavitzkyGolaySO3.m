clearvars; clc; close all; addpath('functions');
set(groot,'defaulttextinterpreter','latex'); set(groot,'defaultAxesTickLabelInterpreter','latex'); set(groot,'defaultLegendInterpreter','latex');
%% ---------------- Savitzky-Golay Filtering on SO(3) ----------------- %%
%% Constants and settings
%User inputs
doSave = false;     %Boolean: set true if you want to save figures
doPlot = false;    %Determine if you want to plot the trajectory in 3D plot
Fs = 1000;         %Sampling frequency fine grid      [Hz]
Fc = 1;            %Signal frequency                  [Hz]
a  = 2;            %Signal amplitude                  [deg]
n  = 13;           %Window size SG-filter             [-]
te = 10;           %Signal length                     [s]
order = 4;         %Savitzky Golay filter order       [-]
alpha = 5;         %Down-sampling rate                [-]
Anoise = 0.02;     %Standard deviation of added noise
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
    g(:,ii) = x0 + x1*a*sin(freq*t1(ii)); %+ x2*a*sin(0.1*freq*t1(ii));
    dg(:,ii) = x1*a*(freq)*cos(freq*t1(ii)); %+x2*a*0.1*freq*cos(0.1*freq*t1(ii));
    ddg(:,ii) = -x1*a*(freq)^2*sin(freq*t1(ii)); %-x2*a*(0.1*freq)^2*sin(freq*t1(ii));
      
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
g_noise(:,tel) = g(:,ii)+Anoise*randn(3,1);
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
   
%plot angular velocity (\omega_x)
figure('rend','painters','pos',[pp{1,1} sizex sizey]);
    ha = tight_subplot(1,1,[.08 .07],[.18 .1],[0.12 0.03]);  %[gap_h gap_w] [lower upper] [left right]
    axes(ha(1));
    plot(t2,omg_FD(1,:)); hold on; grid on
    plot(t3,est_omg(1,:));
    plot(t1,omg(1,:)); 
    xlim([0,2]);
%     ylim([-45,30]);
    xlabel('Time [s]');
    ylabel('Angular velocity [rad/s]');
    legend('$(^F$\boldmath$\tilde{{\omega}}_{F,G})_x$','$(^F$\boldmath$\hat{{\omega}}_{F,G})_x$','$(^F$\boldmath${{\omega}}_{F,G})_x$','location','southeast')
    if doSave
        fig = gcf;
        fig.PaperPositionMode = 'auto';
        fig_pos = fig.PaperPosition;
        fig.PaperSize = [fig_pos(3) fig_pos(4)];
        print(fig,'figures/omg.pdf','-dpdf','-painters')
    end
    
%plot angular velocity (\omega_x)
figure('rend','painters','pos',[pp{1,1} sizex sizey]);
    ha = tight_subplot(1,1,[.08 .07],[.18 .1],[0.12 0.03]);  %[gap_h gap_w] [lower upper] [left right]
    axes(ha(1));
    plot(t2,omg_FD(1,:)); hold on; grid on
    plot(t3,est_omg(1,:));
    plot(t1,omg(1,:)); 
    xlim([0,2]);
%     ylim([-45,30]);
    xlabel('Time [s]');
    ylabel('Angular velocity [rad/s]');
    legend('$(^F$\boldmath$\tilde{{\omega}}_{F,G})_x$','$(^F$\boldmath$\hat{{\omega}}_{F,G})_x$','$(^F$\boldmath${{\omega}}_{F,G})_x$','location','southeast')
    if doSave
        fig = gcf;
        fig.PaperPositionMode = 'auto';
        fig_pos = fig.PaperPosition;
        fig.PaperSize = [fig_pos(3) fig_pos(4)];
        print(fig,'figures/omg.pdf','-dpdf','-painters')
    end    
    
%plot angular acceleration (\dot\omega_x)
figure('rend','painters','pos',[pp{2,1} sizex sizey]);
    ha = tight_subplot(1,1,[.08 .07],[.18 .1],[0.15 0.03]);  %[gap_h gap_w] [lower upper] [left right]
    axes(ha(1));
    plot(t2,domg_FD(1,:)); hold on; grid on
    plot(t3,est_domg(1,:));
    plot(t1,domg(1,:));
%     ylim([-2000,4000]);
    xlim([0,2]);
    xlabel('Time [s]');
    ylabel('Angular acceleration [rad/$s^2$]');
    legend('$(^F$\boldmath$\tilde{\dot{\omega}}_{F,G})_x$',...
        '$(^F$\boldmath$\hat{\dot{\omega}}_{F,G})_x$',...
        '$(^F$\boldmath${\dot{\omega}}_{F,G})_x$',...
        'location','northwest')
    if doSave
        fig = gcf;
        fig.PaperPositionMode = 'auto';
        fig_pos = fig.PaperPosition;
        fig.PaperSize = [fig_pos(3) fig_pos(4)];
        print(fig,'figures/domg.pdf','-dpdf','-painters')
    end
    
%plot angular acceleration (\dot\omega_x) close up
figure('rend','painters','pos',[pp{2,2} sizex sizey]);
    ha = tight_subplot(1,1,[.08 .07],[.18 .1],[0.12 0.03]);  %[gap_h gap_w] [lower upper] [left right]
    axes(ha(1));
    plot(t2,domg_FD(1,:)); hold on; grid on
    plot(t3,est_domg(1,:));
    plot(t1,domg(1,:));
    xlim([0.7,1.4]);
%     ylim([-300,400]);
    xlabel('Time [s]');
    ylabel('Angular acceleration [rad/$s^2$]');
    legend('$(^F$\boldmath$\tilde{\dot{\omega}}_{F,G})_x$','$(^F$\boldmath$\hat{\dot{\omega}}_{F,G})_x$','$(^F$\boldmath${\dot{\omega}}_{F,G})_x$','location','northeast')
    if doSave
        fig = gcf;
        fig.PaperPositionMode = 'auto';
        fig_pos = fig.PaperPosition;
        fig.PaperSize = [fig_pos(3) fig_pos(4)];
        print(fig,'figures/domg_zoom.pdf','-dpdf','-painters')
    end

    
%Plot trajectory in 3D plot
if doPlot
    figure('rend','painters','pos',[pp{1,3} sizex sizey]);
    ha = tight_subplot(1,1,[.08 .07],[.18 .1],[0.12 0.03]);  %[gap_h gap_w] [lower upper] [left right]
    axes(ha(1));
    tel1 = 1;
    plotBool = true;
    for ii = (n*alpha+1):3*alpha:N1
        if plotBool
        %Plot the origin of the world coordinate frame
        tip = [0.3*[1;0;0] 0.3*[0;1;0] 0.3*[0;0;1]];
        plot3([0 tip(1,1)],[0 tip(2,1)],[0 tip(3,1)],'r'); hold on
        plot3([0 tip(1,2)],[0 tip(2,2)],[0 tip(3,2)],'g');
        plot3([0 tip(1,3)],[0 tip(2,3)],[0 tip(3,3)],'b');
        
        %Plot the origin of the robot joint with its unit vectors
        tip = [Fo_G(:,ii)+0.3*R(:,1,ii) Fo_G(:,ii)+0.3*R(:,2,ii) Fo_G(:,ii)+0.3*R(:,3,ii)];
        plot3([Fo_G(1,ii) tip(1,1)],[Fo_G(2,ii) tip(2,1)],[Fo_G(3,ii) tip(3,1)],'r'); hold on
        plot3([Fo_G(1,ii) tip(1,2)],[Fo_G(2,ii) tip(2,2)],[Fo_G(3,ii) tip(3,2)],'g');
        plot3([Fo_G(1,ii) tip(1,3)],[Fo_G(2,ii) tip(2,3)],[Fo_G(3,ii) tip(3,3)],'b');
        
        %Plot the origin of the end effector with its unit vectors
        tip = [Fo_E(:,ii)+0.3*FR_E(:,1,ii) Fo_E(:,ii)+0.3*FR_E(:,2,ii) Fo_E(:,ii)+0.3*FR_E(:,3,ii)];
        plot3([Fo_E(1,ii) tip(1,1)],[Fo_E(2,ii) tip(2,1)],[Fo_E(3,ii) tip(3,1)],'r'); hold on
        plot3([Fo_E(1,ii) tip(1,2)],[Fo_E(2,ii) tip(2,2)],[Fo_E(3,ii) tip(3,2)],'g');
        plot3([Fo_E(1,ii) tip(1,3)],[Fo_E(2,ii) tip(2,3)],[Fo_E(3,ii) tip(3,3)],'b');
        
        %Plot the noisy "measured" coordinate frame
        tip = [Fo_G(:,ii)+0.3*est_R(:,1,tel1) Fo_G(:,ii)+0.3*est_R(:,2,tel1) Fo_G(:,ii)+0.3*est_R(:,3,tel1)];
        plot3([Fo_G(1,ii) tip(1,1)],[Fo_G(2,ii) tip(2,1)],[Fo_G(3,ii) tip(3,1)],'r'); hold on
        plot3([Fo_G(1,ii) tip(1,2)],[Fo_G(2,ii) tip(2,2)],[Fo_G(3,ii) tip(3,2)],'g');
        plot3([Fo_G(1,ii) tip(1,3)],[Fo_G(2,ii) tip(2,3)],[Fo_G(3,ii) tip(3,3)],'b');
        tel1=tel1+3*1;
        
        if tel1 >= length(t3)
            plotBool = false;
        end
        
        text(0,0-0.1,0-0.1,'F');
        text(Fo_E(1,ii),Fo_E(2,ii)-0.1,Fo_E(3,ii)-0.1,'E');
        text(Fo_G(1,1),Fo_G(2,1)-0.1,Fo_G(3,1)-0.1,'G');
        grid on;axis equal;
        axis([-0.1 1.3 -0.5 0.5 -0.1 1.5]);
        xlabel('x [m]');
        ylabel('y [m]');
        zlabel('z [m]');
        drawnow
        hold off;
        end
    end
end    
    
    
    
% Try plotting the rotation

figure('rend','painters','pos',[pp{3,1} 450 200]); 
    ha = tight_subplot(1,3,[.08 -0.2],[-0.05 0],[-0.1 -0.1]);  %[gap_h gap_w] [lower upper] [left right]
    axes(ha(1));    
    [Xding,Yding,Zding] = sphere(20);   % draw the pi-ball
    hSurface = surf(Xding,Yding,Zding);hold on;
    set(hSurface,'FaceColor',[0.9 0.9 0.9],'FaceAlpha',0.3,'FaceLighting','gouraud','EdgeColor','none','LineStyle','none','BackFaceLighting','unlit',...
                 'AmbientStrengt',0.3,'DiffuseStrength',0.7,'SpecularStrength',0.5,'SpecularExponent',5,'SpecularColorReflectance',1);
             
    plot3(squeeze(R(1,1,1)),squeeze(R(2,1,1)),squeeze(R(3,1,1)),'*','MarkerSize',10,'color','k');
    plot3(squeeze(R(1,1,2:1001)),squeeze(R(2,1,2:1001)),squeeze(R(3,1,2:1001)),'color','k','linewidth',2.5);
    plot3(squeeze(est_R(1,1,1:201)),squeeze(est_R(2,1,1:201)),squeeze(est_R(3,1,1:201)),'color',[0.9290 0.6940 0.1250],'linewidth',1.2);
    plot3(squeeze(R_noise(1,1,1:201)),squeeze(R_noise(2,1,1:201)),squeeze(R_noise(3,1,1:201)),'color','b','linewidth',1.2); 

    view(128,31)
    light('Position',[-1 1 1])
    axis equal
    axis off
    text(0,0,-1.5,'x')


    axes(ha(2));    
    [Xding,Yding,Zding] = sphere(20);   % draw the pi-ball
    hSurface = surf(Xding,Yding,Zding);hold on;
    set(hSurface,'FaceColor',[0.9 0.9 0.9],'FaceAlpha',0.3,'FaceLighting','gouraud','EdgeColor','none','LineStyle','none','BackFaceLighting','unlit',...
                 'AmbientStrengt',0.3,'DiffuseStrength',0.7,'SpecularStrength',0.5,'SpecularExponent',5,'SpecularColorReflectance',1);

    plot3(squeeze(R(1,2,1)),squeeze(R(2,2,1)),squeeze(R(3,2,1)),'*','MarkerSize',10,'color','k');
    plot3(squeeze(R(1,2,2:1001)),squeeze(R(2,2,2:1001)),squeeze(R(3,2,2:1001)),'color','k','linewidth',2.5);
    plot3(squeeze(est_R(1,2,1:201)),squeeze(est_R(2,2,1:201)),squeeze(est_R(3,2,1:201)),'color',[0.9290 0.6940 0.1250],'linewidth',1.2);
    plot3(squeeze(R_noise(1,2,1:201)),squeeze(R_noise(2,2,1:201)),squeeze(R_noise(3,2,1:201)),'color','b','linewidth',1.2); 
    
%     plot3(GTR2(1,1),GTR2(2,1),GTR2(3,1),'*','MarkerSize',7,'color','k');
%     g1= plot3(GTR2(1,1:5),GTR2(2,1:5),GTR2(3,1:5),'color','k','linewidth',2);
%     plot3(GTR2(1,5:26),GTR2(2,5:26),GTR2(3,5:26),'color',[0 0 0 0.3],'linewidth',2);
%     plot3(GTR2(1,26:65),GTR2(2,26:65),GTR2(3,26:65),'color','k','linewidth',2);
%         
%     g2 =plot3(R2{ple(2)}(1,1:5),R2{ple(2)}(2,1:5),R2{ple(2)}(3,1:5),'-','color',[0.9290 0.6940 0.1250],'linewidth',1.2);
%     plot3(R2{ple(2)}(1,5:26),R2{ple(2)}(2,5:26),R2{ple(2)}(3,5:26),'-','color',[0.9290 0.6940 0.1250 0.5],'linewidth',1.2);
%     plot3(R2{ple(2)}(1,26:65),R2{ple(2)}(2,26:65),R2{ple(2)}(3,26:65),'-','color',[0.9290 0.6940 0.1250],'linewidth',1.2);
%     
%     g3=plot3(R2{ple(3)}(1,1:5),R2{ple(3)}(2,1:5),R2{ple(3)}(3,1:5),'--','color','b','linewidth',1.2);
%     plot3(R2{ple(3)}(1,5:26),R2{ple(3)}(2,5:26),R2{ple(3)}(3,5:26),'--','color',[0 0 1 0.5],'linewidth',1.2);
%     plot3(R2{ple(3)}(1,26:65),R2{ple(3)}(2,26:65),R2{ple(3)}(3,26:65),'--','color','b','linewidth',1.2);
%     
%     g4=plot3(R2{ple(4)}(1,1:5),R2{ple(4)}(2,1:5),R2{ple(4)}(3,1:5),'--','color','r','linewidth',1.2);
%     plot3(R2{ple(4)}(1,5:26),R2{ple(4)}(2,5:26),R2{ple(4)}(3,5:26),'--','color',[1 0 0 0.5],'linewidth',1.2);
%     plot3(R2{ple(4)}(1,26:65),R2{ple(4)}(2,26:65),R2{ple(4)}(3,26:65),'--','color','r','linewidth',1.2);
    view(128,31)
    light('Position',[-1 1 1])
    axis equal
    axis off
    text(0,0,-1.5,'y')

    axes(ha(3));    
    [Xding,Yding,Zding] = sphere(20);   % draw the pi-ball
    hSurface = surf(Xding,Yding,Zding);hold on;
    set(hSurface,'FaceColor',[0.9 0.9 0.9],'FaceAlpha',0.3,'FaceLighting','gouraud','EdgeColor','none','LineStyle','none','BackFaceLighting','unlit',...
                 'AmbientStrengt',0.3,'DiffuseStrength',0.7,'SpecularStrength',0.5,'SpecularExponent',5,'SpecularColorReflectance',1);
    plot3(squeeze(R(1,3,1)),squeeze(R(2,3,1)),squeeze(R(3,3,1)),'*','MarkerSize',10,'color','k');
    g1 = plot3(squeeze(R(1,3,2:1001)),squeeze(R(2,3,2:1001)),squeeze(R(3,3,2:1001)),'color','k','linewidth',2.5);
    g2 = plot3(squeeze(est_R(1,3,1:201)),squeeze(est_R(2,3,1:201)),squeeze(est_R(3,3,1:201)),'color',[0.9290 0.6940 0.1250],'linewidth',1.2);
    g3 = plot3(squeeze(R_noise(1,3,1:201)),squeeze(R_noise(2,3,1:201)),squeeze(R_noise(3,3,1:201)),'color','b','linewidth',1.2); 
    view(128,31)
    light('Position',[-1 1 1])
    axis equal
    axis off

    L1 = legend([g1 g2 g3],{'$\mathbf{R}$','$\hat{\mathbf{R}}$','$\tilde{\mathbf{R}}$'},'NumColumns',3,'location','northeast');
    L1.Position(2) = 0.88;
    L1.Position(1) = 0.5-(L1.Position(3)/2);
    L1.FontSize = 6;
    text(0,0,-1.5,'z')

    
    if doSave
        fig = gcf;
        fig.PaperPositionMode = 'auto';
        fig_pos = fig.PaperPosition;
        fig.PaperSize = [fig_pos(3) fig_pos(4)];
        print(fig,'figures/Rotation.pdf','-dpdf','-painters')
    end
    
    
%% Plot velocity in 1 plot
% close all;
%plot angular velocity (\omega_x)
figure('rend','painters','pos',[pp{2,3} 2*sizex 0.8*sizey]);
    ha = tight_subplot(1,3,[.05 .04],[.18 .26],[0.06 0.03]);  %[gap_h gap_w] [lower upper] [left right] 
    axes(ha(1));
    g1=plot(t2,omg_FD(1,:)); hold on; grid on
    g2=plot(t3,est_omg(1,:));
    g3=plot(t1,omg(1,:)); 
    xlim([0,tend]);
    xlabel('Time [s]');
    ylabel('Angular velocity [rad/s]');
    t=text(0.5,0.5,'x-component','parent',ha(1),'Fontsize',9); 
    t.Position = [ha(1).XLim(1)+0.5*(abs(ha(1).XLim(1))+abs(ha(1).XLim(2)))-0.5*t.Extent(3) ha(1).YLim(1)+1.1*(abs(ha(1).YLim(1))+abs(ha(1).YLim(2)))];
    
    axes(ha(2));
    plot(t2,omg_FD(2,:)); hold on; grid on
    plot(t3,est_omg(2,:));
    plot(t1,omg(2,:)); 
    xlim([0,tend]);
    xlabel('Time [s]');
    t=text(0.5,0.5,'y-component','parent',ha(2),'Fontsize',9); 
    t.Position = [ha(2).XLim(1)+0.5*(abs(ha(2).XLim(1))+abs(ha(2).XLim(2)))-0.5*t.Extent(3) ha(2).YLim(1)+1.1*(abs(ha(2).YLim(1))+abs(ha(2).YLim(2)))];
    
    axes(ha(3));
    plot(t2,omg_FD(3,:)); hold on; grid on
    plot(t3,est_omg(3,:));
    plot(t1,omg(3,:)); 
    xlim([0,tend]);
    xlabel('Time [s]');
    t=text(0.5,0.5,'z-component','parent',ha(3),'Fontsize',9); 
    t.Position = [ha(3).XLim(1)+0.5*(abs(ha(3).XLim(1))+abs(ha(3).XLim(2)))-0.5*t.Extent(3) ha(3).YLim(1)+1.1*(abs(ha(3).YLim(1))+abs(ha(3).YLim(2)))];
    
    L1 = legend([g3 g1 g2],{'Analytical solution \boldmath${{\omega}}$','Finite differencing \boldmath$\breve{{\omega}}$',...
        'Savitzky-Golay \boldmath$\hat{{\omega}}$'},'NumColumns',3,'location','northeast');
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
    
%% Plot acceleration in 1 plot
% close all;
%plot angular acceleration (\omega_x)
figure('rend','painters','pos',[pp{3,3} 2*sizex 0.8*sizey]);
    ha = tight_subplot(1,3,[.05 .04],[.18 .26],[0.06 0.03]);  %[gap_h gap_w] [lower upper] [left right] 
    axes(ha(1));
    g1=plot(t2,domg_FD(1,:)); hold on; grid on
    g2=plot(t3,est_domg(1,:));
    g3=plot(t1,domg(1,:)); 
    xlim([0,tend]);
    ylim([-600,600]);
    yticks([-600 -400 -200 0 200 400 600])
    yticklabels({'-600','-400','-200','0','200','400','600'})
    xlabel('Time [s]');
    ylabel('Angular acceleration [rad/s$^2$]');
    t=text(0.5,0.5,'x-component','parent',ha(1),'Fontsize',9); 
    t.Position = [ha(1).XLim(1)+0.5*(abs(ha(1).XLim(1))+abs(ha(1).XLim(2)))-0.5*t.Extent(3) ha(1).YLim(1)+1.1*(abs(ha(1).YLim(1))+abs(ha(1).YLim(2)))];
    
    axes(ha(2));
    plot(t2,domg_FD(2,:)); hold on; grid on
    plot(t3,est_domg(2,:));
    plot(t1,domg(2,:)); 
    xlim([0,tend]);
    ylim([-600,600]);
    yticks([-600 -400 -200 0 200 400 600])
    yticklabels({'-600','-400','-200','0','200','400','600'})
    xlabel('Time [s]');
    t=text(0.5,0.5,'y-component','parent',ha(2),'Fontsize',9); 
    t.Position = [ha(2).XLim(1)+0.5*(abs(ha(2).XLim(1))+abs(ha(2).XLim(2)))-0.5*t.Extent(3) ha(2).YLim(1)+1.1*(abs(ha(2).YLim(1))+abs(ha(2).YLim(2)))];
    
    axes(ha(3));
    plot(t2,domg_FD(3,:)); hold on; grid on
    plot(t3,est_domg(3,:));
    plot(t1,domg(3,:)); 
    xlim([0,tend]);
    ylim([-600,600]);
    yticks([-600 -400 -200 0 200 400 600])
    yticklabels({'-600','-400','-200','0','200','400','600'})
    xlabel('Time [s]');
    t=text(0.5,0.5,'z-component','parent',ha(3),'Fontsize',9); 
    t.Position = [ha(3).XLim(1)+0.5*(abs(ha(3).XLim(1))+abs(ha(3).XLim(2)))-0.5*t.Extent(3) ha(3).YLim(1)+1.1*(abs(ha(3).YLim(1))+abs(ha(3).YLim(2)))];
    
    L1 = legend([g3 g1 g2],{'Analytical solution \boldmath$\dot{{\omega}}$','Finite differencing \boldmath$\breve{\dot{{\omega}}}$',...
        'Savitzky-Golay \boldmath$\hat{\dot{{\omega}}}$'},'NumColumns',3,'location','northeast');
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
