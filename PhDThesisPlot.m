%Seperate script to save figures to PhD Thesis format
doSave = false;    %Boolean: set true if you want to save figures


%% Figures
%Check if figures directory exists, if not, it will create one.
if ~isfolder('figures')
    mkdir('figures');
end

%Create a plot grid
sizex = 470;
sizey = 250;
px = (0:7)*(sizex+10)+10;
py = (0:4)*(sizey+40)+45;
for  ii = 1:length(px)
    for jj = 1:length(py)
        pp{jj,ii} = [px(ii) py(jj)];
    end
end

%% Plot the orientation error
figure('rend','painters','pos',[pp{1,1} sizex 1.3*sizey]);
    ha = tight_subplot(1,1,[.08 .07],[.24 .08],[0.09 0.01]);  %[gap_h gap_w] [lower upper] [left right]
    axes(ha(1));
    plot(t2,NeR_meas); hold on; 
    plot(t3,NeR_est);
    xlim([0 2]);
    xlabel('Time [s]');
    ylabel('Orientation error [rad]');
%     L1 = legend('$e_{\widetilde{\mathbf{R}}}$','$e_{\widehat{\mathbf{R}}}$','NumColumns',2);
%     L1.FontSize = 9;
    grid on;
    set(findall(gcf,'-property','FontSize'),'FontSize',10)
    
    if doSave
        fig = gcf;
        fig.PaperPositionMode = 'auto';
        fig_pos = fig.PaperPosition;
        fig.PaperSize = [fig_pos(3) fig_pos(4)];
        print(fig,'figures/THESIS_SG_norm_eR.pdf','-dpdf','-painters')
    end
    
%% Try plotting the rotation
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
        print(fig,'figures/THESIS_SG_Rotation.pdf','-dpdf','-painters')
    end
    
    
%% Plot velocity in 1 plot
figure('rend','painters','pos',[pp{3,1} sizex 2.1*sizey]);
    ha = tight_subplot(3,1,[.03 .05],[.075 .02],[0.13 0.01]);  %[gap_h gap_w] [lower upper] [left right] 
    axes(ha(1));
    g1=plot(t2,omg_FD(1,:),'color',[0 0.4470 0.7410 0.6]); hold on; grid on
    g2=plot(t3,omg_est(1,:),'linewidth',1.5);
    g3=plot(t1,omg(1,:),'linewidth',1.5); 
    xlim([0,te]);
    xticks([0 0.5 1 1.5 2]);
    xticklabels({'','','','',''});
    ylim([-20 20]);
    yticks([-20 -10 0 10 20]);
    yticklabels({'-20','-10','0','10','20'});
    ylabel({'Angular velocity [rad/s]';'x-component'});
    
    axes(ha(2));
    plot(t2,omg_FD(2,:),'color',[0 0.4470 0.7410 0.6]); hold on; grid on
    plot(t3,omg_est(2,:),'linewidth',1.5);
    plot(t1,omg(2,:),'linewidth',1.5); 
    xlim([0,te]);
    xticks([0 0.5 1 1.5 2]);
    xticklabels({'','','','',''});
    ylim([-20 20]);
    yticks([-20 -10 0 10 20]);
    yticklabels({'-20','-10','0','10','20'});
    ylabel({'Angular velocity [rad/s]';'y-component'});
    
    axes(ha(3));
    plot(t2,omg_FD(3,:),'color',[0 0.4470 0.7410 0.6]); hold on; grid on
    plot(t3,omg_est(3,:),'linewidth',1.5);
    plot(t1,omg(3,:),'linewidth',1.5); 
    xlim([0,te]);
    xticks([0 0.5 1 1.5 2]);
    xticklabels({'0','0.5','1','1.5','2'});
    ylim([-20 20]);
    yticks([-20 -10 0 10 20]);
    yticklabels({'-20','-10','0','10','20'});
    ylabel({'Angular velocity [rad/s]';'z-component'});
    xlabel('Time [s]');  
    set(findall(gcf,'-property','FontSize'),'FontSize',10)
    
    if doSave
        fig = gcf;
        fig.PaperPositionMode = 'auto';
        fig_pos = fig.PaperPosition;
        fig.PaperSize = [fig_pos(3) fig_pos(4)];
        print(fig,'figures/THESIS_SG_omg.pdf','-dpdf','-painters')
    end
    
%% Plot acceleration in 1 plot
figure('rend','painters','pos',[pp{3,1} sizex 2.1*sizey]);
    ha = tight_subplot(3,1,[.03 .05],[.075 .02],[0.14 0.01]);  %[gap_h gap_w] [lower upper] [left right] 
    axes(ha(1));
    g1=plot(t2,domg_FD(1,:),'color',[0 0.4470 0.7410 0.6]); hold on; grid on
    g2=plot(t3,domg_est(1,:),'linewidth',1.5);
    g3=plot(t1,domg(1,:),'linewidth',1.5); 
    xlim([0,te]);
    ylim([-200,200]);
    yticks([-200 -150 -100 -50 0 50 100 150 200])
    yticklabels({'-200','-150','-100','-50','0','50','100','150','200'})
    ylabel({'Ang. acceleration [rad/s$^2$]';'x-component'});
    xticks([0 0.5 1 1.5 2]);
    xticklabels({'','','','',''});
    
    axes(ha(2));
    plot(t2,domg_FD(2,:),'color',[0 0.4470 0.7410 0.6]); hold on; grid on
    plot(t3,domg_est(2,:),'linewidth',1.5);
    plot(t1,domg(2,:),'linewidth',1.5); 
    xlim([0,te]);
    ylim([-200,200]);
    yticks([-200 -150 -100 -50 0 50 100 150 200])
    yticklabels({'-200','-150','-100','-50','0','50','100','150','200'})
    ylabel({'Ang. acceleration [rad/s$^2$]';'y-component'});
    xticks([0 0.5 1 1.5 2]);
    xticklabels({'','','','',''});

    axes(ha(3));
    plot(t2,domg_FD(3,:),'color',[0 0.4470 0.7410 0.6]); hold on; grid on
    plot(t3,domg_est(3,:),'linewidth',1.5);
    plot(t1,domg(3,:),'linewidth',1.5); 
    xlim([0,te]);
    ylim([-200,200]);
    yticks([-200 -150 -100 -50 0 50 100 150 200])
    yticklabels({'-200','-150','-100','-50','0','50','100','150','200'})
    ylabel({'Ang. acceleration [rad/s$^2$]';'z-component'});
    xticks([0 0.5 1 1.5 2]);
    xticklabels({'0','0.5','1','1.5','2'});
    xlabel('Time [s]');  
    set(findall(gcf,'-property','FontSize'),'FontSize',10)
    
    if doSave
        fig = gcf;
        fig.PaperPositionMode = 'auto';
        fig_pos = fig.PaperPosition;
        fig.PaperSize = [fig_pos(3) fig_pos(4)];
        print(fig,'figures/THESIS_SG_domg.pdf','-dpdf','-painters')
    end