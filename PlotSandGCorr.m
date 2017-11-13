LoadFromFile = false;

if (0) %LoadFromFile == true)

 Filename = '2017_1006-(15_47_41) - big_stat_fun - steps=4096,musteps=1,t=5,Gamma=(3,0.1),eps=(0,0.16).mat'
 %Filename = '2017_1006-(15_47_05) - big_stat_fun - steps=4096,musteps=1,t=5,Gamma=(3,0.16),eps=(0,0.16).mat'
 Filename = '2017_1006-(15_52_17) - big_stat_fun - steps=4096,musteps=1,t=5,Gamma=(3,0.1),eps=(0,0.16).mat'
 Filename = '2017_1006-(15_57_21) - big_stat_fun - steps=4096,musteps=1,t=5,Gamma=(3,0.2),eps=(0,0.16).mat'
 Filename = '2017_1006-(15_57_21) - big_stat_fun - steps=4096,musteps=1,t=5,Gamma=(3,0.2),eps=(0,0.16).mat'
 %Filename = '2017_1006-(16_02_37)-big_stat_fun-steps=4096,mu(1)=-3.34,musteps=1,t=5,Gamma=(3,0.1),eps=(0,0.16).mat'
 
 Filename =  '2017_1009-(13_54_15)-big_stat_fun-steps=1024,mu(1)=-6,musteps=20,t=5,Gamma=(0.3,0.1),eps=(0,0.16).mat'
 
  MatData = load(Filename);
 %   S_data = MatData.S_data;
 %   G_data = MatData.G_data;
    
    MuRange = MatData.MuRange;
    MuSteps = length(MuRange) %MatData.MuSteps;
    
    Gamma_avg = MatData.Gamma_avg;
    Gamma_std = MatData.Gamma_std;
    
    epsilon = MatData.epsilon;
    epsilon_avg = MatData.epsilon_avg;
    epsilon_std = MatData.epsilon_std;
    
    GammaRange = MatData.GammaRange;
    GammaL = MatData.GammaL;
    GammaR  = MatData.GammaR;
    
 %   G_data = MatData.G_data;
 %   S_data= MatData.S_data;
   
    CorrData = MatData.CorrData;
    
    BigS_data_con = MatData.BigS_data_con;
    BigS_data_des = MatData.BigS_data_des;
    BigG_data_con = MatData.BigG_data_con;
    BigG_data_des = MatData.BigG_data_des;
    
    Steps = MatData.Steps;
end

%% Gamma only run?
if (length(MuRange) == 1 && length(GammaRange) > 1)
  
      f=figure
    imgsize = get(f,'Position');
    f.Units = 'inches';
    f.set('Position',[3,3,5,5]);
    
    subplot(3,1,1)
  %  Line1 = -3.34
    plot(GammaRange,mean(BigG_data_con,3),'b-',GammaRange,mean(BigG_data_des,3),'r--')
    legend({['Without node'],['With node']},'interpreter','latex','Location','best')
    %xlabel('\mu - \mu_0 [eV]')
    ylabel('log_{10}(<G>/G_0)') % [G_0]')
    %title('Conductance')
    title(['Std(\Gamma)=',num2str(Gamma_std),'eV, <\Delta\mu>=',num2str(MuRange),'eV, Std(\Delta\mu)= ',num2str(epsilon_std),'eV'])
       set(gca,'XMinorTick','on','YMinorTick','on')
    
        % rectangle('Position',[Line1,1E-10,.1,1],'FaceColor',[.8 .8 .9],'EdgeColor','w')
         ylim([min(mean(BigG_data_des,3)),1])
         
         
    subplot(3,1,2)

    plot(GammaRange,mean(BigS_data_con,3),'b-',GammaRange,mean(BigS_data_des,3),'r--')
   % legend('No node','With node')
    %xlabel('\mu - \mu_0 [eV]')
    ylabel('<S> [\muV/K]') %,'interpreter','latex')
       set(gca,'XMinorTick','on','YMinorTick','on')
       Smax = max(mean(BigS_data_des,3));

%    rectangle('Position',[Line1,-Smax,.1,2*Smax],'FaceColor',[.8 .8 .9],'EdgeColor','w')
      ylim([-Smax,Smax])
    
    subplot(3,1,3)
   
    plot(GammaRange,CorrData(:,1,1),'b-',GammaRange,CorrData(:,1,2),'r--',GammaRange,zeros(GammaSteps),'k-')
    axis tight
    %line([0 0],[3,3])
   % legend('Without node','With node')
    ylabel('Pearson corr.')
    xlabel('$<\Gamma>$ [eV]','interpreter','latex')
 %      rectangle('Position',[Line1,-1,.1,2],'FaceColor',[.8 .8 .9],'EdgeColor','w')

   
    ylim([-1.0,1.0])          
    

    
    set(gca,'XMinorTick','on','YMinorTick','on')
    set(gcf,'PaperPositionMode','auto')
    
   % mu = MuRange(MuIndex)
    
%  szFilename = ['C:\Users\Justin Bergfield\Desktop\Papers\2016_0927_Statistical Fun\Figures\gamma_scan', ...
%      '_Gammasig_',num2str(Gamma_std),'_epsavg_',num2str(epsilon_avg),'_epssig_',num2str(epsilon_std),'.eps']; %
%   disp(szFilename)
%   
      
 szFilename = ['C:\Users\Justin Bergfield\Google Drive\Papers\2016_0927_Statistical Fun\Figures\gamma_run_mu_',num2str(MuRange(1)),'_epssig_',num2str(epsilon_std),'_Gammasig_',num2str(Gamma_std),'_three_row.eps']; %
  disp(szFilename)

% f.set('Position',[imgsize(1),imgsize(2)+10,4,2.5])

%print(f,['C:\Users\Justin Bergfield\Desktop\Papers\2014_1203_Josh_Paper1\Version 0.5\Figures\EL_g_',num2str(1000.0*hg_1V1),'_kappa_',num2str(1000.0*kappa),'_delta_',num2str((Egap-hwc)*1000),'.eps'],'-depsc2','-r300')

%%
  %  print(f,szFilename,'-depsc2','-r300')

% f.set('Position',[imgsize(1),imgsize(2)+10,4,2.5])

%print(f,['C:\Users\Justin Bergfield\Desktop\Papers\2014_1203_Josh_Paper1\Version 0.5\Figures\EL_g_',num2str(1000.0*hg_1V1),'_kappa_',num2str(1000.0*kappa),'_delta_',num2str((Egap-hwc)*1000),'.eps'],'-depsc2','-r300')

%%
%print(f,szFilename,'-depsc2','-r300')
    
return

%% Mu only run
elseif (length(MuRange) > 1 && length(GammaRange) == 1) 
    
    
    % Let's see if we can't plot the std of the distribution
    figure
    subplot(2,1,1)
        plot(MuRange,std(BigG_data_con,0,3),MuRange,std(BigG_data_des,0,3),MuRange,std(GammaL+GammaR)*ones(length(MuRange)),'k--')
    legend('No node: \sigma_G','Node: \sigma_G','Gamma')
    ylabel('Std. Deviation [G_0]')

    
    subplot(2,1,2)
    plot(MuRange,std(BigS_data_con,0,3),MuRange,std(BigS_data_des,0,3),MuRange,std(epsilon)*ones(length(MuRange)),'k--')
    legend('No node: \sigma_S','Node: \sigma_S','eps')
  %  ylim([0 5])
    ylabel('Std. Deviation [\mu V/K]')
    xlabel('<\Delta\mu> [eV]')
 %   return
    
    %% 3 row version
    f=figure
    
    if (IsMatlab == true)
      imgsize = get(f,'Position');
      f.Units = 'inches';
      f.set('Position',[3,3,5,5]);
    end
    
    subplot(3,1,1)
    DrawRect = false;
    Line1 = 0.0 %-3.34
    semilogy(MuRange,mean(BigG_data_con,3),'b-',MuRange,mean(BigG_data_des,3),'r--')
    legend({['Without node'],['With node']},'interpreter','latex','Location','best')
    %xlabel('\mu - \mu_0 [eV]')
    ylabel('log_{10}(<G>/G_0)') % [G_0]')
    %title('Conductance')
    title(['<\Gamma_\alpha>= ',num2str(Gamma_avg),'eV, Std(\Gamma_\alpha)=',num2str(Gamma_std),'eV, Std(\Delta\mu)= ',num2str(epsilon_std),'eV'])
       set(gca,'XMinorTick','on','YMinorTick','on')
    if (DrawRect)
         rectangle('Position',[Line1,1E-10,.1,1],'FaceColor',[.8 .8 .9],'EdgeColor','w')
    end
         
         ylim([min(mean(BigG_data_des,3)),1])
         
         
    subplot(3,1,2)

    plot(MuRange,mean(BigS_data_con,3),'b-',MuRange,mean(BigS_data_des,3),'r--')
   % legend('No node','With node')
    %xlabel('\mu - \mu_0 [eV]')
    ylabel('<S> [\muV/K]') %,'interpreter','latex')
       set(gca,'XMinorTick','on','YMinorTick','on')
       Smax = max(mean(BigS_data_des,3));
if (DrawRect)
    rectangle('Position',[Line1,-Smax,.1,2*Smax],'FaceColor',[.8 .8 .9],'EdgeColor','w')
end 
    ylim([-Smax,Smax])
    
    subplot(3,1,3)
   
    plot(MuRange,CorrData(1,:,1),'b-',MuRange,CorrData(1,:,2),'r--',MuRange,zeros(MuSteps),'k-')
    axis tight
    %line([0 0],[3,3])
   % legend('Without node','With node')
    ylabel('Pearson corr.')
    xlabel('$<\Delta\mu>$ [eV]','interpreter','latex')
   if (DrawRect)
    rectangle('Position',[Line1,-1,.1,2],'FaceColor',[.8 .8 .9],'EdgeColor','w')
   end
   % ylim([-1.0,1.0])    
   %figure
%        subplot(4,1,4)
%        plot(MuRange,squeeze(mean(Ldata_des(2,:,:,:),4)),MuRange,zeros(length(MuRange)),'k-')
%         
    

    
    set(gca,'XMinorTick','on','YMinorTick','on')
    set(gcf,'PaperPositionMode','auto')
    
   % mu = MuRange(MuIndex)
    
 szFilename = ['C:\Users\Justin Bergfield\Desktop\Papers\2016_0927_Statistical Fun\Figures\', ...
     'Gammaavg_',num2str(Gamma_avg),'_Gammasig_',num2str(Gamma_std),'_epsavg_',num2str(epsilon_avg),'_epssig_',num2str(epsilon_std),'.eps']; %
  disp(szFilename)
  
      
 szFilename = ['C:\Users\Justin Bergfield\Google Drive\Papers\2016_0927_Statistical Fun\Figures\trans_thermo_three_row.eps']; %
  disp(szFilename)

% f.set('Position',[imgsize(1),imgsize(2)+10,4,2.5])

%print(f,['C:\Users\Justin Bergfield\Desktop\Papers\2014_1203_Josh_Paper1\Version 0.5\Figures\EL_g_',num2str(1000.0*hg_1V1),'_kappa_',num2str(1000.0*kappa),'_delta_',num2str((Egap-hwc)*1000),'.eps'],'-depsc2','-r300')

%%
  %  print(f,szFilename,'-depsc2','-r300')

% f.set('Position',[imgsize(1),imgsize(2)+10,4,2.5])

%print(f,['C:\Users\Justin Bergfield\Desktop\Papers\2014_1203_Josh_Paper1\Version 0.5\Figures\EL_g_',num2str(1000.0*hg_1V1),'_kappa_',num2str(1000.0*kappa),'_delta_',num2str((Egap-hwc)*1000),'.eps'],'-depsc2','-r300')

%%
%print(f,szFilename,'-depsc2','-r300')
    
else
    
    disp('here!')
    MuIndex = 1;
  % [A MuIndex] = min(abs(MuRange+.5))
   
  %  if (1)% If we only have a single energy plot the histograms
    
  %%    CorrData
    f=figure;
     f.Units = 'inches';
    f.set('Position',[3,3,2.7,5]);
    
    %set(gcf,'interpreter','latex')
    subplot(2,1,1)
    GUpperBound = max(mean(BigG_data_con(MuIndex,:)),mean(BigG_data_des(MuIndex,:))) + 2.6*max(std(BigG_data_con(MuIndex,:)),std(BigG_data_des(MuIndex,:)));
    Gmax = max(max(BigG_data_con(MuIndex,:)),max(BigG_data_des(MuIndex,:)))
    GPlotRange = [0 GUpperBound];
    
     h2 = histogram(BigG_data_des(MuIndex,:));
    hold on
    h1 = histogram(BigG_data_con(MuIndex,:));
 
    h2.Normalization = 'probability';
    h1.Normalization = 'probability';
    Boxes = 50;
    h1.BinWidth = (GPlotRange(2)-GPlotRange(1))/Boxes;
    h2.BinWidth = (GPlotRange(2)-GPlotRange(1))/Boxes;
    xlim(GPlotRange)
    legend({['1,3BDT,node'],['1,4BDT,no node']},'interpreter','latex','fontsize',8.5,'Location','Best')
    hold off
    xlabel('G $[G_0]$','interpreter','latex')
       ylabel('Normalized Counts')
    %title(['$\mu - \mu_0$=',num2str(MuRange(MuIndex)),'eV'],'interpreter','latex')
     
    set(gca,'XMinorTick','on','YMinorTick','on')
    
    subplot(2,1,2)
       
    h1 = histogram(BigS_data_con(MuIndex,:));
    hold on
    SMax = max(max(abs(BigS_data_des(MuIndex,:))),max(abs(BigS_data_con(MuIndex,:)))); 
    
    Smax = 20; %mean([BigS_data_des BigS_data_con]) + 2.2*std([BigS_data_des, BigS_data_con]);
    Smin = -20; %mean([BigS_data_des, BigS_data_con]) - 2.2*std([BigS_data_des, BigS_data_con]);
    
    SUpperBound = max(abs(mean(BigS_data_des(MuIndex,:))), abs(mean(BigS_data_con(MuIndex,:)))) + 2.5*max(std(BigS_data_des(MuIndex,:)),2*std(BigS_data_con(MuIndex,:)));
   % UpperBound = max(max(BigS_data_des(MuIndex,:)),max(BigS_data_con(MuIndex,:)));
    SPlotRange = [-Smax,Smax];
    h2 = histogram(BigS_data_des(MuIndex,:));
   h1.Normalization = 'probability';
    h1.BinWidth = (SPlotRange(2)-SPlotRange(1))/Boxes;
    h2.Normalization = 'probability';
    h2.BinWidth = (SPlotRange(2)-SPlotRange(1))/Boxes;
    xlim(SPlotRange)
  %  legend({['No node'],['Node']}) %,'interpreter','latex','fontsize',12)
    hold off
    ylabel('Normalized Counts')
    xlabel('S [$\mu V/K$]','interpreter','latex')
    
    

    % write it to pdf
   % f.Units = 'inches';
   % imgsize = get(f,'Position');
   % f.set('Position',[imgsize(1),imgsize(2),4,4]);

 %   legend({['$V_b$=',num2str(Volt(3)),'V'],['$V_b$=',num2str(Volt(2)),'V'],['$V_b$=',num2str(Volt(1)),'V']},'Location','northwest','interpreter','latex','fontsize',12)
  %  xlabel({'Photon Energy [eV]';,''})
  %  ylabel('EL Intensity [arb. units]')
%axis tight


%set(gca,'XTick',[-40:10:40])

   %xlabel('Photon Energy [eV]')
 

   
%title(['$\hbar g_{00}$=',num2str(1000.0*FS.g_matrix(1,1+V)),'meV, $\hbar\kappa$=',num2str(1000.0*kappa),'meV, $\Delta$=',num2str((Egap-hwc)*1000),'meV'],'interpreter','latex','fontweight','bold','fontsize',13)
 %   set(gca, 'FontName', 'Times New Roman','fontsize',12)
    set(gca,'XMinorTick','on','YMinorTick','on')
   % text(1.86,.23,'x25')
   
    set(gcf,'PaperPositionMode','auto')
    
set(gcf, 'Color', 'w');
    mu = MuRange(MuIndex);
    
 szFilename = ['C:\Users\Justin Bergfield\Google Drive\Papers\2016_0927_Statistical Fun\Figures\hist_',num2str(mu),'.png']; %
  disp(szFilename)

  %export_fig szFilename -q101
  
% f.set('Position',[imgsize(1),imgsize(2)+10,4,2.5])

%print(f,['C:\Users\Justin Bergfield\Desktop\Papers\2014_1203_Josh_Paper1\Version 0.5\Figures\EL_g_',num2str(1000.0*hg_1V1),'_kappa_',num2str(1000.0*kappa),'_delta_',num2str((Egap-hwc)*1000),'.eps'],'-depsc2','-r300')

%print(f,szFilename,'-dpng','-r300')

%% 2d histogram
    f=figure;
    set(gca,'XMinorTick','on','YMinorTick','on');
   % text(1.86,.23,'x25')
    f.Units = 'inches';
    f.set('Position',[3,3,3,5]);
     f.set('Position',[3,3,2.7,5]);
    set(gcf,'PaperPositionMode','auto')
    
    subplot(2,1,2)
    
    % We want both paneles on the same axis.
    centers = {linspace(Smin,SMax,256),...
               linspace(0, Gmax, 128)};
           
           
          
     % Destructive data
    [PSG c] = hist3([BigS_data_des(MuIndex,:)',BigG_data_des(MuIndex,:)'],centers);
    %,centers);
    %hist3(data,centers)

    %PSG = PSG ./ sum(sum(PSG));
    imagesc(c{2},c{1},PSG)
    %colorbar
    ylabel('S [$\mu V/K$]','interpreter','latex')
     xlabel('G [$G_0$]','interpreter','latex')
  %  xlabel('G [2e^2/h]','interpreter','latex')
    xlim(GPlotRange);
    ylim(SPlotRange);
   ylim([Smin SMax])
  %  ylim(auto)
    %centers = {Srange,Grange};
    %[PSG c] = hist3(data,centers);
   % title('With Node','interpreter','latex')
    axis ij
    text( mean(GPlotRange),Smin+.1*(Smax-Smin),'1,3-BDT, node','HorizontalAlignment','center','Color','white','FontSize',12,'interpreter','latex')
    text( mean(GPlotRange),Smin+.5*(Smax-Smin),['$\rho_{SG}$=',num2str(CorrData(2),2)],'HorizontalAlignment','center','Color','white','FontSize',12,'interpreter','latex')
    
    
    %title('Node - (meta)BDT')
   %     c=colorbar
   % c.Label.String = 'Counts';
    
    
    subplot(2,1,1)
   % centers = {linspace(SPlotRange(1),SPlotRange(2),64),...
   %            linspace(GPlotRange(1),GPlotRange(2),64)};
    [PSG c] = hist3([BigS_data_con(MuIndex,:)',BigG_data_con(MuIndex,:)'],centers);
    %hist3(data,centers)

   % PSG = PSG ./ sum(sum(PSG));
    imagesc(c{2},c{1},PSG)
   % c=colorbar
   % c.Label.String = 'Counts';
    
    %clabel('Counts')
    xlim(GPlotRange);
    ylim(SPlotRange);
   ylim([Smin SMax])
    ylabel('S [$\mu V/K$]','interpreter','latex')
    xlabel('G [$G_0$]','interpreter','latex')
    %centers = {Srange,Grange};
    %[PSG c] = hist3(data,centers);
  %  title('Without Node','interpreter','latex')
    axis ij
    
   % text(mean(GPlotRange),-30,'Without Node','HorizontalAlignment','center','Color','white','FontSize',12,'interpreter','latex')
        text( mean(GPlotRange),Smin+.1*(Smax-Smin),'1,4-BDT, no node','HorizontalAlignment','center','Color','white','FontSize',12,'interpreter','latex')
            text( mean(GPlotRange),Smin+.5*(Smax-Smin),['$\rho_{SG}$=',num2str(CorrData(1),2)],'HorizontalAlignment','center','Color','white','FontSize',12,'interpreter','latex')
    
  %  title ('No node - (para)BDT')
    
     szFilename = ['C:\Users\Justin Bergfield\Google Drive\Papers\2016_0927_Statistical Fun\Figures\twod_hist_',num2str(mu),'.png']; %
     disp(szFilename)

% f.set('Position',[imgsize(1),imgsize(2)+10,4,2.5])

%print(f,['C:\Users\Justin Bergfield\Desktop\Papers\2014_1203_Josh_Paper1\Version 0.5\Figures\EL_g_',num2str(1000.0*hg_1V1),'_kappa_',num2str(1000.0*kappa),'_delta_',num2str((Egap-hwc)*1000),'.eps'],'-depsc2','-r300')

    % Print the file for publication
%    print(f,szFilename,'-dpng','-r300')
%    print(f,szFilename,'-depsc2','-r600')

end


