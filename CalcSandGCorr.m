
% To complement hte theoretical effort. Calculate the PDF of a simple
% system with thermopower and conductance data.

%% Uncomment this line if we wish to repeat a random sequence
rng('default');


% Use SI units
JtoeV=1.0/1.60217E-19;
kb = 1.380648E-23*JtoeV; % J/K %8.6173e-5
T0 = 300; % K

%q=-1; %1.602176E-19 % Coulomb
h=6.62607E-34*JtoeV; % J*s

IsMatlab = true;

%% Number of members in our data set for each mu
Steps = 2^8


% Benzene dithiol connected to gold 
t = 5.0;


% How larger can gamma get
%figure; x=linspace(0,1,50); plot(x,Gamma_avg+sqrt(Gamma_avg.^2+8*x.^2*log(2)))



% What is the variation over crystal faces / molecular configurations?
% Crystal face shifts? On the order of 0.1eV
epsilon_avg = 0 ; % this value defines the zero of energy. It is the mid-gap energy
epsilon_std = .16%0.16; %1.06; %0.16; % Fig to pm2 e

GammaSteps = 1
GammaRange = linspace(1,1,GammaSteps);
Gamma_std = .05

MuSteps =20;
MuRange = linspace(-6,6,MuSteps);
%MuRange = -1

CorrData = zeros(length(GammaRange),length(MuRange),2);

    DOB = datestr(now);
if (GammaSteps == 1)
   Gamma_avg = GammaRange;
    Filename = [datestr(now,'yyyy_mmdd-(HH_MM_SS)'),'-big_stat_fun-steps=',num2str(Steps),',mu(1)=',num2str(MuRange(1)),',musteps=',num2str(length(MuRange)),',t=',num2str(t),',Gamma=(',num2str(Gamma_avg),',',num2str(Gamma_std),'),eps=(',num2str(epsilon_avg),',',num2str(epsilon_std),').mat']
else
    Filename = [datestr(now,'yyyy_mmdd-(HH_MM_SS)'),'-big_stat_fun-steps=',num2str(Steps),',mu(1)=',num2str(MuRange(1)),',musteps=',num2str(length(MuRange)),',t=',num2str(t),',Gamma=(',num2str(GammaSteps),',',num2str(Gamma_std),'),eps=(',num2str(epsilon_avg),',',num2str(epsilon_std),').mat']
end    
    
    BigS_data_con = zeros(length(GammaRange),length(MuRange),Steps);
    BigS_data_des = zeros(length(GammaRange),length(MuRange),Steps);
    BigG_data_con = zeros(length(GammaRange),length(MuRange),Steps);
    BigG_data_des = zeros(length(GammaRange),length(MuRange),Steps);
    
    Ldata_des = zeros(3,length(GammaRange),length(MuRange),Steps);
    Ldata_con =  zeros(3,length(GammaRange),length(MuRange),Steps);
   tic 
Gamma_count=1;
for Gamma_avg=GammaRange


    % HWHM ~ 1.1774

    Gam_eff = 0.5*(Gamma_avg+sqrt(Gamma_avg.^2+8*Gamma_std^2*log(2)));

    Gam_eff2 =0.5*(Gam_eff+sqrt(Gam_eff.^2 + 8*epsilon_std^2*log(2)));

    Gam_eff3 = Gamma_avg + (2*log(2)/(Gamma_avg^2))*(epsilon_std^2 + Gamma_std^2);

    Gam_eff4 = 0.5*(Gamma_avg + Gamma_std*2.35) + 0.5*sqrt((Gamma_avg + Gamma_std*2.35)^2+8*log(2)*epsilon_std^2);
    % Au 5.47 - 5.31
    % Pt 5.84 - 5.12 (331?)

    %t=5.0 %2.84/2.0; % the gap is 2t. In eV.

    % Build left and right lead-molecule coupling lists (list over all elements
    % of the distribution)

    GammaL = abs(Gamma_avg + Gamma_std*randn(Steps,1));
    GammaR = GammaL; %abs(Gamma_avg + Gamma_std*randn(Steps,1));
    epsilon = epsilon_avg + epsilon_std*randn(Steps,1);

 


    mucount = 1;
    for mu=MuRange

 disp(mu)
        tic
        data_des = zeros(Steps,3);    
        data_con = zeros(Steps,3);
        parfor ii=1:Steps % Gamma Loop (parfor splits across processors)
            GamL = GammaL(ii);
            GamR = GammaR(ii);
            eps = epsilon(ii);

            % 2 site 'stub resonator' (has a node at zero of energy)
            Hmol = [eps, t; t, eps];
            Gamma_L = [GamL,0; 0,0];
            Gamma_R = [GamR,0; 0,0]; 
            [L0,L1,L2] = Lfun(mu,Hmol, Gamma_L, Gamma_R, 200); %-Lfun(mu,Hmol, Gamma_L, Gamma_R, 100);
            data_des(ii,:) = [L0,L1,L2];
            
            % 2 site 'stub resonator' without a node
            Gamma_L = [GamL,0; 0,0];
            Gamma_R = [0,0; 0,GamR]; %[0, 0; 0, Gamma_avg/2.0];
            [L0,L1,L2] = Lfun(mu,Hmol, Gamma_L, Gamma_R, 200);
            data_con(ii,:)=[L0,L1,L2];
        end
        toc
        
        % Calculate G,S,kappa from L0,L1,L2
        G = h*L0;
        S = -1e6*(1/T0)*L1./L0;
        kappa = (1/T0)*(L2 - (L1).^2./L0);

        L0 = data_des(:,1);
        L1 = data_des(:,2);
        L2 = data_des(:,3);
        BigG_data_des(Gamma_count,mucount,:) = h*L0; 
        BigS_data_des(Gamma_count,mucount,:) = -1e6*(1/T0)*L1./L0;

        L0 = data_con(:,1);
        L1 = data_con(:,2);
        L2 = data_con(:,3);
        BigG_data_con(Gamma_count,mucount,:) = h*L0; 
        BigS_data_con(Gamma_count,mucount,:) = -1e6*(1/T0)*L1./L0;
        
       Ldata_des(1,Gamma_count,mucount,:) = L0;
       Ldata_des(2,Gamma_count,mucount,:) = L1;
       Ldata_des(3,Gamma_count,mucount,:) = L2;
       
      
        
        % Put each mu's data into a larger data structure
        
       % BigS_data_con(Gamma_count,mucount,:) = -1e6*(1/T0)*L1./L0;



       % Ldata_con(Gamma_count,mucount,:) = [L0,L1,L2];
       % Ldata_des(Gamma_count,mucount,:)
      %  disp(['<G/G_0>=',num2str(mean(BigG_data_con)),',FWHM(G/G_0)=',num2str(2*sqrt(2*log(2))*std(BigG_data_con)),',<S>=',...
      %      num2str(mean(BigS_data_con)),'uV/K',',FWHM(S)=',num2str(2*sqrt(2*log(2))*std(BigS_data_con))])

        if (IsMatlab == false)
             % For octave
            % Collapse the S and G data
            Sdata = squeeze(BigS_data_con(Gamma_count,mucount,:))
            Gdata = squeeze(BigG_data_con(Gamma_count,mucount,:))
            X = cov(Sdata,Gdata)/(std(Sdata)*std(Gdata));
            CorrData(Gamma_count,mucount,1) = X(1,2); %

            X = cov(BigS_data_des(mucount,:),BigG_data_des(mucount,:))/(std(BigS_data_des(mucount,:))*std(BigG_data_des(mucount,:)));
            CorrData(Gamma_count,mucount,2) = X(1,2);
        else % For matlab
            CorrData(Gamma_count,mucount,1) = corr(squeeze(BigS_data_con(Gamma_count,mucount,:)),squeeze(BigG_data_con(Gamma_count,mucount,:)),'type','Pearson');
            CorrData(Gamma_count,mucount,2) = corr(squeeze(BigS_data_des(Gamma_count,mucount,:)),squeeze(BigG_data_des(Gamma_count,mucount,:)),'type','Pearson');
        end


        
        disp(['mu count=',num2str(mucount), '/',num2str(length(MuRange))])
        mucount = mucount + 1;
    end
 

 
    Gamma_count = Gamma_count + 1
end
toc
   % Save the data
     output=struct('OriginalFilename',Filename,'DOB',DOB,'Steps',Steps, ...
          'MuRange',MuRange,'GammaRange',GammaRange,'GammaL',GammaL,'GammaR',GammaR,'Gamma_avg',Gamma_avg,'Gamma_std',Gamma_std,...
          'epsilon',epsilon,'epsilon_avg',epsilon_avg,'epsilon_std',epsilon_std,'CorrData',CorrData,...
          'BigS_data_des',BigS_data_des,'BigS_data_con',BigS_data_con,'BigG_data_des',BigG_data_des,'BigG_data_con',BigG_data_con,'t',t,'T0',T0);

    save(Filename, '-struct', 'output')
    disp(['CalcSandGCorr7: File "',Filename,'" written'])
    
 % Plot the results   
 PlotSandGCorr
 
 
 %% Extra stuff
%   if 0
%         % Lin's concordance correlation https://en.wikipedia.org/wiki/Concordance_correlation_coefficient
%         rho = corr(data_con(:,1),data_con(:,2));
%         mu_1 = mean(data_con(:,1));
%         mu_2 = mean(data_con(:,2));
%         sigma_1 = var(data_con(:,1));
%         sigma_2 = var(data_con(:,2));
% 
%         rho_c = 2*rho*sqrt(sigma_1)*sqrt(sigma_2) / (sigma_1 + sigma_2 + (mu_1-mu_2)^2);
%         CorrData(mucount,1) = rho_c;
% 
% 
%         rho = corr(data_des(:,1),data_des(:,2));
%         mu_1 = mean(data_des(:,1));
%         mu_2 = mean(data_des(:,2));
%         sigma_1 = var(data_des(:,1));
%         sigma_2 = var(data_des(:,2));
% 
%         rho_c = 2*rho*sqrt(sigma_1)*sqrt(sigma_2) / (sigma_1 + sigma_2 + (mu_1-mu_2)^2);
%         CorrData(mucount,2) = rho_c;
%     end
