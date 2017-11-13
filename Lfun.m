function [L0,L1,L2] = Lfun(Mu, Hmol, GammaI, GammaII, StepsPerEv)

%keyboard
NumberOfSites = size(Hmol,1);

% Sparse identity slow?
IdentyMat=eye(NumberOfSites);

UsefulMatrix = Hmol + 0.5*sqrt(-1).*(GammaI+GammaII);

%global T0
%global h
h= 4.1357E-15;
kT = 8.61738569E-5*300; % in eV
Beta = 1.0/kT;

EnergySteps=3+floor(40.0*kT*StepsPerEv);
Tsig = 0.0;
% If we have a single mu data point
if (length(Mu) == 1)
    EnergyRange = linspace(Mu-20*kT,Mu+20*kT,EnergySteps); 
    Trans12 = zeros(1,EnergySteps);
    %UsefulMatrix = sparse(UsefulMatrix);
    %% Calculate transmission function over energy window
    parfor i=1:EnergySteps
        E=EnergyRange(i);

        %gr=inv(E*IdentyMat-UsefulMatrix);   
        
        % Try without inverse (appears faster/better?)
        gr = (E*IdentyMat-UsefulMatrix)\IdentyMat;
        gr_adj = gr';
        Trans12(i)=real(trace(GammaI*gr*GammaII*gr_adj)) + Tsig;

    end
    
    

%% Calculate L values

    FermiDeriv = -(Beta)./(4.0*(cosh(Beta*(EnergyRange-Mu)./2.0).^2.0));
    L0 = sum((-FermiDeriv).*Trans12);
    L1 = sum((-FermiDeriv).*Trans12.*(EnergyRange-Mu));
    L2 = sum((-FermiDeriv).*Trans12.*((EnergyRange-Mu).^2.0));

else
    
    L0 = zeros(length(Mu),1);
    L1 = zeros(length(Mu),1);
    L2=  zeros(length(Mu),1);

    Trans12 = zeros(1,EnergySteps);

    muIndex = 1;
    for mu=Mu
        EnergyRange = linspace(mu-20*kT,mu+20*kT,EnergySteps); 




        %UsefulMatrix = sparse(UsefulMatrix);
        %% Calculate transmission function over energy window
        parfor i=1:EnergySteps
            E=EnergyRange(i);

            gr=inv(E*IdentyMat-UsefulMatrix);



            % Try without inverse (appears faster/better?)
            %gr = (E*IdentyMat-UsefulMatrix)\IdentyMat;

            Trans12(i)=trace(real(GammaI*gr*GammaII*gr'));

        end



    %% Calculate L values

        FermiDeriv = -(Beta)./(4.0*(cosh(Beta*(EnergyRange-mu)./2.0).^2.0));
        L0(muIndex) = sum((-FermiDeriv).*Trans12);
        L1(muIndex) = sum((-FermiDeriv).*Trans12.*(EnergyRange-mu));
        L2(muIndex) = sum((-FermiDeriv).*Trans12.*((EnergyRange-mu).^2.0));
        muIndex = muIndex + 1;
    end
end

dE = EnergyRange(2)-EnergyRange(1);
L0 = dE*L0/h;
L1 = dE*L1/h;
L2 = dE*L2/h;


   

