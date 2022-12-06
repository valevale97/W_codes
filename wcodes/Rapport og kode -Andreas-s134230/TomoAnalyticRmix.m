%**************************************************************************
% Name:         TomoAnalyticRmix
%
% Purpose:      Contains the function 'TomoAnalyticRmix' which returns 
%               reconstructed 2D fast-ion velocity or (energy,pitch)
%               distribution functions and the weight functions used to
%               obtain them from supplied spectra via a tomographic 
%               inversion procedure based on Tikhonov regularization. The
%               function is capable of using spectra obtained via
%               Collective Thompson Scattering (CTS), Fast-ion D_alpha
%               spectroscopy (FIDA), Neutron Emission Spectrometry (NES),
%               one-step reaction Gamma-ray spectrometry (1stepGRS),
%               two-step reaction Gamma-ray spectrometry (2stepGRS) or any 
%               combination of these.
%               The equations used are explained or derived in the report
%               "Analytic models in velocity-space tomography for fusion
%               plasmas" by Andreas Poulsen.
%               All references refer to this report.
%**************************************************************************

function [WTomo,xalpha] = TomoAnalyticRmix(Svec,Spoint,Sdim,Serror,...
    pspacevec1,pspacevec2,vd,phivec,alpha,Tikhonov,uselsqnonneg,...
    grad_cutoff,method,pspace)
%TomoAnalyticRmix Function used for velocity-space tomography
%
% Output parameters
%------------------
% WTomo     - Matrix with each row containing the weight function
%             corresponding to a particular combination of observation 
%             angle and spectrum X-value (bin). 
%             #rows is equal to length(phivec)*length(Spoint) and #columns
%             is equal to the number of elements in the chosen tomography
%             grid i.e. length(pspacevec1)*length(pspacevec2)
% xalpha    - Matrix with each column containing the reconstructed
%             distribution function obtained for a given value of the
%             regularization parameter alpha.
%             #rows is equal to to the number of elements in the chosen
%             tomography grid i.e. length(pspacevec1)*length(pspacevec2)
%             and #columns is equal to length(alpha)
%
%
% Input parameters
%------------------
% Svec           - Ordered column vector containing the spectrum Y-values 
%                  (signal values) for each observation angle. #values is
%                  equal to sum(Sdim). The values for each observation
%                  angle must be ordered in the same way as the values of
%                  Spoint for that observation angle, and the ordering of
%                  the observation angles must be the same as in phivec.
% Spoint         - Ordered row vector containing the spectrum X-values
%                  (bins) used for every observation angle. The ordering
%                  must be the same as that of the values in Svec.
% Sdim           - Ordered row vector containing the number of measurements
%                  (spectrum values) for each observation angle. The
%                  ordering must be the same as in phivec. Alternatively,
%                  Sdim may be a single value if the numbr of measurements
%                  was the same for all observation angles
% Serror         - Ordered row vector containing the uncertainties in the
%                  spectrum Y-values. The vector may contain a value for
%                  each observation angle or a value for each value in 
%                  Svec. It may also contain a single value.
% pspacevec1     - Ordered row vector containing all of the values of
%                  energy (E) or v_parallel (vpa) to be used for the
%                  tomography grid. [eV] for E and [m/s] for vpa.
% pspacevec2     - Ordered row vector containing all of the values of pitch
%                  (p) or v_perpendicular (vpe) to be used for the
%                  tomography grid. [] for p and [m/s] for vpe.
% vd             - Toroidal drift velocity of the bulk plasma. Only affects
%                  NES, one- and two-step reaction GRS.
% phivec         - Ordered vector containing the used observation
%                  angles in degrees. The values must be ordered in the
%                  same way as in Svec.
% alpha          - Vector containing the desired values of the
%                  regularization parameter to be used for the tomographic
%                  inversion via Tikhonov regularization.
% Tikhonov       - Parameter used to indicate the desired Tikhonov
%                  procedure. The value should be 0 for 0th order Tikhonov,
%                  1 for 1st order Tikhonov or 2 for a mix of 0th and 1st
%                  order Tikhonov.
% uselsgnoneg    - Parameter used to indicate whether or not the output
%                  distribution functions should be constrained to be 
%                  non-negative. The value should be 0 for unconstrained
%                  Tikhonov or 1 for non-negative Tikhonov. Note that using
%                  the non-negativity constraint significantly increases
%                  the computation time.
% grad_cutoff    - The gradient used for the Tikhonov reconstruction is
%                  removed from the regions where the weight function
%                  coverage (sum of all weight functions) is less than
%                  grad_cutoff times the mean of the weight function
%                  coverage. This leads to these regions being disregarded
%                  in the reconstructions.
% method         - Optional string array used to indicate the used method
%                  for each diagnostic. Can be CTS, NESDD, NESDT, 
%                  GRSDpfastp, GRSDpfastD, FIDA or any combination of
%                  these. #elements must be 1 or equal to length(phivec)
% pspace         - Optional string used to indicate the desired plotting
%                  space for the tomographic inversion. Can be Ep or vv.
%
% Fixed constants
Mp = 1.6726e-27; % [kg] Mass of a proton
Mn = 1.6749e-27; % [kg] Mass of a neutron
Mi = 2*Mp; % [kg] Mass of a deuterium ion
M_He = 3*Mp; % [kg] Mass of a helium-3 particle
n_t = 5e19; % [m^-3] The thermal ion density
Qe = 1.6021917e-19; % [C] The elementary charge
c = 3e8; % [m/s] The speed of light

% List of possible methods
methodvec=["CTS", "NESDD", "NESDTfastD", "NESDTfastT", "GRSDpfastp",...
    "GRSDpfastD", "FIDA","2stepGRS"];
methods={'CTS','NESDD','NESDTfastD','NESDTfastT','GRSDpfastp',...
    'GRSDpfastD','FIDA','2stepGRS'};

% The function starts by checking whether or not the variables method and
% pspace were specified and give the user the option to select them if they
% were not. In this case, it is assumed that all diagnostics used the
% chosen method. Note that it is possible to specify method and not pspace,
% but it is not possible to specify pspace and not method.

inputval=1;
inputval2=1;
options.Interpreter = 'tex';
options.Default = 'A';
switch nargin % The number of given input parameters is checked
    case {0,1,2,3,4,5,6,7,8,9,10,11}
        disp('Error: Not enough input arguments.')
        inputval=0;
    case 12
        disp('No diagnostic method or plotting space chosen.')
        [indx, tf]=listdlg('PromptString',...
            'Which diagnostic method was used?',...
            'SelectionMode','single', 'ListString',methods); 
        if ~tf % If no method is chosen, the function stops
            inputval=0;
        else
            method=methodvec(indx);
            pspaceval=questdlg(...
                'Plot in A) (E,p)- or B) (v_{||},v_{\perp})-space?',...
                'Plotting space Menu', 'A', 'B',options);
            % (E,p)-space is the default choice if one presses enter
            if strcmpi(pspaceval,'') 
                % If no plotting space is chosen, the function stops
                inputval=0;
            elseif strcmpi(pspaceval,'A')
                pspace='Ep';
            elseif strcmpi(pspaceval,'B')
                pspace='vv';
            end
        end
    case 13
        disp('No plotting space chosen.')
        pspaceval=questdlg(...
            'Plot in A) (E,p)- or B) (v_{||},v_{\perp})-space?',...
            'Plotting space Menu', 'A', 'B',options);
        % (E,p)-space is the default choice if one presses enter
        if strcmpi(pspaceval,'')
            % If no plotting space is chosen, the function stops
            inputval=0;
        elseif strcmpi(pspaceval,'A')
            pspace='Ep';
        elseif strcmpi(pspaceval,'B')
            pspace='vv';
        end
    case 14
    otherwise
        disp('Error: Too many input arguments.')
        inputval=0;
end

% The function checks whether or not the specified variables method and
% pspace match any of the possible choices and gives the option to select
% them if they do not.

if inputval
    if sum(strcmpi(method,'CTS'))+ sum(strcmpi(method,'NESDD'))+...
            sum(strcmpi(method,'NESDTfastD'))+...
            sum(strcmpi(method,'NESDTfastT'))+...
            sum(strcmpi(method,'GRSDpfastp'))+...
            sum(strcmpi(method,'GRSDpfastD'))+...
            sum(strcmpi(method,'FIDA'))...
            +sum(strcmpi(method,'2stepGRS'))~=length(method)
        tf=1;
        % The function checks each individual string in the method array
        % and gives the option to select replacements for any invalid ones.
        % If a replacement is not chosen, the function stops.
        for i =1:length(method)
            if ~strcmpi(method(i),methods) && tf==1
                fprintf('Diagnostic nr. %s is not a valid choice.\n',...
                    num2str(i,1))
                [indx, tf]=listdlg('PromptString',...
                    'Which diagnostic method was used?','SelectionMode',...
                    'single', 'ListString',methods);
                method(i)=methodvec(indx);
            end
        end
       
        if tf==1 && ~strcmpi(pspace,{'Ep','vv'})
            disp('No valid plotting space chosen')
            pspaceval=questdlg(...
                'Plot in A) (E,p)- or B) (v_{||},v_{\perp})-space?',...
                'Plotting space Menu', 'A', 'B',options);
            % (E,p)-space is the default choice
            if strcmpi(pspaceval,'')
                % If no plotting space is chosen, the function stops
                inputval2=0;
            elseif strcmpi(pspaceval,'A')
                pspace='Ep';
            elseif strcmpi(pspaceval,'B')
                pspace='vv';
            end
        elseif tf==0
            % If the invalid methods are not replaced, the function stops
            inputval2=0; 
        end
    elseif ~strcmpi(pspace,'Ep') && ~strcmpi(pspace,'vv')
        disp('No valid plotting space chosen.')
        pspaceval=questdlg(...
            'Plot in A) (E,p)- or B) (v_{||},v_{\perp})-space?',...
            'Plotting space Menu', 'A', 'B',options);
        % (E,p)-space is the default choice
        if strcmpi(pspaceval,'')
            % If no plotting space is chosen, the function stops
            inputval2=0;
        elseif strcmpi(pspaceval,'A')
            pspace='Ep';
        elseif strcmpi(pspaceval,'B')
            pspace='vv';
        end
    end
end

% If both method and pspace are specified properly, the function continues

if inputval && inputval2
    % If the method variable only contains a single string, all diagnostics
    % are assumed to have used that method
    if length(method)==1
        temp=strings(1,length(phivec));
        temp(1:length(phivec))=method;
        method=temp;
    end
    
    % If Sdim was a single value, it is changed to a vector of length
    % equal to that of phivec where every element has that value
    if length(Sdim)==1
        Sdim=ones(1,length(phivec))*Sdim;
    end
    
    % If Serror was a single value, it is changed to a vector of length
    % equal to that of Svec where every element has that value
    if length(Serror) == 1
        Serror=Serror.*ones(1,length(Svec));
    % If Serror contained a value for each observation angle, it is changed
    % to a vector of length equal to that of Svec where all
    % elements corresponding to a particular observation angle have the
    % value specified for that observation angle
    elseif length(Serror) == length(phivec)
        temp=ones(1,length(Svec));
        num=0;
        for i =1:length(phivec)
            temp((num+1):(num+Sdim(i)))=temp((num+1):(num+Sdim(i)))*...
                Serror(i);
            num=num+Sdim(i);
        end
        Serror=temp;
    end
    
    % The grid points in the tomography grid are generated
    [grid1,grid2]=meshgrid(pspacevec1,pspacevec2);
    % The observation angles are converted to radians
    phivec=phivec.*pi./180;
    Wrow=0;
    Sbin=0;
    % Pre-allocation of the matrix to hold the weight functions for
    % computational efficiency
    Wtomo=zeros(length(Svec),numel(grid1));
    % The weight functions for each individual diagnostic are calculated
    % and stored in the designated matrix
    for i=1:length(phivec)
        phi=phivec(i);
        % Parameters relevant for the considered diagnostic method are
        % defined
        if strcmpi(method(i),'CTS')
            methodval = 'CTS';
            M_f =Mi; % [kg] Mass of the fast ion (assumed to be deuterium)
        elseif strcmpi(method(i),'NESDD') % NES for D-D fusion
            methodval = 'NES';
            reaction = 'NESDD';
            M_f =Mi; % [kg] Mass of the fast ion (deuterium)
            % Coefficient used to calculate the fast-ion center-of-mass
            % energy
            CM_coef = 1/2;
            Q = 3.27e6; % [eV] Energy generated in the D-D fusion reaction
        elseif strcmpi(method(i),'NESDTfastD')
            % NES for D-T fusion with fast D
            methodval = 'NES';
            reaction = 'NESDT';
            M_f =Mi; % [kg] Mass of the fast ion (deuterium)
            M_He =4*Mp; % [kg] Mass of an alpha-particle
            % Coefficient used to calculate the fast-ion center-of-mass energy
            CM_coef = 2/5;
            Q = 17.6e6; % [eV] Energy generated in the D-T fusion reaction
        elseif strcmpi(method(i),'NESDTfastT')
            % NES for D-T fusion with fast T
            methodval = 'NES';
            reaction = 'NESDT';
            M_f =3*Mp; % [kg] Mass of the fast ion (tritium)
            M_He =4*Mp; % [kg] Mass of an alpha-particle
            % Coefficient used to calculate the fast-ion center-of-mass
            % energy
            CM_coef = 3/5;
            Q = 17.6e6; % [eV] Energy generated in the D-T fusion reaction
        elseif strcmpi(method(i),'GRSDpfastp')
            % GRS for D-p fusion with fast p
            methodval = 'GRS';
            reaction = 'GRSDp';
            M_f = Mp; % [kg] Mass of the fast ion (proton)
            M_pr = 3*Mp; % [kg] Mass of a helium-3 particle
            % Coefficient used to calculate the fast-ion center-of-mass
            % energy
            CM_coef = 1/3;
            Q = 5.5e6; % [eV] Energy generated in the D-p fusion reaction
        elseif strcmpi(method(i),'GRSDpfastD')
            % GRS for D-p fusion with fast D
            methodval = 'GRS';
            reaction = 'GRSDp';
            M_f = 2*Mp; % [kg] Mass of the fast ion (deuterium)
            M_pr = 3*Mp; % [kg] Mass of a helium-3 particle
            % Coefficient used to calculate the fast-ion center-of-mass
            % energy
            CM_coef = 2/3;
            Q = 5.5e6; % [eV] Energy generated in the D-p fusion reaction
        elseif strcmpi(method(i),'FIDA')
            methodval = 'FIDA';
            M_f = Mi; % [kg] Mass of the fast ion (deuterium)
            gamma_bar=0; % The phase-shift of the cosine model for the pdf
            B=1.74; % [T] The magnetic field
            E0=1e6; % [eV] The free parameter used in the model for the 
            % amplitude of the cosine model for the pdf
            lambda0=656.1e-9; % [m] The wavelength of the unshifted D_alpha
            % -line
            C_hat =[1,18,16,1681,2304,729,1936,5490,1936,729,2304,1681,...
                16,18,1]/18860; % Constants related to the probabilities of
            % the various Stark lines
            % Proportionality constants for the Stark splitting wavelength
            % shifts
            sl = [-220.2,-165.2,-137.7,-110.2,-82.64,-55.1,-27.56,0,...
                27.57,55.15,82.74,110.3,138.0,165.6,220.9]*1e-18; % [m^2/V]
            % Signs related to the various Stark lines. + for the
            % sigma-lines and - for the pi-lines
            Stark_sign = [-1,1,1,-1,-1,-1,1,1,1,-1,-1,-1,1,1,-1]; 
            % Parameters used for the ad-hoc model for the FIDA intensity
            % function R
            ni = 1e8;     % [m^-3] Ion density
            RTpara = 60e3; % [eV] Parallel plasma temperature
            RTperp = 60e3; % [eV] Perpendicular plasma temperature
            Rvparadrift = 0.4e6; % [m/s] Parallel drift velocity
        elseif strcmpi(method(i),'2stepGRS')
            methodval= '2stepGRS';
            M_f = 4*Mp; % [kg] Mass of the fast ion (alpha-particle)
            Ega0=4.44e6; % [eV]
            Qdot=1.26e6; % [eV]
            n_Be=1e18; % [m^-3] Density of thermal Be
            betaloop=(90.5:0.5:179.5)/180*pi;
            zetaloop=(0:0.5:179.5)/180*pi;
            Esigma=0:0.1:6.1; % [MeV]
            sigma=[0 0 0 0 0 0 0 0 0 0 ...
                0 0 0 0 0 10 30 100 180 230 ...
                220 200 80 70 60 60 80 100 80 60 ...
                60 60 60 60 60 65 70 75 200 320 ...
                320 320 220 240 260 260 250 280 320 340 ...
                340 320 320 320 320 320 320 320 290 250 ...
                250 0]; % [mb]
        end
        if strcmpi(methodval,'CTS') && strcmpi(pspace,'Ep')
            ubroadening=3; % Spectral resolution of the measurements
            % divided by the bin width of the spectra
            % The step sizes in E and p are defined
            dE=pspacevec1(2)-pspacevec1(1);
            dp=pspacevec2(2)-pspacevec2(1);
            % The spectral resolution of the measurements
            dures=ubroadening*(max(Spoint((Sbin+1):(Sbin+Sdim(i))))...
                -min(Spoint((Sbin+1):(Sbin+Sdim(i)))))/(Sdim(i)-1);
            % For each bin in the spectra, the weight function for the
            % chosen tomography grid is calculated using eq. 2.19
            % and stored as a row in a matrix
            for j=(Sbin+1):(Sbin+Sdim(i))
                u=Spoint(j);
                Wrow=Wrow+1;
                % The two gyro-angles (eq. 2.20)
                gamma1 = acos(1./sin(phi).*1./sqrt(1-grid2.^2).*...
                    ((u-dures./2)./sqrt(2.*grid1./Mi*Qe)-grid2.*cos(phi)));
                gamma2 = acos(1./sin(phi).*1./sqrt(1-grid2.^2).*...
                    ((u+dures./2)./sqrt(2.*grid1./Mi*Qe)-grid2.*cos(phi)));
                % The weight function (eq. 2.19) times the step sizes
                wfcEpCTS=real((gamma1-gamma2))./pi./dures.*dE.*dp;
                % The weight function is normalized by the associated
                % uncertainty
                Wtomo(Wrow,:)=reshape(wfcEpCTS,1,numel(grid1))./...
                    Serror(Wrow);
            end
        elseif strcmpi(methodval,'CTS') && strcmpi(pspace,'vv')
            ubroadening=3; % Spectral resolution of the measurements
            % divided by the bin width of the spectra
            % The step sizes in vpa and vpe are defined
            dvpa=pspacevec1(2)-pspacevec1(1);
            dvpe=pspacevec2(2)-pspacevec2(1);
            % The spectral resolution of the measurements
            dures=ubroadening*(max(Spoint((Sbin+1):(Sbin+Sdim(i))))...
                -min(Spoint((Sbin+1):(Sbin+Sdim(i)))))/(Sdim(i)-1);
            % For each bin in the spectra, the weight function for the
            % chosen tomography grid is calculated using eq. 2.19
            % and stored as a row in a matrix
            for j=(Sbin+1):(Sbin+Sdim(i))
                u=Spoint(j);
                Wrow=Wrow+1;
                % The two gyro-angles (eq. 2.20)
                gamma1=acos(1./(sin(phi).*grid2).*(u-dures./2-grid1.*...
                    cos(phi)));
                gamma2=acos(1./(sin(phi).*grid2).*(u+dures./2-grid1.*...
                    cos(phi)));
                % The weight function (eq. 2.19) times the step sizes
                wfcvvCTS=real((gamma1-gamma2))./pi./dures.*dvpa.*dvpe;
                % The weight function is normalized by the associated
                % uncertainty
                Wtomo(Wrow,:)=reshape(wfcvvCTS,1,numel(grid1))./...
                    Serror(Wrow);
            end
        elseif strcmpi(methodval,'NES') && strcmpi(pspace,'Ep')
            Enbroadening=1; % Spectral resolution of the measurements
            % divided by the bin width of the spectra
            % The step sizes in E and p are defined
            dE=pspacevec1(2)-pspacevec1(1);
            dp=pspacevec2(2)-pspacevec2(1);
            % The spectral resolution of the measurements
            dEnres=Enbroadening*(max(Spoint((Sbin+1):(Sbin+Sdim(i))))...
                -min(Spoint((Sbin+1):(Sbin+Sdim(i)))))/(Sdim(i)-1);
            % The values of vpa and vpe corresponding to the chosen
            % (E,p)-grid are calculated (eq. 2.3)
            vpara=grid2.*sqrt(2.*grid1.*Qe./M_f);
            vperp=sqrt(1-grid2.^2).*sqrt(2.*grid1.*Qe./M_f);
            % The center-of-mass energy of the fast ion in [keV]
            E_CM=CM_coef.*M_f.*((vpara-vd).^2+vperp.^2)./Qe/2e3;
            % The rate part of the weight function (eq. 2.51)
            RNES=n_t.*Bosch_Hale_sigma(E_CM,reaction).*10^-31.*...
            sqrt((vpara-vd).^2+vperp.^2);
            RNES(RNES<0)=0; % Negative rates are set to zero
            % For each bin in the spectra, the weight function for the
            % chosen tomography grid is calculated using eq. 2.44
            % and stored as a row in a matrix
            for j=(Sbin+1):(Sbin+Sdim(i))
                En=Spoint(j);
                Wrow=Wrow+1;
                % The two gyro-angles (eq. 2.48)
                gamma1=acos(1./(sqrt(1-grid2.^2).*sin(phi)).*(...
                    1./2.*(M_He+Mn)./sqrt(M_f.*Mn).*sqrt((En-...
                    dEnres./2)./grid1)-1./2.*(M_He-M_f)./sqrt(M_f.*Mn).*...
                    sqrt(grid1./(En-dEnres./2))-1./2.*M_He./sqrt(M_f.*...
                    Mn).*Q./(sqrt(grid1.*(En-dEnres./2)))-...
                    grid2.*cos(phi)));
                gamma2=acos(1./(sqrt(1-grid2.^2).*sin(phi)).*(...
                    1./2.*(M_He+Mn)./sqrt(M_f.*Mn).*sqrt((En+...
                    dEnres./2)./grid1)-1./2.*(M_He-M_f)./sqrt(M_f.*Mn).*...
                    sqrt(grid1./(En+dEnres./2))-1./2.*M_He./sqrt(M_f.*...
                    Mn).*Q./(sqrt(grid1.*(En+dEnres./2)))-...
                    grid2.*cos(phi)));
                % The probability part (eq. 2.45)
                probNES=real((gamma1-gamma2)./pi);
                % The total weight function (eq. 2.44) times the step sizes
                wfcEpNES=probNES.*RNES./dEnres.*dE.*dp;
                % The weight function is normalized by the associated
                % uncertainty
                Wtomo(Wrow,:)=reshape(wfcEpNES,1,numel(grid1))./...
                    Serror(Wrow);
            end
        elseif strcmpi(methodval,'NES') && strcmpi(pspace,'vv')
            Enbroadening=1; % Spectral resolution of the measurements
            % divided by the bin width of the spectra
            % The step sizes in vpa and vpe are defined
            dvpa=pspacevec1(2)-pspacevec1(1);
            dvpe=pspacevec2(2)-pspacevec2(1);
            % The spectral resolution of the measurements
            dEnres=Enbroadening*(max(Spoint((Sbin+1):(Sbin+Sdim(i))))...
                -min(Spoint((Sbin+1):(Sbin+Sdim(i)))))/(Sdim(i)-1);
            % The center-of-mass energy of the fast ion in [keV]
            E_CM=CM_coef.*M_f.*((grid1-vd).^2+grid2.^2)./Qe/2e3;
            % The rate part of the weight function (eq. 2.51)
            RNES=n_t.*Bosch_Hale_sigma(E_CM,reaction).*10^-31.*...
                sqrt((grid1-vd).^2+grid2.^2);
            RNES(RNES<0)=0; % Negative rate are set to zero
            % For each bin in the spectra, the weight function for the
            % chosen tomography grid is calculated using eq. 2.44
            % and stored as a row in a matrix
            for j=(Sbin+1):(Sbin+Sdim(i))
                En=Spoint(j);
                Wrow=Wrow+1;
                % The two gyro-angles (eq. 2.47)
                gamma1=acos(1./(grid2.*sin(phi)).*(0.5.*(M_He+Mn)./M_f.*...
                    sqrt(2.*(En-dEnres./2)./Mn.*Qe)-0.5.*(M_He-...
                    M_f)./Mn.*(grid1.^2+grid2.^2)./sqrt(2.*(En-...
                    dEnres./2)./Mn.*Qe)-M_He./(M_f.*Mn).*Q.*Qe./sqrt(2.*...
                    (En-dEnres./2)./Mn.*Qe)-grid1.*cos(phi)));
                gamma2=acos(1./(grid2.*sin(phi)).*(0.5.*(M_He+Mn)./M_f.*...
                    sqrt(2.*(En+dEnres./2)./Mn.*Qe)-0.5.*(M_He-...
                    M_f)./Mn.*(grid1.^2+grid2.^2)./sqrt(2.*(En+...
                    dEnres./2)./Mn.*Qe)-M_He./(M_f.*Mn).*Q.*Qe./sqrt(2.*...
                    (En+dEnres./2)./Mn.*Qe)-grid1.*cos(phi)));
                % The probability part (eq. 2.45)
                probNES=real((gamma1-gamma2))./pi;
                % The total weight function (eq. 2.44) times the step sizes
                wfcvvNES=probNES.*RNES./dEnres.*dvpa.*dvpe;
                % The weight function is normalized by the associated
                % uncertainty
                Wtomo(Wrow,:)=reshape(wfcvvNES,1,numel(grid1))./...
                    Serror(Wrow);
            end
        elseif strcmpi(methodval,'GRS') && strcmpi(pspace,'Ep')
            Ebroadening=1; % Spectral resolution of the measurements
            % divided by the bin width of the spectra
            % The step sizes in E and p are defined
            dE=pspacevec1(2)-pspacevec1(1);
            dp=pspacevec2(2)-pspacevec2(1);
            % The spectral resolution of the measurements
            dEgares=Ebroadening*(max(Spoint((Sbin+1):(Sbin+Sdim(i))))...
                -min(Spoint((Sbin+1):(Sbin+Sdim(i)))))/(Sdim(i)-1);
            % The values of vpa and vpe corresponding to the chosen
            % (E,p)-grid are calculated (eq. 2.3)
            vpara=grid2.*sqrt(2.*grid1.*Qe./M_f);
            vperp=sqrt(1-grid2.^2).*sqrt(2.*grid1.*Qe./M_f);
            % The center-of-mass energy of the fast ion in [keV]
            E_CM=CM_coef.*M_f.*((vpara-vd).^2+vperp.^2)./Qe/2e3;
            % The rate part of the weight function (eq. 2.51)
            RGRS=n_t.*Bosch_Hale_sigma(E_CM,reaction).*10^-31.*...
            sqrt((vpara-vd).^2+vperp.^2);
            RGRS(RGRS<0)=0; % Negative rates are set to zero
            % For each bin in the spectra, the weight function for the
            % chosen tomography grid is calculated using eq. 2.78
            % and stored as a row in a matrix
            for j=(Sbin+1):(Sbin+Sdim(i))
                Ega=Spoint(j);
                Wrow=Wrow+1;
                % The two gyro-angles (eq. 2.82)
                gamma1=acos(1./(vperp.*sin(phi)).*((M_f-M_pr).*c./...
                    (2.*(Ega-dEgares./2).*Qe).*(vpara.^2+vperp.^2)+...
                    (Ega-dEgares./2).*Qe./(2.*M_f.*c)+M_pr.*c.*...
                    (Ega-dEgares./2-Q)./(M_f.*(Ega-dEgares./2))-vpara.*...
                    cos(phi)));
                gamma2=acos(1./(vperp.*sin(phi)).*((M_f-M_pr).*c./...
                    (2.*(Ega+dEgares./2).*Qe).*(vpara.^2+vperp.^2)+...
                    (Ega+dEgares./2).*Qe./(2.*M_f.*c)+M_pr.*c.*...
                    (Ega+dEgares./2-Q)./(M_f.*(Ega+dEgares./2))-vpara.*...
                    cos(phi)));
                % The probability part (eq. 2.79)
                probGRS=real(gamma1-gamma2)./pi;
                % The total weight function (eq. 2.78) times the step sizes
                wfcEpGRS=probGRS.*RGRS./dEgares.*dE.*dp;
                % The weight function is normalized by the associated
                % uncertainty
                Wtomo(Wrow,:)=reshape(wfcEpGRS,1,numel(grid1))./...
                    Serror(Wrow);
            end
        elseif strcmpi(methodval,'GRS') && strcmpi(pspace,'vv')
            Ebroadening=1; % Spectral resolution of the measurements
            % divided by the bin width of the spectra
            % The step sizes in vpa and vpe are defined
            dvpa=pspacevec1(2)-pspacevec1(1);
            dvpe=pspacevec2(2)-pspacevec2(1);
            % The spectral resolution of the measurements
            dEgares=Ebroadening*(max(Spoint((Sbin+1):(Sbin+Sdim(i))))...
                -min(Spoint((Sbin+1):(Sbin+Sdim(i)))))/(Sdim(i)-1);
            % The center-of-mass energy of the fast ion in [keV]
            E_CM=CM_coef.*M_f.*((grid1-vd).^2+grid2.^2)./Qe/2e3;
            % The rate part of the weight function (eq. 2.51)
            RGRS=n_t.*Bosch_Hale_sigma(E_CM,reaction).*10^-31.*...
                sqrt((grid1-vd).^2+grid2.^2);
            RGRS(RGRS<0)=0; % Negative rates are set to zero
            % For each bin in the spectra, the weight function for the
            % chosen tomography grid is calculated using eq. 2.78
            % and stored as a row in a matrix
            for j=(Sbin+1):(Sbin+Sdim(i))
                Ega=Spoint(j);
                Wrow=Wrow+1;
                % The two gyro-angles (eq. 2.81)
                gamma1=acos(1./(grid2.*sin(phi)).*((M_f-M_pr).*c./...
                    (2.*(Ega-dEgares./2).*Qe).*(grid1.^2+grid2.^2)+...
                    (Ega-dEgares./2).*Qe./(2.*M_f.*c)+M_pr.*c.*...
                    (Ega-dEgares./2-Q)./(M_f.*(Ega-dEgares./2))-grid1.*...
                    cos(phi)));
                gamma2=acos(1./(grid2.*sin(phi)).*((M_f-M_pr).*c./...
                    (2.*(Ega+dEgares./2).*Qe).*(grid1.^2+grid2.^2)+...
                    (Ega+dEgares./2).*Qe./(2.*M_f.*c)+M_pr.*c.*...
                    (Ega+dEgares./2-Q)./(M_f.*(Ega+dEgares./2))-grid1.*...
                    cos(phi)));
                % The probability part of the weight function (eq. 2.79)
                probGRS=real(gamma1-gamma2)./pi;
                % The total weight function (eq. 2.78) times the step sizes
                wfcvvGRS=probGRS.*RGRS./dEgares.*dvpa.*dvpe;
                % The weight function is normalized by the associated
                % uncertainty
                Wtomo(Wrow,:)=reshape(wfcvvGRS,1,numel(grid1))./...
                    Serror(Wrow);
            end
        elseif strcmpi(methodval,'FIDA') && strcmpi(pspace,'Ep')
            lambdabroadening=1; % Spectral resolution of the measurements
            % divided by the bin width of the spectra
            % The step sizes in E and p are defined
            dE=pspacevec1(2)-pspacevec1(1);
            dp=pspacevec2(2)-pspacevec2(1);
            % The spectral resolution of the measurements
            dlambdares=lambdabroadening*(max(Spoint((Sbin+1):(Sbin+...
                Sdim(i))))-min(Spoint((Sbin+1):(Sbin+Sdim(i)))))/...
                (Sdim(i)-1);
            % The ad-hoc model for the FIDA intensity function R
            RFIDA = ni.*sqrt(grid1./(pi.*RTpara.*RTperp.^2.*Qe.^2))...
                .*exp(-(grid2.^2.*grid1.*Qe+0.5.*Mi.*Rvparadrift.^2-...
                Rvparadrift.*grid2.*sqrt(2.*Mi.*grid1.*Qe))/...
                (RTpara.*Qe)-(1-grid2.^2).*grid1./RTperp).*Qe;
            % The values of vpa and vpe corresponding to the chosen
            % (E,p)-grid are calculated (eq. 2.3)
            vpara=grid2.*sqrt(2.*grid1.*Qe./Mi);
            vperp=sqrt(1-grid2.^2).*sqrt(2.*grid1.*Qe./Mi);
            % For each bin in the spectra, the weight function for the
            % chosen tomography grid is calculated using eq. 2.26
            % and stored as a row in a matrix
            for j=(Sbin+1):(Sbin+Sdim(i))
                lambda=Spoint(j);
                Wrow=Wrow+1;
                probFIDA=zeros(size(grid1));
                for l=1:15 % Sum over the 15 Stark lines
                    % The two gyro-angles (eq. 2.31)
                    gamma1=acos(1./(vperp.*sin(phi)).*(c.*((lambda-...
                        dlambdares./2)./(lambda0+sl(l).*vperp.*B)-1)-...
                        vpara.*cos(phi)));
                    gamma2=acos(1./(vperp.*sin(phi)).*(c.*((lambda+...
                        dlambdares./2)./(lambda0+sl(l).*vperp.*B)-1)-...
                        vpara.*cos(phi)));
                    %The probability part of the weight function (eq. 2.37)
                    probFIDA=C_hat(l).*(real(gamma1-gamma2)/pi+...
                        Stark_sign(l).*sin(phi)^2/2.*(real(gamma1-...
                        gamma2)/pi-(sin(2.*real(gamma1))-sin(2.*...
                        real(gamma2)))/(2*pi))+2.*grid1./E0.*...
                        (1-grid2.^2).*cos(gamma_bar).*...
                        (sin(real(gamma1))-sin(real(gamma2))+...
                        Stark_sign(l).*sin(phi)^2/3.*...
                        (sin(real(gamma1)).^3-sin(real(gamma2)).^3)))+...
                        probFIDA;
                end
                % The total weight function (eq. 2.26) times the step sizes
                wfcEpFIDA=probFIDA.*RFIDA.*dE.*dp./dlambdares;
                % The weight function is normalized by the associated
                % uncertainty
                Wtomo(Wrow,:)=reshape(wfcEpFIDA,1,numel(grid1))./...
                    Serror(Wrow);
            end
        elseif strcmpi(methodval,'FIDA') && strcmpi(pspace,'vv')
            lambdabroadening=1; % Spectral resolution of the measurements
            % divided by the bin width of the spectra
            % The step sizes in vpa and vpe are defined
            dvpa=pspacevec1(2)-pspacevec1(1);
            dvpe=pspacevec2(2)-pspacevec2(1);
            % The spectral resolution of the measurements
            dlambdares=lambdabroadening*(max(Spoint((Sbin+1):(Sbin+...
                Sdim(i))))-min(Spoint((Sbin+1):(Sbin+Sdim(i)))))/...
                (Sdim(i)-1);
            % The ad-hoc model for the FIDA intensity function R
            vthpara=sqrt(2*RTpara*Qe/Mi);
            vthperp=sqrt(2*RTperp*Qe/Mi);
            RFIDA = ni./(pi.^(3/2).*vthpara.*vthperp.^2).*...
                exp(-((grid1-Rvparadrift)./vthpara).^2-...
                (grid2./vthperp).^2).*2.*pi.*grid2;
            % The values of E and p corresponding to the chosen
            % (vpa,vpe)-grid are calculated (eq. 2.2)
            E = 1/2.*Mi.*(grid1.^2+grid2.^2)./Qe;
            p = grid1./sqrt(grid1.^2+grid2.^2);
            % For each bin in the spectra, the weight function for the
            % chosen tomography grid is calculated using eq. 2.26
            % and stored as a row in a matrix
            for j=(Sbin+1):(Sbin+Sdim(i))
                lambda=Spoint(j);
                Wrow=Wrow+1;
                probFIDA=zeros(size(grid1));
                for l=1:15 % Sum over the 15 Stark lines
                    % The two gyro-angles (eq. 2.31)
                    gamma1=acos(1./(grid2.*sin(phi)).*(c.*((lambda-...
                        dlambdares./2)./(lambda0+sl(l).*grid2.*B)-1)-...
                        grid1.*cos(phi)));
                    gamma2=acos(1./(grid2.*sin(phi)).*(c.*((lambda+...
                        dlambdares./2)./(lambda0+sl(l).*grid2.*B)-1)-...
                        grid1.*cos(phi)));
                    %The probability part of the weight function (eq. 2.37)
                    probFIDA=C_hat(l).*(real(gamma1-gamma2)/pi+...
                        Stark_sign(l).*sin(phi)^2/2.*(real(gamma1-...
                        gamma2)/pi-(sin(2.*real(gamma1))-...
                        sin(2.*real(gamma2)))/(2*pi))+2.*E./E0.*...
                        (1-p.^2).*cos(gamma_bar).*(sin(real(gamma1))-...
                        sin(real(gamma2))+Stark_sign(l).*sin(phi)^2/3.*...
                        (sin(real(gamma1)).^3-sin(real(gamma2)).^3)))+...
                        probFIDA;
                end
                % The total weight function (eq. 2.26) times the step sizes
                wfcvvFIDA=probFIDA.*RFIDA.*dvpa.*dvpe./dlambdares;
                % The weight function is normalized by the associated
                % uncertainty
                Wtomo(Wrow,:)=reshape(wfcvvFIDA,1,numel(grid1))./...
                    Serror(Wrow);
            end
        elseif strcmpi(methodval,'2stepGRS') && strcmpi(pspace,'Ep')
            Ebroadening=1; % Spectral resolution of the measurements
            % divided by the bin width of the spectra
            % The step sizes in E and p are defined
            dE=pspacevec1(2)-pspacevec1(1);
            dp=pspacevec2(2)-pspacevec2(1);
            % The spectral resolution of the measurements
            dEgares=Ebroadening*(max(Spoint((Sbin+1):(Sbin+Sdim(i))))...
                -min(Spoint((Sbin+1):(Sbin+Sdim(i)))))/(Sdim(i)-1);
            % The values of vpa, vpe and vtot corresponding to the chosen
            % (E,p)-grid are calculated (eq. 2.3)
            vpara=grid2.*sqrt(2.*grid1.*Qe./(4*Mn));
            vperp=sqrt(1-grid2.^2).*sqrt(2.*grid1.*Qe./(4*Mn));
            vtotforwardgrid=sqrt(2.*grid1.*Qe./(4*Mn));
            % The rate part of the weight function
            SIGMA=zeros(size(grid1));
            for k =1:length(pspacevec2)
                for j=1:length(pspacevec1)
                    [~,index]=min(abs(Esigma-grid1(k,j)/1e6));
                    SIGMA(k,j)=sigma(index);
                end
            end
            SIGMAs=medfilt2(SIGMA);
            RGRS=n_Be*SIGMAs.*sqrt((vpara-vd).^2+vperp.^2)*10^(-31);
            % For each bin in the spectra, the weight function for the
            % chosen tomography grid is calculated using eq. 2.88
            % and stored as a row in a matrix
            for j=(Sbin+1):(Sbin+Sdim(i))
                Ega=Spoint(j);
                Wrow=Wrow+1;
                prob2stepGRS=zeros(size(vtotforwardgrid));
                % The integrals are calculated as sums
                for beta=betaloop
                    for zeta=zetaloop
                        % The two values of u_C (eq. 2.90)
                        uC1=((Ega-dEgares/2)/Ega0-1)*c;
                        uC2=((Ega+dEgares/2)/Ega0-1)*c;
                        
                        % The radicands of eq. 2.91 for the two u_C-values
                        radicand1=(sin(beta)^2*cos(zeta)^2*(64*uC1^2*...
                            vtotforwardgrid.^2*cos(beta)^2*(cos(beta)^2+...
                            sin(beta)^2*cos(zeta)^2)-(cos(beta)^2*...
                            (Qdot*Qe/(6*Mn)-vtotforwardgrid.^2)-13*...
                            uC1^2).^2))/(64*uC1^2*cos(beta)^2*...
                            (cos(beta)^2+sin(beta)^2*cos(zeta)^2)^2);
                        radicand2=(sin(beta)^2*cos(zeta)^2*(64*uC2^2*...
                            vtotforwardgrid.^2*cos(beta)^2*(cos(beta)^2+...
                            sin(beta)^2*cos(zeta)^2)-(cos(beta)^2*...
                            (Qdot*Qe/(6*Mn)-vtotforwardgrid.^2)-13*...
                            uC2^2).^2))/(64*uC2^2*cos(beta)^2*...
                            (cos(beta)^2+sin(beta)^2*cos(zeta)^2)^2);
                        
                        % Check for and exclude negative values of the
                        % radicands
                        r1positive=ones(size(radicand1));
                        r1positive(radicand1<0)=0;
                        r2positive=ones(size(radicand2));
                        r2positive(radicand2<0)=0;
                        
                        % The four possible values of u_alpha (eq. 2.91)
                        ua1plus=real((13*uC1^2+cos(beta)^2*...
                            (vtotforwardgrid.^2-Qdot*Qe/(6*Mn)))/(8*...
                            uC1*(cos(beta)^2+sin(beta)^2*cos(zeta)^2))+...
                            sqrt(radicand1));
                        ua2plus=real((13*uC2^2+cos(beta)^2*...
                            (vtotforwardgrid.^2-Qdot*Qe/(6*Mn)))/(8*...
                            uC2*(cos(beta)^2+sin(beta)^2*cos(zeta)^2))+...
                            sqrt(radicand2));
                        ua1minus=real((13*uC1^2+cos(beta)^2*...
                            (vtotforwardgrid.^2-Qdot*Qe/(6*Mn)))/(8*...
                            uC1*(cos(beta)^2+sin(beta)^2*cos(zeta)^2))-...
                            sqrt(radicand1));
                        ua2minus=real((13*uC2^2+cos(beta)^2*...
                            (vtotforwardgrid.^2-Qdot*Qe/(6*Mn)))/(8*...
                            uC2*(cos(beta)^2+sin(beta)^2*cos(zeta)^2))-...
                            sqrt(radicand2));
                        
                        % Check for and exclude solutions that do not keep
                        % eq. 2.89
                        check1plus=Qdot*Qe/(6*Mn)-vtotforwardgrid.^2-...
                            13*uC1^2/cos(beta)^2+8*(ua1plus*uC1+...
                            cos(zeta)*sqrt((vtotforwardgrid.^2-...
                            ua1plus.^2).*uC1^2*tan(beta)^2));
                        check1minus=Qdot*Qe/(6*Mn)-vtotforwardgrid.^2-...
                            13*uC1^2/cos(beta)^2+8*(ua1minus*uC1+...
                            cos(zeta)*sqrt((vtotforwardgrid.^2-...
                            ua1minus.^2).*uC1^2*tan(beta)^2));
                        check2plus=Qdot*Qe/(6*Mn)-vtotforwardgrid.^2-...
                            13*uC2^2/cos(beta)^2+8*(ua2plus*uC2+...
                            cos(zeta)*sqrt((vtotforwardgrid.^2-...
                            ua2plus.^2).*uC2^2*tan(beta)^2));
                        check2minus=Qdot*Qe/(6*Mn)-vtotforwardgrid.^2-...
                            13*uC2^2/cos(beta)^2+8*(ua2minus*uC2+...
                            cos(zeta)*sqrt((vtotforwardgrid.^2-...
                            ua2minus.^2).*uC2^2*tan(beta)^2));
                        % We only exclude solutions that do not keep eq.
                        % 2.89 by more than 9 m^2/s^2
                        checkzero1plus=ones(size(check1plus));
                        checkzero1plus(log10(abs(check1plus)+1)>1)=0;
                        checkzero1minus=ones(size(check1minus));
                        checkzero1minus(log10(abs(check1minus)+1)>1)=0;
                        checkzero2plus=ones(size(check2plus));
                        checkzero2plus(log10(abs(check2plus)+1)>1)=0;
                        checkzero2minus=ones(size(check2minus));
                        checkzero2minus(log10(abs(check2minus)+1)>1)=0;
                        
                        % The checks are combined into acceptance matrices
                        % that define the grid points that are acceptable. 
                        % The other grid points are set to zero.
                        acceptplus=r1positive.*r2positive.*...
                            checkzero1plus.*checkzero2plus;
                        acceptminus=r1positive.*r2positive.*...
                            checkzero1minus.*checkzero2minus;
                        
                        % The four gyro-angles (eq. 2.16)
                        gamma1plus=acos((ua1plus-vpara*cos(phi))./...
                            (vperp*sin(phi)));
                        gamma2plus=acos((ua2plus-vpara*cos(phi))./...
                            (vperp*sin(phi)));
                        gamma1minus=acos((ua1minus-vpara*cos(phi))./...
                            (vperp*sin(phi)));
                        gamma2minus=acos((ua2minus-vpara*cos(phi))./...
                            (vperp*sin(phi)));
                        % The probability part of the weight function (eq.
                        % 2.92) is checked using the acceptance matrices
                        prob=abs(real(gamma2plus-gamma1plus))/pi.*...
                            acceptplus+abs(real(gamma2minus-...
                            gamma1minus))/pi.*acceptminus;
                        prob2stepGRS=prob2stepGRS+prob;
                    end
                end
                prob2stepGRS=prob2stepGRS/(length(betaloop)*...
                    length(zetaloop));
                % The total weight function (eq. 2.88) times the step sizes
                wfcEp2stepGRSforward=RGRS.*prob2stepGRS/dEgares*dE*dp;
                Wtomo(Wrow,:)=reshape(wfcEp2stepGRSforward,1,...
                    numel(grid1))./Serror(Wrow);
            end
        elseif strcmpi(methodval,'2stepGRS') && strcmpi(pspace,'vv')
            Ebroadening=1; % Spectral resolution of the measurements
            % divided by the bin width of the spectra
            % The step sizes in vpa and vpe are defined
            dvpa=pspacevec1(2)-pspacevec1(1);
            dvpe=pspacevec2(2)-pspacevec2(1);
            % The spectral resolution of the measurements
            dEgares=Ebroadening*(max(Spoint((Sbin+1):(Sbin+Sdim(i))))...
                -min(Spoint((Sbin+1):(Sbin+Sdim(i)))))/(Sdim(i)-1);
            % The total velocities and alpha-particle energies are
            % calculated (eq. 2.2)
            vtotforwardgrid=sqrt(grid1.^2+grid2.^2);
            Ealpha=2*Mn*vtotforwardgrid.^2/Qe;
            % The rate part of the weight function
            SIGMA=zeros(size(grid1));
            for k =1:length(pspacevec2)
                for j=1:length(pspacevec1)
                    [~,index]=min(abs(Esigma-Ealpha(k,j)/1e6));
                    SIGMA(k,j)=sigma(index);
                end
            end
            SIGMAs=medfilt2(SIGMA);
            RGRS=n_Be*SIGMAs.*sqrt((grid1-vd).^2+grid2.^2)*10^(-31);
            % For each bin in the spectra, the weight function for the
            % chosen tomography grid is calculated using eq. 2.88
            % and stored as a row in a matrix
            for j=(Sbin+1):(Sbin+Sdim(i))
                Ega=Spoint(j);
                Wrow=Wrow+1;
                prob2stepGRS=zeros(size(vtotforwardgrid));
                % The integrals are calculated as sums
                for beta=betaloop
                    for zeta=zetaloop
                        % The two values of u_C (eq. 2.90)
                        uC1=((Ega-dEgares/2)/Ega0-1)*c;
                        uC2=((Ega+dEgares/2)/Ega0-1)*c;
                        
                        % The radicands of eq. 2.91 for the two u_C-values
                        radicand1=(sin(beta)^2*cos(zeta)^2*(64*uC1^2*...
                            vtotforwardgrid.^2*cos(beta)^2*(cos(beta)^2+...
                            sin(beta)^2*cos(zeta)^2)-(cos(beta)^2*...
                            (Qdot*Qe/(6*Mn)-vtotforwardgrid.^2)-13*...
                            uC1^2).^2))/(64*uC1^2*cos(beta)^2*...
                            (cos(beta)^2+sin(beta)^2*cos(zeta)^2)^2);
                        radicand2=(sin(beta)^2*cos(zeta)^2*(64*uC2^2*...
                            vtotforwardgrid.^2*cos(beta)^2*(cos(beta)^2+...
                            sin(beta)^2*cos(zeta)^2)-(cos(beta)^2*...
                            (Qdot*Qe/(6*Mn)-vtotforwardgrid.^2)-13*...
                            uC2^2).^2))/(64*uC2^2*cos(beta)^2*...
                            (cos(beta)^2+sin(beta)^2*cos(zeta)^2)^2);
                        
                        % Check for and exclude negative values of the
                        % radicands
                        r1positive=ones(size(radicand1));
                        r1positive(radicand1<0)=0;
                        r2positive=ones(size(radicand2));
                        r2positive(radicand2<0)=0;
                        
                        % The four possible values of u_alpha (eq. 2.91)
                        ua1plus=real((13*uC1^2+cos(beta)^2*(...
                            vtotforwardgrid.^2-Qdot*Qe/(6*Mn)))/(8*uC1*...
                            (cos(beta)^2+sin(beta)^2*cos(zeta)^2))+...
                            sqrt(radicand1));
                        ua2plus=real((13*uC2^2+cos(beta)^2*(...
                            vtotforwardgrid.^2-Qdot*Qe/(6*Mn)))/(8*uC2*...
                            (cos(beta)^2+sin(beta)^2*cos(zeta)^2))+...
                            sqrt(radicand2));
                        ua1minus=real((13*uC1^2+cos(beta)^2*(...
                            vtotforwardgrid.^2-Qdot*Qe/(6*Mn)))/(8*uC1*...
                            (cos(beta)^2+sin(beta)^2*cos(zeta)^2))-...
                            sqrt(radicand1));
                        ua2minus=real((13*uC2^2+cos(beta)^2*(...
                            vtotforwardgrid.^2-Qdot*Qe/(6*Mn)))/(8*uC2*...
                            (cos(beta)^2+sin(beta)^2*cos(zeta)^2))-...
                            sqrt(radicand2));
                        
                        % Check for and exclude solutions that do not keep
                        % eq. 2.89
                        check1plus=Qdot*Qe/(6*Mn)-vtotforwardgrid.^2-...
                            13*uC1^2/cos(beta)^2+8*(ua1plus*uC1+...
                            cos(zeta)*sqrt((vtotforwardgrid.^2-...
                            ua1plus.^2).*uC1^2*tan(beta)^2));
                        check1minus=Qdot*Qe/(6*Mn)-vtotforwardgrid.^2-...
                            13*uC1^2/cos(beta)^2+8*(ua1minus*uC1+...
                            cos(zeta)*sqrt((vtotforwardgrid.^2-...
                            ua1minus.^2).*uC1^2*tan(beta)^2));
                        check2plus=Qdot*Qe/(6*Mn)-vtotforwardgrid.^2-...
                            13*uC2^2/cos(beta)^2+8*(ua2plus*uC2+...
                            cos(zeta)*sqrt((vtotforwardgrid.^2-...
                            ua2plus.^2).*uC2^2*tan(beta)^2));
                        check2minus=Qdot*Qe/(6*Mn)-vtotforwardgrid.^2-...
                            13*uC2^2/cos(beta)^2+8*(ua2minus*uC2+...
                            cos(zeta)*sqrt((vtotforwardgrid.^2-...
                            ua2minus.^2).*uC2^2*tan(beta)^2));
                        % We only exclude solutions that do not keep eq.
                        % 2.89 by more than 9 m^2/s^2
                        checkzero1plus=ones(size(check1plus));
                        checkzero1plus(log10(abs(check1plus)+1)>1)=0;
                        checkzero1minus=ones(size(check1minus));
                        checkzero1minus(log10(abs(check1minus)+1)>1)=0;
                        checkzero2plus=ones(size(check2plus));
                        checkzero2plus(log10(abs(check2plus)+1)>1)=0;
                        checkzero2minus=ones(size(check2minus));
                        checkzero2minus(log10(abs(check2minus)+1)>1)=0;
                        
                        % The checks are combined into acceptance matrices
                        % that define the grid points that are acceptable. 
                        % The other grid points are set to zero.
                        acceptplus=r1positive.*r2positive.*...
                            checkzero1plus.*checkzero2plus;
                        acceptminus=r1positive.*r2positive.*...
                            checkzero1minus.*checkzero2minus;
                        
                        % The four gyro-angles (eq. 2.16)
                        gamma1plus=acos((ua1plus-grid1*cos(phi))./...
                            (grid2*sin(phi)));
                        gamma2plus=acos((ua2plus-grid1*cos(phi))./...
                            (grid2*sin(phi)));
                        gamma1minus=acos((ua1minus-grid1*cos(phi))./...
                            (grid2*sin(phi)));
                        gamma2minus=acos((ua2minus-grid1*cos(phi))./...
                            (grid2*sin(phi)));
                        % The probability part of the weight function (eq.
                        % 2.92) is checked using the acceptance matrices
                        prob=abs(real(gamma2plus-gamma1plus))/pi.*...
                            acceptplus+abs(real(gamma2minus-...
                            gamma1minus))/pi.*acceptminus;
                        prob2stepGRS=prob2stepGRS+prob;
                    end
                end
                prob2stepGRS=prob2stepGRS/(length(betaloop)*...
                    length(zetaloop));
                % The total weight function
                wfcvv2stepGRSforward=RGRS.*prob2stepGRS/dEgares*dvpa*dvpe;
                Wtomo(Wrow,:)=reshape(wfcvv2stepGRSforward,1,...
                    numel(grid1))./Serror(Wrow);
            end
        end
        % Variable used to keep track of the number of measurement points
        % per diagnostic
        Sbin=Sbin+Sdim(i);
    end
    % The calculated weight function is used to obtain the
    % reconstructed distribution functions
    if strcmpi(pspace,'Ep')
        xalpha=Tikhonov_reg(Wtomo,Svec,pspacevec1,pspacevec2,M_f,...
            Tikhonov,alpha,uselsqnonneg,grad_cutoff,'Ep');
    elseif strcmpi(pspace,'vv')
        xalpha=Tikhonov_reg(Wtomo,Svec,pspacevec1,pspacevec2,M_f,...
            Tikhonov,alpha,uselsqnonneg,grad_cutoff,'vv');
    end
    WTomo=Wtomo;
end
    
end

function [xalpha] =Tikhonov_reg(W,S,pspacevec1,pspacevec2,M_f,...
    Tikhonov,alpha,uselsqnonneg,gradient_cutoff,pspace)
% Number of points in the grid
n = size(W,2);

%0th order Tikhonov penalty matrix
L0 = eye(n);

% The gradient matrices are determined
if strcmpi(pspace,'Ep')
    [L1,L2] = gradient_matrix(pspacevec1,pspacevec2,M_f,'Ep');
elseif strcmpi(pspace,'vv')
    [L1,L2] = gradient_matrix(pspacevec1,pspacevec2,M_f,'vv');
end

% The regions where the weight function coverage is below the selected
% cutoff are determined
no_weight_coverage = find(sum(W) < gradient_cutoff*mean(sum(W)));

% The gradient is removed from the selected regions with low weight
% function coverage
for i = no_weight_coverage
        
    L1(i,:) = 0;
    L1(i,i) = 1; 
    L2(i,:) = 0;
    L2(i,i) = 1;
end

% The weight functions are rescaled
scaling_factor = 1/max(max(W));
W = W*scaling_factor;

% The penalty operator L^T*L is calculated (eqs. 2.11 and 2.12)
switch Tikhonov
    case 0 % Zeroth-order Tikhonov
        H = L0'*L0;
    case 1 % First-order Tikhonov
        H = L2'*L2 + L1'*L1;
    case 2 % Mix of zeroth- and first-order Tikhonov
        H = L2'*L2 + L1'*L1+1e-12*(L0')*L0;
    otherwise
        error('L must be 0, 1, 2 (0th or 1st order Tikhonov or a mix)')
end

% Pre-allocation for speed
xalpha = zeros(n,length(alpha));

% The reconstructed distribution function is calculated for every chosen 
% value of alpha and saved as a column in a matrix (eq. 2.10)
for i = 1:length(alpha)
    switch Tikhonov
        case 0 % Zeroth-order Tikhonov
            if uselsqnonneg % With the non-negativity constraint
                GalphaL0=zeros(size(L0,2),1);
                WalphaL=double([W; sqrt(alpha(i))*L0]);
                GalphaL=double([S; GalphaL0]);
                xalpha(:,i) = lsqnonneg(WalphaL,GalphaL);
            else % Without the non-negativity constraint
                xalpha(:,i) = (W'*W + alpha(i)*H)\W'*S;
            end
        case 1 % First-order Tikhonov
            if uselsqnonneg % With the non-negativity constraint
                GalphaL1=zeros(2*size(L1,2),1);
                WalphaL=double([W; sqrt(alpha(i))*L1; sqrt(alpha(i))*L2]);
                GalphaL=double([S; GalphaL1]);
                xalpha(:,i) = lsqnonneg(WalphaL,GalphaL);
            else % Without the non-negativity constraint
                xalpha(:,i) = (W'*W + alpha(i)*H)\W'*S;
            end
        case 2 % Mix of zeroth- and first-order Tikhonov
            if uselsqnonneg % With the non-negativity constraint
                GalphaL1=zeros(2*size(L1,2),1);
                GalphaL0=zeros(size(L0,2),1);
                WalphaL=double([W; sqrt(alpha(i))*L1;...
                    sqrt(alpha(i))*L2; 1e-8*sqrt(alpha(i))*L0]);
                GalphaL=double([S; GalphaL1; GalphaL0]);
                xalpha(:,i) = lsqnonneg(WalphaL,GalphaL);
            else % Without the non-negativity constraint
                xalpha(:,i) = (W'*W + alpha(i)*H)\W'*S;
            end
    end
end
% The previous rescaling is undone
xalpha = xalpha*scaling_factor;
end

function [L1, L2] = gradient_matrix(pspacevec1,pspacevec2,M_f,pspace)
% The number of points along each axis
xdim=length(pspacevec1);
ydim=length(pspacevec2);

if strcmpi(pspace,'Ep')
    pspacevec1=pspacevec1*1.6021917e-19; % Conversion to [J]
end

% The step sizes along both axes
dx=pspacevec1(2)-pspacevec1(1);
dy=pspacevec2(2)-pspacevec2(1);

% The gradient matrices are defined according to eqs. 2.11 and 2.12 and the
% definitions in the book "Parameter Estimation and Inverse Problems" by
% Richard Aster, Brian Borchers and Clifford Thurber
Delta_x=diag(ones(1,xdim*ydim))+ diag(-1*ones(1,(xdim-1)*ydim),-ydim);
Delta_x(1:ydim,:) = 0;

if strcmpi(pspace,'Ep')
    Delta_x=Delta_x/dx*sqrt(2*M_f);
    for i=1:xdim
        Delta_x((i-1)*ydim+1:i*ydim,:) = Delta_x((i-1)*ydim+1:i*ydim,:)*...
            sqrt(pspacevec1(i));
    end
elseif strcmpi(pspace,'vv')
    Delta_x=Delta_x/dx;
end

Dy_c0=zeros(xdim*ydim,1);
Dy_c1=zeros(xdim*ydim,1);

for i = 1:xdim
    Dy_c0((i-1)*ydim+1:i*ydim) = [-1*ones(1,ydim-1) 0];
    Dy_c1((i-1)*ydim+1:i*ydim) = [ones(1,ydim-1) 0];
end

Dy_c1 = Dy_c1(1:end-1);
Delta_y = diag(Dy_c0) + diag(Dy_c1,1);

if strcmpi(pspace,'Ep')
    Delta_y = Delta_y/dy*sqrt(M_f/2);
    for i =1:xdim
        for j =1:ydim
            Delta_y((i-1)*ydim + j,:) = Delta_y((i-1)*ydim + j,:)*...
                sqrt(1-pspacevec2(j)^2)/sqrt(pspacevec1(i));
        end
    end
elseif strcmpi(pspace,'vv')
    Delta_y = Delta_y/dy;
end
L1 =Delta_x;
L2 =Delta_y;
end

function sigma = Bosch_Hale_sigma(E,reaction)
% This script calculates the NES-related D-D or D-T or the GRS-related
% D-p fusion cross section using the approximation published by Bosch
% and Hale in NF 1992. E is the energy in the center-of-mass frame
% in units of keV. For a stationary target ion, the center-of-mass
% energy is simply the energy of the fast ion in the lab frame times
% the ratio between its mass and the total mass of the fast and
% stationary ions. The cross section is given in units of mb (milli-barns). 
% 1 mb = 10^-31 m^2
if strcmpi(reaction,'NESDD') % Constants for D-D fusion (Table 2.1)
    B_G = 31.3970;
    A1 = 5.3701e4;
    A2 = 3.3027e2;
    A3 = -1.2706e-1;
    A4 = 2.9327e-5;
    A5 = -2.5151e-9;
    B1 = 0;
    B2 = 0;
    B3 = 0;
    B4 = 0;
elseif strcmpi(reaction,'NESDT') % Constants for D-T fusion (Table 2.1)
    B_G = 34.3827;
    A1 = 6.927e4;
    A2 = 7.454e8;
    A3 = 2.050e6;
    A4 = 5.2002e4;
    A5 = 0;
    B1 = 6.38e1;
    B2 = -9.95e-1;
    B3 = 6.981e-5;
    B4 = 1.728e-4;
elseif strcmpi(reaction,'GRSDp') % Constants for p-D fusion (Table 2.1)
    E = E/1e3;
    B_G = 1.07;
    A1 = 8.09e-4;
    A2 = 1.92e-3;
    A3 = 1.21e-2;
    A4 = -5.26e-3;
    A5 = 6.52e-4;
    B1 = 0;
    B2 = 0;
    B3 = 0;
    B4 = 0;
end

S = (A1 + E.*(A2 + E.*(A3 + E.*(A4 + E.*A5))))./...
    (1+E.*(B1+E.*(B2+E.*(B3+E.*B4)))); % (eq. 2.53)
% The cross section is calculated (eq. 2.52)
sigma = S./(E.*exp(B_G./sqrt(E)));

end
