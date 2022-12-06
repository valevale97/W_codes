% This function is used to calculate a known distribution function,
% generate synthetic spectra based on this distribution function and then
% use the function TomoAnalyticRmix to reconstruct distribution functions
% from these synthetic spectra. The known distribution function and the
% reconstructed ones are all plotted along with the synthetic spectra with
% and without noise as well as the spectra obtained from the
% reconstructions. The calculated quality factor Q for each reconstruction
% is also plotted as a function of the regularization parameter.
% The equations used are explained or derived in the report
% "Analytic models in velocity-space tomography for fusion plasmas" by
% Andreas Poulsen. All references refer to this report.

clear;

dist =1; %1 for Bi-Maxwellian, 2 for NBI and [1 2] for a sum of these
pspace =2; %1 for (E,p)-space and 2 for (vpa,vpe)-space

%Choose the used method for each diagnostic. 1 for CTS, 2 for NES of D-D
%fusion, 3 for NES of D-T fusion with fast D, 4 for NES of D-T fusion with
%fast T, 5 for one-step reaction GRS of D-p fusion with fast p, 6 for
%one-step reaction GRS of D-p fusion with fast D, 7 for FIDA and 8 for
%two-step reaction GRS of Be-alpha fusion with fast alpha. #elements must 
%be 1 or equal to length(phivec)
method =2;

%observation angles, Sdim(i) measurements probing the target velocity space
%for view i,
phivec=[10 40 85];

%number of points in each spectrum: length must be equal to length(phivec)
% or 1, total number of measurements: sum(Sdim) or length(phivec)*Sdim if
% length(Sdim)=1
Sdim=100;

%The gradient used for the Tikhonov reconstruction is removed from the
%regions where the weight function coverage (sum of all weight functions)
%is less than gradient_cutoff times the mean of the weight function
%coverage. This leads to these regions being disregarded in the
%reconstructions
gradient_cutoff=0.00;

%natural constants
Mn   = 1.6749e-27; % [kg] Mass of a neutron
Mp   = 1.6726e-27; % [kg] Mass of a proton
Qe   = 1.6021917e-19; % [C] Elementary charge
c    = 3e8; % [m/s] The speed of light

%deuterium plasma
Mi=2*Mp; % [kg] Mass of a deuterium ion
M_He=3*Mp; % [kg] Mass of a helium-3 particle
n_t = 5e19; % [m^-3] The thermal ion density
vd=0; % [m/s] Toroidal drift of the bulk plasma

%plasma parameters used for the bi-Maxwellian distribution, 
%isotropic means Tpara=Tperp
Tpara=200e3; % [eV] Parallel plasma temperature
Tperp=1.5e6; % [eV] Perpendicular plasma temperature
ni=1e19;     % [m^-3] Ion density
vparadrift=0; % [m/s] Parallel drift velocity

%Parameters used for the NBI distribution
nibeam=1e19; % [m^-3] Beam ion density
Mibeam=2*Mp; % [kg] Mass of the beam ions
Ecritbeam=30e3; % [eV] 30 keV for 60 keV injection,
% 400 keV for 1 MeV injection
Ebirthbeamfull=60e3; %keV
vbirthwidthbeam=3e4; %3e4 for 60 keV injection, 2e5 for 1 MeV injection
pbirthbeam=0.5; % Pitch of the injection
pitchwidthbeam=0.15;

%noiselevel of the signal
noiselevel=0.03;

%0 does 0th order Tikhonov, 1-1st order, 2-mix of 0th and 1st order
Tikhonov=1;
%0 does unconstrained Tikhonov, 1 does non-negative Tikhonov
uselsqnonneg=0;

% Values of the regularization parameter alpha
switch Tikhonov
    case 0
        alpha = logspace(-2,2,9)';
    case 1
        alpha= logspace(8,13,9)';
    case 2
        alpha = logspace(8,11,9)';
    otherwise
        error('L must be 0, 1, 2 (0th or 1st order Tikhonov or a mix)')
end

%%%%parameters for (E,p)-space
%boundaries of the target (E,p)-space
Emin=5e2; % 0.1e6
Emax=80e3;
pmin=-0.99;
pmax=0.99;

%resolution of the tomography grid in (E,p)
Edimtomo=30;
pdimtomo=30;

%resolution of the forward model grid in (E,p)
Edimforward=40;
pdimforward=40;

%%%%parameters for (vpa,vpe)-space
%boundaries of the target (v_para,v_perp)-space
vparamax=10e6;
vparamin=-10e6;
vperpmin=1e3;
vperpmax=30e6;

%resolution of the tomography grid in (v_para,v_perp)
vparadimtomo=30;
vperpdimtomo=30;

%resolution of the forward model grid in (v_para,v_perp)
vparadimforward=40;
vperpdimforward=40;

%The possible methods
methodvec=["CTS", "NESDD", "NESDTfastD", "NESDTfastT", "GRSDpfastp",...
    "GRSDpfastD", "FIDA", "2stepGRS"];

%If method is a single value, it is changed to a vector of length equal
%to that of phivec where every element has that value
if length(method)==1
    method=ones(1,length(phivec))*method;
end
%The same is done to Sdim
if length(Sdim)==1
    Sdim=ones(1,length(phivec))*Sdim;
end

%A string array denoting the used methods is created
methods=strings(1,length(method));
for i=1:length(method)
    methods(i)=methodvec(method(i));
end

% The forward and tomography grids are defined based on the chosen plotting
% space via vectors that contain the values to be used in the grids. The
% grids are then used to calculate the values of the chosen 2D distribution
% on both the forward and the tomography grid. The axes and their labels
% are also defined

if pspace == 1 % (E,p)-space
    dEforward=(Emax-Emin)/(Edimforward-1);
    dpforward=(pmax-pmin)/(pdimforward-1);

    dEtomo=(Emax-Emin)/(Edimtomo-1);
    dptomo=(pmax-pmin)/(pdimtomo-1);

    Eforward = Emin:dEforward:Emax;
    pforward = pmin:dpforward:pmax;
    [Eforwardgrid,pforwardgrid] = meshgrid(Eforward,pforward);
    forwardgrid1=Eforwardgrid;
    forwardgrid2=pforwardgrid;

    Etomo = Emin:dEtomo:Emax;
    ptomo = pmin:dptomo:pmax;
    [Etomogrid,ptomogrid] = meshgrid(Etomo,ptomo);
    tomogrid1=Etomogrid;
    tomogrid2=ptomogrid;
        
    distsimforward=zeros(size(Eforwardgrid));
    distsimtomo=zeros(size(Etomogrid));
    % The model distribution can be bi-Maxwellian, slowing-down or a sum of
    % these
    for j =1:length(dist)
        if dist(j) ==1 % bi-Maxwellian distribution in (E,p)-space
            distsimforward=ni.*sqrt(Eforwardgrid./(pi.*Tpara.*Tperp.^2.*...
                Qe.^2)).*exp(-(pforwardgrid.^2.*Eforwardgrid.*Qe+0.5.*...
                Mi.*vparadrift.^2-vparadrift.*pforwardgrid.*sqrt(2.*...
                Mi.*Eforwardgrid.*Qe))/(Tpara.*Qe)-(1-pforwardgrid.^2).*...
                Eforwardgrid./Tperp).*Qe+distsimforward; % eq. 2.94
            
            distsimtomo=ni.*sqrt(Etomogrid./(pi.*Tpara.*Tperp.^2.*...
                Qe.^2)).*exp(-(ptomogrid.^2.*Etomogrid.*Qe+0.5.*Mi.*...
                vparadrift.^2-vparadrift.*ptomogrid.*sqrt(2.*Mi.*...
                Etomogrid.*Qe))/(Tpara.*Qe)-(1-ptomogrid.^2).*...
                Etomogrid./Tperp).*Qe+distsimtomo; % eq. 2.94
        elseif dist(j) ==2 % Slowing-down distribution in (E,p)-space
            % Supplied code
            vcritbeam=sqrt(2*Ecritbeam*Qe/Mibeam);
            vbirthbeamfull=sqrt(2*Ebirthbeamfull*Qe/Mibeam);
            vbirthbeamhalf=sqrt(2*Ebirthbeamfull/2*Qe/Mibeam);
            vbirthbeamthird=sqrt(2*Ebirthbeamfull/3*Qe/Mibeam);
            pitchkernel=exp(-(pforward-pbirthbeam).^2/pitchwidthbeam^2);
            intkernel=trapz(pforward,pitchkernel);
            pitchkernel=pitchkernel/intkernel;
            
            vmaxE=sqrt(2*Emax*Qe/Mi);
            vminE=sqrt(2*Emin*Qe/Mi);
            v=vminE:(vmaxE/200):vmaxE;
            [Vforward,Pforward2]=meshgrid(v,pforward);
            
            VRATIOBCFULL=((Vforward.^3/vbirthbeamfull^3).*(...
                vbirthbeamfull^3+vcritbeam^3)./(Vforward.^3+...
                vcritbeam^3)).^(1/6);
            VRATIOBCHALF=((Vforward.^3/vbirthbeamhalf^3).*(...
                vbirthbeamhalf^3+vcritbeam^3)./(Vforward.^3+...
                vcritbeam^3)).^(1/6);
            VRATIOBCTHIRD=((Vforward.^3/vbirthbeamthird^3).*(...
                vbirthbeamthird^3+vcritbeam^3)./(Vforward.^3+...
                vcritbeam^3)).^(1/6);
            SUMLEGENDREFULL=zeros(size(VRATIOBCFULL));
            SUMLEGENDREHALF=zeros(size(VRATIOBCFULL));
            SUMLEGENDRETHIRD=zeros(size(VRATIOBCFULL));
            
            for n=0:15
                legendrep=legendreP(n,Pforward2);
                legendrepbirth=trapz(pforward,pitchkernel.*...
                    legendreP(n,pforward));
                SUMLEGENDREFULL=SUMLEGENDREFULL+(n+0.5)*legendrep.*...
                    legendrepbirth.*VRATIOBCFULL.^(n*(n+1)).*...
                    erfc((Vforward-vbirthbeamfull)/vbirthwidthbeam)/2;
                SUMLEGENDREHALF=SUMLEGENDREHALF+(n+0.5)*legendrep.*...
                    legendrepbirth.*VRATIOBCHALF.^(n*(n+1)).*...
                    erfc((Vforward-vbirthbeamhalf)/vbirthwidthbeam)/2;
                SUMLEGENDRETHIRD=SUMLEGENDRETHIRD+(n+0.5)*legendrep.*...
                    legendrepbirth.*VRATIOBCTHIRD.^(n*(n+1)).*...
                    erfc((Vforward-vbirthbeamthird)/vbirthwidthbeam)/2;
            end
            SUMLEGENDREFULL(SUMLEGENDREFULL<0)=0;
            SUMLEGENDREHALF(SUMLEGENDREHALF<0)=0;
            SUMLEGENDRETHIRD(SUMLEGENDRETHIRD<0)=0;
            fvpbeam=Vforward.^2./(Vforward.^3+vcritbeam^3).*(0.5*...
                SUMLEGENDREFULL+0.3*SUMLEGENDREHALF+0.2*SUMLEGENDRETHIRD);
            
            v_on_E_grid = sqrt(2*Eforwardgrid*Qe/Mibeam);
            fEpbeam=interp2(v,pforward,fvpbeam./(Mibeam*Vforward),...
                v_on_E_grid,pforwardgrid,'linear',0);
            fEpbeam=nibeam/trapz(pforward,trapz(Eforward,fEpbeam'))*...
                fEpbeam;
            distsimforward=fEpbeam+distsimforward;
            
            pitchkernel=exp(-(ptomo-pbirthbeam).^2/pitchwidthbeam^2);
            intkernel=trapz(ptomo,pitchkernel);
            pitchkernel=pitchkernel/intkernel;
            [Vtomo,Ptomo2]=meshgrid(v,ptomo);
            
            VRATIOBCFULL=((Vtomo.^3/vbirthbeamfull^3).*(...
                vbirthbeamfull^3+vcritbeam^3)./(Vtomo.^3+...
                vcritbeam^3)).^(1/6);
            VRATIOBCHALF=((Vtomo.^3/vbirthbeamhalf^3).*(...
                vbirthbeamhalf^3+vcritbeam^3)./(Vtomo.^3+...
                vcritbeam^3)).^(1/6);
            VRATIOBCTHIRD=((Vtomo.^3/vbirthbeamthird^3).*(...
                vbirthbeamthird^3+vcritbeam^3)./(Vtomo.^3+...
                vcritbeam^3)).^(1/6);
            SUMLEGENDREFULL=zeros(size(VRATIOBCFULL));
            SUMLEGENDREHALF=zeros(size(VRATIOBCFULL));
            SUMLEGENDRETHIRD=zeros(size(VRATIOBCFULL));
            for n=0:15
                legendrep=legendreP(n,Ptomo2);
                legendrepbirth=trapz(ptomo,pitchkernel.*...
                    legendreP(n,ptomo));
                SUMLEGENDREFULL=SUMLEGENDREFULL+(n+0.5)*legendrep.*...
                    legendrepbirth.*VRATIOBCFULL.^(n*(n+1)).*...
                    erfc((Vtomo-vbirthbeamfull)/vbirthwidthbeam)/2;
                SUMLEGENDREHALF=SUMLEGENDREHALF+(n+0.5)*legendrep.*...
                    legendrepbirth.*VRATIOBCHALF.^(n*(n+1)).*...
                    erfc((Vtomo-vbirthbeamhalf)/vbirthwidthbeam)/2;
                SUMLEGENDRETHIRD=SUMLEGENDRETHIRD+(n+0.5)*legendrep.*...
                    legendrepbirth.*VRATIOBCTHIRD.^(n*(n+1)).*...
                    erfc((Vtomo-vbirthbeamthird)/vbirthwidthbeam)/2;
            end
            SUMLEGENDREFULL(SUMLEGENDREFULL<0)=0;
            SUMLEGENDREHALF(SUMLEGENDREHALF<0)=0;
            SUMLEGENDRETHIRD(SUMLEGENDRETHIRD<0)=0;
            fvpbeam=Vtomo.^2./(Vtomo.^3+vcritbeam^3).*(0.5*...
                SUMLEGENDREFULL+0.3*SUMLEGENDREHALF+0.2*SUMLEGENDRETHIRD);
            
            v_on_E_grid = sqrt(2*Etomogrid*Qe/Mibeam);
            fEpbeam=interp2(v,ptomo,fvpbeam./(Mibeam*Vtomo),v_on_E_grid,...
                ptomogrid,'linear',0);
            fEpbeam=nibeam/trapz(ptomo,trapz(Etomo,fEpbeam'))*fEpbeam;
            distsimtomo=fEpbeam+distsimtomo;

        end
    end
    
    % The labels and axes are defined
    if Emax>=1e6
        axes = [Emin/1e6 Emax/1e6 pmin pmax];
        xlab ='Energy [MeV]';
        xaxisforward = Eforward/1e6;
        xaxistomo = Etomo/1e6;
    elseif Emax<1e6
        axes = [0 Emax/1e3 pmin pmax];
        xlab ='Energy [keV]';
        xaxisforward = Eforward/1e3;
        xaxistomo = Etomo/1e3;
    end
    yaxisforward = pforward;
    yaxistomo = ptomo;
    ylab ='Pitch [-]';
elseif pspace == 2 % (vpa,vpe)-space
    dvparaforward=(vparamax-vparamin)/(vparadimforward-1);
    dvperpforward=(vperpmax-vperpmin)/(vperpdimforward-1);
        
    dvparatomo=(vparamax-vparamin)/(vparadimtomo-1);
    dvperptomo=(vperpmax-vperpmin)/(vperpdimtomo-1);

    vparaforward = vparamin:dvparaforward:vparamax;
    vperpforward = vperpmin:dvperpforward:vperpmax;
    [vparaforwardgrid,vperpforwardgrid] = meshgrid(vparaforward,...
        vperpforward);
    forwardgrid1=vparaforwardgrid;
    forwardgrid2=vperpforwardgrid;

    vparatomo = vparamin:dvparatomo:vparamax;
    vperptomo = vperpmin:dvperptomo:vperpmax;
    [vparatomogrid,vperptomogrid] = meshgrid(vparatomo,vperptomo);
    tomogrid1=vparatomogrid;
    tomogrid2=vperptomogrid;
    
    distsimforward=zeros(size(vparaforwardgrid));
    distsimtomo=zeros(size(vparatomogrid));
    % The model distribution can be bi-Maxwellian, slowing-down or a sum of
    % these
    for j=1:length(dist)
        if dist(j)==1 % bi-Maxwellian distribution in (vpa,vpe)-space
            %define thermal velocities
            vthpara=sqrt(2*Tpara*Qe/Mi);
            vthperp=sqrt(2*Tperp*Qe/Mi);
        
            fvpavpe3DbiMaxDriftforward=ni./(pi.^(3/2).*vthpara.*...
                vthperp.^2).*exp(-((vparaforwardgrid-vparadrift)./...
                vthpara).^2-(vperpforwardgrid./vthperp).^2);
        
            fvpavpe3DbiMaxDrifttomo=ni./(pi.^(3/2).*vthpara.*...
                vthperp.^2).*exp(-((vparatomogrid-vparadrift)./...
                vthpara).^2-(vperptomogrid./vthperp).^2); 
        
            distsimforward=fvpavpe3DbiMaxDriftforward.*2.*pi.*...
                vperpforwardgrid+distsimforward; % eq. 2.93
        
            distsimtomo=fvpavpe3DbiMaxDrifttomo.*2.*pi.*vperptomogrid+...
                distsimtomo; % eq. 2.93
        elseif dist(j)==2 % Slowing-down distribution in (vpa,vpe)-space
            % Supplied code
            pmin=-0.99;
            pmax=0.99;
            Emax=0.5*Mi*(vparamax^2+vperpmax^2)/Qe;
            Emin=0.5*Mi*vperpmin^2/Qe;
            pdimforward=vperpdimforward;
            Edimforward=vparadimforward;
            pdimtomo=vperpdimtomo;
            Edimtomo=vparadimtomo;
            
            dEforward=(Emax-Emin)/(Edimforward-1);
            dpforward=(pmax-pmin)/(pdimforward-1);

            dEtomo=(Emax-Emin)/(Edimtomo-1);
            dptomo=(pmax-pmin)/(pdimtomo-1);

            Eforward = Emin:dEforward:Emax;
            pforward = pmin:dpforward:pmax;
            [Eforwardgrid,pforwardgrid] = meshgrid(Eforward,pforward);

            Etomo = Emin:dEtomo:Emax;
            ptomo = pmin:dptomo:pmax;
            [Etomogrid,ptomogrid] = meshgrid(Etomo,ptomo);
            
            vcritbeam=sqrt(2*Ecritbeam*Qe/Mibeam);
            vbirthbeamfull=sqrt(2*Ebirthbeamfull*Qe/Mibeam);
            vbirthbeamhalf=sqrt(2*Ebirthbeamfull/2*Qe/Mibeam);
            vbirthbeamthird=sqrt(2*Ebirthbeamfull/3*Qe/Mibeam);
            pitchkernel=exp(-(pforward-pbirthbeam).^2/pitchwidthbeam^2);
            intkernel=trapz(pforward,pitchkernel);
            pitchkernel=pitchkernel/intkernel;
            
            vmax=sqrt(vparamax^2+vperpmax^2);
            vmin=vperpmin;
            v=vmin:(vmax/200):vmax;
            [Vforward,Pforward2]=meshgrid(v,pforward);
            
            VRATIOBCFULL=((Vforward.^3/vbirthbeamfull^3).*(...
                vbirthbeamfull^3+vcritbeam^3)./(Vforward.^3+...
                vcritbeam^3)).^(1/6);
            VRATIOBCHALF=((Vforward.^3/vbirthbeamhalf^3).*(...
                vbirthbeamhalf^3+vcritbeam^3)./(Vforward.^3+...
                vcritbeam^3)).^(1/6);
            VRATIOBCTHIRD=((Vforward.^3/vbirthbeamthird^3).*(...
                vbirthbeamthird^3+vcritbeam^3)./(Vforward.^3+...
                vcritbeam^3)).^(1/6);
            SUMLEGENDREFULL=zeros(size(VRATIOBCFULL));
            SUMLEGENDREHALF=zeros(size(VRATIOBCFULL));
            SUMLEGENDRETHIRD=zeros(size(VRATIOBCFULL));
            
            for n=0:15
                legendrep=legendreP(n,Pforward2);
                legendrepbirth=trapz(pforward,pitchkernel.*...
                    legendreP(n,pforward));
                SUMLEGENDREFULL=SUMLEGENDREFULL+(n+0.5)*legendrep.*...
                    legendrepbirth.*VRATIOBCFULL.^(n*(n+1)).*...
                    erfc((Vforward-vbirthbeamfull)/vbirthwidthbeam)/2;
                SUMLEGENDREHALF=SUMLEGENDREHALF+(n+0.5)*legendrep.*...
                    legendrepbirth.*VRATIOBCHALF.^(n*(n+1)).*...
                    erfc((Vforward-vbirthbeamhalf)/vbirthwidthbeam)/2;
                SUMLEGENDRETHIRD=SUMLEGENDRETHIRD+(n+0.5)*legendrep.*...
                    legendrepbirth.*VRATIOBCTHIRD.^(n*(n+1)).*...
                    erfc((Vforward-vbirthbeamthird)/vbirthwidthbeam)/2;
            end
            SUMLEGENDREFULL(SUMLEGENDREFULL<0)=0;
            SUMLEGENDREHALF(SUMLEGENDREHALF<0)=0;
            SUMLEGENDRETHIRD(SUMLEGENDRETHIRD<0)=0;
            fvpbeam=Vforward.^2./(Vforward.^3+vcritbeam^3).*(0.5*...
                SUMLEGENDREFULL+0.3*SUMLEGENDREHALF+0.2*SUMLEGENDRETHIRD);
            v_on_E_grid = sqrt(2*Eforwardgrid*Qe/Mibeam);
            fEpbeam=interp2(v,pforward,fvpbeam./(Mibeam*Vforward),...
                v_on_E_grid,pforwardgrid,'linear',0);
            fEpbeam=nibeam/trapz(pforward,trapz(Eforward,fEpbeam'))*...
                fEpbeam;
            JacobianEpvpavpe=Mi*vperpforwardgrid./...
                sqrt(vparaforwardgrid.^2+vperpforwardgrid.^2)/Qe;
            E_on_vpavpeMESH=0.5*Mi*(vparaforwardgrid.^2+...
                vperpforwardgrid.^2)/Qe;
            p_on_vpavpeMESH=vparaforwardgrid./sqrt(vparaforwardgrid.^2+...
                vperpforwardgrid.^2);
            fvpavpe2DtransformedNBI=interp2(Eforward,pforward,fEpbeam,...
                E_on_vpavpeMESH,p_on_vpavpeMESH,'*spline');
            fvpavpe2DtransformedNBI=fvpavpe2DtransformedNBI.*...
                JacobianEpvpavpe;
            distsimforward=distsimforward+fvpavpe2DtransformedNBI;
            
            pitchkernel=exp(-(ptomo-pbirthbeam).^2/pitchwidthbeam^2);
            intkernel=trapz(ptomo,pitchkernel);
            pitchkernel=pitchkernel/intkernel;
            [Vtomo,Ptomo2]=meshgrid(v,ptomo);
            
            VRATIOBCFULL=((Vtomo.^3/vbirthbeamfull^3).*(...
                vbirthbeamfull^3+vcritbeam^3)./(Vtomo.^3+...
                vcritbeam^3)).^(1/6);
            VRATIOBCHALF=((Vtomo.^3/vbirthbeamhalf^3).*(...
                vbirthbeamhalf^3+vcritbeam^3)./(Vtomo.^3+...
                vcritbeam^3)).^(1/6);
            VRATIOBCTHIRD=((Vtomo.^3/vbirthbeamthird^3).*(...
                vbirthbeamthird^3+vcritbeam^3)./(Vtomo.^3+...
                vcritbeam^3)).^(1/6);
            SUMLEGENDREFULL=zeros(size(VRATIOBCFULL));
            SUMLEGENDREHALF=zeros(size(VRATIOBCFULL));
            SUMLEGENDRETHIRD=zeros(size(VRATIOBCFULL));
            for n=0:15
                legendrep=legendreP(n,Ptomo2);
                legendrepbirth=trapz(ptomo,pitchkernel.*...
                    legendreP(n,ptomo));
                SUMLEGENDREFULL=SUMLEGENDREFULL+(n+0.5)*legendrep.*...
                    legendrepbirth.*VRATIOBCFULL.^(n*(n+1)).*...
                    erfc((Vtomo-vbirthbeamfull)/vbirthwidthbeam)/2;
                SUMLEGENDREHALF=SUMLEGENDREHALF+(n+0.5)*legendrep.*...
                    legendrepbirth.*VRATIOBCHALF.^(n*(n+1)).*...
                    erfc((Vtomo-vbirthbeamhalf)/vbirthwidthbeam)/2;
                SUMLEGENDRETHIRD=SUMLEGENDRETHIRD+(n+0.5)*legendrep.*...
                    legendrepbirth.*VRATIOBCTHIRD.^(n*(n+1)).*...
                    erfc((Vtomo-vbirthbeamthird)/vbirthwidthbeam)/2;
            end
            SUMLEGENDREFULL(SUMLEGENDREFULL<0)=0;
            SUMLEGENDREHALF(SUMLEGENDREHALF<0)=0;
            SUMLEGENDRETHIRD(SUMLEGENDRETHIRD<0)=0;
            fvpbeam=Vtomo.^2./(Vtomo.^3+vcritbeam^3).*(0.5*...
                SUMLEGENDREFULL+0.3*SUMLEGENDREHALF+0.2*SUMLEGENDRETHIRD);
            
            v_on_E_grid = sqrt(2*Etomogrid*Qe/Mibeam);
            fEpbeam=interp2(v,ptomo,fvpbeam./(Mibeam*Vtomo),v_on_E_grid,...
                ptomogrid,'linear',0);
            fEpbeam=nibeam/trapz(ptomo,trapz(Etomo,fEpbeam'))*fEpbeam;
            
            JacobianEpvpavpe=Mi*vperptomogrid./...
                sqrt(vparatomogrid.^2+vperptomogrid.^2)/Qe;
            E_on_vpavpeMESH=0.5*Mi*(vparatomogrid.^2+...
                vperptomogrid.^2)/Qe;
            p_on_vpavpeMESH=vparaforwardgrid./sqrt(vparatomogrid.^2+...
                vperptomogrid.^2);
            fvpavpe2DtransformedNBI=interp2(Etomo,ptomo,fEpbeam,...
                E_on_vpavpeMESH,p_on_vpavpeMESH,'*spline');
            fvpavpe2DtransformedNBI=fvpavpe2DtransformedNBI.*...
                JacobianEpvpavpe;
            distsimtomo=distsimtomo+fvpavpe2DtransformedNBI;
        end
    end
    
    % The axes and labels are defined
    axes=[vparamin vparamax vperpmin vperpmax]/1e6;
    xlab='v_{||} [10^6 m/s]';
    ylab='v_{\perp} [10^6 m/s]';
    xaxisforward=vparaforward/1e6;
    yaxisforward=vperpforward/1e6;
    xaxistomo=vparatomo/1e6;
    yaxistomo=vperptomo/1e6;
end

% Pre-allocation for speed
%Matrix to hold the weight functions
Wforward=zeros(sum(Sdim),numel(forwardgrid1)); 
Spoint=zeros(1,sum(Sdim)); %Vector to hold the spectrum bin values
Wrow=0;
Sbin=0;

% The weight functions for each individual diagnostic are calculated
% and stored in the designated matrix

for i=1:length(phivec)
    phi=phivec(i)*pi/180;
    
    % Parameters that differ for different methods
    if method(i) == 2 %NES for D-D fusion
        M_f=Mi; % [kg] Mass of the fast ion (deuterium)
        Q = 3.27e6; % [eV] Energy released in the D-D fusion reaction
        CM_coef = 1/2; % Coefficient for the fast-ion center-of-mass energy
        reaction = 'NESDD';
    elseif method(i) == 3 %NES for D-T fusion with fast D
        M_f =Mi; % [kg] Mass of the fast ion (deuterium)
        Q = 17.6e6; % [eV] Energy released in the D-T fusion reaction
        CM_coef = 2/5; % Coefficient for the fast-ion center-of-mass energy
        M_He = 4*Mp; % [kg] Mass of an alpha-particle
        reaction = 'NESDT';
    elseif method(i) == 4 %NES for D-T fusion with fast T
        M_f =3*Mp; % [kg] Mass of the fast ion (deuterium)
        Q = 17.6e6; % [eV] Energy released in the D-T fusion reaction
        CM_coef = 3/5; % Coefficient for the fast-ion center-of-mass energy
        M_He = 4*Mp; % [kg] Mass of an alpha-particle
        reaction = 'NESDT';
    elseif method(i) == 5 %One-step reaction GRS for D-p fusion with fast p
        Q = 5.5e6; % [eV] Energy released in the D-p fusion reaction
        CM_coef = 1/3; % Coefficient for the fast-ion center-of-mass energy
        M_f = Mp; % [kg] Mass of the fast ion (proton)
        M_pr = 3*Mp; % [kg] Mass of a helium-3 particle
        reaction = 'GRSDp';
    elseif method(i) == 6 %One-step reaction GRS for D-p fusion with fast D
        Q = 5.5e6; % [eV] Energy released in the D-p fusion reaction
        CM_coef = 2/3; % Coefficient for the fast-ion center-of-mass energy
        M_f = 2*Mp; % [kg] Mass of the fast ion (deuterium)
        M_pr = 3*Mp; % [kg] Mass of a helium-3 particle
        reaction = 'GRSDp';
    elseif method(i) == 7 % FIDA
        gamma_bar=0; % The phase-shift of the cosine model for the pdf
        B=1.74; % [T] The magnetic field
        E0=1e6; % [eV] The free parameter used in the model for the 
        % amplitude of the cosine model for the pdf
        lambda0=656.1e-9; % [m] The wavelength of the unshifted
        % D_alpha-line
        C_hat =[1,18,16,1681,2304,729,1936,5490,1936,729,2304,1681,...
            16,18,1]/18860; % Constants related to the probabilities of
        % the various Stark lines
        % Proportionality constants for the Stark splitting wavelength
        % shifts
        sl = [-220.2,-165.2,-137.7,-110.2,-82.64,-55.1,-27.56,0,27.57,...
            55.15,82.74,110.3,138.0,165.6,220.9]*1e-18; % [m^2/V]
        % Signs related to the various Stark lines. + for the
        % sigma-lines and - for the pi-lines
        Stark_sign = [-1,1,1,-1,-1,-1,1,1,1,-1,-1,-1,1,1,-1]; 
        % Parameters used for the ad-hoc model for the FIDA intensity
        % function R
        Rni=1e8; % [m^-3] Ion density
        RTpara=60e3; % [eV] Parallel plasma temperature
        RTperp=60e3; % [eV] Perpendicular plasma temperature
        Rvparadrift=0.4e6; % [m/s] Parallel drift velocity
    elseif method(i) == 8 % Two-step reaction GRS
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

    if method(i) ==1 && pspace ==1 % CTS in (E,p)-space
        ubroadening =3; % Spectral resolution of the measurements divided
        % by the bin width of the spectra
        % The value of the upper/lower limit of the spectra bins
        umax=sqrt(2*Emax*Qe/Mi); % [m/s] (eq. 2.21)
        du=(2*umax)/(Sdim(i)-1); % The bin width of the spectra
        dures=ubroadening*du; % Spectral resolution of the measurements
        uvec=(-umax:du:umax); % The bin values of the spectra
        % The calculated bin values are stored in the designated vector
        Spoint((Sbin+1):(Sbin+Sdim(i)))=uvec;
        % For each bin in the spectra, the weight function for the chosen
        % forward grid is calculated using eq. 2.19 and stored as a row in
        % a matrix
        for u=uvec
            Wrow=Wrow+1;
            % The two gyro-angles (eq. 2.20)
            gamma1forward = acos(1./sin(phi).*1./sqrt(1-...
                pforwardgrid.^2).*((u-dures./2)./sqrt(2.*...
                Eforwardgrid./Mi*Qe)-pforwardgrid.*cos(phi)));
            gamma2forward = acos(1./sin(phi).*1./sqrt(1-...
                pforwardgrid.^2).*((u+dures./2)./sqrt(2.*...
                Eforwardgrid./Mi*Qe)-pforwardgrid.*cos(phi)));
            % The weight function (eq. 2.19) times the step sizes
            wfcEpCTSforward=real((gamma1forward-gamma2forward))./pi./...
                dures.*dEforward.*dpforward; 
            Wforward(Wrow,:)=reshape(wfcEpCTSforward,1,...
                numel(Eforwardgrid));
        end
    elseif method(i) ==1 && pspace ==2 % CTS in (vpa,vpe)-space
        ubroadening =3; % Spectral resolution of the measurements divided
        % by the bin width of the spectra
        % The value of the upper/lower limit of the spectra bins
        if phi <= pi/2
            umin=vparamin*cos(phi)-vperpmax*sin(phi); % [m/s] (eq. 2.22)
            umax=vparamax*cos(phi)+vperpmax*sin(phi); % [m/s] (eq. 2.23)
        elseif phi > pi/2
            umin=vparamax*cos(phi)-vperpmax*sin(phi); % [m/s] (eq. 2.24)
            umax=vparamin*cos(phi)+vperpmax*sin(phi); % [m/s] (eq. 2.25)
        end
        du=(umax-umin)/(Sdim(i)-1); % The bin width of the spectra
        dures=ubroadening*du; % The spectral resolution of the measurements
        uvec=(umin:du:umax); % The bin values of the spectra
        % The calculated bin values are stored in the designated vector
        Spoint((Sbin+1):(Sbin+Sdim(i)))=uvec;
        % For each bin in the spectra, the weight function for the chosen
        % forward grid is calculated using eq. 2.19 and stored as a row in
        % a matrix
        for u=uvec
            Wrow=Wrow+1;
            % The two gyro-angles (eq. 2.16)
            gamma1forward = acos(1./(sin(phi).*vperpforwardgrid).*...
                (u-dures./2-vparaforwardgrid.*cos(phi)));
            gamma2forward = acos(1./(sin(phi).*vperpforwardgrid).*...
                (u+dures./2-vparaforwardgrid.*cos(phi)));
            % The weight function (eq. 2.19) times the step sizes
            wfcvvCTSforward=real((gamma1forward-gamma2forward))./pi./...
                dures.*dvparaforward.*dvperpforward; 
            Wforward(Wrow,:)=reshape(wfcvvCTSforward,1,...
                numel(vparaforwardgrid));
        end
    elseif (method(i)==2 || method(i)==3 || method(i)==4) && pspace ==1
        %NES in (E,p)-space
        Enbroadening=1; % Spectral resolution of the measurements divided
        % by the bin width of the spectra
        % The values of the upper/lower limits of the spectra bins
        % The lower limit of the spectra bins (eq. 2.62)
        Enmin=M_He*(M_He-M_f)*Q/((M_He+Mn)*(M_He-M_f)+Mn*M_f); % [eV]
        % The upper limit of the spectra bins (eqs. 2.70-2.73)
        A=4*(M_He+Mn)^2;
        B=-(8*M_He*(M_He+Mn)*Q+8*(M_He+Mn)*(M_He-M_f)*Emax+16*Mn*M_f*Emax);
        C=4*M_He^2*Q^2+8*M_He*(M_He-M_f)*Q*Emax+4*(M_He-M_f)^2*Emax^2;
        Enmax=(-B+sqrt(B^2-4*A*C))/(2*A); % [eV] 
        dEn=(Enmax-Enmin)/(Sdim(i)-1); % The bin width of the spectra
        dEnres=Enbroadening*dEn; % Spectral resolution of the measurements
        Envec=(Enmin:dEn:Enmax); % The bin values of the spectra
        % The calculated bin values are stored in the designated vector
        Spoint((Sbin+1):(Sbin+Sdim(i)))=Envec;
        % The values of vpa and vpe corresponding to the chosen
        % (E,p)-grid are calculated (eq. 2.3)
        vpara=pforwardgrid.*sqrt(2.*Eforwardgrid.*Qe./M_f);
        vperp=sqrt(1-pforwardgrid.^2).*sqrt(2.*Eforwardgrid.*Qe./M_f);
        % The center-of-mass energy of the fast ion in [keV]
        E_CM=CM_coef.*M_f.*((vpara-vd).^2+vperp.^2)./Qe/2e3;
        % The rate part of the weight function (eq. 2.51)
        RNES=n_t.*Bosch_Hale_sigma(E_CM,reaction).*10^-31.*...
            sqrt((vpara-vd).^2+vperp.^2);
        RNES(RNES<0)=0; % Negative rates are set to zero
        % For each bin in the spectra, the weight function for the chosen
        % forward grid is calculated using eq. 2.44 and stored as a row in
        % a matrix
        for En=Envec
            Wrow=Wrow+1;
            % The two gyro-angles (eq. 2.48)
            gamma1forward=acos(1./(sqrt(1-pforwardgrid.^2).*sin(phi)).*(...
                1./2.*(M_He+Mn)./sqrt(M_f.*Mn).*sqrt((En-dEnres./2)./...
                Eforwardgrid)-1./2.*(M_He-M_f)./sqrt(M_f.*Mn).*sqrt(...
                Eforwardgrid./(En-dEnres./2))-1./2.*M_He./sqrt(M_f.*...
                Mn).*Q./(sqrt(Eforwardgrid.*(En-dEnres./2)))-...
                pforwardgrid.*cos(phi)));
            gamma2forward=acos(1./(sqrt(1-pforwardgrid.^2).*sin(phi)).*(...
                1./2.*(M_He+Mn)./sqrt(M_f.*Mn).*sqrt((En+dEnres./2)./...
                Eforwardgrid)-1./2.*(M_He-M_f)./sqrt(M_f.*Mn).*sqrt(...
                Eforwardgrid./(En+dEnres./2))-1./2.*M_He./sqrt(M_f.*...
                Mn).*Q./(sqrt(Eforwardgrid.*(En+dEnres./2)))-...
                pforwardgrid.*cos(phi)));
            % The probability part of the weight function (eq. 2.45)
            probNES=real((gamma1forward-gamma2forward))./pi;
            % The total weight function (eq. 2.44) times the step sizes
            wfcEpNESforward=probNES.*RNES./dEnres.*...
                dEforward.*dpforward; 
            Wforward(Wrow,:)=reshape(wfcEpNESforward,1,...
                numel(Eforwardgrid));
        end
    elseif (method(i)==2 || method(i)==3 || method(i)==4) && pspace ==2
        % NES in (vpa,vpe)-space
        Enbroadening=1; % Spectral resolution of the measurements divided
        % by the bin width of the spectra
        % The values of the upper/lower limits of the spectra bins
        % The lower limit of the spectra bins (eq. 2.62)
        Enmin=M_He*(M_He-M_f)*Q/((M_He+Mn)*(M_He-M_f)+Mn*M_f); % [eV]
        % The maximum possible fast ion energy for the considered
        % velocity-space region
        EmaxW=0.5*M_f*(max(abs(vparatomo))^2+vperpmax^2)/Qe;
        % The upper limit of the spectra bins (eq. 2.70-2.73)
        A=4*(M_He+Mn)^2;
        B=-(8*M_He*(M_He+Mn)*Q+8*(M_He+Mn)*(M_He-M_f)*EmaxW+...
            16*Mn*M_f*EmaxW);
        C=4*M_He^2*Q^2+8*M_He*(M_He-M_f)*Q*EmaxW+4*(M_He-M_f)^2*EmaxW^2;
        Enmax=(-B+sqrt(B^2-4*A*C))/(2*A); % [eV]
        dEn=(Enmax-Enmin)/(Sdim(i)-1); % The bin width of the spectra
        dEnres=Enbroadening*dEn; % Spectral resolution of the measurements
        Envec=(Enmin:dEn:Enmax); % The bin values of the spectra
        % The calculated bin values are stored in the designated vector
        Spoint((Sbin+1):(Sbin+Sdim(i)))=Envec;
        % The center-of-mass energy of the fast ion in [keV]
        E_CM = CM_coef.*M_f.*((vparaforwardgrid-vd).^2+...
                vperpforwardgrid.^2)./Qe/2e3;
        % The rate part of the weight function (eq. 2.51)
        RNES=n_t.*Bosch_Hale_sigma(E_CM,reaction).*10^-31.*...
                sqrt((vparaforwardgrid-vd).^2+vperpforwardgrid.^2);
        RNES(RNES<0)=0; % Negative rates are set to zero
        % For each bin in the spectra, the weight function for the chosen
        % forward grid is calculated using eq. 2.44 and stored as a row in
        % a matrix
        for En=Envec
            Wrow=Wrow+1;
            % The two gyro-angles (eq. 2.47)
            gamma1forward=acos(1./(vperpforwardgrid.*sin(phi)).*...
                (0.5.*(M_He+Mn)./M_f.*sqrt(2.*(En-dEnres./2)./Mn.*Qe)-...
                0.5.*(M_He-M_f)./Mn.*(vparaforwardgrid.^2+...
                vperpforwardgrid.^2)./sqrt(2.*(En-dEnres./2)./Mn.*Qe)-...
                M_He./(M_f.*Mn).*Q.*Qe./sqrt(2.*(En-dEnres./2)./Mn.*Qe)-...
                vparaforwardgrid.*cos(phi)));
            gamma2forward=acos(1./(vperpforwardgrid.*sin(phi)).*...
                (0.5.*(M_He+Mn)./M_f.*sqrt(2.*(En+dEnres./2)./Mn.*Qe)-...
                0.5.*(M_He-M_f)./Mn.*(vparaforwardgrid.^2+...
                vperpforwardgrid.^2)./sqrt(2.*(En+dEnres./2)./Mn.*Qe)-...
                M_He./(M_f.*Mn).*Q.*Qe./sqrt(2.*(En+dEnres./2)./Mn.*Qe)-...
                vparaforwardgrid.*cos(phi)));
            % The probability part of the weight function (eq. 2.45)
            probNES=real(gamma1forward-gamma2forward)./pi;
            % The total weight function (eq. 2.44) times the step sizes
            wfcvvNESforward=probNES.*RNES./dEnres.*...
                dvparaforward.*dvperpforward;
            Wforward(Wrow,:)=reshape(wfcvvNESforward,1,...
                numel(vparaforwardgrid));
        end
    elseif (method(i)==5 || method(i)==6) && pspace==1
        % One-step reaction GRS in (E,p)-space
        Ebroadening=1; % Spectral resolution of the measurements divided by
        % the bin width of the spectra
        % The values of the upper/lower limits of the spectra bins
        % The lower limit of the spectra bins (eq. 2.85)
        Egamin=(sqrt(2*(M_pr-M_f)*c^2*Q*Qe+(M_pr-M_f)^2*c^4)-(M_pr-M_f)*...
            c^2)/Qe; % [eV] 
        syms x
        % The upper limit of the spectra bins (eq. 2.64 solution)
        Egamax=double(solve(Emax*Qe==0.5*M_f*(sqrt(M_pr/M_f*(x^2/...
            ((M_pr-M_f)^2*c^2)+2*(x-Q*Qe)/(M_pr-M_f)))-...
            x/((M_pr-M_f)*c))^2,x))/Qe; % [eV]
        dEga=(Egamax-Egamin)/(Sdim(i)-1); % The bin width of the spectra
        dEgares=Ebroadening*dEga; % Spectral resolution of the measurements
        Egavec=(Egamin:dEga:Egamax); % The bin values of the spectra
        % The calculated bin values are stored in the designated vector
        Spoint((Sbin+1):(Sbin+Sdim(i)))=Egavec;
        % The values of vpa and vpe corresponding to the chosen
        % (E,p)-grid are calculated (eq. 2.3)
        vpara=pforwardgrid.*sqrt(2.*Eforwardgrid.*Qe./M_f);
        vperp=sqrt(1-pforwardgrid.^2).*sqrt(2.*Eforwardgrid.*Qe./M_f);
        % The center-of-mass energy of the fast ion in [keV]
        E_CM=CM_coef.*M_f.*((vpara-vd).^2+vperp.^2)./Qe/2e3;
        % The rate part of the weight function (eq. 2.51)
        RGRS=n_t.*Bosch_Hale_sigma(E_CM,reaction).*10^-31.*...
            sqrt((vpara-vd).^2+vperp.^2);
        RGRS(RGRS<0)=0; % Negative rates are set to zero
        % For each bin in the spectra, the weight function for the chosen
        % forward grid is calculated using eq. 2.78 and stored as a row in
        % a matrix
        for Ega=Egavec
            Wrow=Wrow+1;
            % The two gyro-angles (eq. 2.82)
            gamma1forward=acos(1./(vperp.*sin(phi)).*((M_f-M_pr).*c./...
                (2.*(Ega-dEgares./2).*Qe).*(vpara.^2+vperp.^2)+...
                (Ega-dEgares./2).*Qe./(2.*M_f.*c)+M_pr.*c.*...
                (Ega-dEgares./2-Q)./(M_f.*(Ega-dEgares./2))-...
                vpara.*cos(phi)));
            gamma2forward=acos(1./(vperp.*sin(phi)).*((M_f-M_pr).*c./...
                (2.*(Ega+dEgares./2).*Qe).*(vpara.^2+vperp.^2)+...
                (Ega+dEgares./2).*Qe./(2.*M_f.*c)+M_pr.*c.*...
                (Ega+dEgares./2-Q)./(M_f.*(Ega+dEgares./2))-...
                vpara.*cos(phi)));
            % The probability part of the weight function (eq. 2.79)
            probGRS=real(gamma1forward-gamma2forward)./pi;
            % The total weight function (eq. 2.78) times the step sizes
            wfcEpGRSforward=probGRS.*RGRS./dEgares.*...
                dEforward.*dpforward;
            Wforward(Wrow,:)=reshape(wfcEpGRSforward,1,...
                numel(Eforwardgrid));
        end
    elseif (method(i)==5 || method(i)==6) && pspace==2
        % One-step reaction GRS in (vpa,vpe)-space
        Ebroadening=1; % Spectral resolution of the measurements divided by
        % the bin width of the spectra
        % The values of the upper/lower limits of the spectra bins
        % The lower limit of the spectra bins (eq. 2.85)
        Egamin=(sqrt(2*(M_pr-M_f)*c^2*Q*Qe+(M_pr-M_f)^2*c^4)-...
            (M_pr-M_f)*c^2)/Qe; % [eV]
        syms x
        % The largest possible fast ion energy for the considered region of
        % velocity-space
        EmaxW=0.5*M_f*(max(abs(vparatomo))^2+vperpmax^2);
        % The upper limit of the spectra bins (eq. 2.64 solution)
        Egamax=double(solve(EmaxW==0.5*M_f*(sqrt(M_pr/M_f*(x^2/...
            ((M_pr-M_f)^2*c^2)+2*(x-Q*Qe)/(M_pr-M_f)))-x/((M_pr-...
            M_f)*c))^2,x))/Qe; % [eV]
        dEga=(Egamax-Egamin)/(Sdim(i)-1); % The bin width of the spectra
        dEgares=Ebroadening*dEga; % Spectral resolution of the measurements
        Egavec=(Egamin:dEga:Egamax); % The bin values of the spectra
        % The calculated bin values are stored in the designated vector
        Spoint((Sbin+1):(Sbin+Sdim(i)))=Egavec;
        % The center-of-mass energy of the fast ion in [keV]
        E_CM = CM_coef.*M_f.*((vparaforwardgrid-vd).^2+...
            vperpforwardgrid.^2)./Qe/2e3;
        % The rate part of the weight function (eq. 2.51)
        RGRS=n_t.*Bosch_Hale_sigma(E_CM,reaction).*10^-31.*...
            sqrt((vparaforwardgrid-vd).^2+vperpforwardgrid.^2);
        RGRS(RGRS<0)=0; % Negative rates are set to zero
        % For each bin in the spectra, the weight function for the chosen
        % forward grid is calculated using eq. 2.78 and stored as a row in
        % a matrix
        for Ega=Egavec
            Wrow=Wrow+1;
            % The two gyro-angles (eq. 2.81)
            gamma1forward=acos(1./(vperpforwardgrid.*sin(phi)).*((M_f-...
                M_pr).*c./(2.*(Ega-dEgares./2).*Qe).*(...
                vparaforwardgrid.^2+vperpforwardgrid.^2)+(Ega-dEgares./...
                2).*Qe./(2.*M_f.*c)+M_pr.*c.*(Ega-dEgares./2-Q)./(M_f.*...
                (Ega-dEgares./2))-vparaforwardgrid.*cos(phi)));
            gamma2forward=acos(1./(vperpforwardgrid.*sin(phi)).*((M_f-...
                M_pr).*c./(2.*(Ega+dEgares./2).*Qe).*(...
                vparaforwardgrid.^2+vperpforwardgrid.^2)+(Ega+dEgares./...
                2).*Qe./(2.*M_f.*c)+M_pr.*c.*(Ega+dEgares./2-Q)./(M_f.*...
                (Ega+dEgares./2))-vparaforwardgrid.*cos(phi)));
            % The probability part of the weight function (eq. 2.79)
            probGRS=real(gamma1forward-gamma2forward)./pi;
            % The total weight function (eq. 2.78) times the step sizes
            wfcvvGRSforward=probGRS.*RGRS./dEgares.*...
                dvparaforward.*dvperpforward;
            Wforward(Wrow,:)=reshape(wfcvvGRSforward,1,...
                numel(vparaforwardgrid));
        end
    elseif method(i)==7 && pspace==1 % FIDA in (E,p)-space
        lambdabroadening=1; % Spectral resolution of the measurements 
        % divided by the bin width of the spectra
        % The values of vpa and vpe corresponding to the chosen (E,p)-grid
        % are calculated (eq. 2.3)
        vpara=pforwardgrid.*sqrt(2.*Eforwardgrid.*Qe./Mi);
        vperp=sqrt(1-pforwardgrid.^2).*sqrt(2.*Eforwardgrid.*Qe./Mi);
        % The values of the upper/lower limits of the spectra bins
        if phi <= pi/2
            % The lower limit (eq. 2.38)
            lambdamin = (lambda0+sl(1)*max(max(vperp))*B)*(1+1/c*(min(...
                min(vpara))*cos(phi)-max(max(vperp))*sin(phi))); % [m]
            % The upper limit (eq. 2.39)
            lambdamax = (lambda0+sl(15)*max(max(vperp)*B))*(1+1/c*(max(...
                max(vpara))*cos(phi)+max(max(vperp))*sin(phi))); % [m]
        elseif phi > pi/2
            % The lower limit (eq. 2.40)
            lambdamin = (lambda0+sl(1)*max(max(vperp))*B)*(1+1/c*(max(...
                max(vpara))*cos(phi)-max(max(vperp))*sin(phi))); % [m]
            % The upper limit (eq. 2.41)
            lambdamax = (lambda0+sl(15)*max(max(vperp))*B)*(1+1/c*(min(...
                min(vpara))*cos(phi)+max(max(vperp))*sin(phi))); % [m]
        end
        % The bin width of the spectra
        dlambda=(lambdamax-lambdamin)/(Sdim(i)-1);
        % Spectral resolution of the measurements
        dlambdares=lambdabroadening*dlambda;
        % The bin values of the spectra
        lambdavec=(lambdamin:dlambda:lambdamax);
        % The calculated bin values are stored in the designated vector
        Spoint((Sbin+1):(Sbin+Sdim(i)))=lambdavec;
        % The ad-hoc model for the FIDA intensity function R
        RFIDA = Rni.*sqrt(Eforwardgrid./(pi.*RTpara.*RTperp.^2.*Qe.^...
            2)).*exp(-(pforwardgrid.^2.*Eforwardgrid.*Qe+0.5.*Mi.*...
            Rvparadrift.^2-Rvparadrift.*pforwardgrid.*sqrt(2.*Mi.*...
            Eforwardgrid.*Qe))/(RTpara.*Qe)-(1-pforwardgrid.^2).*...
            Eforwardgrid./RTperp).*Qe;
        % For each bin in the spectra, the weight function for the chosen
        % forward grid is calculated using eq. 2.26 and stored as a row in
        % a matrix
        for lambda=lambdavec
            Wrow=Wrow+1;
            probFIDA=zeros(size(Eforwardgrid));
            for l=1:15 % sum over the 15 Stark lines
                % The two gyro-angles (eq. 2.31)
                gamma1=acos(1./(vperp.*sin(phi)).*(c.*((lambda-...
                    dlambdares./2)./(lambda0+sl(l).*vperp.*B)-1)-...
                    vpara.*cos(phi)));
                gamma2=acos(1./(vperp.*sin(phi)).*(c.*((lambda+...
                    dlambdares./2)./(lambda0+sl(l).*vperp.*B)-1)-...
                    vpara.*cos(phi)));
                % The probability part of the weight function (eq. 2.37)
                probFIDA=C_hat(l).*(real(gamma1-gamma2)/pi+...
                    Stark_sign(l).*sin(phi)^2/2.*(real(gamma1-gamma2)/...
                    pi-(sin(2.*real(gamma1))-sin(2.*real(gamma2)))/...
                    (2*pi))+2.*Eforwardgrid./E0.*(1-pforwardgrid.^2).*...
                    cos(gamma_bar).*(sin(real(gamma1))-sin(real(gamma2))...
                    +Stark_sign(l).*sin(phi)^2/3.*(sin(real(gamma1)).^3-...
                    sin(real(gamma2)).^3)))+probFIDA;
            end
            % The total weight function (eq. 2.26) times the step sizes
            wfcEpFIDAforward=probFIDA.*RFIDA.*dEforward.*dpforward./...
                dlambdares;
            Wforward(Wrow,:)=reshape(wfcEpFIDAforward,1,...
                numel(Eforwardgrid));
        end
    elseif method(i) ==7 && pspace ==2 % FIDA in (vpa,vpe)-space
        lambdabroadening=1; % Spectral resolution of the measurements
        % divided by the bin width of the spectra
        % The value of vpe0 corresponding to the chosen E0-value
        % is calculated
        vperp0=sqrt(2*E0/Mi*Qe);
        % The values of the upper/lower limits of the spectra bins
        if phi <= pi/2
            % The lower limit (eq. 2.38)
            lambdamin = (lambda0+sl(1)*vperpmax*B)*(1+1/c*(vparamin*...
                cos(phi)-vperpmax*sin(phi))); % [m]
            % The upper limit (eq. 2.39)
            lambdamax = (lambda0+sl(15)*vperpmax*B)*(1+1/c*(vparamax*...
                cos(phi)+vperpmax*sin(phi))); % [m] 
        elseif phi > pi/2
            % The lower limit (eq. 2.40)
            lambdamin = (lambda0+sl(1)*vperpmax*B)*(1+1/c*(vparamax*...
                cos(phi)-vperpmax*sin(phi))); % [m]
            % The upper limit (eq. 2.41)
            lambdamax = (lambda0+sl(15)*vperpmax*B)*(1+1/c*(vparamin*...
                cos(phi)+vperpmax*sin(phi))); % [m] 
        end
        % The bin width of the spectra
        dlambda=(lambdamax-lambdamin)/(Sdim(i)-1);
        % Spectral resolution of the measurements
        dlambdares=lambdabroadening*dlambda;
        % The bin values of the spectra
        lambdavec=(lambdamin:dlambda:lambdamax);
        % The calculated bin values are stored in the designated vector
        Spoint((Sbin+1):(Sbin+Sdim(i)))=lambdavec;
        % The ad-hoc model for the FIDA intensity function R
        vthpara=sqrt(2*RTpara*Qe/Mi);
        vthperp=sqrt(2*RTperp*Qe/Mi);
        RFIDA = Rni./(pi.^(3/2).*vthpara.*vthperp.^2).*...
            exp(-((vparaforwardgrid-Rvparadrift)./vthpara).^2-...
            (vperpforwardgrid./vthperp).^2).*2.*pi.*vperpforwardgrid;
        % For each bin in the spectra, the weight function for the chosen
        % forward grid is calculated using eq. 2.26 and stored as a row in
        % a matrix
        for lambda=lambdavec
            Wrow=Wrow+1;
            probFIDA=zeros(size(vparaforwardgrid));
            for l=1:15 % Sum over the 15 Stark lines
                % The two gyro-angles (eq. 2.31)
                gamma1=acos(1./(vperpforwardgrid.*sin(phi)).*(c.*...
                    ((lambda-dlambdares./2)./(lambda0+sl(l).*...
                    vperpforwardgrid.*B)-1)-vparaforwardgrid.*cos(phi)));
                gamma2=acos(1./(vperpforwardgrid.*sin(phi)).*(c.*...
                    ((lambda+dlambdares./2)./(lambda0+sl(l).*...
                    vperpforwardgrid.*B)-1)-vparaforwardgrid.*cos(phi)));
                % The probability part of the weight function (eq. 2.37)
                probFIDA=C_hat(l).*(real(gamma1-gamma2)/pi+...
                    Stark_sign(l).*sin(phi)^2/2.*(real(gamma1-gamma2)/...
                    pi-(sin(2.*real(gamma1))-sin(2.*real(gamma2)))/...
                    (2*pi))+2.*vperpforwardgrid.^2./vperp0^2.*...
                    cos(gamma_bar).*(sin(real(gamma1))-sin(real(gamma2))...
                    +Stark_sign(l).*sin(phi)^2/3.*(sin(real(gamma1)).^3-...
                    sin(real(gamma2)).^3)))+probFIDA;
            end
            % The total weight function (eq. 2.26) times the step sizes
            wfcvvFIDAforward=probFIDA.*RFIDA.*dvparaforward.*...
                dvperpforward./dlambdares;
            Wforward(Wrow,:)=reshape(wfcvvFIDAforward,1,...
                numel(vparaforwardgrid));
        end
    elseif method(i) ==8 && pspace ==1
        % Two-step reaction GRS in (E,p)-space
        Ebroadening=1; % Spectral resolution of the measurements divided by
        % the bin width of the spectra
        % The values of the upper/lower limits of the spectra bins
        Egamin=Ega0-60e3; % [eV] The lower limit of the spectra bins 
        Egamax=Ega0+60e3; % [eV] The upper limit of the spectra bins 
        dEga=(Egamax-Egamin)/(Sdim(i)-1); % The bin width of the spectra
        dEgares=Ebroadening*dEga; % Spectral resolution of the measurements
        Egavec=(Egamin:dEga:Egamax); % The bin values of the spectra
        % The calculated bin values are stored in the designated vector
        Spoint((Sbin+1):(Sbin+Sdim(i)))=Egavec;
        % The values of vpa, vpe and vtot corresponding to the chosen
        % (E,p)-grid are calculated (eq. 2.3)
        vpara=pforwardgrid.*sqrt(2.*Eforwardgrid.*Qe./(4*Mn));
        vperp=sqrt(1-pforwardgrid.^2).*sqrt(2.*Eforwardgrid.*Qe./(4*Mn));
        vtotforwardgrid=sqrt(2.*Eforwardgrid.*Qe./(4*Mn));
        % The rate part of the weight function (supplied code)
        SIGMA=zeros(size(Eforwardgrid));
        for k =1:length(pforward)
            for j=1:length(Eforward)
                [~,index]=min(abs(Esigma-Eforwardgrid(k,j)/1e6));
                SIGMA(k,j)=sigma(index);
            end
        end
        SIGMAs=medfilt2(SIGMA);
        RGRS=n_Be*SIGMAs.*sqrt((vpara-vd).^2+vperp.^2)*10^(-31);
        % For each bin in the spectra, the weight function for the chosen
        % forward grid is calculated using eq. 2.88 and stored as a row in
        % a matrix
        for Ega=Egavec
            Wrow=Wrow+1;
            prob2stepGRS=zeros(size(vtotforwardgrid));
            % The integrals are calculated as sums
            for beta=betaloop
                for zeta=zetaloop
                    % The two values of u_C (eq. 2.90)
                    uC1=((Ega-dEga/2)/Ega0-1)*c;
                    uC2=((Ega+dEga/2)/Ega0-1)*c;
                    % The radicands of eq. 2.91 for the two values of u_C
                    radicand1=(sin(beta)^2*cos(zeta)^2*(64*uC1^2*...
                        vtotforwardgrid.^2*cos(beta)^2*(cos(beta)^2+...
                        sin(beta)^2*cos(zeta)^2)-(cos(beta)^2*(Qdot*Qe/...
                        (6*Mn)-vtotforwardgrid.^2)-13*uC1^2).^2))/(64*...
                        uC1^2*cos(beta)^2*(cos(beta)^2+sin(beta)^2*...
                        cos(zeta)^2)^2);
                    radicand2=(sin(beta)^2*cos(zeta)^2*(64*uC2^2*...
                        vtotforwardgrid.^2*cos(beta)^2*(cos(beta)^2+...
                        sin(beta)^2*cos(zeta)^2)-(cos(beta)^2*(Qdot*Qe/...
                        (6*Mn)-vtotforwardgrid.^2)-13*uC2^2).^2))/(64*...
                        uC2^2*cos(beta)^2*(cos(beta)^2+sin(beta)^2*...
                        cos(zeta)^2)^2);
                    
                    %Check for and exclude negative values of the radicands
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
                    
                    % Check for and exclude solutions that do not keep eq.
                    % 2.89
                    check1plus=Qdot*Qe/(6*Mn)-vtotforwardgrid.^2-13*...
                        uC1^2/cos(beta)^2+8*(ua1plus*uC1+cos(zeta)*...
                        sqrt((vtotforwardgrid.^2-ua1plus.^2).*uC1^2*...
                        tan(beta)^2));
                    check1minus=Qdot*Qe/(6*Mn)-vtotforwardgrid.^2-13*...
                        uC1^2/cos(beta)^2+8*(ua1minus*uC1+cos(zeta)*...
                        sqrt((vtotforwardgrid.^2-ua1minus.^2).*uC1^2*...
                        tan(beta)^2));
                    check2plus=Qdot*Qe/(6*Mn)-vtotforwardgrid.^2-13*...
                        uC2^2/cos(beta)^2+8*(ua2plus*uC2+cos(zeta)*...
                        sqrt((vtotforwardgrid.^2-ua2plus.^2).*uC2^2*...
                        tan(beta)^2));
                    check2minus=Qdot*Qe/(6*Mn)-vtotforwardgrid.^2-13*...
                        uC2^2/cos(beta)^2+8*(ua2minus*uC2+cos(zeta)*...
                        sqrt((vtotforwardgrid.^2-ua2minus.^2).*uC2^2*...
                        tan(beta)^2));
                    % We only exclude solutions that do not keep eq. 2.89
                    % by more than 9 m^2/s^2
                    checkzero1plus=ones(size(check1plus));
                    checkzero1plus(log10(abs(check1plus)+1)>1)=0;
                    checkzero1minus=ones(size(check1minus));
                    checkzero1minus(log10(abs(check1minus)+1)>1)=0;
                    checkzero2plus=ones(size(check2plus));
                    checkzero2plus(log10(abs(check2plus)+1)>1)=0;
                    checkzero2minus=ones(size(check2minus));
                    checkzero2minus(log10(abs(check2minus)+1)>1)=0;
                    
                    % The checks are combined into acceptance matrices that
                    % define the grid points that are acceptable. The other
                    % grid points are set to zero.
                    acceptplus=r1positive.*r2positive.*checkzero1plus.*...
                        checkzero2plus;
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
                    % The probability part of the weight function
                    % (eq. 2.92) is checked using the acceptance matrices
                    prob=abs(real(gamma2plus-gamma1plus))/pi.*...
                        acceptplus+abs(real(gamma2minus-gamma1minus))/...
                        pi.*acceptminus;
                    prob2stepGRS=prob2stepGRS+prob;
                end
            end
            prob2stepGRS=prob2stepGRS/(length(betaloop)*length(zetaloop));
            % The total weight function (eq. 2.88) times the step sizes
            wfcEp2stepGRSforward=RGRS.*prob2stepGRS/dEga*dEforward*...
                dpforward;
            Wforward(Wrow,:)=reshape(wfcEp2stepGRSforward,1,...
                numel(Eforwardgrid));
        end
    elseif method(i) ==8 && pspace ==2
        % Two-step reaction GRS in (vpa,vpe)-space
        Ebroadening=1; % Spectral resolution of the measurements divided by
        % the bin width of the spectra
        % The values of the upper/lower limits of the spectra bins
        Egamin=Ega0-60e3; % [eV] The lower limit of the spectra bins
        Egamax=Ega0+60e3; % [eV] The upper limit of the spectra bins
        dEga=(Egamax-Egamin)/(Sdim(i)-1); % The bin width of the spectra
        dEgares=Ebroadening*dEga; % Spectral resolution of the measurements
        Egavec=(Egamin:dEga:Egamax); % The bin values of the spectra
        % The calculated bin values are stored in the designated vector
        Spoint((Sbin+1):(Sbin+Sdim(i)))=Egavec;
        % The total velocities and alpha-particle energies corresponding to
        % the chosen grid are calculated (eq. 2.2)
        vtotforwardgrid=sqrt(vparaforwardgrid.^2+vperpforwardgrid.^2);
        Ealpha=2*Mn*vtotforwardgrid.^2/Qe;
        % The rate part of the weight function (supplied code)
        SIGMA=zeros(size(vparaforwardgrid));
        for k =1:length(vperpforward)
            for j=1:length(vparaforward)
                [~,index]=min(abs(Esigma-Ealpha(k,j)/1e6));
                SIGMA(k,j)=sigma(index);
            end
        end
        SIGMAs=medfilt2(SIGMA);
        RGRS=n_Be*SIGMAs.*sqrt((vparaforwardgrid-vd).^2+...
            vperpforwardgrid.^2)*10^(-31);
        % For each bin in the spectra, the weight function for the chosen
        % forward grid is calculated using eq. 2.88 and stored as a row in
        % a matrix
        for Ega=Egavec
            Wrow=Wrow+1;
            prob2stepGRS=zeros(size(vtotforwardgrid));
            % The integrals are calculated as sums
            for beta=betaloop
                for zeta=zetaloop
                    % The two values of u_C (eq. 2.90)
                    uC1=((Ega-dEga/2)/Ega0-1)*c;
                    uC2=((Ega+dEga/2)/Ega0-1)*c;
                    
                    % The radicands of eq. 2.91 for the two values of u_C
                    radicand1=(sin(beta)^2*cos(zeta)^2*(64*uC1^2*...
                        vtotforwardgrid.^2*cos(beta)^2*(cos(beta)^2+...
                        sin(beta)^2*cos(zeta)^2)-(cos(beta)^2*(Qdot*Qe/...
                        (6*Mn)-vtotforwardgrid.^2)-13*uC1^2).^2))/(64*...
                        uC1^2*cos(beta)^2*(cos(beta)^2+sin(beta)^2*...
                        cos(zeta)^2)^2);
                    radicand2=(sin(beta)^2*cos(zeta)^2*(64*uC2^2*...
                        vtotforwardgrid.^2*cos(beta)^2*(cos(beta)^2+...
                        sin(beta)^2*cos(zeta)^2)-(cos(beta)^2*(Qdot*Qe/...
                        (6*Mn)-vtotforwardgrid.^2)-13*uC2^2).^2))/(64*...
                        uC2^2*cos(beta)^2*(cos(beta)^2+sin(beta)^2*...
                        cos(zeta)^2)^2);
                    
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
                    
                    % Check for and exclude solutions that do not keep eq.
                    % 2.89
                    check1plus=Qdot*Qe/(6*Mn)-vtotforwardgrid.^2-13*...
                        uC1^2/cos(beta)^2+8*(ua1plus*uC1+cos(zeta)*...
                        sqrt((vtotforwardgrid.^2-ua1plus.^2).*uC1^2*...
                        tan(beta)^2));
                    check1minus=Qdot*Qe/(6*Mn)-vtotforwardgrid.^2-13*...
                        uC1^2/cos(beta)^2+8*(ua1minus*uC1+cos(zeta)*...
                        sqrt((vtotforwardgrid.^2-ua1minus.^2).*uC1^2*...
                        tan(beta)^2));
                    check2plus=Qdot*Qe/(6*Mn)-vtotforwardgrid.^2-13*...
                        uC2^2/cos(beta)^2+8*(ua2plus*uC2+cos(zeta)*...
                        sqrt((vtotforwardgrid.^2-ua2plus.^2).*uC2^2*...
                        tan(beta)^2));
                    check2minus=Qdot*Qe/(6*Mn)-vtotforwardgrid.^2-13*...
                        uC2^2/cos(beta)^2+8*(ua2minus*uC2+cos(zeta)*...
                        sqrt((vtotforwardgrid.^2-ua2minus.^2).*uC2^2*...
                        tan(beta)^2));
                    % We only exclude solutions that do not keep eq. 2.89
                    % by more than 9 m^2/s^2
                    checkzero1plus=ones(size(check1plus));
                    checkzero1plus(log10(abs(check1plus)+1)>1)=0;
                    checkzero1minus=ones(size(check1minus));
                    checkzero1minus(log10(abs(check1minus)+1)>1)=0;
                    checkzero2plus=ones(size(check2plus));
                    checkzero2plus(log10(abs(check2plus)+1)>1)=0;
                    checkzero2minus=ones(size(check2minus));
                    checkzero2minus(log10(abs(check2minus)+1)>1)=0;
                    
                    % The checks are combined into acceptance matrices that
                    % define the grid points that are acceptable. The other
                    % grid points are set to zero.
                    acceptplus=r1positive.*r2positive.*checkzero1plus.*...
                        checkzero2plus;
                    acceptminus=r1positive.*r2positive.*...
                        checkzero1minus.*checkzero2minus;
                    
                    % The four gyro-angles (eq. 2.16)
                    gamma1plus=acos((ua1plus-vparaforwardgrid*...
                        cos(phi))./(vperpforwardgrid*sin(phi)));
                    gamma2plus=acos((ua2plus-vparaforwardgrid*...
                        cos(phi))./(vperpforwardgrid*sin(phi)));
                    gamma1minus=acos((ua1minus-vparaforwardgrid*...
                        cos(phi))./(vperpforwardgrid*sin(phi)));
                    gamma2minus=acos((ua2minus-vparaforwardgrid*...
                        cos(phi))./(vperpforwardgrid*sin(phi)));
                    % The probability part of the weight function 
                    % (eq. 2.92) is checked using the acceptance matrices
                    prob=abs(real(gamma2plus-gamma1plus))/pi.*...
                        acceptplus+abs(real(gamma2minus-gamma1minus))/...
                        pi.*acceptminus;
                    prob2stepGRS=prob2stepGRS+prob;
                end
            end
            prob2stepGRS=prob2stepGRS/(length(betaloop)*length(zetaloop));
            % The total weight function (eq. 2.88) times the step sizes
            wfcvv2stepGRSforward=RGRS.*prob2stepGRS/dEga*dvparaforward*...
                dvperpforward;
            Wforward(Wrow,:)=reshape(wfcvv2stepGRSforward,1,...
                numel(vparaforwardgrid));
        end
    end
    % Variable used to keep track of the number of measurement points per
    % diagnostic
    Sbin=Sbin+Sdim(i);
end
% Temporary synthetic spectra are calculated from the forward weight
% functions and the known distribution (eq. 2.7)
Stemp=Wforward*reshape(distsimforward,numel(forwardgrid1),1);

% Pre-allocation for speed
Serror=zeros(1,length(phivec));
% For each observation angle, the uncertainty of the spectrum values is
% estimated as 1/10 of the maximum value in the temporary spectrum for
% that observation angle, and the weight functions for that observation
% angle are normalized by this uncertainty
num=0;
for i=1:length(phivec)
    Serror(i)=0.1*max(Stemp((1+num):(num+Sdim(i))));
    Wforward((1+num):(num+Sdim(i)),:)=Wforward((1+num):(num+Sdim(i)),:)/...
        Serror(i);
    num=num+Sdim(i);
end

% Synthetic spectra are calculated from the normalized forward weight
% functions and the known distribution (eq. 2.7)
Sforward=Wforward*reshape(distsimforward,numel(forwardgrid1),1);
% Noise is added to the synthetic spectra
Sforwardnoise=Sforward+max(Sforward).*randn(sum(Sdim),1)...
    .*noiselevel;
% The normalized tomography weight functions and the reconstructed
% distribution functions are obtained by calling the function 
% TomoAnalyticRmix.
switch pspace
    case 1
        [Wtomo,xalpha]=TomoAnalyticRmix(Sforwardnoise,Spoint,Sdim,...
            Serror,Etomo,ptomo,vd,phivec,alpha,Tikhonov,uselsqnonneg,...
            gradient_cutoff,methods,'Ep');
    case 2
        [Wtomo,xalpha]=TomoAnalyticRmix(Sforwardnoise,Spoint,Sdim,...
            Serror,vparatomo,vperptomo,vd,phivec,alpha,Tikhonov,...
            uselsqnonneg,gradient_cutoff,methods,'vv');
end

% Spectra are calculated from the normalized tomography weight
% functions and the known distribution (eq. 2.7)
Stomo=Wtomo*reshape(distsimtomo,numel(tomogrid1),1);

% Pre-allocation for speed
balpha=zeros(sum(Sdim),length(alpha));
% A set of spectra is calculated from the normalized tomography weight
% functions and the various reconstructed distribution functions
for ialpha=1:length(alpha)
    balpha(:,ialpha)=Wtomo*xalpha(:,ialpha); % (eq. 2.7)
end

% The dimensions of the tomography grid
ydim=size(distsimtomo,1);
xdim=size(distsimtomo,2);
% Pre-allocation for speed
Qvec=zeros(1,length(alpha));
% The Q-factor for each reconstruction is calculated (eq. 2.13)
for i=1:length(alpha)
    Qvec(i)=sum(sum((distsimtomo-reshape(xalpha(:,i),ydim,xdim)).^2))/...
        sum(sum(distsimtomo.^2));
end

load('colormap_rel_diff.mat')

% All of the non-temporary spectra are plotted in a single figure
figure(10);clf;hold on;
plot(Sforward)
plot(Sforwardnoise)
plot(Stomo)
for ialpha=1:length(alpha)
    plot(balpha(:,ialpha))
end

% The chosen known distribution is plotted on the forward grid
fig1=figure(20); clf; hold on; box on;set(gca, 'Layer','top');
set(gcf,'DefaultLineLineWidth',2)
imagesc(xaxisforward,yaxisforward,distsimforward)
% caxis([0 280])
axis equal
axis(axes)
colorbar
colormap(flipud(hot))
set(gca,'Fontsize',16)
xlabel(xlab)
ylabel(ylab)

% The chosen known distribution is plotted on the tomography grid
fig2=figure(30); clf; hold on; box on;set(gca, 'Layer','top');
set(gcf,'DefaultLineLineWidth',2)
imagesc(xaxistomo,yaxistomo,distsimtomo)
% caxis([0 280])
%axis equal
axis(axes)
%pbaspect([1 3 1])
colorbar
colormap(flipud(hot))
set(gca,'Fontsize',16)
xlabel(xlab)
ylabel(ylab)

% The various reconstructed distribution functions are plotted in a set of
% subplots
figure(40); clf; hold on
for i=1:length(alpha)
    subplot(ceil(sqrt(length(alpha))),ceil(sqrt(length(alpha))),i)
    imagesc(xaxistomo,yaxistomo,reshape(xalpha(:,i),ydim,xdim))
    colorbar
    caxis([-max(abs(xalpha(:,i))) max(abs(xalpha(:,i)))])
    
    % caxis([0 120])
    %axis equal
    axis(axes)
    set(gca,'ydir','normal')
    xlabel(xlab)
    ylabel(ylab)
    set(gca,'Fontsize',13)
end
set(gcf,'colormap',colormap_rel_diff)

% The Q-factors are plotted as a function of the alpha-values
figure(50); clf;
loglog(alpha,Qvec)
xlabel('Regularization parameter, \alpha')
ylabel('Quality factor, Q')

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
