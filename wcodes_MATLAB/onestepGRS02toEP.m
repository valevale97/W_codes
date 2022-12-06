%one-step reactions
% constants
Qe=1.6e-19;
mp=1.67e-27;
c=3e8;
h=6.6e-34;
Ep=mp*c^2/Qe;


mprod=3*mp;
mf=mp;
Q=5.5e6*Qe;

Etyp=0.5*mf*(2e7)^2;
EtypMeV=Etyp/1e6/Qe;%MeV;
Momtyp=mf*1e7;

if 0
    mprod=5*mp;
    mf=2*mp;
    Q=16.85e6*Qe;
    %    Q=16.3e6*Qe;
    Qstar=[12:0.001:22];
    Qwidth=0.685;
    Lorentz=1/pi*0.5*Qwidth./((Q/1e6/Qe-Qstar).^2+(0.5*Qwidth)^2);
    trapz(Qstar,Lorentz)
    figure; clf; hold on;
    plot(Qstar, Lorentz)
end

if 0
    mprod=4*mp;
    mf=2*mp;
    Q=23.8e6*Qe;
end

if 0
    mprod=4*mp;
    mf=mp;
    Q=19.7e6*Qe;
end

if 0
    mprod=1.0002*mp;
    mf=mp;
    Q=2e3*Qe;
end

nt=2e19;




Qvelocity=sqrt(Q/(6*mp));
Eproton=mp*c^2/Qe/1e6;
Egammamin=(sqrt(((mprod-mf)*c^2)^2+2*(mprod-mf)*c^2*Q)...
    -(mprod-mf)*c^2)/1e6/Qe;



%EalphaQvelocity=0.5*4*Mn*Qvelocity^2/Qe/1e6;

%grid
vpara=[-20:0.1:20]*1e6;
vperp=[0.2:0.1:20]*1e6;
gamma=[0:0.1:360]/180*pi;

[VPARA,VPERP]=meshgrid(vpara,vperp);
VMAG=sqrt(VPARA.^2+VPERP.^2);

savefigs=0;
fs=28;
ax=[-15.5 15.5 -0.5 15.5];

EFASTMeV=0.5*mf*VMAG.^2/Qe/1e6;

betagamow=1.07;
A0=8.09e-4;
A1=1.92e-3;
A2=1.21e-2;
A3=-5.26e-3;
A4=6.52e-4;
Sastro=zeros(size(EFASTMeV));
Sastro=A0+A1*EFASTMeV+A2*EFASTMeV.^2+A3*EFASTMeV.^3+A4*EFASTMeV.^4;
sigma=Sastro./EFASTMeV.*exp(-betagamow./sqrt(EFASTMeV))*10^-31;
R=nt*sigma.*VMAG;
% figure(5); clf; hold on; box on;
% plot(log10(EFASTMeV(:,101)),log10(sigma(:,101)))

if 0
fig=figure(10); clf; hold on; box on;
set(gcf,'DefaultLineLineWidth',4)
[dummy,H]=contourf(vpara/1e6,vperp/1e6, log10(R),40);
for i = 1:length(H) set(H(i),'Edgecolor','none');end
axis equal
axis(ax)
%caxis([-8 -2])
set(gca,'fontsize',fs,'fontweight','bold')
set(gca,'Position', [0.22 0.2 0.6 0.65])
%title(strcat('n, \phi=',num2str(phi/pi*180),', E\gamma=',num2str(Egamma1/1e6/Qe)),'fontsize',fs,'fontweight','bold')
xlabel('v_{||} [10^6 m/s]','fontsize',fs,'fontweight','bold')
ylabel('v_{\perp} [10^6 m/s]','fontsize',fs,'fontweight','bold')
set(gca,'ytick',[0 10 20])
set(gca,'xtick',[-20 -10 0 10 20])
%         set(gca,'yticklabel',[3 2 1 0 1 2 3])
%caxis([-14 -7])
handlecolorbar=colorbar;%('horiz');
colormapInverseHotBWwfc;
set(handlecolorbar,'fontsize',fs,'fontweight','bold');
set(handlecolorbar,'Position',[0.85 0.3 0.03 0.43]);
if savefigs
    saveas(fig,strcat('figs\R.jpg'));
    saveas(fig,strcat('figs\R.epsc'));
end
end
Egammarange=[5:0.01:7]*1e6*Qe;
sumprobfcanalytic=zeros(size(VPARA));
sumwfcanalytic=zeros(size(VPARA));
sumwfcanalyticnormalized=zeros(size(VPARA));
%intgration limits

%for Egamma1=[5.55 5.65 5.8 6]*1e6*Qe;
for Egamma1=[5.7 5.85 6.05 6.3]*1e6*Qe;
    %for Egamma1=[5.55 6]*1e6*Qe;
    %for Egamma1=[23.85 24.3]*1e6*Qe;
    %for Egamma1=[19.75 19.9 20.1 20.4]*1e6*Qe;
    
    %for Egamma1=[19.75 20.2]*1e6*Qe;
    
    %for Egamma1=1.02*Q
%    for Egamma1=[2.01]*1e3*Qe;
    %for Egamma1=[17 17.4 17.8]*1e6*Qe;
    %for Egamma1=[17.2]*1e6*Qe
    %   for Egamma2=[Egamma1+0.02*1e6*Qe Egamma1+0.03*1e6*Qe Egamma1+0.05*1e6*Qe...
    %       Egamma1+0.1*1e6*Qe Egamma1+0.2*1e6*Qe Egamma1+0.3*1e6*Qe];
    
   for Egamma2=Egamma1+0.01*1e6*Qe;
  %  for Egamma2=1.001*Egamma1    
        [dummy,Epickindex]=min(abs(Egammarange-(Egamma1+Egamma2)/2));
        
        for phi=[90]/180*pi
            
            if 0
                probfcnum=zeros(length(vpara),length(vperp),length(Egammarange));
                for ipara=1:length(vpara)
                    %    vpara(ipara)
                    for iperp=1:length(vperp)
                        %           Ealpha=0.5*4*Mn*(vpara(ipara)^2+vperp(iperp)^2)/Qe/1e6;
                        uspectrum=vpara(ipara)*cos(phi)+vperp(iperp)*sin(phi)*cos(gamma);
                        vmag=sqrt(vpara(ipara)^2+vperp(iperp)^2);
                        minusphalf=mf*c*uspectrum-mHe3*c^2;
                        q=-(2*mHe3*c^2*Q+mf*(mHe3-mf)*c^2*vmag^2);
                        Egamma=(minusphalf+sqrt(minusphalf.^2-q));
                        histEgamma=histc(Egamma,Egammarange);
                        probfcnum(ipara,iperp,:)=histEgamma/sum(histEgamma);
                    end
                end
                
                figure(20); clf; hold on;
                plot(gamma,Egamma)
                
                for iwfc=Epickindex
                    probfcnumpick=squeeze(probfcnum(:,:,iwfc));
                    wfcnumpick=probfcnumpick'.*R;
                    
                    fig=figure(50); clf; hold on; box on;
                    set(gcf,'DefaultLineLineWidth',4)
                    [dummy,H]=contourf(vpara/1e6,vperp/1e6, probfcnumpick',40);
                    for i = 1:length(H) set(H(i),'Edgecolor','none');end
                    axis equal
                    axis(ax)
                    caxis([0 1])
                    set(gca,'fontsize',fs,'fontweight','bold')
                    set(gca,'Position', [0.22 0.2 0.6 0.65])
                    title(strcat('n, \phi=',num2str(phi/pi*180),', E\gamma=',num2str(Egamma1/1e6/Qe)),'fontsize',fs,'fontweight','bold')
                    xlabel('v_{||} [10^6 m/s]','fontsize',fs,'fontweight','bold')
                    ylabel('v_{\perp} [10^6 m/s]','fontsize',fs,'fontweight','bold')
                    set(gca,'ytick',[0 10 20])
                    set(gca,'xtick',[-20 -10 0 10 20])
                    %         set(gca,'yticklabel',[3 2 1 0 1 2 3])
                    %         caxis([-3 0])
                    handlecolorbar=colorbar;%('horiz');
                    colormapInverseHotBWwfc;
                    set(handlecolorbar,'fontsize',fs,'fontweight','bold');
                    set(handlecolorbar,'Position',[0.85 0.3 0.03 0.43]);
                    % set(handlecolorbar,'zticklabel',[-3 -2 -1 0])
                    %             if savefigs
                    saveas(fig,strcat('figs\probfcEgam1_',num2str(Egamma1/1e6/Qe),'Egam2_',num2str(Egamma2/Qe/1e6),'p',num2str(phi/pi*180),'n.jpg'));
                    saveas(fig,strcat('figs\probfcEgam1_',num2str(Egamma1/1e6/Qe),'Egam2_',num2str(Egamma2/1e6/Qe),'p',num2str(phi/pi*180),'n.epsc'));
                end
                
                for ipara=1:length(vpara)
                    for iperp=1:length(vperp)
                        if wfcnumpick(iperp,ipara)~=0
                            wfcnumpick(iperp,ipara)=log10(wfcnumpick(iperp,ipara));
                        else
                            wfcanumpick(iperp,ipara)=NaN;
                        end
                    end
                end
                
                
                
                fig=figure(51); clf; hold on; box on;
                set(gcf,'DefaultLineLineWidth',4)
                [dummy,H]=contourf(vpara/1e6,vperp/1e6, wfcnumpick,40);
                for i = 1:length(H) set(H(i),'Edgecolor','none');end
                axis equal
                axis(ax)
                %caxis([0 1])
                set(gca,'fontsize',fs,'fontweight','bold')
                set(gca,'Position', [0.22 0.2 0.6 0.65])
                title(strcat('n, \phi=',num2str(phi/pi*180),', E\gamma=',num2str(Egamma1/1e6/Qe)),'fontsize',fs,'fontweight','bold')
                xlabel('v_{||} [10^6 m/s]','fontsize',fs,'fontweight','bold')
                ylabel('v_{\perp} [10^6 m/s]','fontsize',fs,'fontweight','bold')
                set(gca,'ytick',[0 10 20])
                set(gca,'xtick',[-20 -10 0 10 20])
                %         set(gca,'yticklabel',[3 2 1 0 1 2 3])
                %caxis([-14 -7])
                handlecolorbar=colorbar;%('horiz');
                colormap(flipud(bone));
                set(handlecolorbar,'fontsize',fs,'fontweight','bold');
                set(handlecolorbar,'Position',[0.85 0.3 0.03 0.43]);
                % set(handlecolorbar,'zticklabel',[-3 -2 -1 0])
                if savefigs
                    saveas(fig,strcat('figs\wfcEgam1_',num2str(Egamma1/1e6/Qe),'Egam2_',num2str(Egamma2/Qe/1e6),'p',num2str(phi/pi*180),'n.jpg'));
                    saveas(fig,strcat('figs\wfcEgam1_',num2str(Egamma1/1e6/Qe),'Egam2_',num2str(Egamma2/1e6/Qe),'p',num2str(phi/pi*180),'n.epsc'));
                end
                
            end
            
            u1=(mf-mprod)*c/(2*Egamma1)*VMAG.^2+Egamma1/(2*mf*c)...
                +mprod*c*(Egamma1-Q)/(mf*Egamma1);
            u2=(mf-mprod)*c/(2*Egamma2)*VMAG.^2+Egamma2/(2*mf*c)...
                +mprod*c*(Egamma2-Q)/(mf*Egamma2);
            
            gamma1=real(acos((u1-VPARA*cos(phi))./(VPERP*sin(phi))));
            gamma2=real(acos((u2-VPARA*cos(phi))./(VPERP*sin(phi))));
            probfcanalytic=(gamma1-gamma2)/pi;
            wfcanalytic=probfcanalytic.*R;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%normalize
            wfcanalyticnormalized=wfcanalytic/max(max(wfcanalytic));
            sumwfcanalyticnormalized=sumwfcanalyticnormalized+wfcanalyticnormalized;
            
            
            %plasma parameters, isotropic means Tpara=Tperp
            Tpara=15; %keV
            Tperp=150; %keV
            ne=6e15;     %1/m^3
            
            %define thermal velocities
            vthpara=sqrt(2*Tpara*Qe*1000/mf);
            vthperp=sqrt(2*Tperp*Qe*1000/mf);
            
            %bi-Maxwellian distribution function
            fvpavpe3DbiMax=ne/(pi^(3/2)*vthpara*vthperp^2)*exp(-(VPARA/vthpara).^2-(VPERP/vthperp).^2);
            
            %fvpavpe in 2D
            fvpavpe2DbiMax=fvpavpe3DbiMax.*(2*pi*VPERP);
            
            sumprobfcanalytic=sumprobfcanalytic+probfcanalytic;
            sumwfcanalytic=sumprobfcanalytic.*R;
            sumorigin=fvpavpe2DbiMax.*sumwfcanalytic;
            
            for ipara=1:length(vpara)
                for iperp=1:length(vperp)
                    if wfcanalytic(iperp,ipara)~=0
                        wfcanalytic(iperp,ipara)=log10(wfcanalytic(iperp,ipara));
                    else
                        wfcanalytic(iperp,ipara)=NaN;
                    end
                end
            end
            
            for ipara=1:length(vpara)
                for iperp=1:length(vperp)
                    if probfcanalytic(iperp,ipara)~=0
                        probfcanalytic(iperp,ipara)=log10(probfcanalytic(iperp,ipara));
                    else
                        probfcanalytic(iperp,ipara)=NaN;
                    end
                end
            end
            %
                       
            figure(1000); clf; hold on; box on;
            set(gcf,'DefaultLineLineWidth',4)
            [dummy,H]=contourf(vpara/1e6,vperp/1e6,probfcanalytic,40);
            for i = 1:length(H) set(H(i),'Edgecolor','none');end
            axis equal
            axis(ax)
            %caxis([-4 0])
            set(gca,'fontsize',fs,'fontweight','bold')
            set(gca,'Position', [0.22 0.2 0.6 0.65])
            %title(strcat('a, \phi=',num2str(phi/pi*180),', E\gamma=',num2str(Egamma1/1e6/Qe)),'fontsize',fs,'fontweight','bold')
            xlabel('v_{||} [10^6 m/s]','fontsize',fs,'fontweight','bold')
            ylabel('v_{\perp} [10^6 m/s]','fontsize',fs,'fontweight','bold')
            set(gca,'ytick',[0 10 20])
            set(gca,'xtick',[-20 -10 0 10 20])
            %        set(gca,'yticklabel',[3 2 1 0 1 2 3])
            %caxis([0 0.1])
            handlecolorbar=colorbar;%('horiz');
            colormap(flipud(bone));
            set(handlecolorbar,'fontsize',fs,'fontweight','bold');
            set(handlecolorbar,'Position',[0.85 0.3 0.03 0.43]);
            
            if savefigs
                saveas(fig,strcat('figs\probfcEgam1_',num2str(Egamma1/1e6/Qe),'Egam2_',num2str(Egamma2/Qe/1e6),'p',num2str(phi/pi*180),'Q',num2str(Q/Qe/1e6),'a.jpg'));
                saveas(fig,strcat('figs\probfcEgam1_',num2str(Egamma1/1e6/Qe),'Egam2_',num2str(Egamma2/1e6/Qe),'p',num2str(phi/pi*180),'Q',num2str(Q/Qe/1e6),'a.epsc'));
            end
               
            if 0
            
            fig=figure(1001); clf; hold on; box on;
            set(gcf,'DefaultLineLineWidth',4)
            [dummy,H]=contourf(vpara/1e6,vperp/1e6,wfcanalytic,40);
            for i = 1:length(H) set(H(i),'Edgecolor','none');end
            axis equal
            axis(ax)
            %caxis([-14 -7])
            set(gca,'fontsize',fs,'fontweight','bold')
            set(gca,'Position', [0.22 0.2 0.6 0.65])
            %title(strcat('a, \phi=',num2str(phi/pi*180),', E\gamma=',num2str(Egamma1/1e6/Qe)),'fontsize',fs,'fontweight','bold')
            xlabel('v_{||} [10^6 m/s]','fontsize',fs,'fontweight','bold')
            ylabel('v_{\perp} [10^6 m/s]','fontsize',fs,'fontweight','bold')
            set(gca,'ytick',[0 10 20])
            set(gca,'xtick',[-20 -10 0 10 20])
            %        set(gca,'yticklabel',[3 2 1 0 1 2 3])
            %caxis([0 1])
            handlecolorbar=colorbar;%('horiz');
            colormapInverseHotBWwfc;
            set(handlecolorbar,'fontsize',fs,'fontweight','bold');
            set(handlecolorbar,'Position',[0.85 0.3 0.03 0.43]);
            
            if savefigs
                saveas(fig,strcat('figs\wfcEgam1_',num2str(Egamma1/1e6/Qe),'Egam2_',num2str(Egamma2/Qe/1e6),'p',num2str(phi/pi*180),'Q',num2str(Q/Qe/1e6),'a.jpg'));
                saveas(fig,strcat('figs\wfcEgam1_',num2str(Egamma1/1e6/Qe),'Egam2_',num2str(Egamma2/1e6/Qe),'p',num2str(phi/pi*180),'Q',num2str(Q/Qe/1e6),'a.epsc'));
            end
            
            
            
            
            isocontours=[1e-11 1e-10 1e-9 1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1 10];
            
            
            fig=figure(2000);clf; hold on; box on;
            set(gcf,'DefaultLineLineWidth',4)
            contour(vpara/1e6,vperp/1e6,fvpavpe3DbiMax,isocontours);
            contour(vpara/1e6,vperp/1e6,fvpavpe3DbiMax,20);
            axis('equal')
            axis(ax)
            %text(1, 2.6,'f [s^3/m^6]','fontsize',16,'fontweight','bold')
            set(gca,'fontsize',16,'fontweight','bold')
            set(gca,'Position', [0.22 0.2 0.6 0.65])
            xlabel('v_{||} [10^6 m/s]','fontsize',16,'fontweight','bold')
            ylabel('v_{\perp} [10^6 m/s]','fontsize',16,'fontweight','bold')
            %   set(gca,'ytick',[-3 -2 -1 0 1 2 3])
            %   set(gca,'xtick',[-3 -2 -1 0 1 2 3])
            %   set(gca,'yticklabel',[3 2 1 0 1 2 3])
            handlecolorbar=colorbar;%('horiz');
            set(handlecolorbar,'fontsize',16,'fontweight','bold');
            set(handlecolorbar,'Position',[0.85 0.3 0.03 0.43]);
            
            if savefigs
                saveas(fig,strcat('figs\biMax3D.jpg'));
                saveas(fig,strcat('figs\biMax3D.epsc'));
            end
            
            fig=figure(2001);clf; hold on; box on;
            set(gcf,'DefaultLineLineWidth',4)
            contour(vpara/1e6,vperp/1e6,fvpavpe2DbiMax,isocontours*1e6);
            contour(vpara/1e6,vperp/1e6,fvpavpe2DbiMax,20);
            axis('equal')
            axis(ax)
            %text(1, 2.6,'f [s^2/m^5]','fontsize',16,'fontweight','bold')
            set(gca,'fontsize',16,'fontweight','bold')
            set(gca,'Position', [0.22 0.2 0.6 0.65])
            xlabel('v_{||} [10^6 m/s]','fontsize',16,'fontweight','bold')
            ylabel('v_{\perp} [10^6 m/s]','fontsize',16,'fontweight','bold')
            %   set(gca,'ytick',[-3 -2 -1 0 1 2 3])
            %   set(gca,'xtick',[-3 -2 -1 0 1 2 3])
            %   set(gca,'yticklabel',[3 2 1 0 1 2 3])
            handlecolorbar=colorbar;%('horiz');
            set(handlecolorbar,'fontsize',16,'fontweight','bold');
            set(handlecolorbar,'Position',[0.85 0.3 0.03 0.43]);
            if savefigs
                saveas(fig,strcat('figs\biMax2D.jpg'));
                saveas(fig,strcat('figs\biMax2D.epsc'));
            end
            
            end
            if 0
                r1=sqrt((mprod*Egamma1^2+2*mprod*(mprod-mf)*c^2*(Egamma1-Q))/(mf*(mprod-mf)^2*c^2));
                
                vcenterpara1=-cos(phi)*Egamma1/((mprod-mf)*c);
                vcenterperpplus1=sin(phi)*Egamma1/((mprod-mf)*c);
                vcenterperpminus1=-sin(phi)*Egamma1/((mprod-mf)*c);
                vcenter1=Egamma1/((mprod-mf)*c);
                
                rmax1=r1-vcenter1;
                rmin1=sqrt(r1^2-(vcenter1*sin(phi))^2)-vcenter1*abs(cos(phi));
                
                
                r2=sqrt((mprod*Egamma2^2+2*mprod*(mprod-mf)*c^2*(Egamma2-Q))/(mf*(mprod-mf)^2*c^2));
                vcenterpara2=-cos(phi)*Egamma2/((mprod-mf)*c);
                vcenterperpplus2=sin(phi)*Egamma2/((mprod-mf)*c);
                vcenterperpminus2=-sin(phi)*Egamma2/((mprod-mf)*c);
                
                vcenter2=Egamma2/((mprod-mf)*c);
                rmax2=r2+vcenter2;
                rmin2=sqrt(r2^2-(vcenter2*sin(phi))^2)-vcenter2*abs(cos(phi));
                
                vparacirc=[-15:0.01:15]*1e6;
                
                vmin1=sqrt(rmin1^2-vparacirc.^2);
                vmax1=sqrt(rmax1^2-vparacirc.^2);
                vmin2=sqrt(rmin2^2-vparacirc.^2);
                vmax2=sqrt(rmax2^2-vparacirc.^2);
                
                
                vcircleperppp1=vcenterperpminus1+sqrt(r1^2-(vparacirc-vcenterpara1).^2);
                vcircleperppm1=vcenterperpminus1-sqrt(r1^2-(vparacirc-vcenterpara1).^2);
                vcircleperpmp1=vcenterperpplus1+sqrt(r1^2-(vparacirc-vcenterpara1).^2);
                vcircleperpmm1=vcenterperpplus1-sqrt(r1^2-(vparacirc-vcenterpara1).^2);
                
                vcircleperppp2=vcenterperpminus2+sqrt(r2^2-(vparacirc-vcenterpara2).^2);
                vcircleperppm2=vcenterperpminus2-sqrt(r2^2-(vparacirc-vcenterpara2).^2);
                vcircleperpmp2=vcenterperpplus2+sqrt(r2^2-(vparacirc-vcenterpara2).^2);
                vcircleperpmm2=vcenterperpplus2-sqrt(r2^2-(vparacirc-vcenterpara2).^2);
                
                
                for i=1:length(vcircleperppp1)
                    if ~isreal(vcircleperppp1(i)) vcircleperppp1(i)=NaN; end
                    if ~isreal(vcircleperppm1(i)) vcircleperppm1(i)=NaN; end
                    if ~isreal(vcircleperpmp1(i)) vcircleperpmp1(i)=NaN; end
                    if ~isreal(vcircleperpmm1(i)) vcircleperpmm1(i)=NaN; end
                    if ~isreal(vcircleperppp2(i)) vcircleperppp2(i)=NaN; end
                    if ~isreal(vcircleperppm2(i)) vcircleperppm2(i)=NaN; end
                    if ~isreal(vcircleperpmp2(i)) vcircleperpmp2(i)=NaN; end
                    if ~isreal(vcircleperpmm2(i)) vcircleperpmm2(i)=NaN; end
                    if ~isreal(vmin1(i)) vmin1(i)=NaN; end
                    if ~isreal(vmax1(i)) vmax1(i)=NaN; end
                    if ~isreal(vmin2(i)) vmin2(i)=NaN; end
                    if ~isreal(vmax2(i)) vmax2(i)=NaN; end
                end
                
                fig=figure(1011); hold on;
                plot(vparacirc/1e6, vcircleperppp1/1e6,'--k')
                plot(vparacirc/1e6, vmax1/1e6,'--g')
                % plot(vpara/1e6, vmin1/1e6,'--m')
                plot(vparacirc/1e6, vmax2/1e6,'--g')
                %plot(vpara/1e6, vmin2/1e6,'--c')
                % plot(vpara/1e6, vcircleperppm1/1e6,'--g')
                % plot(vpara/1e6, vcircleperpmp1/1e6,'--r')
                % plot(vpara/1e6, vcircleperpmm1/1e6,'--y')
                % plot(vpara/1e6, vcircleperppp2/1e6,'k')
                % plot(vpara/1e6, vcircleperppm2/1e6,'g')
                plot(vparacirc/1e6, vcircleperpmp2/1e6,'--k')
                %plot(vpara/1e6, vcircleperpmm2/1e6,'y')
                
                if savefigs
                    saveas(fig,strcat('figs\wfcBoundEgam1_',num2str(Egamma1/1e6/Qe),'Egam2_',num2str(Egamma2/Qe/1e6),'p',num2str(phi/pi*180),'a.jpg'));
                    saveas(fig,strcat('figs\wfcBoundEgam1_',num2str(Egamma1/1e6/Qe),'Egam2_',num2str(Egamma2/1e6/Qe),'p',num2str(phi/pi*180),'a.epsc'));
                end
            end
            
            if 0
                sumprobfcanalytic=zeros(size(probfcanalytic));
                
                
                for iq=1:length(Qstar)
                    u1=(mf-mprod)*c/(2*Egamma1)*VMAG.^2+Egamma1/(2*mf*c)...
                        +mprod*c*(Egamma1-Qstar(iq)*1e6*Qe)/(mf*Egamma1);
                    u2=(mf-mprod)*c/(2*Egamma2)*VMAG.^2+Egamma2/(2*mf*c)...
                        +mprod*c*(Egamma2-Qstar(iq)*1e6*Qe)/(mf*Egamma2);
                    
                    gamma1=real(acos((u1-VPARA*cos(phi))./(VPERP*sin(phi))));
                    gamma2=real(acos((u2-VPARA*cos(phi))./(VPERP*sin(phi))));
                    probfcanalytic=abs(gamma1-gamma2)/pi;
                    sumprobfcanalytic=sumprobfcanalytic+Lorentz(iq)*(Qstar(2)-Qstar(1))*probfcanalytic;
                end
                
                
                
                fig=figure(4000); clf; hold on; box on;
                set(gcf,'DefaultLineLineWidth',4)
                [dummy,H]=contourf(vpara/1e6,vperp/1e6,log10(sumprobfcanalytic),40);
                for i = 1:length(H) set(H(i),'Edgecolor','none');end
                axis equal
                axis(ax)
                caxis([-5 -3])
                set(gca,'fontsize',fs,'fontweight','bold')
                set(gca,'Position', [0.22 0.2 0.6 0.65])
                %title(strcat('a, \phi=',num2str(phi/pi*180),', E\gamma=',num2str(Egamma1/1e6/Qe)),'fontsize',fs,'fontweight','bold')
                xlabel('v_{||} [10^6 m/s]','fontsize',fs,'fontweight','bold')
                ylabel('v_{\perp} [10^6 m/s]','fontsize',fs,'fontweight','bold')
                set(gca,'ytick',[0 10 20])
                set(gca,'xtick',[-20 -10 0 10 20])
                %        set(gca,'yticklabel',[3 2 1 0 1 2 3])
                %caxis([0 0.1])
                handlecolorbar=colorbar;%('horiz');
                colormapInverseHotBWwfc;
                set(handlecolorbar,'fontsize',fs,'fontweight','bold');
                set(handlecolorbar,'Position',[0.85 0.3 0.03 0.43]);
                if savefigs
                    saveas(fig,strcat('figs\probEgam1_',num2str(Egamma1/1e6/Qe),'Egam2_',num2str(Egamma2/Qe/1e6),'p',num2str(phi/pi*180),'Qwidth',num2str(Qwidth),'a.jpg'));
                    saveas(fig,strcat('figs\probEgam1_',num2str(Egamma1/1e6/Qe),'Egam2_',num2str(Egamma2/1e6/Qe),'p',num2str(phi/pi*180),'Qwidth',num2str(Qwidth),'a.epsc'));
                end
            end
        end
    end
    end

if 0
for ipara=1:length(vpara)
    for iperp=1:length(vperp)
        if sumwfcanalytic(iperp,ipara)~=0
            sumwfcanalytic(iperp,ipara)=log10(sumwfcanalytic(iperp,ipara));
        else
            sumwfcanalytic(iperp,ipara)=NaN;
        end
    end
end

for ipara=1:length(vpara)
    for iperp=1:length(vperp)
        if sumprobfcanalytic(iperp,ipara)~=0
            sumprobfcanalytic(iperp,ipara)=log10(sumprobfcanalytic(iperp,ipara));
        else
            sumprobfcanalytic(iperp,ipara)=NaN;
        end
    end
end

for ipara=1:length(vpara)
    for iperp=1:length(vperp)
        if sumorigin(iperp,ipara)~=0
            sumorigin(iperp,ipara)=log10(sumorigin(iperp,ipara));
        else
            sumorigin(iperp,ipara)=NaN;
        end
    end
end
fig=figure(1010); clf; hold on; box on;
set(gcf,'DefaultLineLineWidth',4)
[dummy,H]=contourf(vpara/1e6,vperp/1e6,sumprobfcanalytic,40);
for i = 1:length(H) set(H(i),'Edgecolor','none');end
axis equal
axis(ax)
%caxis([-4 0])
set(gca,'fontsize',fs,'fontweight','bold')
set(gca,'Position', [0.22 0.2 0.6 0.65])
%title(strcat('a, \phi=',num2str(phi/pi*180),', E\gamma=',num2str(Egamma1/1e6/Qe)),'fontsize',fs,'fontweight','bold')
xlabel('v_{||} [10^6 m/s]','fontsize',fs,'fontweight','bold')
ylabel('v_{\perp} [10^6 m/s]','fontsize',fs,'fontweight','bold')
set(gca,'ytick',[0 10 20])
set(gca,'xtick',[-20 -10 0 10 20])
%        set(gca,'yticklabel',[3 2 1 0 1 2 3])
caxis([-3 0])
handlecolorbar=colorbar;%('horiz');
colormapInverseHotBWwfc;
set(handlecolorbar,'fontsize',fs,'fontweight','bold');
set(handlecolorbar,'Position',[0.85 0.3 0.03 0.43]);

if savefigs
    saveas(fig,strcat('figs\sumprobfc_p',num2str(phi/pi*180),'Q',num2str(Q/Qe/1e6),'a.jpg'));
    saveas(fig,strcat('figs\sumprobfc_p',num2str(phi/pi*180),'Q',num2str(Q/Qe/1e6),'a.epsc'));
end

fig=figure(1011); clf; hold on; box on;
set(gcf,'DefaultLineLineWidth',4)
[dummy,H]=contourf(vpara/1e6,vperp/1e6,sumwfcanalytic,40);
for i = 1:length(H) set(H(i),'Edgecolor','none');end
axis equal
axis(ax)
%caxis([-11 -8])
set(gca,'fontsize',fs,'fontweight','bold')
set(gca,'Position', [0.22 0.2 0.6 0.65])
%title(strcat('a, \phi=',num2str(phi/pi*180),', E\gamma=',num2str(Egamma1/1e6/Qe)),'fontsize',fs,'fontweight','bold')
xlabel('v_{||} [10^6 m/s]','fontsize',fs,'fontweight','bold')
ylabel('v_{\perp} [10^6 m/s]','fontsize',fs,'fontweight','bold')
set(gca,'ytick',[0 10 20])
set(gca,'xtick',[-20 -10 0 10 20])
%        set(gca,'yticklabel',[3 2 1 0 1 2 3])
%caxis([0 1])
handlecolorbar=colorbar;%('horiz');
colormap(flipud(bone));
set(handlecolorbar,'fontsize',fs,'fontweight','bold');
set(handlecolorbar,'Position',[0.85 0.3 0.03 0.43]);
% set(handlecolorbar,'zticklabel',[-3 -2 -1 0])

if savefigs
    saveas(fig,strcat('figs\sumwfc_p',num2str(phi/pi*180),'Q',num2str(Q/Qe/1e6),'a.jpg'));
    saveas(fig,strcat('figs\sumwfc_p',num2str(phi/pi*180),'Q',num2str(Q/Qe/1e6),'a.epsc'));
end

%%%%%%FST2018
    fs=22;
    deltaenergy=10e3;
E_coarse=1:deltaenergy:1600e3;
deltapitch=0.01;
p_coarse=-0.99:deltapitch:0.99;
[energyMESH, pitchMESH] = meshgrid(E_coarse, p_coarse);
v_par_on_E_grid = -pitchMESH.*sqrt(2*energyMESH*Qe/mf);
v_perp_on_E_grid = sqrt(1-pitchMESH.^2).*sqrt(2*energyMESH*Qe/mf);
WF_temp=sumwfcanalyticnormalized;
         WFEp= interp2(vpara,vperp,WF_temp,v_par_on_E_grid,v_perp_on_E_grid,'linear',0);
        fig=figure(10+i); clf; hold on; box on;set(gca, 'Layer','top');
        set(gcf,'DefaultLineLineWidth',4)
                [dummy,H]=contourf(E_coarse/1e6,p_coarse, WFEp,40);
                for j = 1:length(H) set(H(j),'Edgecolor','none');end
    % figure; imagesc(vpara/1e6,vperp/1e6, WF_temp);
               colormap(flipud(bone))
    %    colorbar;
   
        %set(gca,'Fontsize',20)
        xlabel('E [MeV]')
        ylabel('p [-]')
               %caxis([-8 -2])
        set(gca,'fontsize',fs)
        
       saveas(gcf,['figs\oneGRSDpgammaHe3.png'])
       saveas(gcf,['figs\oneGRSDpgammaHe3.eps'],'epsc')     
%%%%%%%%%%

fig=figure(1012); clf; hold on; box on;
set(gcf,'DefaultLineLineWidth',4)
[dummy,H]=contourf(vpara/1e6,vperp/1e6,sumorigin,40);
for i = 1:length(H) set(H(i),'Edgecolor','none');end
axis equal
axis(ax)
caxis([-15 -5])
set(gca,'fontsize',fs,'fontweight','bold')
set(gca,'Position', [0.22 0.2 0.6 0.65])
%title(strcat('a, \phi=',num2str(phi/pi*180),', E\gamma=',num2str(Egamma1/1e6/Qe)),'fontsize',fs,'fontweight','bold')
xlabel('v_{||} [10^6 m/s]','fontsize',fs,'fontweight','bold')
ylabel('v_{\perp} [10^6 m/s]','fontsize',fs,'fontweight','bold')
set(gca,'ytick',[0 10 20])
set(gca,'xtick',[-20 -10 0 10 20])
%        set(gca,'yticklabel',[3 2 1 0 1 2 3])
%caxis([0 1])
handlecolorbar=colorbar;%('horiz');
colormapInverseHotBWwfc;
set(handlecolorbar,'fontsize',fs,'fontweight','bold');
set(handlecolorbar,'Position',[0.85 0.3 0.03 0.43]);
% set(handlecolorbar,'zticklabel',[-3 -2 -1 0])

[dummy,H]=contour(vpara/1e6,vperp/1e6,sumprobfcanalytic,[1e-6],'k');
%caxis([min(min(sumorigin)) max(max(sumorigin))])
%for i = 1:length(H) set(H(i),'Edgecolor','none');end
if savefigs
    saveas(fig,strcat('figs\origin_p',num2str(phi/pi*180),'Q',num2str(Q/Qe/1e6),'a.jpg'));
    saveas(fig,strcat('figs\origin_p',num2str(phi/pi*180),'Q',num2str(Q/Qe/1e6),'a.epsc'));
end

Egamma1scan=[5.5 6.5]*1e6*Qe;
Egamma2scan=Egamma1scan+0.0001*1e6*Qe;

r1scan=sqrt((mprod*Egamma1scan.^2+2*mprod*(mprod-mf)*c^2*(Egamma1scan-Q))/(mf*(mprod-mf)^2*c^2));
vcenter1scan=Egamma1scan/((mprod-mf)*c);
rmin1scan=r1scan-vcenter1scan;
r2scan=sqrt((mprod*Egamma2scan.^2+2*mprod*(mprod-mf)*c^2*(Egamma2scan-Q))/(mf*(mprod-mf)^2*c^2));
vcenter2scan=Egamma2scan/((mprod-mf)*c);
rmax2scan=r2scan+vcenter2scan;
Emin=0.5*mf*rmin1scan.^2/1e6/Qe;
Emax=0.5*mf*rmax2scan.^2/1e6/Qe;

Egammamid=(Egamma1scan+Egamma2scan)/2/1e6/Qe;

fig=figure(3000);clf;box on;hold on;
set(gcf,'DefaultLineLineWidth',4)
set(gca,'fontsize',28,'fontweight','bold')
xlabel('E_\gamma [MeV]','fontsize',28,'fontweight','bold')
ylabel('E_p [MeV]','fontsize',28,'fontweight','bold')
set(gca,'xlim',[5.5 6.5])
%set(gca,'ylim',[-0.1 1.8])
%set(gca,'xtick',[5.6 5.8 6 6.2 6.4])
%set(gca,'ytick',[0.2:0.2:1.8])
area(Egammamid,Emax,'facecolor','g')
area(Egammamid,Emin,'facecolor','w')
plot(Egammamid,Emax,'k')
plot(Egammamid,Emin,'k')
if savefigs
    saveas(fig,strcat('figs\energylimits.jpg'));
    saveas(fig,strcat('figs\energylimits.epsc'));
end
end