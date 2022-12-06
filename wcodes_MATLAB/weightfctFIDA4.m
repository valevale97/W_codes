clear;
fs=28;
%natconst;
Qe=1.602e-19;
Mp=1.673e-27;

lambdabalmer=656.1e-9;
lambda=[-257.56,-193.17,-128.78,-64.39,0.000001,64.39,128.78,193.17,257.56];
stark_intens=[1.d0,   18.d0,   16.d0, 1681.d0, 2304.d0, 729.d0, 1936.d0, 5490.d0, 1936.d0, 729.d0...
        2304.d0, 1681.d0, 16.d0, 18.d0, 1.d0];
stark_pi=[1,0,0,1,1,1,0,0,0,1,1,1,0,0,1];
stark_sigma=1-stark_pi;
stark_wavel =[-2.20200d-06,-1.65200d-06,-1.37700d-06,-1.10200d-06, -8.26400d-07,-5.51000d-07,-2.75600d-07,...
        0.00000d0, 2.75700d-07, 5.51500d-07, 8.27400d-07, 1.10300d-06, 1.38000d-06, 1.65600d-06, 2.20900d-06]*1e-10;


    
    
lambdaexp=659.1e-9;

%for lambdaexp=[662.1 659.1 653.1 650.1]*1e-9;
%for phi=[20 40 60 80];
phi=40;
gamma=0:360;

stark_intens_prob=stark_intens/sum(stark_intens);

for j=1:length(stark_intens)
    for i=1:length(gamma);
        stark_intens_prob_gyro(j,i)=stark_intens_prob(j)*(1-stark_pi(j)*(sin(phi/180*pi)*sin(gamma(i)/180*pi))^2 ...
            +stark_sigma(j)*(sin(phi/180*pi)*sin(gamma(i)/180*pi))^2);
    end
end

stark_intens_prob_gyro_sum=sum(stark_intens_prob_gyro,1);


figure(433);clf;hold on;
for i=1:length(stark_intens)
plot(gamma, log10(stark_intens_prob_gyro(i,:)))
end
plot(gamma, log10(stark_intens_prob_gyro_sum),'r')


savefigs=0;

%construct 2D distribution function as a function of parallel and perpendicular velocity. 
dv=3.5e6/200;

%vperp = [0.01e6:dv:3.5e6];
vperpt = [dv:dv:3.5e6];
vparat  = [-3.5e06:dv:3.5e6];
[VPARAT,VPERPT]=meshgrid(vparat,vperpt);



[rows,columns]=size(VPERPT);
ax = [-3.6 3.6 -0.1 3.6];
dlambda=.1e-9;
%for dlambda=[0.1 0.2 0.5 1]*1e-9

du=3e8*(dlambda/lambdabalmer);
%du=1e4;

% %Stark splitting
Efield=1e5;

Bfield=1.79;
vcrossBfield=VPERPT*Bfield;


for i=1:length(stark_wavel)  
    %for i=6:9
    %    lambdastark(:,:,i)=lambdabalmer*(1 + 1.e-6*vcrossBfield*lambdabalmer*lambda(i));
    lambdastark(:,:,i)=lambdabalmer + vcrossBfield*stark_wavel(i);
end

ustark=3e8*(lambdaexp./lambdastark-1);


intennorm=stark_intens/sum(stark_intens);

wvstark=zeros(size(VPARAT));
for i=1:15
    gamma1=acos((ustark(:,:,i)-du/2-cos(phi/180*pi).*(-VPARAT))./(sin(phi/180*pi).*VPERPT));
    gamma2=acos((ustark(:,:,i)+du/2-cos(phi/180*pi).*(-VPARAT))./(sin(phi/180*pi).*VPERPT));
    wv=real( (gamma1-gamma2)/pi - stark_pi(i)*(sin(phi/180*pi))^2/2*((gamma1-gamma2)/pi - (sin(2*gamma1)-sin(2*gamma2))/(2*pi))...
        + stark_sigma(i)*(sin(phi/180*pi))^2/2*((gamma1-gamma2)/pi - (sin(2*gamma1)-sin(2*gamma2))/(2*pi)));
    wvstark=wvstark+stark_intens_prob(i)*wv;
end

gamma1=acos((ustark(:,:,8)-du/2-cos(phi/180*pi).*(-VPARAT))./(sin(phi/180*pi).*VPERPT));
gamma2=acos((ustark(:,:,8)+du/2-cos(phi/180*pi).*(-VPARAT))./(sin(phi/180*pi).*VPERPT));
wv=real((gamma1-gamma2))/pi;%/du*dv^2;

energy0=5.53e6;
phaseshift=0;
vperp0=sqrt(energy0*2*Qe/(2*Mp));
amplitudeV=VPERPT/vperp0;
wvcosine=wv+2*amplitudeV.^2*cos(phaseshift/180*pi).*(real(sin(gamma1)-sin(gamma2)));
       
fig=figure(200);clf; hold on; box on;
set(gcf,'DefaultLineLineWidth',4)
[dummy,H]=contourf(vparat/1e6,vperpt/1e6,real(log10(wvstark)),200);
for i = 1:length(H) set(H(i),'Edgecolor','none');end
axis('equal')
axis(ax)
set(gca,'fontsize',fs,'fontweight','bold')
xlabel('v_{||} [10^6 m/s]','fontsize',fs,'fontweight','bold')
ylabel('v_{\perp} [10^6 m/s]','fontsize',fs,'fontweight','bold')
set(gca,'ytick',[-3 -2 -1 0 1 2 3])
set(gca,'xtick',[-3 -2 -1 0 1 2 3])
set(gca,'yticklabel',[3 2 1 0 1 2 3])
caxis([-3 0])
handlecolorbar=colorbar;%('horiz');
colormapInverseHotBWwfc;
set(handlecolorbar,'fontsize',fs,'fontweight','bold'); 
set(handlecolorbar,'Position',[0.85 0.3 0.03 0.43]);
%set(handlecolorbar,'zticklabel',[-3 -2 -1 0])
 
 

fig=figure(210);clf; hold on; box on;
set(gcf,'DefaultLineLineWidth',4)
[dummy,H]=contourf(vparat/1e6,vperpt/1e6,log10(wv),200);
for i = 1:length(H) set(H(i),'Edgecolor','none');end
axis('equal')
axis(ax)
%set(gca,'Position', [0.2 0.2 0.6 0.75])
set(gca,'fontsize',fs,'fontweight','bold')
xlabel('v_{||} [10^6 m/s]','fontsize',fs,'fontweight','bold')
ylabel('v_{\perp} [10^6 m/s]','fontsize',fs,'fontweight','bold')
set(gca,'ytick',[-3 -2 -1 0 1 2 3])
set(gca,'xtick',[-3 -2 -1 0 1 2 3])
set(gca,'yticklabel',[3 2 1 0 1 2 3])
caxis([-3 0])
handlecolorbar=colorbar;%('horiz');
colormapInverseHotBWwfc;
set(handlecolorbar,'fontsize',fs,'fontweight','bold'); 
set(handlecolorbar,'Position',[0.85 0.3 0.03 0.43]);
%set(handlecolorbar,'zticklabel',[-3 -2 -1 0])

    if 1
fig=figure(220);clf; hold on; box on;
set(gcf,'DefaultLineLineWidth',4)
[dummy,H]=contourf(vparat/1e6,vperpt/1e6,log10(wvcosine),200);
for i = 1:length(H) set(H(i),'Edgecolor','none');end
axis('equal')
axis(ax)
%set(gca,'Position', [0.2 0.2 0.6 0.75])
set(gca,'fontsize',fs,'fontweight','bold')
xlabel('v_{||} [10^6 m/s]','fontsize',fs,'fontweight','bold')
ylabel('v_{\perp} [10^6 m/s]','fontsize',fs,'fontweight','bold')
set(gca,'ytick',[-3 -2 -1 0 1 2 3])
set(gca,'xtick',[-3 -2 -1 0 1 2 3])
set(gca,'yticklabel',[3 2 1 0 1 2 3])
caxis([-3 0])
handlecolorbar=colorbar;%('horiz');
colormapInverseHotBWwfc;
set(handlecolorbar,'fontsize',fs,'fontweight','bold'); 
set(handlecolorbar,'Position',[0.85 0.3 0.03 0.43]);
%set(handlecolorbar,'zticklabel',[-3 -2 -1 0])



fig=figure(230);clf; hold on; box on;
set(gcf,'DefaultLineLineWidth',4)
[dummy,H]=contourf(vparat/1e6,vperpt/1e6,amplitudeV,50);
for i = 1:length(H) set(H(i),'Edgecolor','none');end
axis('equal')
axis(ax)
%set(gca,'Position', [0.2 0.2 0.6 0.75])
set(gca,'fontsize',fs,'fontweight','bold')
xlabel('v_{||} [10^6 m/s]','fontsize',fs,'fontweight','bold')
ylabel('v_{\perp} [10^6 m/s]','fontsize',fs,'fontweight','bold')
set(gca,'ytick',[-3 -2 -1 0 1 2 3])
set(gca,'xtick',[-3 -2 -1 0 1 2 3])
set(gca,'yticklabel',[3 2 1 0 1 2 3])
%caxis([-3 0])
handlecolorbar=colorbar;%('horiz');
colormapInverseHotBWwfc;
set(handlecolorbar,'fontsize',fs,'fontweight','bold'); 
set(handlecolorbar,'Position',[0.85 0.3 0.03 0.43]);
%set(handlecolorbar,'zticklabel',[-3 -2 -1 0])
if savefigs
    saveas(fig,strcat('C:\CTS\MSAL\investigations\weightfct\FIDA\figs\wfcFIDAnoStarkp',num2str(phi),'l',num2str(lambdaexp*1e9),'dl',num2str(dlambda*1e9),'.epsc'));
    saveas(fig,strcat('C:\CTS\MSAL\investigations\weightfct\FIDA\figs\wfcFIDAnoStarkp',num2str(phi),'l',num2str(lambdaexp*1e9),'dl',num2str(dlambda*1e9),'.jpg'));
end  
end


if 1
% energy = 2e3:2e2:140e3;
% pitch = -0.998:0.004:0.998;
energy = 2e3:5e2:140e3;
pitch = -0.998:0.01:0.998;


m = 2*1.67e-27;

[ENERGY,PITCH] = meshgrid(energy,pitch);

VPERPEP=sqrt((1-PITCH.^2)*2.*ENERGY*Qe/(2*Mp));

vcrossBfieldEP=VPERPEP*Bfield;


for i=1:15  
    %for i=6:9
    lambdastarkEP(:,:,i)=lambdabalmer+stark_wavel(i)*vcrossBfieldEP;
end

ustarkEP=3e8*(lambdaexp./lambdastarkEP-1);

% for jp= 1: length(Energy)
%     vpartt(jp,:) =  -Pitch'.*sqrt(2*Energy(jp)*Qe/mD);
%     vperptt(jp,:) = sqrt((1-Pitch'.^2)*2*Energy(jp)*Qe/mD);
% end

wEPstark=zeros(size(ENERGY));
for i=1:15
    gamma1 = acos(1./sin(phi.*pi./180).*1./sqrt(1-PITCH.^2).*((ustarkEP(:,:,i)-du./2)./sqrt(2.*ENERGY./m*Qe)-PITCH.*cos(phi.*pi./180)));
    gamma2 = acos(1./sin(phi.*pi./180).*1./sqrt(1-PITCH.^2).*((ustarkEP(:,:,i)+du./2)./sqrt(2.*ENERGY./m*Qe)-PITCH.*cos(phi.*pi./180)));       
    wEP=real( (gamma1-gamma2)/pi - stark_pi(i)*(sin(phi/180*pi))^2/2*((gamma1-gamma2)/pi - (sin(2*gamma1)-sin(2*gamma2))/(2*pi))...
        + stark_sigma(i)*(sin(phi/180*pi))^2/2*((gamma1-gamma2)/pi - (sin(2*gamma1)-sin(2*gamma2))/(2*pi)));
    wEPstark=wEPstark+stark_intens_prob(i)*wEP;
end

gamma1 = acos(1./sin(phi.*pi./180).*1./sqrt(1-PITCH.^2).*((ustarkEP(:,:,5)-du./2)./sqrt(2.*ENERGY./m*Qe)-PITCH.*cos(phi.*pi./180)));
gamma2 = acos(1./sin(phi.*pi./180).*1./sqrt(1-PITCH.^2).*((ustarkEP(:,:,5)+du./2)./sqrt(2.*ENERGY./m*Qe)-PITCH.*cos(phi.*pi./180)));       
wEP=real((gamma1-gamma2))/pi;%/du*dv^2;

amplitudeEP=sqrt(ENERGY./energy0.*(1-PITCH.^2));
wEPcosine=wEP+2*amplitudeEP*cos(phaseshift/180*pi).*(real(sin(gamma1)-sin(gamma2)));




axEP = [0 142 -1.02 1.02];

if 1 
    fig=figure(300);clf; hold on; box on;
    set(gcf,'DefaultLineLineWidth',4)
    [dummy,H]=contourf(energy/1000,pitch,wEPstark,50);
    for i = 1:length(H) set(H(i),'Edgecolor','none');end
    %axis('equal')
    axis(axEP)
    set(gca,'fontsize',fs,'fontweight','bold')
    xlabel('Energy [keV]','fontsize',fs,'fontweight','bold')
    ylabel('Pitch [-]','fontsize',fs,'fontweight','bold')
    %set(gca,'ytick',[-3 -2 -1 0 1 2 3])
    %set(gca,'xtick',[-3 -2 -1 0 1 2 3])
    %set(gca,'yticklabel',[3 2 1 0 1 2 3])
    handlecolorbar=colorbar;%('horiz');
    colormapInverseHotBWwfc;
    set(handlecolorbar,'fontsize',fs,'fontweight','bold'); 
    set(handlecolorbar,'Position',[0.85 0.3 0.03 0.4]);
    set(gca,'Position', [0.150 0.150 0.65 0.80])
    if savefigs
        saveas(fig,strcat('C:\CTS\MSAL\investigations\weightfct\FIDA\figs\wfcEPFIDAStarkp',num2str(phi),'l',num2str(lambdaexp*1e9),'dl',num2str(dlambda*1e9),'.epsc'));
        saveas(fig,strcat('C:\CTS\MSAL\investigations\weightfct\FIDA\figs\wfcEPFIDAStarkp',num2str(phi),'l',num2str(lambdaexp*1e9),'dl',num2str(dlambda*1e9),'.jpg'));
        
    end  
end

fig=figure(310);clf; hold on; box on;
set(gcf,'DefaultLineLineWidth',4)
%[dummy,H]=contourf(energy/1000,pitch,log10(wEP),100);
[dummy,H]=contourf(energy/1000,pitch,wEP,100);
for i = 1:length(H) set(H(i),'Edgecolor','none');end
axis(axEP)
set(gca,'fontsize',fs,'fontweight','bold')
xlabel('Energy [keV]','fontsize',fs,'fontweight','bold')
ylabel('Pitch [-]','fontsize',fs,'fontweight','bold')
%set(gca,'ytick',[-3 -2 -1 0 1 2 3])
%set(gca,'xtick',[-3 -2 -1 0 1 2 3])
%set(gca,'yticklabel',[3 2 1 0 1 2 3])
%caxis([-3 0])
handlecolorbar=colorbar;%('horiz');
%colormapInverseHotBWwfc;
colormap(flipud(pink))
set(handlecolorbar,'fontsize',fs,'fontweight','bold'); 
set(handlecolorbar,'Position',[0.85 0.2 0.03 0.7]);
set(gca,'Position', [0.180 0.150 0.62 0.80])
if savefigs
    saveas(fig,strcat('C:\CTS\MSAL\investigations\weightfct\FIDA\figs\wfcEPFIDAnoStarkp',num2str(phi),'l',num2str(lambdaexp*1e9),'dl',num2str(dlambda*1e9),'.epsc'));
    saveas(fig,strcat('C:\CTS\MSAL\investigations\weightfct\FIDA\figs\wfcEPFIDAnoStarkp',num2str(phi),'l',num2str(lambdaexp*1e9),'dl',num2str(dlambda*1e9),'.jpg'));
    
end  


if 0
fig=figure(320);clf; hold on; box on;
set(gcf,'DefaultLineLineWidth',4)
[dummy,H]=contourf(energy/1000,pitch,log10(wEPcosine),100);
for i = 1:length(H) set(H(i),'Edgecolor','none');end
axis(axEP)
set(gca,'fontsize',fs,'fontweight','bold')
xlabel('Energy [keV]','fontsize',fs,'fontweight','bold')
ylabel('Pitch [-]','fontsize',fs,'fontweight','bold')
%set(gca,'ytick',[-3 -2 -1 0 1 2 3])
%set(gca,'xtick',[-3 -2 -1 0 1 2 3])
%set(gca,'yticklabel',[3 2 1 0 1 2 3])
caxis([-3 0])
handlecolorbar=colorbar;%('horiz');
colormapInverseHotBWwfc;
set(handlecolorbar,'fontsize',fs,'fontweight','bold'); 
set(handlecolorbar,'Position',[0.85 0.2 0.03 0.7]);
set(gca,'Position', [0.180 0.150 0.62 0.80])
if savefigs
    saveas(fig,strcat('C:\CTS\MSAL\investigations\weightfct\FIDA\figs\wfcEPFIDAcosine',num2str(phi),'l',num2str(lambdaexp*1e9),'dl',num2str(dlambda*1e9),'.epsc'));
    saveas(fig,strcat('C:\CTS\MSAL\investigations\weightfct\FIDA\figs\wfcEPFIDAcosine',num2str(phi),'l',num2str(lambdaexp*1e9),'dl',num2str(dlambda*1e9),'.jpg'));
    
end  

fig=figure(330);clf; hold on; box on;
set(gcf,'DefaultLineLineWidth',4)
[dummy,H]=contourf(energy/1000,pitch,amplitudeEP,50);
for i = 1:length(H) set(H(i),'Edgecolor','none');end
axis(axEP)
set(gca,'fontsize',fs,'fontweight','bold')
xlabel('Energy [keV]','fontsize',fs,'fontweight','bold')
ylabel('Pitch [-]','fontsize',fs,'fontweight','bold')
%set(gca,'ytick',[-3 -2 -1 0 1 2 3])
%set(gca,'xtick',[-3 -2 -1 0 1 2 3])
%set(gca,'yticklabel',[3 2 1 0 1 2 3])
%caxis([-3 0])
handlecolorbar=colorbar;%('horiz');
colormapInverseHotBWwfc;
set(handlecolorbar,'fontsize',fs,'fontweight','bold'); 
set(handlecolorbar,'Position',[0.85 0.2 0.03 0.7]);
set(gca,'Position', [0.180 0.150 0.62 0.80])
if savefigs
    saveas(fig,strcat('C:\CTS\MSAL\investigations\weightfct\FIDA\figs\wfcEPFIDamplitude',num2str(phi),'l',num2str(lambdaexp*1e9),'dl',num2str(dlambda*1e9),'.epsc'));
    saveas(fig,strcat('C:\CTS\MSAL\investigations\weightfct\FIDA\figs\wfcEPFIDamplitude',num2str(phi),'l',num2str(lambdaexp*1e9),'dl',num2str(dlambda*1e9),'.jpg'));
    
end  

end
end
if 0
    for i=1:9
        figure(700+i);clf;hold on;
        [dummy,H]=contourf(vparat/1e6,vperpt/1e6,lambdastark(:,:,i)*1e9,20); 
        %for i = 1:length(H) set(H(i),'Edgecolor','none');end
        %[dummy,H]=contour(vpara/1e6,vperp/1e6,fVvpavpe2pivperp,10,'k'); 
        %for i = 1:length(H) set(H(i),'linewidth',1,'color',[0.8 0.8 0.8]);end
        axis('equal')
        axis(ax)
        set(gca,'fontsize',fs,'fontweight','bold')
        xlabel('v_{||} [10^6 m/s]','fontsize',fs,'fontweight','bold')
        ylabel('v_{\perp} [10^6 m/s]','fontsize',fs,'fontweight','bold')
        %   set(gca,'ytick',[-3 -2 -1 0 1 2 3])
        %   set(gca,'xtick',[-3 -2 -1 0 1 2 3])
        %   set(gca,'yticklabel',[3 2 1 0 1 2 3])
        handlecolorbar=colorbar;%('horiz');
        %caxis([0 max(max(fVvpavpe2pivperp))])
        %colormap(jet);
        set(handlecolorbar,'fontsize',fs,'fontweight','bold'); 
        %set(handlecolorbar,'Position',[0.85 0.3 0.03 0.43]);
        
        figure(800+i);clf;hold on;
        [dummy,H]=contourf(vparat/1e6,vperpt/1e6,ustark(:,:,i),20); 
        %for i = 1:length(H) set(H(i),'Edgecolor','none');end
        %[dummy,H]=contour(vpara/1e6,vperp/1e6,fVvpavpe2pivperp,10,'k'); 
        %for i = 1:length(H) set(H(i),'linewidth',1,'color',[0.8 0.8 0.8]);end
        axis('equal')
        axis(ax)
        set(gca,'fontsize',fs,'fontweight','bold')
        xlabel('v_{||} [10^6 m/s]','fontsize',fs,'fontweight','bold')
        ylabel('v_{\perp} [10^6 m/s]','fontsize',fs,'fontweight','bold')
        %   set(gca,'ytick',[-3 -2 -1 0 1 2 3])
        %   set(gca,'xtick',[-3 -2 -1 0 1 2 3])
        %   set(gca,'yticklabel',[3 2 1 0 1 2 3])
        handlecolorbar=colorbar;%('horiz');
        %caxis([0 max(max(fVvpavpe2pivperp))])
        %colormap(jet);
        set(handlecolorbar,'fontsize',fs,'fontweight','bold'); 
        %set(handlecolorbar,'Position',[0.85 0.3 0.03 0.43]);
    end
    
end

%end
