% This script is made for plotting analytical neutron spectrometry weight
% functions in both velocity space and energy-pitch space.

clear all

m_n = 1.67e-27;
m_D = 2*m_n;
m_He = 3*m_n;
Q = 3.27e6*1.6e-19;

phi_degrees = 10;
phi = phi_degrees*pi/180;

E_n = linspace(0, 10e6, 100);
E_n = E_n*1.6e-19;

dE = E_n(2)-E_n(1);

dv = 3.5/250;

v_parallel = (-2.5:dv:2.5)*1e7;
v_perpendicular = (dv:dv:5+dv)*1e7;

%v_parallel = (-5:dv:5)*1e7;
%v_perpendicular = (dv:dv:10+dv)*1e7;

[V_PARALLEL,V_PERPENDICULAR] = meshgrid(v_parallel,v_perpendicular);

dE_d = 1/5;

E_d = (dE_d:dE_d:100)*5e3*1.6e-19;

dp = 0.002;

pitch = -1+0.0001:dp:1-0.0001;

[E_D,PITCH] = meshgrid(E_d,pitch);

V_PAR_NY = PITCH.*sqrt(2.*E_D./m_D);
V_PERP_NY = sqrt(1-PITCH.^2).*sqrt(2.*E_D./m_D);

E_counter = 1;

[sizeA,sizeB] = size(V_PARALLEL);

wf_neutron_array = zeros(sizeA,sizeB,length(E_n));
wf_neutron_array_E_p = zeros(length(pitch),length(E_d),length(E_n));

for E = E_n;
    
    if E > dE/2
        
    gamma_1 = acos( (m_He + m_n)/m_D*sqrt(E - dE/2)./V_PERPENDICULAR*1/sqrt(2*m_n)*1/sin(phi)...
        - 1/2*(m_He-m_D)/sqrt(2*m_n)*1/sqrt(E - dE/2)*(V_PARALLEL.^2 + V_PERPENDICULAR.^2)./V_PERPENDICULAR*1/sin(phi)...
        - m_He/m_D*Q./V_PERPENDICULAR*1/sqrt(2*m_n)*1/sqrt(E - dE/2)*1/sin(phi)...
        - V_PARALLEL./V_PERPENDICULAR*cot(phi) );
    
    gamma_2 = acos( (m_He + m_n)/m_D*sqrt(E + dE/2)./V_PERPENDICULAR*1/sqrt(2*m_n)*1/sin(phi)...
        - 1/2*(m_He-m_D)/sqrt(2*m_n)*1/sqrt(E + dE/2)*(V_PARALLEL.^2 + V_PERPENDICULAR.^2)./V_PERPENDICULAR*1/sin(phi)...
        - m_He/m_D*Q./V_PERPENDICULAR*1/sqrt(2*m_n)*1/sqrt(E + dE/2)*1/sin(phi)...
        - V_PARALLEL./V_PERPENDICULAR*cot(phi) );    
    
    gamma_1_E_p = acos(1./2.*(m_He + m_n)./sqrt(m_D.*m_n).*sqrt((E-dE/2)./E_D).*1./sqrt(1-PITCH.^2).*1./sin(phi)...
        - 1/2.*(m_He - m_D)./sqrt(m_D.*m_n).*sqrt(E_D./(E-dE/2))*1./sqrt(1-PITCH.^2).*1./sin(phi)...
        - m_He./sqrt(m_D.*m_n).*Q./(2.*sqrt(E_D.*(E-dE./2))).*1./sqrt(1-PITCH.^2)*1/sin(phi) - PITCH./sqrt(1-PITCH.^2).*cot(phi) );
    
    gamma_2_E_p = acos(1./2.*(m_He + m_n)./sqrt(m_D.*m_n).*sqrt((E+dE/2)./E_D).*1./sqrt(1-PITCH.^2).*1./sin(phi)...
        - 1/2.*(m_He - m_D)./sqrt(m_D.*m_n).*sqrt(E_D./(E+dE/2))*1./sqrt(1-PITCH.^2).*1./sin(phi)...
        - m_He./sqrt(m_D.*m_n).*Q./(2.*sqrt(E_D.*(E+dE./2))).*1./sqrt(1-PITCH.^2)*1/sin(phi) - PITCH./sqrt(1-PITCH.^2).*cot(phi) );
    
    weight_function_E_p = real((gamma_1_E_p - gamma_2_E_p)/pi);
    %transfer_matrix_neutrons_E_p(:,E_counter) = reshape(weight_function_E_p,length(E_d)*length(pitch),1);
    wf_neutron_array_E_p(:,:,E_counter) = weight_function_E_p;
    
    %gamma_1 = acos(sqrt(E-1/2*dE)*(m_alpha+m_n)./(sqrt(2*m_n)*m_D*V_PERPENDICULAR)...
    %    -1/sqrt(E-1/2*dE)*m_alpha/m_D*Q./V_PERPENDICULAR./sqrt(2*m_n)...
    %    - 1/sqrt(E-1/2*dE)*(m_alpha-m_D)*(V_PARALLEL.^2 + V_PERPENDICULAR.^2)./(2*sqrt(2*m_n)*V_PERPENDICULAR));
    %gamma_2 = acos(sqrt(E+1/2*dE)*(m_alpha+m_n)./(sqrt(2*m_n)*m_D*V_PERPENDICULAR)...
    %    -1/sqrt(E+1/2*dE)*m_alpha/m_D*Q./V_PERPENDICULAR./sqrt(2*m_n)...
    %    - 1/sqrt(E+1/2*dE)*(m_alpha-m_D)*(V_PARALLEL.^2 + V_PERPENDICULAR.^2)./(2*sqrt(2*m_n)*V_PERPENDICULAR));
    
    %gamma_1 = acos(-1./(2.*m_D.*(v_n-dv/2).*V_PERPENDICULAR).*(m_D./m_n.*(m_alpha-m_D).*(V_PARALLEL.^2 + V_PERPENDICULAR.^2)...
    %    + 2.*m_alpha./m_n.*Q - (v_n-dv/2).^2).*(m_alpha+m_n));
    %gamma_2 = acos(-1./(2.*m_D.*(v_n+dv/2).*V_PERPENDICULAR).*(m_D./m_n.*(m_alpha-m_D).*(V_PARALLEL.^2 + V_PERPENDICULAR.^2)...
    %    + 2.*m_alpha./m_n.*Q - (v_n+dv/2).^2).*(m_alpha+m_n));
    
    weight_test = (gamma_1 - gamma_2)/pi;
    weight_function = real((gamma_1 - gamma_2)/pi);
    %transfer_matrix_neutrons(:,E_counter) = reshape(weight_function,sizeA*sizeB,1);
    wf_neutron_array(:,:,E_counter) = weight_function;
    end
    E_counter = E_counter + 1;
end

load('MyColormap.mat')

if 0

frame_counter = 1;

pause on

figure

for j = 17:100
    
    subplot(1,3,1)
    contourf(v_parallel,v_perpendicular,wf_neutron_array(:,:,j),50)
    set(gcf,'Colormap',mycolormap)
    axis equal
    axis tight
    colorbar
    caxis([0 1])
    
    subplot(1,3,2)
    contourf(E_d/1.6e-19,pitch,wf_neutron_array_E_p(:,:,j),50)
    set(gcf,'Colormap',mycolormap)
    axis square
    colorbar
    caxis([0 1])
    
    wf_neutron_E_p_test = interp2(v_parallel,v_perpendicular,wf_neutron_array(:,:,j),V_PAR_NY,V_PERP_NY);
    
    subplot(1,3,3)
    contourf(E_d/1.6e-19,pitch,wf_neutron_E_p_test,50)
    set(gcf,'Colormap',mycolormap)
    axis square
    colorbar
    caxis([0 1])
    
    suptitle(['E_n = ' num2str(E_n(j)/1.6e-19/1e6) ' MeV'])
    
    Frame(frame_counter) = getframe(gcf);
    
    frame_counter = frame_counter + 1;
    
    pause(0.2)
        
end
end


if 1
    
    v_n = sqrt(2*E_n/m_n);
    
    j = [21 25 31];
    
    %j = [17 21 25 31 36];
    %j = [16 17 18 19 20];
    
    %E_n_plots = char('1_6', '2', '2_4', '3', '3_5');
    E_n_plots = char('2', '2_4', '3');
    E_n_plot_counter = 1;
        
    for l = j
        
        figure
        %hold on
        %log_wf = log10(wf_neutron_array(:,:,l));
        %log_wf(isinf(log_wf)) = -3;
        lscale = logspace(-3, 0);
        %contourf(v_parallel/1e6,v_perpendicular/1e6,log_wf,log10(lscale))
        %contourf(v_parallel/1e6,v_perpendicular/1e6,wf_neutron_array(:,:,l),lscale)
        contourf(v_parallel/1e6,v_perpendicular/1e6,wf_neutron_array(:,:,l),50)
        axis equal
        axis tight
        colorbar
        caxis([0 1])
        %caxis([1e-3 1])
        set(gcf,'Colormap',mycolormap);
        set(gca,'FontSize',28)
        xlabel('v_{||} [10^6 m/s]')
        ylabel('v_\perp [10^6 m/s]')
        
        if 0
        v_par_0 = -cos(phi)/(m_He - m_D)*m_n*v_n(l);
        v_perp_0_1 = sin(phi)/(m_He - m_D)*m_n*v_n(l);
        
        v_perp_0_2 = -sin(phi)/(m_He - m_D)*m_n*v_n(l);
        
        r_1 = sqrt(3/m_n*(2*E_n(l) + dE/2 - Q));
        r_2 = sqrt(3/m_n*(2*E_n(l) - dE/2 - Q));
        
        circ_eq_1 = (V_PARALLEL - v_par_0).^2 + (V_PERPENDICULAR - v_perp_0_1).^2 - r_1^2;
        circ_eq_2 = (V_PARALLEL - v_par_0).^2 + (V_PERPENDICULAR - v_perp_0_2).^2 - r_2^2;
                
        wf_edge_1 = zeros(size(V_PARALLEL));
        wf_edge_2 = zeros(size(V_PARALLEL));
                
        wf_edge_1(abs(circ_eq_1)< 1e13) = 1;
        wf_edge_2(abs(circ_eq_2)< 1e13) = 1;
                
        contour(v_parallel/1e6,v_perpendicular/1e6,wf_edge_1,50,'LineWidth',5);
        contour(v_parallel/1e6,v_perpendicular/1e6,wf_edge_2,50,'LineWidth',5);
            
        end
        %hold off
        
        %saveas(gcf,strcat('results\Article_2\Extra_wf_plots\wf_v_space_phi_', num2str(phi_degrees), '_En_', E_n_plots(E_n_plot_counter,:), '_MeV.png'))
        %saveas(gcf,strcat('results\Article_2\Extra_wf_plots\wf_v_space_phi_', num2str(phi_degrees), '_En_', E_n_plots(E_n_plot_counter,:), '_MeV.epsc'))
        
        %saveas(gcf,strcat('results\Article_2\wf_plots\wf_v_space_phi_', num2str(phi_degrees), '_En_', E_n_plots(E_n_plot_counter,:), '_MeV.png'))
        %saveas(gcf,strcat('results\Article_2\wf_plots\wf_v_space_phi_', num2str(phi_degrees), '_En_', E_n_plots(E_n_plot_counter,:), '_MeV.epsc'))
        
        
        figure
        contourf(E_d/1.6e-19/1e3,pitch,wf_neutron_array_E_p(:,:,l),20)
        axis square
        axis tight
        ylim([-1 1])
        xlim([0 500])
        set(gca,'Xtick',[0,250,500])
        caxis([0 1])
        %caxis([1e-3 1])
        set(gcf,'Colormap',mycolormap);
        set(gca,'FontSize',28)
        xlabel('E [keV]')
        ylabel('pitch [-]')
        
        
        %saveas(gcf,strcat('results\Article_2\Extra_wf_plots\wf_E_p_phi_', num2str(phi_degrees), '_En_', E_n_plots(E_n_plot_counter,:), '_MeV.png'))
        %saveas(gcf,strcat('results\Article_2\Extra_wf_plots\wf_E_p_phi_', num2str(phi_degrees), '_En_', E_n_plots(E_n_plot_counter,:), '_MeV.epsc'))
        
        %saveas(gcf,strcat('results\Article_2\wf_plots\wf_E_p_phi_', num2str(phi_degrees), '_En_', E_n_plots(E_n_plot_counter,:), '_MeV_zoom.png'))
        %saveas(gcf,strcat('results\Article_2\wf_plots\wf_E_p_phi_', num2str(phi_degrees), '_En_', E_n_plots(E_n_plot_counter,:), '_MeV_zoom.epsc'))
        
        
        E_n_plot_counter = E_n_plot_counter + 1;
        
    end
end
        
        
        
        
        
        
    
    