% This script calculates the neutron rate part of the weight functions
% analytically and compares them with numerically calculated ones

clear all

m_n = 1.67e-27;
m_D = 2*m_n;

energy = (0:100:10000)*1e3*1.6e-19;
pitch_analytic = -1:0.02:1;

[ENERGY,PITCH] = meshgrid(energy,pitch_analytic);

v_par_on_Ep_grid = PITCH.*sqrt(2*ENERGY/m_D);
v_perp_on_Ep_grid = sqrt(1-PITCH.^2).*sqrt(2*ENERGY/m_D);

Thermal_drift = 1; % 1 = yes, 0 = no

switch Thermal_drift
    
    case 0
        
        v_rel = sqrt(v_par_on_Ep_grid.^2 + v_perp_on_Ep_grid.^2);
        E_cm = 1/2 * m_D/2 * v_rel.^2;
        
    case 1
        
        v_par_thermal_drift = 2.1e5;
        v_rel = sqrt((v_par_on_Ep_grid - v_par_thermal_drift).^2 + v_perp_on_Ep_grid.^2);
        
        E_cm = 1/2 * m_D/2 * v_rel.^2;
        
    otherwise
        
        error('Must choose 0 or 1 for Thermal_drift')
        
end

E_cm_keV = E_cm/1e3/1.6e-19;

cross_section_B_H = Bosch_Hale_cross_section(E_cm_keV);
cross_section_B_H = cross_section_B_H*1e-31;

load('mycolormap.mat')

n_t = 5e19;

figure
[~,h] = contourf(energy/1.6e-19/1e3,pitch_analytic,n_t*cross_section_B_H.*v_rel,20)
%xlim([20 100])
%colorbar
set(gca,'fontsize',20)
set(gcf,'colormap',mycolormap)
set(h,'Edgecolor','none')
xlabel('E [keV]')
ylabel('pitch [-]')
axis square
%caxis([0 3]*1e-4)


switch Thermal_drift
    
    case 0
        
        %saveas(gcf,'results\Article_2\rate_function.png')
        %saveas(gcf,'results\Article_2\rate_function.epsc')
        saveas(gcf,'results\Article_2\rate_function.eps','epsc')
        
    case 1 
        %saveas(gcf,'results\Article_2\rate_function_with_drift.png')
        %saveas(gcf,'results\Article_2\rate_function_with_drift.epsc')
        saveas(gcf,'results\Article_2\rate_function_with_drift.eps','epsc')
end


if 0

%load('C:\Local Documents\MATLAB_local\WeightFunctions\Neutron_weight_functions\Neutron_weight_function_19.mat')

%load('C:\Local Documents\MATLAB_local\WeightFunctions\Neutron_weight_functions\Neutron_weight_function_27.mat')

load('C:\Local Documents\MATLAB_local\WeightFunctions\Neutron_weight_functions\Neutron_weight_function_32.mat')

R_numeric = zeros(size(pitch,1),size(E_ion(2:end),1));

for a = 1:length(pitch)
    for b = 1:length(E_ion(2:end))
        R_numeric(a,b) = sum(wf(:,a,b));
    end
end

figure
contourf(E_ion(2:end),pitch,R_numeric,20)
colorbar

figure
contourf(E_ion(2:end),pitch,R_numeric/5e19,20)
colorbar

E_ion_cm = E_ion(2:end)/2;

[E_ION_CM,PITCH_NUM] = meshgrid(E_ion_cm,pitch);

cross_section_B_H_num = Bosch_Hale_cross_section(E_ION_CM);

cross_section_B_H_num(isnan(cross_section_B_H_num)) = 0;

figure
contourf(E_ion(2:end),pitch,cross_section_B_H_num);
colorbar

R_numeric = R_numeric./cross_section_B_H_num;
R_numeric = R_numeric/5e19;

figure
contourf(E_ion(2:end),pitch,R_numeric,20)
colorbar

end