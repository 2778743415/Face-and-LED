clear; close all; clc;

% --- 1. 定义光源和系统参数 ---
peak_wavelength_red_nm = 635;
fwhm_red_nm = 20;
intensity_max_red = 1.0;

peak_wavelength_blue_nm = 470;
fwhm_blue_nm = 25;
intensity_max_blue = 0.8;

attenuation_factor_led_surface = 0.7; % LED出射面本身的衰减

% LED 物理尺寸参数
led_length = 1.0; 
led_width = 1.0;  
led_height = 0.4; % LED 器件本身的高度
source_distance = 0.6; 

% 光源几何参数
source_cuboid_dim = 0.06; 
src_half_dim = source_cuboid_dim / 2;

% --- 2. 定义组织层参数 ---
% Z 坐标定义 (cm)，假设 LED 的顶面 (出射面) 在 z = led_height
z_led_emitting_surface = led_height; % = 0.4 cm

z_air_start_abs = z_led_emitting_surface; % 空气层从LED表面开始
z_air_end_abs = 1.4;       % 空气层结束 Z = 1.4 cm
z_skin_start_abs = 1.4;    % 皮肤层开始 Z = 1.4 cm
z_skin_end_abs = 1.8;      % 皮肤层结束 Z = 1.8 cm
z_muscle_start_abs = 1.8;  % 肌肉层开始 Z = 1.8 cm
z_muscle_end_abs = 2.0;    % 肌肉层结束 Z = 2.0 cm
z_vessel_start_abs = 2.0;  % 血管层开始 Z = 2.0 cm
z_vessel_end_abs = 2.2;    % 血管层结束 Z = 2.2 cm

% 计算各层厚度 (cm)
thickness_air = z_air_end_abs - z_air_start_abs;
thickness_skin = z_skin_end_abs - z_skin_start_abs;  
thickness_muscle = z_muscle_end_abs - z_muscle_start_abs; 
thickness_vessel = z_vessel_end_abs - z_vessel_start_abs;

if thickness_air < 0 || thickness_skin < 0 || thickness_muscle < 0 || thickness_vessel < 0
    error('层厚度计算为负值，请检查Z坐标定义。');
end

% 定义波长数组 (nm)
wavelengths_nm = linspace(400, 750, 351); 

% --- 辅助函数：获取组织的吸收系数 mu_a (cm^-1) ---
function mu_a = get_mu_a_skin(lambda_nm_array)
    mu_a = zeros(size(lambda_nm_array));
    mu_a_blue_skin = 4.0; 
    mu_a_red_skin = 0.8;  
    for i = 1:length(lambda_nm_array)
        lambda = lambda_nm_array(i);
        if lambda <= 470
            mu_a(i) = mu_a_blue_skin * (1 + (470-lambda)/100 * 0.5) ; 
        elseif lambda > 470 && lambda < 635
            mu_a(i) = mu_a_blue_skin - (mu_a_blue_skin - mu_a_red_skin) * (lambda - 470) / (635 - 470);
        else 
            mu_a(i) = mu_a_red_skin * (1 - (lambda-635)/100 * 0.3); 
        end
        if mu_a(i) < 0.05, mu_a(i) = 0.05; end 
    end
end

function mu_a = get_mu_a_muscle(lambda_nm_array)
    mu_a = zeros(size(lambda_nm_array));
    mu_a_blue_muscle = 8.0;  
    mu_a_red_muscle = 1.5;   
    for i = 1:length(lambda_nm_array)
        lambda = lambda_nm_array(i);
        if lambda <= 470
            mu_a(i) = mu_a_blue_muscle * (1 + (470-lambda)/100 * 0.6);
        elseif lambda > 470 && lambda < 635
            mu_a(i) = mu_a_blue_muscle - (mu_a_blue_muscle - mu_a_red_muscle) * (lambda - 470) / (635 - 470);
        else 
            mu_a(i) = mu_a_red_muscle * (1 - (lambda-635)/100 * 0.4);
        end
         if mu_a(i) < 0.1, mu_a(i) = 0.1; end
    end
end


function mu_a = get_mu_a_vessel(lambda_nm_array)
    mu_a = zeros(size(lambda_nm_array));
    % 示例: 简化血液吸收模型 (血红蛋白吸收峰在蓝绿光区和近红外Q谷)
    mu_a_blue_blood = 50.0; % cm^-1 at 420nm (Soret band 附近，高吸收) - 简化
    mu_a_green_blood = 30.0; % cm^-1 at 550nm (Q band 附近) - 简化
    mu_a_red_blood = 10.0;  % cm^-1 at 635nm (相对较低) - 简化
    mu_a_nir_blood = 5.0;   % cm^-1 at 750nm (更低) - 简化

    for i = 1:length(lambda_nm_array)
        lambda = lambda_nm_array(i);
        if lambda <= 450
            mu_a(i) = mu_a_blue_blood * exp(-(lambda-420)^2/(2*30^2)); % Soret band 附近高斯近似
        elseif lambda > 450 && lambda <= 590
             mu_a(i) = mu_a_green_blood * exp(-(lambda-550)^2/(2*40^2)); % Q band 附近高斯近似
        elseif lambda > 590 && lambda < 700
             mu_a(i) = mu_a_red_blood + (mu_a_green_blood - mu_a_red_blood) * (1 - (lambda - 590)/(700-590)) * 0.5; % 简化过渡
        else % lambda >= 700
            mu_a(i) = mu_a_red_blood - (mu_a_red_blood - mu_a_nir_blood) * (lambda - 700)/(750-700); % 向NIR下降
        end
        if mu_a(i) < 1.0, mu_a(i) = 1.0; end % 血液吸收通常较高
    end
end


mu_a_air_val = 0.001; 

% --- 3. 生成LED初始光谱 ---
sigma_red = fwhm_red_nm / (2 * sqrt(2 * log(2)));
spectrum_red_source = intensity_max_red * exp(-(wavelengths_nm - peak_wavelength_red_nm).^2 / (2 * sigma_red^2));
sigma_blue = fwhm_blue_nm / (2 * sqrt(2 * log(2)));
spectrum_blue_source = intensity_max_blue * exp(-(wavelengths_nm - peak_wavelength_blue_nm).^2 / (2 * sigma_blue^2));

spectrum_combined_source = spectrum_red_source + spectrum_blue_source;
spectrum_after_led_surface = spectrum_combined_source * attenuation_factor_led_surface; 

% --- 4. 计算光谱在各层中的衰减 ---
mu_a_air = ones(size(wavelengths_nm)) * mu_a_air_val; 
transmittance_air = exp(-mu_a_air .* thickness_air);
spectrum_after_air = spectrum_after_led_surface .* transmittance_air;

mu_a_skin_values = get_mu_a_skin(wavelengths_nm);
transmittance_skin = exp(-mu_a_skin_values .* thickness_skin);
spectrum_after_skin = spectrum_after_air .* transmittance_skin;

mu_a_muscle_values = get_mu_a_muscle(wavelengths_nm);
transmittance_muscle = exp(-mu_a_muscle_values .* thickness_muscle);
spectrum_entering_vessel = spectrum_after_skin .* transmittance_muscle; % 进入血管的光谱

mu_a_vessel_values = get_mu_a_vessel(wavelengths_nm);
transmittance_vessel = exp(-mu_a_vessel_values .* thickness_vessel);
spectrum_exiting_vessel = spectrum_entering_vessel .* transmittance_vessel; % 穿出血管的光谱


figure(1); 
set(gcf, 'Name', 'LED光谱在多层组织中的衰减分析', 'Position', [50, 50, 800, 600]);
plot(wavelengths_nm, spectrum_after_led_surface, 'k-', 'LineWidth', 2, 'DisplayName', 'LED出射');
hold on;
if thickness_air > 1e-3 
    plot(wavelengths_nm, spectrum_after_air, 'c--', 'LineWidth', 1.5, 'DisplayName', sprintf('空气后 (%.2fcm)', thickness_air));
end
plot(wavelengths_nm, spectrum_after_skin, 'g-.', 'LineWidth', 1.5, 'DisplayName', sprintf('皮肤后 (%.2fcm)', thickness_skin));
plot(wavelengths_nm, spectrum_entering_vessel, 'm:', 'LineWidth', 1.5, 'DisplayName', sprintf('肌肉后/入血管 (%.2fcm)', thickness_muscle));
plot(wavelengths_nm, spectrum_exiting_vessel, 'r-', 'LineWidth', 2, 'DisplayName', sprintf('血管后 (%.2fcm)', thickness_vessel));
title('LED光谱在多层组织中的衰减', 'FontSize', 16, 'FontWeight', 'bold');
xlabel('波长 (nm)', 'FontSize', 14);
ylabel('相对光强 (任意单位)', 'FontSize', 14);
grid on; axis tight; 
current_ylim = ylim; ylim([0, current_ylim(2) * 1.1]);
legend('show', 'Location', 'northeast', 'FontSize', 10); legend('boxoff');
hold off;


fprintf('--- LED光在多层组织中衰减分析结果 ---\n');
total_intensity_led_surface = trapz(wavelengths_nm, spectrum_after_led_surface);
total_intensity_after_air = trapz(wavelengths_nm, spectrum_after_air);
total_intensity_after_skin = trapz(wavelengths_nm, spectrum_after_skin);
total_intensity_entering_vessel = trapz(wavelengths_nm, spectrum_entering_vessel);
total_intensity_exiting_vessel = trapz(wavelengths_nm, spectrum_exiting_vessel);

fprintf('LED出射总光强 (相对值): %.4f\n', total_intensity_led_surface);
fprintf('穿过空气层 (厚度 %.2fcm) 后总光强: %.4f (衰减 %.2f%%)\n', thickness_air, total_intensity_after_air, (1-total_intensity_after_air/total_intensity_led_surface)*100);
fprintf('穿过皮肤层 (厚度 %.2fcm) 后总光强: %.4f (相对空气后衰减 %.2f%%)\n', thickness_skin, total_intensity_after_skin, (1-total_intensity_after_skin/total_intensity_after_air)*100);
fprintf('进入血管层 (穿过肌肉厚度 %.2fcm) 总光强: %.4f (相对皮肤后衰减 %.2f%%)\n',thickness_muscle, total_intensity_entering_vessel, (1-total_intensity_entering_vessel/total_intensity_after_skin)*100);
fprintf('穿出血管层 (厚度 %.2fcm) 总光强: %.4f (相对入血管衰减 %.2f%%)\n', thickness_vessel, total_intensity_exiting_vessel, (1-total_intensity_exiting_vessel/total_intensity_entering_vessel)*100);
fprintf('总衰减 (从LED出射到达血管后): %.2f%%\n', (1-total_intensity_exiting_vessel/total_intensity_led_surface)*100);

[~, idx_red_peak] = min(abs(wavelengths_nm - peak_wavelength_red_nm));
[~, idx_blue_peak] = min(abs(wavelengths_nm - peak_wavelength_blue_nm));

intensity_red_entering_vessel = spectrum_entering_vessel(idx_red_peak);
intensity_blue_entering_vessel = spectrum_entering_vessel(idx_blue_peak);
intensity_red_exiting_vessel = spectrum_exiting_vessel(idx_red_peak);
intensity_blue_exiting_vessel = spectrum_exiting_vessel(idx_blue_peak);
fprintf('在 %.0fnm (红光峰值) 进入血管的光强: %.4e\n', peak_wavelength_red_nm, intensity_red_entering_vessel);
fprintf('在 %.0fnm (蓝光峰值) 进入血管的光强: %.4e\n', peak_wavelength_blue_nm, intensity_blue_entering_vessel);
fprintf('在 %.0fnm (红光峰值) 穿出血管的光强: %.4e\n', peak_wavelength_red_nm, intensity_red_exiting_vessel);
fprintf('在 %.0fnm (蓝光峰值) 穿出血管的光强: %.4e\n', peak_wavelength_blue_nm, intensity_blue_exiting_vessel);

figure(2); 
set(gcf, 'Name', 'LED发光及组织层3D模拟图 (分层粒子云与体积渲染)', 'Position', [900, 50, 850, 750]);
clf; 


hl_led = led_length / 2; hw_led = led_width / 2; hh_led_device = led_height; 
vertices_led_body = [-hl_led,-hw_led,0; hl_led,-hw_led,0; hl_led,hw_led,0; -hl_led,hw_led,0; -hl_led,-hw_led,hh_led_device; hl_led,-hw_led,hh_led_device; hl_led,hw_led,hh_led_device; -hl_led,hw_led,hh_led_device];
faces_led_body = [1,2,3,4; 5,6,7,8; 1,2,6,5; 2,3,7,6; 3,4,8,7; 4,1,5,8];
patch('Vertices', vertices_led_body, 'Faces', faces_led_body, 'FaceColor', [0.75,0.75,0.85], 'EdgeColor','k', 'FaceAlpha',0.15, 'LineWidth',1); 
hold on;

src_center_z_offset = src_half_dim; 
src_red_center_pos = [-source_distance/2, 0, src_center_z_offset]; 
src_blue_center_pos = [source_distance/2, 0, src_center_z_offset];
color_red_source = [1,0,0]; color_blue_source = [0,0,1];
src_faces_cuboid = [1,2,3,4; 5,6,7,8; 1,2,6,5; 2,3,7,6; 3,4,8,7; 4,1,5,8];
v_red_src = [src_red_center_pos(1)-src_half_dim, src_red_center_pos(2)-src_half_dim, src_red_center_pos(3)-src_half_dim; src_red_center_pos(1)+src_half_dim, src_red_center_pos(2)-src_half_dim, src_red_center_pos(3)-src_half_dim; src_red_center_pos(1)+src_half_dim, src_red_center_pos(2)+src_half_dim, src_red_center_pos(3)-src_half_dim; src_red_center_pos(1)-src_half_dim, src_red_center_pos(2)+src_half_dim, src_red_center_pos(3)-src_half_dim; src_red_center_pos(1)-src_half_dim, src_red_center_pos(2)-src_half_dim, src_red_center_pos(3)+src_half_dim; src_red_center_pos(1)+src_half_dim, src_red_center_pos(2)-src_half_dim, src_red_center_pos(3)+src_half_dim; src_red_center_pos(1)+src_half_dim, src_red_center_pos(2)+src_half_dim, src_red_center_pos(3)+src_half_dim; src_red_center_pos(1)-src_half_dim, src_red_center_pos(2)+src_half_dim, src_red_center_pos(3)+src_half_dim];
h_patch_red_src = patch('Vertices',v_red_src,'Faces',src_faces_cuboid,'FaceColor',color_red_source,'EdgeColor','k','LineWidth',0.5,'DisplayName','红光光源');
v_blue_src = [src_blue_center_pos(1)-src_half_dim, src_blue_center_pos(2)-src_half_dim, src_blue_center_pos(3)-src_half_dim; src_blue_center_pos(1)+src_half_dim, src_blue_center_pos(2)-src_half_dim, src_blue_center_pos(3)-src_half_dim; src_blue_center_pos(1)+src_half_dim, src_blue_center_pos(2)+src_half_dim, src_blue_center_pos(3)-src_half_dim; src_blue_center_pos(1)-src_half_dim, src_blue_center_pos(2)+src_half_dim, src_blue_center_pos(3)-src_half_dim; src_blue_center_pos(1)-src_half_dim, src_blue_center_pos(2)-src_half_dim, src_blue_center_pos(3)+src_half_dim; src_blue_center_pos(1)+src_half_dim, src_blue_center_pos(2)-src_half_dim, src_blue_center_pos(3)+src_half_dim; src_blue_center_pos(1)+src_half_dim, src_blue_center_pos(2)+src_half_dim, src_blue_center_pos(3)+src_half_dim; src_blue_center_pos(1)-src_half_dim, src_blue_center_pos(2)+src_half_dim, src_blue_center_pos(3)+src_half_dim];
h_patch_blue_src = patch('Vertices',v_blue_src,'Faces',src_faces_cuboid,'FaceColor',color_blue_source,'EdgeColor','k','LineWidth',0.5,'DisplayName','蓝光光源');

% 7.3 绘制组织层体积
layer_plane_xy_half_size = max(led_length, led_width) * 0.9; 

function [V, F] = create_layer_volume(z_bottom, z_top, xy_half_size)
    V = [
        -xy_half_size, -xy_half_size, z_bottom; xy_half_size, -xy_half_size, z_bottom;
         xy_half_size,  xy_half_size, z_bottom; -xy_half_size,  xy_half_size, z_bottom;
        -xy_half_size, -xy_half_size, z_top;    xy_half_size, -xy_half_size, z_top;
         xy_half_size,  xy_half_size, z_top;   -xy_half_size,  xy_half_size, z_top;
    ];
    F = [1,2,3,4; 5,6,7,8; 1,2,6,5; 2,3,7,6; 3,4,8,7; 4,1,5,8]; 
end

color_air_volume = [0.8, 0.9, 1.0]; 
color_skin_volume = [0.96, 0.87, 0.80];    
color_muscle_volume = [0.82, 0.60, 0.55];  
color_vessel_volume = [0.80, 0.25, 0.25];  % Vessel layer color

h_air_vol = gobjects(1);
if thickness_air > 1e-4
    [v_air, f_air] = create_layer_volume(z_air_start_abs, z_air_end_abs, layer_plane_xy_half_size);
    h_air_vol = patch('Vertices', v_air, 'Faces', f_air, 'FaceColor', color_air_volume, 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'DisplayName', sprintf('空气层 (Z=%.2f-%.2fcm)',z_air_start_abs,z_air_end_abs));
end

h_skin_vol = gobjects(1);
if thickness_skin > 1e-4
    [v_skin, f_skin] = create_layer_volume(z_skin_start_abs, z_skin_end_abs, layer_plane_xy_half_size);
    h_skin_vol = patch('Vertices', v_skin, 'Faces', f_skin, 'FaceColor', color_skin_volume, 'FaceAlpha', 0.25, 'EdgeColor', 'none', 'DisplayName', sprintf('皮肤层体积 (Z=%.2f-%.2fcm)',z_skin_start_abs,z_skin_end_abs));
end

h_muscle_vol = gobjects(1);
if thickness_muscle > 1e-4
    [v_muscle, f_muscle] = create_layer_volume(z_muscle_start_abs, z_muscle_end_abs, layer_plane_xy_half_size);
    h_muscle_vol = patch('Vertices', v_muscle, 'Faces', f_muscle, 'FaceColor', color_muscle_volume, 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', sprintf('肌肉层体积 (Z=%.2f-%.2fcm)',z_muscle_start_abs,z_muscle_end_abs));
end

h_vessel_vol = gobjects(1);
if thickness_vessel > 1e-4
    [v_vessel, f_vessel] = create_layer_volume(z_vessel_start_abs, z_vessel_end_abs, layer_plane_xy_half_size);
    h_vessel_vol = patch('Vertices', v_vessel, 'Faces', f_vessel, 'FaceColor', color_vessel_volume, 'FaceAlpha', 0.35, 'EdgeColor', 'none', 'DisplayName', sprintf('血管层体积 (Z=%.2f-%.2fcm)',z_vessel_start_abs,z_vessel_end_abs));
end


% 7.4 模拟和绘制分层粒子云
max_initial_particles_per_source = 1500; 
angle_spread_cone = pi/3.0;     

particle_size_base = 12; 
particle_alpha_base = 0.5;        
marker_type = 'o';                 

legend_handles = [h_patch_red_src, h_patch_blue_src]; 
legend_labels = {get(h_patch_red_src, 'DisplayName'), get(h_patch_blue_src, 'DisplayName')}; 

if isgraphics(h_air_vol), legend_handles(end+1)=h_air_vol; legend_labels{end+1}=get(h_air_vol,'DisplayName'); end
if isgraphics(h_skin_vol), legend_handles(end+1)=h_skin_vol; legend_labels{end+1}=get(h_skin_vol,'DisplayName'); end
if isgraphics(h_muscle_vol), legend_handles(end+1)=h_muscle_vol; legend_labels{end+1}=get(h_muscle_vol,'DisplayName'); end
if isgraphics(h_vessel_vol), legend_handles(end+1)=h_vessel_vol; legend_labels{end+1}=get(h_vessel_vol,'DisplayName'); end % Add vessel volume to legend


plotted_legend_particles_air_red = false; plotted_legend_particles_skin_red = false; plotted_legend_particles_muscle_red = false; plotted_legend_particles_vessel_red = false;
plotted_legend_particles_air_blue = false; plotted_legend_particles_skin_blue = false; plotted_legend_particles_muscle_blue = false; plotted_legend_particles_vessel_blue = false;

sources_center_pos_for_emission = {src_red_center_pos, src_blue_center_pos}; 
base_colors = {color_red_source, color_blue_source};
source_names = {'红光', '蓝光'};

function [px, py, pz, h_scatter, legend_flag_updated] = generate_particle_cloud( ...
    num_particles, start_z_layer, layer_thickness_current, emission_center_xy, ...
    angle_spread, particle_color, particle_size, particle_alpha, marker, ...
    display_name_base, source_name_current, legend_flag_current, led_emission_surface_z_coord) 
    
    px = zeros(num_particles, 1);
    py = zeros(num_particles, 1);
    pz = zeros(num_particles, 1);
    h_scatter = gobjects(1); 
    legend_flag_updated = legend_flag_current;

    if num_particles == 0 || layer_thickness_current <= 1e-5 
        return; 
    end

    for k = 1:num_particles
        phi = 2 * pi * rand(); 
        theta = angle_spread * sqrt(rand()); 
        
        particle_depth_in_layer = layer_thickness_current * rand();
        current_pz = start_z_layer + particle_depth_in_layer;
        
        dist_from_led_surface_z = current_pz - led_emission_surface_z_coord; 
        if dist_from_led_surface_z < 1e-4, dist_from_led_surface_z = 1e-4; end 

        radial_spread_at_pz = dist_from_led_surface_z * tan(theta);
        
        px(k) = emission_center_xy(1) + radial_spread_at_pz * cos(phi);
        py(k) = emission_center_xy(2) + radial_spread_at_pz * sin(phi);
        pz(k) = current_pz;
    end
    
    current_disp_name = '';
    if ~legend_flag_current
        current_disp_name = sprintf('%s: %s', source_name_current, display_name_base);
        legend_flag_updated = true;
    end
    
    h_scatter = scatter3(px, py, pz, particle_size, particle_color, marker, ...
        'filled', 'MarkerFaceAlpha', particle_alpha, 'MarkerEdgeAlpha', particle_alpha*0.7, ...
        'DisplayName', current_disp_name);
end

for src_idx = 1:length(sources_center_pos_for_emission)
    current_source_emission_origin = sources_center_pos_for_emission{src_idx}; 
    current_base_color = base_colors{src_idx};
    current_source_name = source_names{src_idx};
    emission_center_xy_on_led_surface = [current_source_emission_origin(1), current_source_emission_origin(2)];


    num_p_air = round(max_initial_particles_per_source * (total_intensity_after_air / total_intensity_led_surface));
    current_legend_flag_air = ifthen(src_idx==1, plotted_legend_particles_air_red, plotted_legend_particles_air_blue);
    if thickness_air > 1e-4 
        [~,~,~,h_scatter_air, new_flag_air] = generate_particle_cloud(num_p_air, z_air_start_abs, thickness_air, emission_center_xy_on_led_surface, ...
                                               angle_spread_cone, current_base_color, particle_size_base, particle_alpha_base, marker_type, ...
                                               '空气中粒子', current_source_name, current_legend_flag_air, z_led_emitting_surface); 
        if src_idx==1, plotted_legend_particles_air_red = new_flag_air; else, plotted_legend_particles_air_blue = new_flag_air; end
        if ~isempty(get(h_scatter_air,'DisplayName')) && isgraphics(h_scatter_air), legend_handles(end+1)=h_scatter_air; legend_labels{end+1}=get(h_scatter_air,'DisplayName'); end
    end

    num_p_skin = round(max_initial_particles_per_source * (total_intensity_after_skin / total_intensity_led_surface));
    current_legend_flag_skin = ifthen(src_idx==1, plotted_legend_particles_skin_red, plotted_legend_particles_skin_blue);
     if thickness_skin > 1e-4
        [~,~,~,h_scatter_skin, new_flag_skin] = generate_particle_cloud(num_p_skin, z_skin_start_abs, thickness_skin, emission_center_xy_on_led_surface, ...
                                                 angle_spread_cone, current_base_color, particle_size_base*0.9, particle_alpha_base*0.85, marker_type, ...
                                                 '皮肤中粒子', current_source_name, current_legend_flag_skin, z_led_emitting_surface); 
        if src_idx==1, plotted_legend_particles_skin_red = new_flag_skin; else, plotted_legend_particles_skin_blue = new_flag_skin; end
        if ~isempty(get(h_scatter_skin,'DisplayName')) && isgraphics(h_scatter_skin), legend_handles(end+1)=h_scatter_skin; legend_labels{end+1}=get(h_scatter_skin,'DisplayName'); end
     end
    
    num_p_muscle = round(max_initial_particles_per_source * (total_intensity_entering_vessel / total_intensity_led_surface));
    current_legend_flag_muscle = ifthen(src_idx==1, plotted_legend_particles_muscle_red, plotted_legend_particles_muscle_blue);
    if thickness_muscle > 1e-4
        [~,~,~,h_scatter_muscle, new_flag_muscle] = generate_particle_cloud(num_p_muscle, z_muscle_start_abs, thickness_muscle, emission_center_xy_on_led_surface, ...
                                                   angle_spread_cone, current_base_color, particle_size_base*0.8, particle_alpha_base*0.7, marker_type, ...
                                                   '肌肉中粒子', current_source_name, current_legend_flag_muscle, z_led_emitting_surface); 
        if src_idx==1, plotted_legend_particles_muscle_red = new_flag_muscle; else, plotted_legend_particles_muscle_blue = new_flag_muscle; end
        if ~isempty(get(h_scatter_muscle,'DisplayName')) && isgraphics(h_scatter_muscle), legend_handles(end+1)=h_scatter_muscle; legend_labels{end+1}=get(h_scatter_muscle,'DisplayName'); end
    end
    
    num_p_vessel = round(max_initial_particles_per_source * (total_intensity_exiting_vessel / total_intensity_led_surface)); % Particles remaining after vessel
    current_legend_flag_vessel = ifthen(src_idx==1, plotted_legend_particles_vessel_red, plotted_legend_particles_vessel_blue);
    if thickness_vessel > 1e-4
        [~,~,~,h_scatter_vessel, new_flag_vessel] = generate_particle_cloud(num_p_vessel, z_vessel_start_abs, thickness_vessel, emission_center_xy_on_led_surface, ...
                                                   angle_spread_cone, current_base_color, particle_size_base*0.7, particle_alpha_base*0.6, marker_type, ...
                                                   '血管中粒子', current_source_name, current_legend_flag_vessel, z_led_emitting_surface); 
        if src_idx==1, plotted_legend_particles_vessel_red = new_flag_vessel; else, plotted_legend_particles_vessel_blue = new_flag_vessel; end
        if ~isempty(get(h_scatter_vessel,'DisplayName')) && isgraphics(h_scatter_vessel), legend_handles(end+1)=h_scatter_vessel; legend_labels{end+1}=get(h_scatter_vessel,'DisplayName'); end
    end
end

text_z_offset = 0.05; 
intensity_initial_red_led_exit = spectrum_after_led_surface(idx_red_peak); 
intensity_initial_blue_led_exit = spectrum_after_led_surface(idx_blue_peak);

intensity_skin_red_exit = spectrum_after_skin(idx_red_peak);
intensity_skin_blue_exit = spectrum_after_skin(idx_blue_peak);
text(-layer_plane_xy_half_size*0.9, layer_plane_xy_half_size*0.9, z_skin_end_abs + text_z_offset, ... 
    sprintf('皮肤出口: 红光剩%.1f%%\n             蓝光剩%.1f%%', ...
    (intensity_skin_red_exit/intensity_initial_red_led_exit)*100, (intensity_skin_blue_exit/intensity_initial_blue_led_exit)*100), ...
    'Color', 'k', 'FontSize', 8, 'BackgroundColor', [1,1,1,0.6], 'EdgeColor', 'k', 'HorizontalAlignment', 'left');

intensity_vessel_red_entry = spectrum_entering_vessel(idx_red_peak);
intensity_vessel_blue_entry = spectrum_entering_vessel(idx_blue_peak);
text(-layer_plane_xy_half_size*0.9, layer_plane_xy_half_size*0.9, z_vessel_start_abs + text_z_offset, ... 
    sprintf('血管入口: 红光剩%.1f%%\n             蓝光剩%.1f%%', ...
    (intensity_vessel_red_entry/intensity_initial_red_led_exit)*100, (intensity_vessel_blue_entry/intensity_initial_blue_led_exit)*100), ...
    'Color', 'k', 'FontSize', 8, 'BackgroundColor', [1,1,1,0.6], 'EdgeColor', 'k', 'HorizontalAlignment', 'left');
    
intensity_vessel_red_exit_val = spectrum_exiting_vessel(idx_red_peak);
intensity_vessel_blue_exit_val = spectrum_exiting_vessel(idx_blue_peak);
text(-layer_plane_xy_half_size*0.9, layer_plane_xy_half_size*0.9, z_vessel_end_abs + text_z_offset, ... 
    sprintf('血管出口: 红光剩%.1f%%\n             蓝光剩%.1f%% (相对LED出射)', ...
    (intensity_vessel_red_exit_val/intensity_initial_red_led_exit)*100, (intensity_vessel_blue_exit_val/intensity_initial_blue_led_exit)*100), ...
    'Color', 'k', 'FontSize', 8, 'BackgroundColor', [1,1,1,0.6], 'EdgeColor', 'k', 'HorizontalAlignment', 'left');


title('LED发光及在组织层中衰减的3D示意图 (分层粒子云与体积渲染)', 'FontSize', 16, 'FontWeight','bold'); 
xlabel('X (cm)', 'FontSize', 13); ylabel('Y (cm)', 'FontSize', 13); zlabel('Z (cm, 深度)', 'FontSize', 13);
axis equal; grid on; view(35, 30); rotate3d on; 
max_z_coord_plot = z_vessel_end_abs + thickness_vessel * 0.5; % Adjusted Z limit
max_xy_coord_plot = layer_plane_xy_half_size * 1.1;
xlim([-max_xy_coord_plot max_xy_coord_plot]); ylim([-max_xy_coord_plot max_xy_coord_plot]); zlim([0 max_z_coord_plot]);

if ~isempty(legend_handles)
    valid_handles = []; unique_display_names = {}; temp_labels = {};
    for i = 1:length(legend_handles)
        if isgraphics(legend_handles(i)) && ~isempty(get(legend_handles(i), 'DisplayName')) 
            current_disp_name = get(legend_handles(i), 'DisplayName');
            is_already_in_legend = false;
            for k_leg = 1:length(unique_display_names) 
                if strcmp(current_disp_name, unique_display_names{k_leg})
                    is_already_in_legend = true;
                    break;
                end
            end
            if ~is_already_in_legend
                valid_handles(end+1) = legend_handles(i);
                unique_display_names{end+1} = current_disp_name;
                temp_labels{end+1} = current_disp_name; 
            end
        end
    end
    if ~isempty(valid_handles)
        legend(valid_handles, temp_labels, 'Location', 'northeastoutside', 'FontSize', 8); 
        legend('boxoff');
    end
end

camlight('headlight'); lighting gouraud; 
hold off;

disp(' ');
disp('多层组织衰减分析和3D可视化完成。');
disp('图形窗口1: 光谱衰减图');
disp('图形窗口2: 3D LED发光及组织层示意图 (分层粒子云与体积渲染)');

function result = ifthen(condition, true_val, false_val)
    if condition
        result = true_val;
    else
        result = false_val;
    end
end
