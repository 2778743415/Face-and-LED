clear; close all; clc;

peak_wavelength_red_nm = 635;
fwhm_red_nm = 20;
intensity_max_red = 1.0;

peak_wavelength_blue_nm = 470;
fwhm_blue_nm = 25;
intensity_max_blue = 0.8;

attenuation_factor = 0.7;

led_length = 1.0;
led_width = 1.0;
led_height = 0.4;
source_distance = 0.6;

source_cuboid_dim = 0.06;
src_half_dim = source_cuboid_dim / 2;

wavelengths_nm = linspace(380, 780, 1000);
sigma_red = fwhm_red_nm / (2 * sqrt(2 * log(2)));
spectrum_red_unattenuated = intensity_max_red * exp(-(wavelengths_nm - peak_wavelength_red_nm).^2 / (2 * sigma_red^2));
sigma_blue = fwhm_blue_nm / (2 * sqrt(2 * log(2)));
spectrum_blue_unattenuated = intensity_max_blue * exp(-(wavelengths_nm - peak_wavelength_blue_nm).^2 / (2 * sigma_blue^2));
spectrum_combined_unattenuated = spectrum_red_unattenuated + spectrum_blue_unattenuated;
spectrum_final_in_air = spectrum_combined_unattenuated * attenuation_factor;

figure(1); 
set(gcf, 'Name', '双光源LED在空气中的光谱分析', 'Position', [50, 100, 800, 500]);
plot(wavelengths_nm, spectrum_red_unattenuated, 'r--', 'LineWidth', 1.5, 'DisplayName', sprintf('红光LED (内部, 峰值: %.0fnm)', peak_wavelength_red_nm));
hold on;
plot(wavelengths_nm, spectrum_blue_unattenuated, 'b--', 'LineWidth', 1.5, 'DisplayName', sprintf('蓝光LED (内部, 峰值: %.0fnm)', peak_wavelength_blue_nm));
plot(wavelengths_nm, spectrum_combined_unattenuated, 'm-', 'LineWidth', 2, 'DisplayName', '叠加光谱 (LED内部, 未衰减)');
plot(wavelengths_nm, spectrum_final_in_air, 'k-', 'LineWidth', 2.5, 'DisplayName', sprintf('最终光谱 (在空气中, 衰减后 x%.2f)', attenuation_factor));
title('双光源 LED 在空气中的组合光谱分析 (考虑透出衰减)', 'FontSize', 16, 'FontWeight', 'bold');
xlabel('波长 (nm)', 'FontSize', 14);
ylabel('相对光强 (任意单位)', 'FontSize', 14);
grid on; axis tight; current_ylim = ylim; ylim([0, current_ylim(2) * 1.1]);
legend('show', 'Location', 'northeast', 'FontSize', 10); legend('boxoff');
hold off;

fprintf('--- 双光源 LED 最终在空气中光谱的分析结果 ---\n');
[peaks_final, locs_final_idx] = findpeaks(spectrum_final_in_air, 'MinPeakHeight', 0.05 * max(spectrum_final_in_air), 'MinPeakDistance', round(numel(wavelengths_nm)*0.05) );
locs_final_wl = wavelengths_nm(locs_final_idx);
fprintf('最终在空气中光谱的峰值信息:\n');
if isempty(peaks_final)
    fprintf('  未检测到明显峰值。\n');
else
    for i = 1:length(peaks_final), fprintf('  峰值 %d: 波长 %.2f nm, 相对强度 %.3f\n', i, locs_final_wl(i), peaks_final(i)); end
end
centroid_wavelength_final = sum(wavelengths_nm .* spectrum_final_in_air) / sum(spectrum_final_in_air);
fprintf('最终在空气中光谱的质心波长: %.2f nm\n', centroid_wavelength_final);
total_intensity_final = trapz(wavelengths_nm, spectrum_final_in_air);
fprintf('最终在空气中光谱的总相对光强 (积分面积): %.2f\n', total_intensity_final);

figure(2); 
set(gcf, 'Name', 'LED 发光3D模拟图 (粒子云可视化)', 'Position', [900, 100, 800, 700]);
clf; 

hl = led_length / 2; hw = led_width / 2; hh = led_height;
vertices_led = [
    -hl, -hw, 0; hl, -hw, 0; hl, hw, 0; -hl, hw, 0;
    -hl, -hw, hh; hl, -hw, hh; hl, hw, hh; -hl, hw, hh
];
faces_led = [
    1, 2, 3, 4; 5, 6, 7, 8; 1, 2, 6, 5; 
    2, 3, 7, 6; 3, 4, 8, 7; 4, 1, 5, 8
];
patch('Vertices', vertices_led, 'Faces', faces_led, ...
      'FaceColor', [0.75, 0.75, 0.85], 'EdgeColor', 'k', ...
      'FaceAlpha', 0.20, 'LineWidth', 1); 
hold on;

src_red_center_pos = [-source_distance/2, 0, src_half_dim]; 
src_blue_center_pos = [source_distance/2, 0, src_half_dim];

color_red_source = [1, 0, 0]; color_blue_source = [0, 0, 1];

src_faces = [1,2,3,4; 5,6,7,8; 1,2,6,5; 2,3,7,6; 3,4,8,7; 4,1,5,8];

v_red_src = [
    src_red_center_pos(1)-src_half_dim, src_red_center_pos(2)-src_half_dim, src_red_center_pos(3)-src_half_dim;
    src_red_center_pos(1)+src_half_dim, src_red_center_pos(2)-src_half_dim, src_red_center_pos(3)-src_half_dim;
    src_red_center_pos(1)+src_half_dim, src_red_center_pos(2)+src_half_dim, src_red_center_pos(3)-src_half_dim;
    src_red_center_pos(1)-src_half_dim, src_red_center_pos(2)+src_half_dim, src_red_center_pos(3)-src_half_dim;
    src_red_center_pos(1)-src_half_dim, src_red_center_pos(2)-src_half_dim, src_red_center_pos(3)+src_half_dim;
    src_red_center_pos(1)+src_half_dim, src_red_center_pos(2)-src_half_dim, src_red_center_pos(3)+src_half_dim;
    src_red_center_pos(1)+src_half_dim, src_red_center_pos(2)+src_half_dim, src_red_center_pos(3)+src_half_dim;
    src_red_center_pos(1)-src_half_dim, src_red_center_pos(2)+src_half_dim, src_red_center_pos(3)+src_half_dim;
];
h_patch_red_src = patch('Vertices', v_red_src, 'Faces', src_faces, 'FaceColor', color_red_source, 'EdgeColor', 'k', 'LineWidth', 0.5, 'DisplayName', '红光光源');

v_blue_src = [
    src_blue_center_pos(1)-src_half_dim, src_blue_center_pos(2)-src_half_dim, src_blue_center_pos(3)-src_half_dim;
    src_blue_center_pos(1)+src_half_dim, src_blue_center_pos(2)-src_half_dim, src_blue_center_pos(3)-src_half_dim;
    src_blue_center_pos(1)+src_half_dim, src_blue_center_pos(2)+src_half_dim, src_blue_center_pos(3)-src_half_dim;
    src_blue_center_pos(1)-src_half_dim, src_blue_center_pos(2)+src_half_dim, src_blue_center_pos(3)-src_half_dim;
    src_blue_center_pos(1)-src_half_dim, src_blue_center_pos(2)-src_half_dim, src_blue_center_pos(3)+src_half_dim;
    src_blue_center_pos(1)+src_half_dim, src_blue_center_pos(2)-src_half_dim, src_blue_center_pos(3)+src_half_dim;
    src_blue_center_pos(1)+src_half_dim, src_blue_center_pos(2)+src_half_dim, src_blue_center_pos(3)+src_half_dim;
    src_blue_center_pos(1)-src_half_dim, src_blue_center_pos(2)+src_half_dim, src_blue_center_pos(3)+src_half_dim;
];
h_patch_blue_src = patch('Vertices', v_blue_src, 'Faces', src_faces, 'FaceColor', color_blue_source, 'EdgeColor', 'k', 'LineWidth', 0.5, 'DisplayName', '蓝光光源');

num_particles_top_cone = 700;       
num_particles_per_side_group = 1200; 
angle_spread_top_cone = pi/3.0;     
angle_spread_side = pi/2.8;         

particle_max_dist_top = led_height * 1.8;  
particle_max_dist_side = led_width * 1.1;  

particle_size_top_attenuated = 10; 
particle_size_side_unattenuated = 18;
particle_alpha_top = 0.4;          
particle_alpha_side = 0.6;         
marker_type = 'o';                 

legend_handles = [h_patch_red_src, h_patch_blue_src]; 
legend_labels = {get(h_patch_red_src, 'DisplayName'), get(h_patch_blue_src, 'DisplayName')}; 

plotted_legend_top_red = false; plotted_legend_side_red = false;
plotted_legend_top_blue = false; plotted_legend_side_blue = false;

sources_center_pos = {src_red_center_pos, src_blue_center_pos}; 
base_colors = {color_red_source, color_blue_source};
source_names = {'红光', '蓝光'};

for src_idx = 1:length(sources_center_pos)
    current_source_center = sources_center_pos{src_idx};
    current_base_color = base_colors{src_idx};
    current_source_name = source_names{src_idx};

    P_top_emission_center_on_led = [current_source_center(1), current_source_center(2), hh]; 

    top_particles_x = zeros(num_particles_top_cone, 1);
    top_particles_y = zeros(num_particles_top_cone, 1);
    top_particles_z = zeros(num_particles_top_cone, 1);
    for k_part = 1:num_particles_top_cone
        phi = 2 * pi * rand(); 
        theta = angle_spread_top_cone * sqrt(rand()); 
        dist = particle_max_dist_top * rand()^(1/3); 
        
        dx = sin(theta) * cos(phi);
        dy = sin(theta) * sin(phi);
        dz = cos(theta);
        if dz < 0.01, dz = 0.01 + rand()*0.1; end 
        
        dir_vec = [dx, dy, dz] / norm([dx,dy,dz]);
        
        top_particles_x(k_part) = P_top_emission_center_on_led(1) + dist * dir_vec(1);
        top_particles_y(k_part) = P_top_emission_center_on_led(2) + dist * dir_vec(2);
        top_particles_z(k_part) = P_top_emission_center_on_led(3) + dist * dir_vec(3);
    end
    
    lh_top_particles_disp_name = '';
    if src_idx == 1 && ~plotted_legend_top_red
        lh_top_particles_disp_name = sprintf('%s: 顶部粒子云 (衰减)', current_source_name); 
        plotted_legend_top_red = true;
    elseif src_idx == 2 && ~plotted_legend_top_blue
        lh_top_particles_disp_name = sprintf('%s: 顶部粒子云 (衰减)', current_source_name); 
        plotted_legend_top_blue = true;
    end
    h_top_particles = scatter3(top_particles_x, top_particles_y, top_particles_z, ...
        particle_size_top_attenuated, current_base_color, marker_type, ...
        'filled', 'MarkerFaceAlpha', particle_alpha_top, 'MarkerEdgeAlpha', particle_alpha_top*0.7, ...
        'DisplayName', lh_top_particles_disp_name);
    if ~isempty(lh_top_particles_disp_name) && ishandle(h_top_particles)
        legend_handles(end+1) = h_top_particles;
        legend_labels{end+1} = lh_top_particles_disp_name;
    end

    side_emission_origins = { 
        [hl, current_source_center(2), current_source_center(3) + (rand-0.5)*hh*0.3], 
        [-hl, current_source_center(2), current_source_center(3) + (rand-0.5)*hh*0.3],
        [current_source_center(1), hw, current_source_center(3) + (rand-0.5)*hh*0.3], 
        [current_source_center(1), -hw, current_source_center(3) + (rand-0.5)*hh*0.3] 
    };
    side_main_normals = {[1,0,0], [-1,0,0], [0,1,0], [0,-1,0]}; 

    lh_side_particles_disp_name = '';
     if src_idx == 1 && ~plotted_legend_side_red
        lh_side_particles_disp_name = sprintf('%s: 侧面粒子云 (未衰减)', current_source_name); 
        plotted_legend_side_red = true;
    elseif src_idx == 2 && ~plotted_legend_side_blue
        lh_side_particles_disp_name = sprintf('%s: 侧面粒子云 (未衰减)', current_source_name); 
        plotted_legend_side_blue = true;
    end

    all_side_particles_x = []; all_side_particles_y = []; all_side_particles_z = [];

    for j_side = 1:length(side_emission_origins)
        P_side_emit_origin = side_emission_origins{j_side}; 
        main_normal_side = side_main_normals{j_side};
        
        num_particles_this_side_segment = round(num_particles_per_side_group / length(side_emission_origins));
        side_particles_x_temp = zeros(num_particles_this_side_segment, 1);
        side_particles_y_temp = zeros(num_particles_this_side_segment, 1);
        side_particles_z_temp = zeros(num_particles_this_side_segment, 1);

        for k_part_side = 1:num_particles_this_side_segment
            rand_vec = randn(1,3); 
            ortho_vec = cross(main_normal_side, rand_vec); 
            if norm(ortho_vec) < 1e-6, ortho_vec = cross(main_normal_side, [rand() rand() rand()]+1e-3); end 
            ortho_vec = ortho_vec / norm(ortho_vec);
            
            angle = angle_spread_side * rand(); 
            
            perturb_x = (rand - 0.5) * 2 * sin(angle_spread_side/2);
            perturb_y = (rand - 0.5) * 2 * sin(angle_spread_side/2);
            perturb_z = (rand - 0.5) * 2 * sin(angle_spread_side/2);
            
            dir_vec = main_normal_side + [perturb_x, perturb_y, perturb_z]*0.9; 
            dir_vec = dir_vec / norm(dir_vec);
            
            if dot(dir_vec, main_normal_side) < cos(angle_spread_side * 1.15) 
                phi_side = 2 * pi * rand();
                theta_side = angle_spread_side * sqrt(rand());
                
                if abs(main_normal_side(3)) < 0.95
                    axis1 = cross(main_normal_side, [0,0,1]); axis1 = axis1/norm(axis1);
                else
                    axis1 = cross(main_normal_side, [1,0,0]); axis1 = axis1/norm(axis1);
                end
                axis2 = cross(main_normal_side, axis1); axis2 = axis2/norm(axis2);
                
                dir_vec_local_z = cos(theta_side);
                dir_vec_local_x = sin(theta_side) * cos(phi_side);
                dir_vec_local_y = sin(theta_side) * sin(phi_side);
                
                dir_vec = main_normal_side * dir_vec_local_z + axis1 * dir_vec_local_x + axis2 * dir_vec_local_y;
            end
            dir_vec = dir_vec / norm(dir_vec);
            if dot(dir_vec, main_normal_side) < 0.05 
                dir_vec = main_normal_side; 
            end

            dist_side = particle_max_dist_side * rand()^(1/3);
            
            side_particles_x_temp(k_part_side) = P_side_emit_origin(1) + dist_side * dir_vec(1);
            side_particles_y_temp(k_part_side) = P_side_emit_origin(2) + dist_side * dir_vec(2);
            side_particles_z_temp(k_part_side) = P_side_emit_origin(3) + dist_side * dir_vec(3);
        end
        all_side_particles_x = [all_side_particles_x; side_particles_x_temp];
        all_side_particles_y = [all_side_particles_y; side_particles_y_temp];
        all_side_particles_z = [all_side_particles_z; side_particles_z_temp];
    end
    
    if ~isempty(all_side_particles_x) 
        h_side_particles = scatter3(all_side_particles_x, all_side_particles_y, all_side_particles_z, ...
            particle_size_side_unattenuated, current_base_color, marker_type, ...
            'filled', 'MarkerFaceAlpha', particle_alpha_side, 'MarkerEdgeAlpha', particle_alpha_side*0.7, ...
            'DisplayName', lh_side_particles_disp_name); 
        if ~isempty(lh_side_particles_disp_name) && ishandle(h_side_particles)
            legend_handles(end+1) = h_side_particles;
            legend_labels{end+1} = lh_side_particles_disp_name;
        end
    end
end

title('LED 模型及其在空气中发光的3D模拟图 (粒子云)', 'FontSize', 16, 'FontWeight','bold'); 
xlabel('X 轴 (单位)', 'FontSize', 13);
ylabel('Y 轴 (单位)', 'FontSize', 13);
zlabel('Z 轴 (高度, 单位)', 'FontSize', 13);
axis equal; grid on; view(35, 30); rotate3d on; 

if ~isempty(legend_handles)
    valid_handles = []; unique_display_names = {}; temp_labels = {};
    for i = 1:length(legend_handles)
        if isgraphics(legend_handles(i)) && ~isempty(get(legend_handles(i), 'DisplayName')) 
            current_disp_name = get(legend_handles(i), 'DisplayName');
            if ~ismember(current_disp_name, unique_display_names)
                valid_handles(end+1) = legend_handles(i);
                unique_display_names{end+1} = current_disp_name;
                temp_labels{end+1} = current_disp_name; 
            end
        end
    end
    if ~isempty(valid_handles)
        legend(valid_handles, temp_labels, 'Location', 'northeastoutside', 'FontSize', 10);
        legend('boxoff');
    end
end

camlight('headlight'); lighting gouraud; 
hold off;

disp(' ');
disp('双光源频谱分析和3D粒子云可视化完成。');
disp('图形窗口1: 光谱图');
disp('图形窗口2: 3D LED发光模拟图 (粒子云可视化)');
