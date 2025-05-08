
function stl_face_gui_mcmatlab_modified
    % 创建主界面
    fig = figure('Name', 'STL眼镜模型与人脸曲面可视化 (蒙特卡洛 - 批量LED组合与多角度旋转)', ...
                 'NumberTitle', 'off', ...
                 'Position', [100 100 850 550]);

    % --- UI Controls ---
    uicontrol('Style', 'pushbutton', 'String', '加载眼镜与人脸', ...
              'FontSize', 12, 'Position', [20 490 150 30], ...
              'Callback', @load_models_callback);
    uicontrol('Style', 'pushbutton', 'String', '改变皮肤颜色', ...
              'FontSize', 12, 'Position', [20 450 150 30], ...
              'Callback', @change_surface_color_callback);
    uicontrol('Style', 'pushbutton', 'String', '输入光源模块中心 (手动)', ...
              'FontSize', 10, 'Position', [20 410 200 30], ...
              'Callback', @input_points_and_set_lights_callback);
    uicontrol('Style', 'pushbutton', 'String', '输出当前光强数据 (手动MC)', ...
              'FontSize', 10, 'Position', [20 370 200 30], ...
              'Callback', @export_light_intensity_callback);
    uicontrol('Style', 'pushbutton', 'String', '测试光源模块 (手动)', ...
              'FontSize', 10, 'Position', [20 330 200 30], ...
              'Callback', @test_point_callback);
    
    uicontrol('Style', 'pushbutton', 'String', '运行所有预定义LED组合 (多角度)', ...
              'FontSize', 11, 'FontWeight', 'bold', 'Position', [20 280 210 40], ...
              'BackgroundColor', [0.8 1 0.8], ... 
              'Callback', @run_all_combinations_multi_angle_callback);
    uicontrol('Style', 'text', 'String', '表面漫反射系数:', ...
              'FontSize', 10, 'Position', [20 240 100 20], 'HorizontalAlignment', 'right');
    h_diffuse_slider = uicontrol('Style', 'slider', ...
                             'Position', [125 240 70 20], 'Min', 0, 'Max', 1, 'Value', 0.9, ...
                             'Callback', @diffuse_intensity_slider_callback);
    h_diffuse_value_display = uicontrol('Style', 'text', 'String', '0.90', ...
                                      'FontSize', 10, 'Position', [200 240 40 20], 'HorizontalAlignment', 'left');
    
    uicontrol('Style', 'pushbutton', 'String', '输出所有历史模块中心 (手动输入)', ...
              'FontSize', 9, 'Position', [20 200 200 30], ...
              'Callback', @export_all_history_callback);
    
    ax = axes('Parent', fig, 'Units', 'pixels', 'Position', [250 50 580 480]);
    
    h_surf_skin = []; 
    h_surf_glasses = []; 
    points_data = []; 
    h_all_points_markers = []; 
    h_custom_lights = []; 
    h_rectangle_plots = []; 
    all_input_points_history = {}; 
    current_light_color = [1 0 0]; 

    num_rays_per_light_source_param = 200; 
    MAX_RAY_BOUNCES = 2; 
    RAY_ENERGY_THRESHOLD = 1e-4; 
    EPSILON_SHIFT = 1e-5; 
    TOP_FILTER_ATTENUATION_COEFFICIENT = 0.7; 
    ROTATION_ANGLE_DEGREES = -30; 
    
    MATERIAL_SKIN = 1;
    MATERIAL_LED_OPAQUE_SIDE = 2;
    MATERIAL_LED_FILTER_NON_ATTENUATING = 3; 
    MATERIAL_LED_FILTER_ATTENUATING = 4;     
    
    g_skin_vertices = []; 
    g_skin_triangle_indices = []; 
    g_led_boxes_vertices = []; 
    g_led_boxes_triangle_indices = []; 
    g_led_boxes_triangle_materials = []; 
    
    function load_models_callback(~, ~)
        glassesFilePath = 'C:\\Users\\WC\\Desktop\\new\\glasses.stl'; 

        cla(ax); 
        if ~isempty(h_surf_skin) && isgraphics(h_surf_skin)
            delete(h_surf_skin); 
            h_surf_skin = []; 
        end
        if ~isempty(h_surf_glasses) && isgraphics(h_surf_glasses)
            delete(h_surf_glasses); 
            h_surf_glasses = []; 
        end

        g_skin_vertices = []; 
        g_skin_triangle_indices = []; 
        g_led_boxes_vertices = []; 
        g_led_boxes_triangle_indices = []; 
        g_led_boxes_triangle_materials = [];
        set_lights_and_markers([]); 
        all_input_points_history = {}; 
        
        axes(ax); 
        hold on;  

        try
            model_glasses = stlread(glassesFilePath);
            h_surf_glasses = trisurf(model_glasses.ConnectivityList, ...
                                   model_glasses.Points(:,1), model_glasses.Points(:,2), model_glasses.Points(:,3), ...
                                   'Parent', ax, 'FaceColor', [0.3,0.3,0.4], 'EdgeColor', 'none', 'FaceAlpha', 0.35);
            disp(['成功加载眼镜模型: ' glassesFilePath]);
        catch ME
            errordlg(['读取 STL 文件 (glasses.stl) 出错: ' ME.message], 'STL 读取错误');
        end

        v_center_face = [0,-80,0];      
        x_r = linspace(-70,70,50);      
        z_r = linspace(-60,60,50);      
        [Xf, Zf] = meshgrid(x_r, z_r);  
        a_el=0.008; b_el=0.004;         
        Yf = v_center_face(2) + (a_el*(Xf-v_center_face(1)).^2 + b_el*(Zf-v_center_face(3)).^2)-10; 
        
        h_surf_skin = surf(Xf, Yf, Zf, 'Parent', ax); 
        
        % --- 从生成的表面提取顶点和面片数据 (MODIFIED) ---
        if ~isempty(h_surf_skin) && isgraphics(h_surf_skin)
            % 使用 surf2patch 转换为三角面片
            fv_skin = surf2patch(h_surf_skin, 'triangles');
            g_skin_vertices = fv_skin.vertices;
            g_skin_triangle_indices = fv_skin.faces;

            % 检查提取是否成功
            if isempty(g_skin_vertices) || isempty(g_skin_triangle_indices)
                errordlg('未能从生成的模拟人脸曲面提取顶点/面片数据。光照模拟可能无法进行。', '表面数据提取错误');
                % 清理以确保后续检查能正确报告“模型未加载”
                g_skin_vertices = [];
                g_skin_triangle_indices = [];
                if isgraphics(h_surf_skin), delete(h_surf_skin); h_surf_skin = []; end
                % 不直接 return，允许后续的显示设置，但光照模拟会失败
            end
        else
            errordlg('未能生成模拟人脸曲面对象。光照模拟无法进行。', '表面生成错误');
            g_skin_vertices = []; % 确保这些为空
            g_skin_triangle_indices = [];
            % h_surf_skin 此时应为空或无效
        end
        % --- 提取结束 ---
        
        init_diff = 0.9; 
        if ~isempty(h_surf_skin) && isgraphics(h_surf_skin) % 仅当皮肤表面有效时设置其属性
            set(h_surf_skin,'FaceColor',[1,0.8,0.6],'EdgeColor','none','FaceAlpha',0.9,...
                            'SpecularStrength',0,'SpecularExponent',10); 
        end
        
        set(h_diffuse_slider,'Value',init_diff); 
        set(h_diffuse_value_display,'String',sprintf('%.2f',init_diff));
        update_surface_lighting(); 
        
        daspect(ax,[1 1 1]);     
        lighting(ax,'gouraud');  
        xlabel(ax,'X'); ylabel(ax,'Y'); zlabel(ax,'Z');
        title(ax,'眼镜模型与模拟人脸曲面'); 
        view(ax,3); 
        hold off; 
        
        if ~isempty(g_skin_vertices) && ~isempty(g_skin_triangle_indices)
            helpdlg('已自动加载眼镜模型并生成模拟人脸曲面。光照模拟将针对人脸曲面进行。', '加载与生成完成');
        else
            helpdlg('已自动加载眼镜模型，但生成或提取模拟人脸曲面数据时遇到问题。请检查错误提示。', '加载完成但有问题');
        end
    end

    function change_surface_color_callback(~, ~)
        if ~isempty(h_surf_skin) && isgraphics(h_surf_skin) 
            curr_color = get(h_surf_skin, 'FaceColor');
            new_color = uisetcolor(curr_color, '选择皮肤基础颜色 (预览)');
            if ~isequal(new_color, 0) 
                set(h_surf_skin, 'FaceColor', new_color);
            end
        else
            warndlg('请先加载模型以创建皮肤曲面 (点击 "加载眼镜与人脸")', '提示');
        end
    end

    function input_points_and_set_lights_callback(~, ~)
        prompt_num = {'请输入光源模块的数量 (手动测试):'}; 
        dlg_title_num = '输入光源模块数量 (手动)';
        num_pts_str = inputdlg(prompt_num, dlg_title_num, [1 40], {'1'});
        if isempty(num_pts_str), return; end 
        num_cuboids = str2double(num_pts_str{1}); 
        
        if ~(isscalar(num_cuboids) && isfinite(num_cuboids) && num_cuboids >= 0 && num_cuboids == round(num_cuboids))
            errordlg('请输入一个有效的非负整数作为模块数量。', '输入错误'); return;
        end
        if num_cuboids == 0 
            set_lights_and_markers([]);
            helpdlg('已清除所有手动设置的光源模块。', '操作完成');
            return;
        end

        if ~isvalid(ax) 
            errordlg('绘图区域无效。', '内部错误'); return;
        end
        
        temp_centers = zeros(num_cuboids, 3); 
        valid_cnt = 0;
        for i_src = 1:num_cuboids
            prompt_coords = {sprintf('模块中心 X%d:', i_src), sprintf('模块中心 Y%d:', i_src), sprintf('模块中心 Z%d:', i_src)};
            dlg_title_coords = sprintf('输入光源模块 %d 的中心坐标 (手动)', i_src);
            default_coords_str = {'0','-70','0'}; 
            coords_str = inputdlg(prompt_coords, dlg_title_coords, [1 35; 1 35; 1 35], default_coords_str);
            
            if isempty(coords_str) 
                warning('GUI:PointInputCancelled', '模块 %d 输入取消。', i_src); 
                continue; 
            end
            center_c = str2double(coords_str); 
            
            if any(isnan(center_c)) || numel(center_c) ~= 3
                warning('GUI:InvalidPointCoords', '模块 %d 中心坐标无效 (非数字或数量不足)。跳过此模块。', i_src); 
                continue; 
            end
            valid_cnt = valid_cnt + 1; 
            temp_centers(valid_cnt, :) = center_c'; 
        end
        
        if valid_cnt == 0
            warndlg('没有成功输入任何光源模块中心。', '提示'); 
            set_lights_and_markers([]); 
            return; 
        end
        
        curr_centers = temp_centers(1:valid_cnt, :);
        set_lights_and_markers(curr_centers); 
        
        for k_h = 1:size(curr_centers, 1)
            current_center_to_add = curr_centers(k_h, :);
            is_dup = false; 
            tol = 1e-6; 
            for hist_idx = 1:length(all_input_points_history)
                if isequal(size(all_input_points_history{hist_idx}), size(current_center_to_add)) && ...
                   all(abs(all_input_points_history{hist_idx} - current_center_to_add) < tol)
                    is_dup = true; 
                    break;
                end
            end
            if ~is_dup
                all_input_points_history{end+1} = current_center_to_add;
            end
        end
        
        msg_str = sprintf('成功手动设置 %d 个光源模块。', valid_cnt);
        if valid_cnt < num_cuboids
            msg_str = [msg_str sprintf(' (原计划 %d 个模块, %d 个被跳过或取消)。', num_cuboids, num_cuboids - valid_cnt)];
        end
        helpdlg(msg_str, '手动设置完成');
    end

    function run_all_combinations_multi_angle_callback(~, ~)
        if isempty(g_skin_vertices) || isempty(g_skin_triangle_indices)
            warndlg('模拟人脸曲面数据未加载或无效。请先点击 "加载眼镜与人脸"。', '错误：无皮肤数据'); return;
        end
        if ~isgraphics(h_surf_skin) 
             warndlg('模拟人脸预览曲面无效。请重新加载模型。', '错误：皮肤预览无效'); return;
        end

        led_groups_config = { ...
            {{[38,-84,23]}, {[-38,-84,23]}}, 
            {{[55,-84,-15]}, {[-55,-84,-15]}}, 
            {{[27,-83,-18]}, {[-27,-83,-18]}}, 
            {{[61,-85,19]}, {[-61,-85,19]}}, 
            {{[20,-86,22]}, {[-20,-86,22]}}, 
            {{[12,-95,2]}, {[-12,-95,2]}},   
            {{[68,-80,2]}, {[-68,-80,2]}}    
        };
        num_led_groups = length(led_groups_config);
        total_combinations = 2^num_led_groups; 
        
        angles_to_process_magnitudes = [30, 40, 50, 60, 70]; 
        base_output_dir = 'C:\\Users\\WC\\Desktop\\test_output_matlab'; 
        if ~exist(base_output_dir, 'dir'), mkdir(base_output_dir); end
        
        num_angles = length(angles_to_process_magnitudes);
        total_runs = num_angles * total_combinations; 
        run_counter = 0;
        
        h_main_prog = waitbar(0, sprintf('总进度: 0/%d 角度-组合', total_runs), ...
                               'Name', '处理所有LED组合与角度', 'CreateCancelBtn', 'setappdata(gcbf,''canceling_all'',1)');
        setappdata(h_main_prog,'canceling_all',0); 

        for angle_idx = 1:num_angles
            current_angle_magnitude = angles_to_process_magnitudes(angle_idx);
            ROTATION_ANGLE_DEGREES = -current_angle_magnitude; 
            
            angle_specific_output_dir = fullfile(base_output_dir, sprintf('angle_%d_deg', current_angle_magnitude));
            if ~exist(angle_specific_output_dir, 'dir'), mkdir(angle_specific_output_dir); end
            
            disp(sprintf('--- 开始处理角度: %d degrees (全局旋转变量 ROTATION_ANGLE_DEGREES = %.1f) ---', current_angle_magnitude, ROTATION_ANGLE_DEGREES));
            
            for i_combo = 0:(total_combinations - 1) 
                run_counter = run_counter + 1;
                if getappdata(h_main_prog,'canceling_all'), break; end 
                
                waitbar(run_counter / total_runs, h_main_prog, ...
                        sprintf('角度 %d° (%d/%d), 组合 %d/%d. 总进度: %d/%d', ...
                        current_angle_magnitude, angle_idx, num_angles, ...
                        i_combo, total_combinations - 1, run_counter, total_runs));
                
                active_module_centers_for_combo = []; 
                active_group_names_str = ''; 
                
                for k_group = 1:num_led_groups
                    if bitget(i_combo, k_group) 
                        current_group_modules = led_groups_config{k_group};
                        for m_idx = 1:length(current_group_modules)
                            active_module_centers_for_combo = [active_module_centers_for_combo; current_group_modules{m_idx}{:}];
                        end
                        if isempty(active_group_names_str)
                            active_group_names_str = sprintf('G%d',k_group);
                        else
                            active_group_names_str = [active_group_names_str, sprintf('-G%d',k_group)];
                        end
                    end
                end
                
                set_lights_and_markers(active_module_centers_for_combo); 
                
                if isempty(active_module_centers_for_combo) 
                    if i_combo == 0
                        active_group_names_str = 'all_off';
                        disp(sprintf('组合 %d (%s), 角度 %d°: 所有LED关闭，将生成零强度文件。', i_combo, active_group_names_str, current_angle_magnitude));
                    else 
                         active_group_names_str = sprintf('empty_combo_%d', i_combo); 
                         warning('组合 %d (非 all_off) 角度 %d° 意外地没有活动光源点。将生成零强度文件。', i_combo, current_angle_magnitude);
                    end
                else
                    disp(sprintf('组合 %d (%s), 角度 %d°: 开始蒙特卡洛计算...', i_combo, active_group_names_str, current_angle_magnitude));
                end
                
                combo_filename = sprintf('combo_%s_angle%d.txt', active_group_names_str, current_angle_magnitude);
                combo_filepath = fullfile(angle_specific_output_dir, combo_filename);
                
                V_scene_combined = [g_skin_vertices; g_led_boxes_vertices]; 
                
                num_skin_tri = size(g_skin_triangle_indices,1);
                led_tri_offset = size(g_skin_vertices,1); 
                
                adj_led_box_tri_idx = [];
                if ~isempty(g_led_boxes_triangle_indices) 
                    adj_led_box_tri_idx = g_led_boxes_triangle_indices + led_tri_offset;
                end
                
                T_scene_idx = [g_skin_triangle_indices; adj_led_box_tri_idx];
                
                M_scene_materials = [repmat(MATERIAL_SKIN,num_skin_tri,1); g_led_boxes_triangle_materials]; 
                                
                curr_diff_coeff = get(h_diffuse_slider,'Value'); 
                skin_base_mc_color = get(h_surf_skin,'FaceColor'); 
                
                LightPos_MC_combo = points_data; 
                AccumulatedEnergy_SkinVertex = zeros(size(g_skin_vertices,1),3); 
                
                if ~isempty(LightPos_MC_combo)
                    BaseLightPwr_src = repmat(current_light_color,size(LightPos_MC_combo,1),1); 
                    num_eff_lights_combo = size(LightPos_MC_combo,1);
                    
                    h_mc_inner = waitbar(0,sprintf('角度 %d°, 组合 %s: 光源...',current_angle_magnitude, active_group_names_str), ...
                                         'Name','单次蒙特卡洛模拟','CreateCancelBtn','setappdata(gcbf,"canceling_mc_inner",1)');
                    setappdata(h_mc_inner,'canceling_mc_inner',0); 
                    
                    try
                        for j_eff = 1:num_eff_lights_combo 
                            if getappdata(h_main_prog,'canceling_all') || getappdata(h_mc_inner,'canceling_mc_inner'), break; end
                            
                            waitbar(j_eff/num_eff_lights_combo, h_mc_inner, ...
                                    sprintf('角度 %d°, 组合 %s: 光源 %d/%d',current_angle_magnitude,active_group_names_str,j_eff,num_eff_lights_combo));
                            
                            light_start_pos = LightPos_MC_combo(j_eff,:);
                            initial_energy_per_ray_rgb = BaseLightPwr_src(j_eff,:)/num_rays_per_light_source_param; 
                            
                            for r_idx = 1:num_rays_per_light_source_param 
                                if getappdata(h_main_prog,'canceling_all') || getappdata(h_mc_inner,'canceling_mc_inner'), break; end
                                
                                current_ray_origin = light_start_pos;
                                phi_rand = 2*pi*rand(); 
                                costheta_rand = 2*rand()-1; 
                                sintheta_rand = sqrt(max(0,1-costheta_rand^2));
                                current_ray_direction = [sintheta_rand*cos(phi_rand), sintheta_rand*sin(phi_rand), costheta_rand];
                                current_ray_energy_rgb = initial_energy_per_ray_rgb;
                                
                                for bounce_n = 1:MAX_RAY_BOUNCES
                                    if sum(current_ray_energy_rgb(:)) < RAY_ENERGY_THRESHOLD, break; end 
                                    
                                    min_t_intersection = inf; 
                                    hit_found = false; 
                                    hit_triangle_global_idx = -1; 
                                    hit_point_coords = []; 
                                    hit_barycentric_coords = [];
                                    
                                    for tri_g_idx = 1:size(T_scene_idx,1) 
                                        vertex_indices_global = T_scene_idx(tri_g_idx,:);
                                        p0 = V_scene_combined(vertex_indices_global(1),:);
                                        p1 = V_scene_combined(vertex_indices_global(2),:);
                                        p2 = V_scene_combined(vertex_indices_global(3),:);
                                        
                                        edge1 = p1-p0; 
                                        edge2 = p2-p0;
                                        h_vec = cross(current_ray_direction, edge2);
                                        a_det = dot(edge1, h_vec);
                                        
                                        if(abs(a_det) < EPSILON_SHIFT), continue; end 
                                        
                                        f_inv_det = 1.0/a_det;
                                        s_vec = current_ray_origin - p0;
                                        u_bary = f_inv_det * dot(s_vec, h_vec);
                                        
                                        if(u_bary < 0.0 || u_bary > 1.0), continue; end 
                                        
                                        q_vec = cross(s_vec, edge1);
                                        v_bary = f_inv_det * dot(current_ray_direction, q_vec);
                                        
                                        if(v_bary < 0.0 || u_bary + v_bary > 1.0), continue; end 
                                        
                                        t_intersection = f_inv_det * dot(edge2, q_vec); 
                                        
                                        if(t_intersection > EPSILON_SHIFT) 
                                            if t_intersection < min_t_intersection
                                                min_t_intersection = t_intersection;
                                                hit_found = true;
                                                hit_triangle_global_idx = tri_g_idx;
                                                hit_point_coords = current_ray_origin + t_intersection * current_ray_direction;
                                                hit_barycentric_coords = [u_bary, v_bary]; 
                                            end
                                        end
                                    end 
                                    
                                    if ~hit_found, break; end 
                                    
                                    hit_material_id = M_scene_materials(hit_triangle_global_idx);
                                    hit_vertex_global_indices = T_scene_idx(hit_triangle_global_idx,:);
                                    
                                    p0_hit = V_scene_combined(hit_vertex_global_indices(1),:);
                                    p1_hit = V_scene_combined(hit_vertex_global_indices(2),:);
                                    p2_hit = V_scene_combined(hit_vertex_global_indices(3),:);
                                    face_normal = cross(p1_hit-p0_hit, p2_hit-p0_hit);
                                    face_normal = face_normal / (norm(face_normal) + eps); 
                                    
                                    if hit_material_id == MATERIAL_SKIN
                                        effective_albedo_rgb = skin_base_mc_color .* curr_diff_coeff; 
                                        effective_albedo_rgb = min(max(effective_albedo_rgb,0),1); 
                                        
                                        reflected_energy_rgb = current_ray_energy_rgb .* effective_albedo_rgb; 
                                        
                                        u_bc = hit_barycentric_coords(1); 
                                        v_bc = hit_barycentric_coords(2);
                                        w_bc = 1 - u_bc - v_bc; 
                                        
                                        AccumulatedEnergy_SkinVertex(hit_vertex_global_indices(1),:) = AccumulatedEnergy_SkinVertex(hit_vertex_global_indices(1),:) + reflected_energy_rgb * w_bc;
                                        AccumulatedEnergy_SkinVertex(hit_vertex_global_indices(2),:) = AccumulatedEnergy_SkinVertex(hit_vertex_global_indices(2),:) + reflected_energy_rgb * u_bc; 
                                        AccumulatedEnergy_SkinVertex(hit_vertex_global_indices(3),:) = AccumulatedEnergy_SkinVertex(hit_vertex_global_indices(3),:) + reflected_energy_rgb * v_bc; 
                                        
                                        break; 
                                        
                                    elseif hit_material_id == MATERIAL_LED_OPAQUE_SIDE
                                        break; 
                                        
                                    elseif hit_material_id == MATERIAL_LED_FILTER_NON_ATTENUATING
                                        current_ray_origin = hit_point_coords + current_ray_direction * EPSILON_SHIFT;
                                        
                                    elseif hit_material_id == MATERIAL_LED_FILTER_ATTENUATING
                                        current_ray_energy_rgb = current_ray_energy_rgb * TOP_FILTER_ATTENUATION_COEFFICIENT;
                                        current_ray_origin = hit_point_coords + current_ray_direction * EPSILON_SHIFT;
                                    else
                                        warning('MC_Trace:UnknownMaterial', 'Ray hit unknown material ID: %d. Stopping ray.', hit_material_id);
                                        break; 
                                    end
                                end 
                                if getappdata(h_main_prog,'canceling_all') || getappdata(h_mc_inner,'canceling_mc_inner'), break; end 
                            end 
                            if getappdata(h_main_prog,'canceling_all') || getappdata(h_mc_inner,'canceling_mc_inner'), break; end 
                        end 
                    catch ME_mc_inner_run
                        if ishandle(h_mc_inner), delete(h_mc_inner); end 
                        rethrow(ME_mc_inner_run); 
                    end
                    
                    if ishandle(h_mc_inner), delete(h_mc_inner); end 
                    if getappdata(h_main_prog,'canceling_all'), disp('全组合处理被用户取消 (MC内部)。'); break; end 
                end 
                
                ExportData = zeros(size(g_skin_vertices,1), 4); 
                max_observed_intensity_skin = 0;
                
                if any(AccumulatedEnergy_SkinVertex(:) > 0)
                    scalar_intensity_skin = sum(AccumulatedEnergy_SkinVertex, 2); 
                    max_observed_intensity_skin = max(scalar_intensity_skin);
                end
                
                for i_sv = 1:size(g_skin_vertices,1) 
                    intensity_scalar = sum(AccumulatedEnergy_SkinVertex(i_sv,:));
                    intensity_scalar_out = 0; 
                    if max_observed_intensity_skin > 1e-9 
                        intensity_scalar_out = 3 * (intensity_scalar / max_observed_intensity_skin);
                        intensity_scalar_out = min(3, max(0, intensity_scalar_out)); 
                    end
                    ExportData(i_sv,:) = [g_skin_vertices(i_sv,:), intensity_scalar_out];
                end
                
                try
                    writematrix(ExportData, combo_filepath); 
                    disp(['数据已保存到: ' combo_filepath]);
                catch ME_save
                    errordlg(['保存组合 ' combo_filename ' 数据出错: ' ME_save.message],'文件保存错误');
                end
                
                if getappdata(h_main_prog,'canceling_all'), break; end 
            end 
            if getappdata(h_main_prog,'canceling_all'), break; end 
        end 
        
        if ishandle(h_main_prog), delete(h_main_prog); end 
        
        if getappdata(fig,'canceling_all') 
            disp('全组合与多角度处理被用户取消。');
            rmappdata(fig,'canceling_all'); 
        else
            helpdlg(sprintf('所有 %d 个角度-组合的光强数据已处理并尝试保存。',total_runs),'批量处理完成');
        end
    end

    function export_light_intensity_callback(~, ~) 
        if isempty(g_skin_vertices) || isempty(g_skin_triangle_indices)
            warndlg('模拟人脸数据未加载或无效。请先点击 "加载眼镜与人脸"。','错误：无皮肤数据'); return;
        end
        if isempty(points_data) 
            warndlg('没有设置光源。请先通过 "输入光源模块中心" 或 "测试光源模块" 设置光源。','错误：无光源'); return;
        end
        if ~isgraphics(h_surf_skin) 
            warndlg('模拟人脸预览曲面无效。请重新加载。','错误：皮肤预览无效'); return;
        end
        
        if ispc, default_desktop_path = fullfile(getenv('USERPROFILE'),'Desktop'); else, default_desktop_path = fullfile(getenv('HOME'),'Desktop'); end
        [manual_fname, manual_pname] = uiputfile({'*.txt';'*.csv'}, ...
                                                 '保存当前光强数据 (手动MC)', ...
                                                 fullfile(default_desktop_path,'face_vertex_mc_manual_output.txt')); 
        if isequal(manual_fname,0) || isequal(manual_pname,0) 
            warndlg('手动导出操作已取消。','提示'); return;
        end
        manual_filepath = fullfile(manual_pname, manual_fname);
        
        disp('开始蒙特卡洛计算 (手动导出)...');
        h_prog_manual_mc = waitbar(0,'蒙特卡洛计算中 (手动)...','Name','手动光强计算','CreateCancelBtn','setappdata(gcbf,"canceling_manual_mc",1)');
        setappdata(h_prog_manual_mc,'canceling_manual_mc',0); 
        
        V_scene_comb_manual = [g_skin_vertices; g_led_boxes_vertices];
        num_skin_T_manual = size(g_skin_triangle_indices,1);
        led_T_offset_manual = size(g_skin_vertices,1);
        
        adj_led_T_idx_manual = [];
        if ~isempty(g_led_boxes_triangle_indices)
            adj_led_T_idx_manual = g_led_boxes_triangle_indices + led_T_offset_manual;
        end
        T_scene_comb_idx_manual = [g_skin_triangle_indices; adj_led_T_idx_manual];
        M_scene_comb_mat_manual = [repmat(MATERIAL_SKIN, num_skin_T_manual,1); g_led_boxes_triangle_materials];
        
        curr_diff_coeff_slider_val = get(h_diffuse_slider,'Value'); 
        skin_base_mc_color_manual = get(h_surf_skin,'FaceColor'); 
        
        LightPos_MC_manual = points_data; 
        BaseLightPwr_manual = repmat(current_light_color, size(LightPos_MC_manual,1), 1);
        Accum_E_SkinVert_manual = zeros(size(g_skin_vertices,1), 3); 
        num_tot_eff_lights_manual = size(LightPos_MC_manual,1);
        
        if isempty(LightPos_MC_manual) 
            if ishandle(h_prog_manual_mc), delete(h_prog_manual_mc); end
            warndlg('手动导出: 无光源激活。将生成零强度文件。','提示');
        end

        try 
            if ~isempty(LightPos_MC_manual)
                for j_eff_s = 1:num_tot_eff_lights_manual 
                    if getappdata(h_prog_manual_mc,'canceling_manual_mc'), break; end 
                    
                    waitbar(j_eff_s/num_tot_eff_lights_manual, h_prog_manual_mc, ...
                            sprintf('手动MC: 光源 %d/%d...',j_eff_s,num_tot_eff_lights_manual));
                    
                    light_start_s = LightPos_MC_manual(j_eff_s,:); 
                    init_E_ray_s = BaseLightPwr_manual(j_eff_s,:)/num_rays_per_light_source_param;
                    
                    for r_idx_s = 1:num_rays_per_light_source_param 
                        if getappdata(h_prog_manual_mc,'canceling_manual_mc'), break; end 
                        
                        curr_O_s = light_start_s; 
                        phi_s=2*pi*rand(); cost_s=2*rand()-1; sint_s=sqrt(max(0,1-cost_s^2));
                        curr_D_s=[sint_s*cos(phi_s),sint_s*sin(phi_s),cost_s]; 
                        curr_E_s = init_E_ray_s;
                        
                        for bounce_s = 1:MAX_RAY_BOUNCES 
                            if sum(curr_E_s(:)) < RAY_ENERGY_THRESHOLD, break; end
                            
                            min_t_s=inf; hit_f_s=false; hit_tri_idx_s=-1; hit_pt_s=[]; hit_bary_s=[];
                            for tri_g_s = 1:size(T_scene_comb_idx_manual,1) 
                                v_g_s=T_scene_comb_idx_manual(tri_g_s,:);
                                p0s=V_scene_comb_manual(v_g_s(1),:); p1s=V_scene_comb_manual(v_g_s(2),:); p2s=V_scene_comb_manual(v_g_s(3),:);
                                e1s=p1s-p0s; e2s=p2s-p0s; hs=cross(curr_D_s,e2s); as=dot(e1s,hs);
                                if(abs(as)<EPSILON_SHIFT), continue; end
                                fs=1.0/as; ss=curr_O_s-p0s; us=fs*dot(ss,hs);
                                if(us<0.0||us>1.0), continue; end
                                qs=cross(ss,e1s); vs=fs*dot(curr_D_s,qs);
                                if(vs<0.0||us+vs>1.0), continue; end
                                ts=fs*dot(e2s,qs);
                                if(ts>EPSILON_SHIFT)
                                    if ts<min_t_s
                                        min_t_s=ts; hit_f_s=true; hit_tri_idx_s=tri_g_s; 
                                        hit_pt_s=curr_O_s+ts*curr_D_s; hit_bary_s=[us,vs];
                                    end
                                end
                            end 
                            
                            if ~hit_f_s, break; end 
                            
                            hit_m_s = M_scene_comb_mat_manual(hit_tri_idx_s); 
                            v_hit_s_global_indices = T_scene_comb_idx_manual(hit_tri_idx_s,:);
                            
                            if hit_m_s == MATERIAL_SKIN
                                eff_as = skin_base_mc_color_manual .* curr_diff_coeff_slider_val;
                                eff_as = min(max(eff_as,0),1);
                                refl_es = curr_E_s .* eff_as; 
                                bus=hit_bary_s(1); bvs=hit_bary_s(2); bws=1-bus-bvs;
                                
                                Accum_E_SkinVert_manual(v_hit_s_global_indices(1),:) = Accum_E_SkinVert_manual(v_hit_s_global_indices(1),:) + refl_es*bws;
                                Accum_E_SkinVert_manual(v_hit_s_global_indices(2),:) = Accum_E_SkinVert_manual(v_hit_s_global_indices(2),:) + refl_es*bus;
                                Accum_E_SkinVert_manual(v_hit_s_global_indices(3),:) = Accum_E_SkinVert_manual(v_hit_s_global_indices(3),:) + refl_es*bvs;
                                break; 
                            elseif hit_m_s == MATERIAL_LED_OPAQUE_SIDE
                                break; 
                            elseif hit_m_s == MATERIAL_LED_FILTER_NON_ATTENUATING
                                curr_O_s = hit_pt_s + curr_D_s*EPSILON_SHIFT;
                            elseif hit_m_s == MATERIAL_LED_FILTER_ATTENUATING
                                curr_E_s = curr_E_s * TOP_FILTER_ATTENUATION_COEFFICIENT;
                                curr_O_s = hit_pt_s + curr_D_s*EPSILON_SHIFT;
                            else
                                warning('MC_Manual:UnknownMaterial', '手动MC中光线击中未知材质 ID: %d。停止光线。', hit_m_s);
                                break; 
                            end
                        end 
                        if getappdata(h_prog_manual_mc,'canceling_manual_mc'), break; end
                    end 
                    if getappdata(h_prog_manual_mc,'canceling_manual_mc'), break; end
                end 
            end 
        catch ME_mc_manual_run
            if ishandle(h_prog_manual_mc), delete(h_prog_manual_mc); end
            rethrow(ME_mc_manual_run);
        end
        
        if ishandle(h_prog_manual_mc), delete(h_prog_manual_mc); end 
        if getappdata(fig,'canceling_manual_mc') 
            disp('手动蒙特卡洛计算被用户取消。');
            rmappdata(fig,'canceling_manual_mc'); 
            return;
        end
        
        ExportData_s = zeros(size(g_skin_vertices,1),4);
        max_obs_s = 0;
        if any(Accum_E_SkinVert_manual(:) > 0)
            sc_ints_s = sum(Accum_E_SkinVert_manual,2);
            max_obs_s = max(sc_ints_s);
        end
        
        for i_s = 1:size(g_skin_vertices,1)
            int_sc_s = sum(Accum_E_SkinVert_manual(i_s,:));
            int_sc_out_s = 0;
            if max_obs_s > 1e-9
                int_sc_out_s = 3 * (int_sc_s / max_obs_s);
                int_sc_out_s = min(3, max(0, int_sc_out_s));
            end
            ExportData_s(i_s,:) = [g_skin_vertices(i_s,:), int_sc_out_s];
        end
        
        try
            writematrix(ExportData_s, manual_filepath);
            helpdlg(sprintf('手动导出的光强数据已保存到:\n%s',manual_filepath),'手动保存成功');
        catch ME_manual_save
            errordlg(['手动保存数据时出错: ' ME_manual_save.message],'手动保存错误');
        end
        disp('手动蒙特卡洛计算与导出完成。');
    end

    function test_point_callback(~, ~)
        if ~isvalid(ax) 
            errordlg('绘图区域无效。', '内部错误'); return;
        end

        test_cuboid_center_coords = [38,-82,23]; 
        set_lights_and_markers(test_cuboid_center_coords); 
        
        is_dup_t = false; 
        tol_t = 1e-6;
        for hist_idx_t = 1:length(all_input_points_history)
             if isequal(size(all_input_points_history{hist_idx_t}), size(test_cuboid_center_coords)) && ...
                all(abs(all_input_points_history{hist_idx_t} - test_cuboid_center_coords) < tol_t)
                is_dup_t = true; 
                break;
             end
        end
        if ~is_dup_t
            all_input_points_history{end+1} = test_cuboid_center_coords;
        end
        
        helpdlg(sprintf('测试旋转光源模块中心 (%g,%g,%g) 已设置。',test_cuboid_center_coords),'手动测试点设置');
    end

    function export_all_history_callback(~, ~)
        if isempty(all_input_points_history)
            warndlg('无历史手动输入模块中心点可导出。','提示'); return;
        end 
        
        if ispc, desk_p_h = fullfile(getenv('USERPROFILE'),'Desktop'); else, desk_p_h = fullfile(getenv('HOME'),'Desktop'); end
        [fname_h,pname_h] = uiputfile({'*.txt';'*.csv'}, ...
                                      '保存历史手动模块中心坐标', ...
                                      fullfile(desk_p_h,'historical_manual_module_centers.txt')); 
        if isequal(fname_h,0) || isequal(pname_h,0), return; end 
        
        fpath_h = fullfile(pname_h,fname_h); 
        fid_h = fopen(fpath_h,'w');
        if fid_h == -1
            errordlg('无法打开文件进行写入。请检查权限或路径。','文件写入错误'); return;
        end
        
        try
            fprintf(fid_h,'# Manually Input Rotated Light Module Center Coordinates (X,Y,Z)\n'); 
            for k_h_e = 1:length(all_input_points_history)
                if ~isempty(all_input_points_history{k_h_e})
                    fprintf(fid_h,'%g,%g,%g\n',all_input_points_history{k_h_e});
                end
            end
            fclose(fid_h); 
            helpdlg(sprintf('历史手动模块中心坐标已保存到:\n%s',fpath_h),'保存成功');
        catch ME_h
            if fid_h ~= -1, fclose(fid_h); end 
            errordlg(['保存历史坐标时出错: ' ME_h.message],'保存错误');
        end
    end
     
    function set_lights_and_markers(center_coords_matrix) 
        if ~isvalid(ax), return; end 

        if ~isempty(h_custom_lights) && all(isgraphics(h_custom_lights(isvalid(h_custom_lights))))
            delete(h_custom_lights(isgraphics(h_custom_lights))); 
        end
        h_custom_lights=[];
        
        if ~isempty(h_all_points_markers) && all(isgraphics(h_all_points_markers(isvalid(h_all_points_markers))))
            delete(h_all_points_markers(isgraphics(h_all_points_markers)));
        end
        h_all_points_markers=[];
        
        if ~isempty(h_rectangle_plots) && all(isgraphics(h_rectangle_plots(isvalid(h_rectangle_plots))))
            delete(h_rectangle_plots(isgraphics(h_rectangle_plots)));
        end
        h_rectangle_plots=[];
        
        points_data=[]; 
        g_led_boxes_vertices=[]; 
        g_led_boxes_triangle_indices=[]; 
        g_led_boxes_triangle_materials=[];
        
        axes(ax); 
        hold on;  

        rect_dims = [1, 0.4, 1]; 
        Hx = rect_dims(1)/2; Hy = rect_dims(2)/2; Hz = rect_dims(3)/2; 
        
        light_separation_within_module = 0.6; 
        light_offset_from_center_x = light_separation_within_module/2; 
        
        rot_ang_rad = ROTATION_ANGLE_DEGREES * pi/180; 
        R_x_mat = [1, 0, 0; 
                   0, cos(rot_ang_rad), -sin(rot_ang_rad); 
                   0, sin(rot_ang_rad), cos(rot_ang_rad)];
        
        vertex_offset_for_global_indices = 0; 
        
        if ~isempty(center_coords_matrix)
            for i_c = 1:size(center_coords_matrix,1) 
                center_crd = center_coords_matrix(i_c,:); 
                Cx_c = center_crd(1); Cy_c = center_crd(2); Cz_c = center_crd(3);
                
                v_box_initial_local = [ 
                    -Hx, -Hy, -Hz; +Hx, -Hy, -Hz; +Hx, +Hy, -Hz; -Hx, +Hy, -Hz; 
                    -Hx, -Hy, +Hz; +Hx, -Hy, +Hz; +Hx, +Hy, +Hz; -Hx, +Hy, +Hz  
                ];
                
                v_box_centered_at_origin_for_rotation = v_box_initial_local; 
                
                v_rotated_local = (R_x_mat * v_box_centered_at_origin_for_rotation')';
                
                v_box_final_world = v_rotated_local + repmat(center_crd, 8, 1);
                
                g_led_boxes_vertices = [g_led_boxes_vertices; v_box_final_world]; 
                
                local_tris_std = [ 
                    1,2,6; 1,6,5; 
                    4,8,7; 4,7,3; 
                    1,5,8; 1,8,4; 
                    2,3,7; 2,7,6; 
                    1,4,3; 1,3,2; 
                    5,6,7; 5,7,8  
                ];
                tris_mats_box_std = [ 
                    MATERIAL_LED_FILTER_NON_ATTENUATING; MATERIAL_LED_FILTER_NON_ATTENUATING;
                    MATERIAL_LED_FILTER_ATTENUATING; MATERIAL_LED_FILTER_ATTENUATING;
                    MATERIAL_LED_OPAQUE_SIDE; MATERIAL_LED_OPAQUE_SIDE;
                    MATERIAL_LED_OPAQUE_SIDE; MATERIAL_LED_OPAQUE_SIDE;
                    MATERIAL_LED_OPAQUE_SIDE; MATERIAL_LED_OPAQUE_SIDE;
                    MATERIAL_LED_OPAQUE_SIDE; MATERIAL_LED_OPAQUE_SIDE
                ];

                current_box_tris_global_idx = local_tris_std + vertex_offset_for_global_indices;
                g_led_boxes_triangle_indices = [g_led_boxes_triangle_indices; current_box_tris_global_idx];
                g_led_boxes_triangle_materials = [g_led_boxes_triangle_materials; tris_mats_box_std];
                
                vertex_offset_for_global_indices = vertex_offset_for_global_indices + 8; 
                
                edges_vis_idx = [1,2;2,3;3,4;4,1; 5,6;6,7;7,8;8,5; 1,5;2,6;3,7;4,8]; 
                for i_e = 1:size(edges_vis_idx,1)
                    v1_idx = edges_vis_idx(i_e,1); v2_idx = edges_vis_idx(i_e,2);
                    h_r_edge = plot3(ax, [v_box_final_world(v1_idx,1), v_box_final_world(v2_idx,1)], ...
                                          [v_box_final_world(v1_idx,2), v_box_final_world(v2_idx,2)], ...
                                          [v_box_final_world(v1_idx,3), v_box_final_world(v2_idx,3)], ...
                                          'm-', 'LineWidth', 0.5); 
                    h_rectangle_plots = [h_rectangle_plots, h_r_edge];
                end
                
                light_p1_local_unrotated = [-light_offset_from_center_x, -Hy, 0]; 
                light_p2_local_unrotated = [+light_offset_from_center_x, -Hy, 0]; 

                light_p1_local_rotated = (R_x_mat * light_p1_local_unrotated')';
                light_p2_local_rotated = (R_x_mat * light_p2_local_unrotated')';

                light_p1_world = light_p1_local_rotated + center_crd;
                light_p2_world = light_p2_local_rotated + center_crd;
                
                points_data = [points_data; light_p1_world; light_p2_world]; 
                
                for eff_l_p_cell = {light_p1_world, light_p2_world}
                    pos_for_light_obj = eff_l_p_cell{:};
                    if isvalid(ax) 
                        h_l_prev = light('Parent',ax,'Position',pos_for_light_obj,'Style','local','Color',current_light_color); 
                        h_custom_lights = [h_custom_lights, h_l_prev];
                        
                        h_m_prev = plot3(ax, pos_for_light_obj(1), pos_for_light_obj(2), pos_for_light_obj(3), ...
                                         'o','MarkerFaceColor','r','MarkerEdgeColor','k','MarkerSize',5); 
                        h_all_points_markers = [h_all_points_markers, h_m_prev];
                    end
                end
            end 
        end 
        
        hold off; 
        drawnow;  
    end

    function update_surface_lighting() 
        if ~isempty(h_surf_skin) && isgraphics(h_surf_skin) 
            diff_val = get(h_diffuse_slider,'Value');
            set(h_surf_skin,'DiffuseStrength',max(0,min(1,diff_val))); 
            
            if isempty(h_custom_lights) 
                if isempty(findobj(ax,'Type','light')) 
                    camlight(ax,'headlight'); 
                end
            end
            drawnow; 
        end
    end

    function diffuse_intensity_slider_callback(~,~)
        val = get(h_diffuse_slider,'Value');
        set(h_diffuse_value_display,'String',sprintf('%.2f',val));
        update_surface_lighting();
    end

    load_models_callback(fig, []); 

end
