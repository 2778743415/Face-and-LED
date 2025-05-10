function custom_thermal_simulation_cli() 

    disp('=======================================================================');
    disp('  光子路线与热传导模型 (命令行版本)');
    disp('  模拟LED在金属眼镜上的热量传递到人脸');
    disp('=======================================================================');

    original_path_at_start = path;
    original_dir_at_start = pwd;

    try
       
        glasses_stl_path_default = 'C:\Users\WC\Desktop\new\glasses.stl'; 

        glasses_stl_path = input(sprintf('请输入眼镜STL文件路径 [%s]: ', glasses_stl_path_default), 's');
        if isempty(glasses_stl_path)
            glasses_stl_path = glasses_stl_path_default;
        end
        if ~exist(glasses_stl_path, 'file')
            error('眼镜STL文件未找到: "%s"', glasses_stl_path);
        end

             disp('--- 定义模拟参数 (单位使用 cm, g, s, W, K/°C) ---');
        
        total_photons_to_launch = str2double(input('输入总的光子发射数量 [1e6]: ', 's'));
        if isnan(total_photons_to_launch) || total_photons_to_launch <= 0, total_photons_to_launch = 1e6; end
        
        optical_source_power_W = str2double(input('总的光学光源功率 (W) [0.1]: ', 's'));
        if isnan(optical_source_power_W), optical_source_power_W = 0.1; end
        
        ambient_temp_C = 25;         
        sim_duration_s = str2double(input('热模拟总时长 (s) [60]: ', 's'));
        if isnan(sim_duration_s), sim_duration_s = 60; end
        num_time_steps_output = 20; 

        disp('-- 空气光学属性 --');
        air_mua_1_cm = str2double(input('空气 mua (1/cm) [0.0001]: ','s')); if isnan(air_mua_1_cm), air_mua_1_cm = 0.0001; end
        air_mus_1_cm = str2double(input('空气 mus (1/cm) [0.0001]: ','s')); if isnan(air_mus_1_cm), air_mus_1_cm = 0.0001; end
        air_g = str2double(input('空气 g [1.0]: ','s')); if isnan(air_g), air_g = 1.0; end
        air_n = str2double(input('空气 n [1.0]: ','s')); if isnan(air_n), air_n = 1.0; end

        disp('-- 皮肤光学属性 --');
        skin_mua_1_cm = str2double(input('皮肤 mua (1/cm) [0.2]: ','s')); if isnan(skin_mua_1_cm), skin_mua_1_cm = 0.2; end
        skin_mus_1_cm = str2double(input('皮肤 mus (1/cm) [10.0]: ','s')); if isnan(skin_mus_1_cm), skin_mus_1_cm = 10.0; end
        skin_g = str2double(input('皮肤 g [0.85]: ','s')); if isnan(skin_g), skin_g = 0.85; end
        skin_n = str2double(input('皮肤 n [1.4]: ','s')); if isnan(skin_n), skin_n = 1.4; end

        disp('-- 眼镜金属光学属性 --');
        glasses_mua_1_cm = str2double(input('眼镜 mua (1/cm) [1000]: ','s')); if isnan(glasses_mua_1_cm), glasses_mua_1_cm = 1000; end 
        glasses_mus_1_cm = str2double(input('眼镜 mus (1/cm) [100]: ','s')); if isnan(glasses_mus_1_cm), glasses_mus_1_cm = 100; end 
        glasses_g = str2double(input('眼镜 g [0.5]: ','s')); if isnan(glasses_g), glasses_g = 0.5; end
        glasses_n = str2double(input('眼镜 n [2.5]: ','s')); if isnan(glasses_n), glasses_n = 2.5; end 

        disp('-- LED模块材料光学属性 --');
        led_mua_1_cm = str2double(input('LED mua (1/cm) [5.0]: ','s')); if isnan(led_mua_1_cm), led_mua_1_cm = 5.0; end 
        led_mus_1_cm = str2double(input('LED mus (1/cm) [5.0]: ','s')); if isnan(led_mus_1_cm), led_mus_1_cm = 5.0; end
        led_g = str2double(input('LED g [0.7]: ','s')); if isnan(led_g), led_g = 0.7; end
        led_n = str2double(input('LED n [1.5]: ','s')); if isnan(led_n), led_n = 1.5; end

        optical_props(1) = struct('mua', air_mua_1_cm, 'mus', air_mus_1_cm, 'g', air_g, 'n', air_n, 'name', 'Air'); 
        optical_props(2) = struct('mua', skin_mua_1_cm, 'mus', skin_mus_1_cm, 'g', skin_g, 'n', skin_n, 'name', 'Skin'); 
        optical_props(3) = struct('mua', glasses_mua_1_cm, 'mus', glasses_mus_1_cm, 'g', glasses_g, 'n', glasses_n, 'name', 'Glasses'); 
        optical_props(4) = struct('mua', led_mua_1_cm, 'mus', led_mus_1_cm, 'g', led_g, 'n', led_n, 'name', 'LED'); 
        disp('光学属性结构体已创建。');

        disp('-- 皮肤热物性 --');
        skin_TC_WmK = str2double(input('皮肤热导率 (W/mK) [0.5]: ', 's')); if isnan(skin_TC_WmK), skin_TC_WmK = 0.5; end
        skin_Rho_kgm3 = str2double(input('皮肤密度 (kg/m^3) [1050]: ', 's')); if isnan(skin_Rho_kgm3), skin_Rho_kgm3 = 1050; end
        skin_Cp_JkgK = str2double(input('皮肤比热容 (J/kgK) [3500]: ', 's')); if isnan(skin_Cp_JkgK), skin_Cp_JkgK = 3500; end

        disp('-- 眼镜金属热物性 --');
        glasses_TC_WmK = str2double(input('眼镜热导率 (W/mK) [200]: ', 's')); if isnan(glasses_TC_WmK), glasses_TC_WmK = 200; end
        glasses_Rho_kgm3 = str2double(input('眼镜密度 (kg/m^3) [2700]: ', 's')); if isnan(glasses_Rho_kgm3), glasses_Rho_kgm3 = 2700; end
        glasses_Cp_JkgK = str2double(input('眼镜比热容 (J/kgK) [900]: ', 's')); if isnan(glasses_Cp_JkgK), glasses_Cp_JkgK = 900; end

        disp('-- LED模块材料热物性 (简化) --');
        led_TC_WmK = str2double(input('LED热导率 (W/mK) [1.0]: ', 's')); if isnan(led_TC_WmK), led_TC_WmK = 1.0; end
        led_Rho_kgm3 = str2double(input('LED密度 (kg/m^3) [1200]: ', 's')); if isnan(led_Rho_kgm3), led_Rho_kgm3 = 1200; end
        led_Cp_JkgK = str2double(input('LED比热容 (J/kgK) [1000]: ', 's')); if isnan(led_Cp_JkgK), led_Cp_JkgK = 1000; end

        air_TC_WmK = 0.026; air_Rho_kgm3 = 1.2; air_Cp_JkgK = 1005;

        convection_h_Wm2K = str2double(input('对流换热系数 (空气 W/m^2K) [10]: ', 's'));
        if isnan(convection_h_Wm2K), convection_h_Wm2K = 10; end

        voxel_size_mm = str2double(input('体素大小 (mm) [5.0]: ', 's')); 
        if isnan(voxel_size_mm) || voxel_size_mm <= 0, voxel_size_mm = 5.0; end
        voxel_size_cm = voxel_size_mm / 10;
        fprintf('!!! 注意: 使用的体素大小为 %.2f mm (%.3f cm).\n', voxel_size_mm, voxel_size_cm);

        active_module_centers_mm = [-38, -84, 23; 
                                     38, -84, 23]; 
        rect_dims_led_mm = [1, 0.4, 1];      
        disp('自动生成LED模块，中心坐标 (mm):'); 
        disp(active_module_centers_mm);
        disp(['每个LED模块尺寸 (mm): ', mat2str(rect_dims_led_mm)]);
        
        skin_TC_WcmK    = skin_TC_WmK / 100; skin_Rho_gcm3   = skin_Rho_kgm3 / 1000; skin_Cp_JgK     = skin_Cp_JkgK / 1000;
        skin_VHC_Jcm3K  = skin_Rho_gcm3 * skin_Cp_JgK; skin_alpha_cm2s = skin_TC_WcmK / skin_VHC_Jcm3K; 
        glasses_TC_WcmK   = glasses_TC_WmK / 100; glasses_Rho_gcm3  = glasses_Rho_kgm3 / 1000; glasses_Cp_JgK    = glasses_Cp_JkgK / 1000;
        glasses_VHC_Jcm3K = glasses_Rho_gcm3 * glasses_Cp_JgK; glasses_alpha_cm2s= glasses_TC_WcmK / glasses_VHC_Jcm3K;
        led_TC_WcmK     = led_TC_WmK / 100; led_Rho_gcm3    = led_Rho_kgm3 / 1000; led_Cp_JgK      = led_Cp_JkgK / 1000;
        led_VHC_Jcm3K   = led_Rho_gcm3 * led_Cp_JgK; led_alpha_cm2s  = led_TC_WcmK / led_VHC_Jcm3K;
        air_TC_WcmK     = air_TC_WmK / 100; air_Rho_gcm3    = air_Rho_kgm3 / 1000; air_Cp_JgK      = air_Cp_JkgK / 1000;
        air_VHC_Jcm3K   = air_Rho_gcm3 * air_Cp_JgK; air_alpha_cm2s  = air_TC_WcmK / air_VHC_Jcm3K;
        convection_h_Wcm2K = convection_h_Wm2K / (100*100);
        
        disp('--- 正在加载模型和定义几何形状 (转换为cm) ---');
        model_glasses_stl = stlread(glasses_stl_path);
        g_glasses_vertices_cm = model_glasses_stl.Points / 10; 
        g_glasses_faces = model_glasses_stl.ConnectivityList;
        disp(['成功加载眼镜模型并转换为cm: ' glasses_stl_path]);

        v_center_face_cm = [0,-8,0]; 
        x_r_cm = linspace(-7,7,50); z_r_cm = linspace(-6,6,50);      
        [Xf_cm, Zf_cm] = meshgrid(x_r_cm, z_r_cm); a_el=0.08; b_el=0.04; 
        Yf_cm = v_center_face_cm(2) + (a_el*(Xf_cm-v_center_face_cm(1)).^2 + b_el*(Zf_cm-v_center_face_cm(3)).^2)-1; 
        
        face_surf_X = Xf_cm; face_surf_Y = Yf_cm; face_surf_Z = Zf_cm;

        temp_fig_skin = figure('Visible','off'); temp_ax_skin = axes(temp_fig_skin);
        h_surf_skin_temp = surf(temp_ax_skin, face_surf_X, face_surf_Y, face_surf_Z); 
        fv_skin = surf2patch(h_surf_skin_temp, 'triangles');
        g_skin_vertices_cm = fv_skin.vertices; 
        g_skin_triangle_indices = fv_skin.faces;
        close(temp_fig_skin); 
        if isempty(g_skin_vertices_cm) || isempty(g_skin_triangle_indices), error('未能提取人脸几何数据。'); end
        disp('模拟人脸曲面已生成 (cm)。');

        active_module_centers_cm = active_module_centers_mm / 10;
        rect_dims_led_cm = rect_dims_led_mm / 10;
        g_all_led_box_geometries_cm = {};
        
        points_data_mm_temp = []; 
        
        ROTATION_ANGLE_DEGREES = -30; 
        rot_ang_rad = ROTATION_ANGLE_DEGREES * pi/180; 
        R_x_mat = [1,0,0;0,cos(rot_ang_rad),-sin(rot_ang_rad);0,sin(rot_ang_rad),cos(rot_ang_rad)];

        if ~isempty(active_module_centers_cm) 
            for i_c = 1:size(active_module_centers_cm,1) 
                center_crd_cm = active_module_centers_cm(i_c,:); 
                Hx_cm = rect_dims_led_cm(1)/2; Hy_cm = rect_dims_led_cm(2)/2; Hz_cm = rect_dims_led_cm(3)/2;
                v_box_initial_local_cm = [-Hx_cm,-Hy_cm,-Hz_cm; +Hx_cm,-Hy_cm,-Hz_cm; +Hx_cm,+Hy_cm,-Hz_cm; -Hx_cm,+Hy_cm,-Hz_cm; ...
                                       -Hx_cm,-Hy_cm,+Hz_cm; +Hx_cm,-Hy_cm,+Hz_cm; +Hx_cm,+Hy_cm,+Hz_cm; -Hx_cm,+Hy_cm,+Hz_cm];
                
                v_rotated_local_cm = (R_x_mat * v_box_initial_local_cm')'; 
                v_box_final_world_current_led_cm = v_rotated_local_cm + repmat(center_crd_cm,8,1); 

                local_tris_std_faces = [1,2,6;1,6,5;  4,8,7;4,7,3;  1,5,8;1,8,4;  2,3,7;2,7,6;  1,4,3;1,3,2;  5,6,7;5,7,8]; 
                g_all_led_box_geometries_cm{i_c}.vertices = v_box_final_world_current_led_cm; 
                g_all_led_box_geometries_cm{i_c}.faces = local_tris_std_faces; 
                
                light_point_local_cm = [0, -Hy_cm, 0]; 
                light_point_rotated_local_cm = (R_x_mat * light_point_local_cm')';
                light_point_world_cm = light_point_rotated_local_cm + center_crd_cm;
                points_data_mm_temp = [points_data_mm_temp; light_point_world_cm * 10]; 
            end
        end
        points_data_cm = points_data_mm_temp / 10; 
        disp(['已定义 ' num2str(size(active_module_centers_cm,1)) ' 个LED模块 (cm)。']);


        disp('--- 可视化导入的STL和生成的人脸曲面 ---');
        h_geom_fig = figure('Name', '几何模型预览 (眼镜与人脸)');
        ax_geom = axes(h_geom_fig);
        hold(ax_geom, 'on');
        if ~isempty(g_glasses_vertices_cm)
            trisurf(g_glasses_faces, ...
                    g_glasses_vertices_cm(:,1), g_glasses_vertices_cm(:,2), g_glasses_vertices_cm(:,3), ...
                    'Parent', ax_geom, 'FaceColor', [0.6 0.6 0.7], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
        end
        if ~isempty(face_surf_X) 
             surf(ax_geom, face_surf_X, face_surf_Y, face_surf_Z, ...
                 'FaceColor', [1,0.8,0.6], 'EdgeColor', 'none', 'FaceAlpha', 0.9);
        end
        if ~isempty(g_all_led_box_geometries_cm)
            for i_led_vis = 1:length(g_all_led_box_geometries_cm)
                led_geom_vis = g_all_led_box_geometries_cm{i_led_vis};
                patch('Parent', ax_geom, 'Vertices', led_geom_vis.vertices, 'Faces', led_geom_vis.faces, ...
                      'FaceColor', 'red', 'EdgeColor', 'black', 'FaceAlpha', 0.7);
            end
        end
        hold(ax_geom, 'off'); axis(ax_geom, 'equal');
        xlabel(ax_geom, 'X (cm)'); ylabel(ax_geom, 'Y (cm)'); zlabel(ax_geom, 'Z (cm)');
        title(ax_geom, '几何模型预览'); view(ax_geom, 3); grid(ax_geom, 'on');
        light('Parent',ax_geom,'Position',[100 100 100],'Style','infinite'); 
        lighting(ax_geom, 'gouraud'); rotate3d(ax_geom, 'on'); 
        disp('几何预览图已生成。请检查图形窗口。');
        
        disp('--- 创建体素网格 ---');
        all_geom_vertices_cm = g_skin_vertices_cm;
        if ~isempty(g_glasses_vertices_cm), all_geom_vertices_cm = [all_geom_vertices_cm; g_glasses_vertices_cm]; end
        if ~isempty(g_all_led_box_geometries_cm) 
            for i_led = 1:length(g_all_led_box_geometries_cm)
                if ~isempty(g_all_led_box_geometries_cm{i_led}) 
                    all_geom_vertices_cm = [all_geom_vertices_cm; g_all_led_box_geometries_cm{i_led}.vertices];
                end
            end
        end
                    
        min_coords_cm = min(all_geom_vertices_cm, [], 1); 
        max_coords_cm = max(all_geom_vertices_cm, [], 1);
        padding_cm = 1.0; 
        domain_min_cm = min_coords_cm - padding_cm; 
        domain_max_cm = max_coords_cm + padding_cm;
        
        Lx_cm = domain_max_cm(1) - domain_min_cm(1); 
        Ly_cm = domain_max_cm(2) - domain_min_cm(2); 
        Lz_cm = domain_max_cm(3) - domain_min_cm(3); 

        nx = max(10, round(Lx_cm / voxel_size_cm)); 
        ny = max(10, round(Ly_cm / voxel_size_cm)); 
        nz = max(10, round(Lz_cm / voxel_size_cm));
        
        dx = Lx_cm / nx; dy = Ly_cm / ny; dz = Lz_cm / nz; 

        fprintf('  仿真域 (cm): Lx=%.2f, Ly=%.2f, Lz=%.2f\n', Lx_cm, Ly_cm, Lz_cm);
        fprintf('  网格点数: nx=%d, ny=%d, nz=%d (总体素: %d)\n', nx, ny, nz, nx*ny*nz);
        fprintf('  体素尺寸 (cm): dx=%.3f, dy=%.3f, dz=%.3f\n', dx, dy, dz);

        if nx*ny*nz > 3e6 
            warning('警告: 网格点数过大 (%d)，可能导致内存不足或计算缓慢。', nx*ny*nz);
            cont = input('是否继续? (y/n) [n]: ', 's');
            if ~strcmpi(cont, 'y'), disp('仿真中止。'); return; end
        end

        x_coords_centers = domain_min_cm(1) + (0.5:1:nx-0.5)*dx; 
        y_coords_centers = domain_min_cm(2) + (0.5:1:ny-0.5)*dy;
        z_coords_centers = domain_min_cm(3) + (0.5:1:nz-0.5)*dz;

        disp('--- 为体素分配材料属性 ---');
        material_map = zeros(nx, ny, nz, 'uint8'); 
        
        [X_grid_cm, Y_grid_cm, Z_grid_cm] = meshgrid(x_coords_centers, y_coords_centers, z_coords_centers);
        grid_points_cm = [X_grid_cm(:), Y_grid_cm(:), Z_grid_cm(:)];
        disp('体素中心点网格已创建。开始进行inpolyhedron判断...');

        inpolyhedron_exists = (exist('inpolyhedron','file') == 2);
       
        if inpolyhedron_exists && ~isempty(g_all_led_box_geometries_cm)
            disp('  分配LED材料...');
            for i_led = 1:length(g_all_led_box_geometries_cm)
                if isempty(g_all_led_box_geometries_cm{i_led}), continue; end
                led_V_cm = g_all_led_box_geometries_cm{i_led}.vertices;
                led_F    = g_all_led_box_geometries_cm{i_led}.faces;
                min_led = min(led_V_cm, [], 1); max_led = max(led_V_cm, [], 1);
                idx_in_bb_logical = ( grid_points_cm(:,1) >= min_led(1) & grid_points_cm(:,1) <= max_led(1) & ...
                                      grid_points_cm(:,2) >= min_led(2) & grid_points_cm(:,2) <= max_led(2) & ...
                                      grid_points_cm(:,3) >= min_led(3) & grid_points_cm(:,3) <= max_led(3) );
                idx_in_bb = find(idx_in_bb_logical);
                if ~isempty(idx_in_bb)
                    in_led_relative_to_bb = inpolyhedron(struct('faces',led_F,'vertices',led_V_cm), grid_points_cm(idx_in_bb,:));
                    material_map(idx_in_bb(in_led_relative_to_bb)) = 3; 
                end
                fprintf('    LED模块 %d 材料分配完成。\n', i_led);
            end
        elseif ~isempty(g_all_led_box_geometries_cm) && ~inpolyhedron_exists
             disp('由于inpolyhedron缺失，跳过LED材料分配。');
        end
        
        if inpolyhedron_exists && ~isempty(g_glasses_vertices_cm) 
            disp('  分配眼镜材料...');
            min_gla = min(g_glasses_vertices_cm,[],1); max_gla = max(g_glasses_vertices_cm,[],1);
            points_to_check_indices_logical = (material_map(:) == 0 & ... 
                                           grid_points_cm(:,1) >= min_gla(1) & grid_points_cm(:,1) <= max_gla(1) & ...
                                           grid_points_cm(:,2) >= min_gla(2) & grid_points_cm(:,2) <= max_gla(2) & ...
                                           grid_points_cm(:,3) >= min_gla(3) & grid_points_cm(:,3) <= max_gla(3) );
            points_to_check_indices = find(points_to_check_indices_logical);
            if ~isempty(points_to_check_indices)
                in_glasses = inpolyhedron(struct('faces',g_glasses_faces,'vertices',g_glasses_vertices_cm), grid_points_cm(points_to_check_indices,:));
                material_map(points_to_check_indices(in_glasses)) = 2; 
            end
            disp('  眼镜材料分配完成。');
        elseif ~isempty(g_glasses_vertices_cm) && ~inpolyhedron_exists
            disp('由于inpolyhedron缺失，跳过眼镜材料分配。');
        end

        if inpolyhedron_exists && ~isempty(g_skin_vertices_cm)
            disp('  分配皮肤材料...');
            min_skin = min(g_skin_vertices_cm,[],1); max_skin = max(g_skin_vertices_cm,[],1);
            points_to_check_indices_logical = (material_map(:) == 0 & ... 
                                           grid_points_cm(:,1) >= min_skin(1) & grid_points_cm(:,1) <= max_skin(1) & ...
                                           grid_points_cm(:,2) >= min_skin(2) & grid_points_cm(:,2) <= max_skin(2) & ...
                                           grid_points_cm(:,3) >= min_skin(3) & grid_points_cm(:,3) <= max_skin(3) );
            points_to_check_indices = find(points_to_check_indices_logical);
            if ~isempty(points_to_check_indices)
                in_skin = inpolyhedron(struct('faces',g_skin_triangle_indices,'vertices',g_skin_vertices_cm), grid_points_cm(points_to_check_indices,:));
                material_map(points_to_check_indices(in_skin)) = 1; 
            end
            disp('  皮肤材料分配完成。');
        elseif ~isempty(g_skin_vertices_cm) && ~inpolyhedron_exists
             disp('由于inpolyhedron缺失，跳过皮肤材料分配。');
        end
        material_map = reshape(material_map, nx, ny, nz); 

        fprintf('体素材料分布: Air=%d, Skin=%d, Glasses=%d, LED=%d\n', ...
                sum(material_map(:)==0), sum(material_map(:)==1), ...
                sum(material_map(:)==2), sum(material_map(:)==3));

        K_map   = zeros(nx,ny,nz); VHC_map = zeros(nx,ny,nz); Alpha_map = zeros(nx,ny,nz); 
        K_map(material_map == 0) = air_TC_WcmK; VHC_map(material_map == 0) = air_VHC_Jcm3K; Alpha_map(material_map == 0) = air_alpha_cm2s;
        K_map(material_map == 1) = skin_TC_WcmK; VHC_map(material_map == 1) = skin_VHC_Jcm3K; Alpha_map(material_map == 1) = skin_alpha_cm2s;
        K_map(material_map == 2) = glasses_TC_WcmK; VHC_map(material_map == 2) = glasses_VHC_Jcm3K; Alpha_map(material_map == 2) = glasses_alpha_cm2s;
        K_map(material_map == 3) = led_TC_WcmK; VHC_map(material_map == 3) = led_VHC_Jcm3K; Alpha_map(material_map == 3) = led_alpha_cm2s;
        
        if any(VHC_map(:) <= 1e-9), error('体积热容为零或过小，请检查材料分配和单位。'); end 

        disp('--- 开始光学蒙特卡洛模拟 ---');
        AbsorbedEnergyGrid_J_cm3 = zeros(nx,ny,nz); 
        photon_weight_threshold = 1e-4;
        max_steps = 1000; 
        num_photons_to_plot = 100; 
        photon_paths = cell(min(total_photons_to_launch, num_photons_to_plot), 1);

        if ~exist('total_photons_to_launch','var') || isempty(total_photons_to_launch)
            error('变量 total_photons_to_launch 未定义或为空。');
        end
        energy_per_photon = optical_source_power_W * sim_duration_s / total_photons_to_launch; 

        if ~exist('points_data_cm', 'var') || isempty(points_data_cm)
            error('光学MC错误: 光源位置数据 (points_data_cm) 未定义或为空。');
        end
        disp('当前工作区变量 (光学MC循环前):');
        whos points_data_cm domain_min_cm dx dy dz material_map optical_props AbsorbedEnergyGrid_J_cm3 total_photons_to_launch

        for p_idx = 1:total_photons_to_launch
            if mod(p_idx, round(total_photons_to_launch/20)) == 0 || p_idx == 1 
                fprintf('  模拟光子: %d / %d\n', p_idx, total_photons_to_launch);
            end

            current_launch_point_idx = mod(p_idx-1, size(points_data_cm,1)) + 1; 
            pos = points_data_cm(current_launch_point_idx, :); 
            
            initial_dir_local = [0; -1; 0]; 
            ROTATION_ANGLE_DEGREES_local = -30; 
            rot_ang_rad_local = ROTATION_ANGLE_DEGREES_local * pi/180;
            R_x_mat_local = [1,0,0;0,cos(rot_ang_rad_local),-sin(rot_ang_rad_local);0,sin(rot_ang_rad_local),cos(rot_ang_rad_local)];
            dir = (R_x_mat_local * initial_dir_local)'; 
            dir = dir / norm(dir); 

            weight = 1.0; 
            path_segments = pos; 

            for step_idx = 1:max_steps
                ix = max(1, min(nx, floor((pos(1) - domain_min_cm(1))/dx) + 1));
                iy = max(1, min(ny, floor((pos(2) - domain_min_cm(2))/dy) + 1));
                iz = max(1, min(nz, floor((pos(3) - domain_min_cm(3))/dz) + 1));

                current_material_id = material_map(ix,iy,iz);
                if current_material_id == 0, current_material_id = 1; end 
                
                mat_idx_for_optical_props = current_material_id + 1; 
                if mat_idx_for_optical_props > length(optical_props)
                    warning('光学MC: 无效的 material_id (%d) 用于 optical_props。默认为空气。', current_material_id);
                    mat_idx_for_optical_props = 1; 
                end

                mua = optical_props(mat_idx_for_optical_props).mua; 
                mus = optical_props(mat_idx_for_optical_props).mus;
                g   = optical_props(mat_idx_for_optical_props).g;
                mut = mua + mus;

                if mut <= 1e-6 
                    s = Inf; 
                else
                    s = -log(rand()) / mut; 
                end
                
                pos_new = pos + s * dir;
                
                if pos_new(1) < domain_min_cm(1) || pos_new(1) > domain_max_cm(1) || ...
                   pos_new(2) < domain_min_cm(2) || pos_new(2) > domain_max_cm(2) || ...
                   pos_new(3) < domain_min_cm(3) || pos_new(3) > domain_max_cm(3)
                    break; 
                end
                
                d_weight_abs = weight * (mua / mut); 
                if mut <= 1e-6, d_weight_abs = 0; end 
                
                AbsorbedEnergyGrid_J_cm3(ix,iy,iz) = AbsorbedEnergyGrid_J_cm3(ix,iy,iz) + d_weight_abs; 
                weight = weight * (mus / mut); 
                                                
                pos = pos_new; 
                if p_idx <= num_photons_to_plot
                    path_segments = [path_segments; pos]; %#ok<AGROW>
                end

                if mus > 1e-6 && weight > 0 
                    if abs(g) < 1e-3 
                        costheta = 2*rand() - 1;
                    else
                        costheta = (1/(2*g)) * (1 + g^2 - ((1-g^2)/(1-g+2*g*rand()))^2);
                    end
                    sintheta = sqrt(max(0,1 - costheta^2)); 
                    phi = 2*pi*rand();
                    
                    ux = sintheta * cos(phi); uy = sintheta * sin(phi); uz = costheta;
                    
                    if abs(dir(3)) > 0.99999 
                        dir_new_x = ux; dir_new_y = uy; dir_new_z = uz * sign(dir(3));
                    else
                        sqrt_term = sqrt(max(0,1 - dir(3)^2)); 
                        if sqrt_term < 1e-9 
                           dir_new_x = ux; dir_new_y = uy; dir_new_z = uz * sign(dir(3)); 
                        else
                            dir_new_x = (ux * dir(1) * dir(3) - uy * dir(2)) / sqrt_term + dir(1) * uz;
                            dir_new_y = (ux * dir(2) * dir(3) + uy * dir(1)) / sqrt_term + dir(2) * uz;
                            dir_new_z = -ux * sqrt_term + dir(3) * uz;
                        end
                    end
                    dir = [dir_new_x, dir_new_y, dir_new_z];
                    dir = dir / (norm(dir) + eps); 
                else 
                    if weight <= 0, break; end 
                end
                
                if weight < photon_weight_threshold
                    m_rr = 10; 
                    if rand() > (1/m_rr)
                        break; 
                    else
                        weight = weight * m_rr; 
                    end
                end
            end 
            if p_idx <= num_photons_to_plot
                photon_paths{p_idx} = path_segments;
            end
        end 
        disp('光学蒙特卡洛模拟完成。');

        Q_gen_map_W_cm3 = (AbsorbedEnergyGrid_J_cm3 / total_photons_to_launch) * optical_source_power_W / (dx*dy*dz);
        
        fprintf('  光学吸收产生的最大热源: %.4e W/cm^3\n', max(Q_gen_map_W_cm3(:)));
        if max(Q_gen_map_W_cm3(:)) < 1e-6
            disp('警告: 光学吸收产生的热源非常小，温度变化可能不明显。');
        end

        if num_photons_to_plot > 0 && exist('h_geom_fig','var') && ishandle(h_geom_fig) && isvalid(h_geom_fig)
            disp('--- 可视化部分光子路径 ---');
            figure(h_geom_fig); 
            hold(ax_geom, 'on');
            for i_path = 1:min(total_photons_to_launch, num_photons_to_plot)
                if ~isempty(photon_paths{i_path}) && size(photon_paths{i_path},1) > 1
                    plot3(ax_geom, photon_paths{i_path}(:,1), photon_paths{i_path}(:,2), photon_paths{i_path}(:,3), '-c', 'LineWidth', 0.5);
                end
            end
            hold(ax_geom, 'off');
            title(ax_geom, '几何模型预览与部分光子路径');
            disp('部分光子路径已添加到几何预览图。');
        end


        disp('--- 开始有限差分热传导模拟 ---');
        alpha_max = max(Alpha_map(:));
        if alpha_max <= 1e-9, error('最大热扩散率为零或过小，无法计算dt。'); end 
        dt_crit = min([dx^2, dy^2, dz^2]) / (6 * alpha_max); 
        dt = 0.9 * dt_crit; 
        actual_num_time_steps = ceil(sim_duration_s / dt);
        
        fprintf('  临界时间步长 dt_crit: %.4e s\n', dt_crit);
        fprintf('  使用时间步长 dt: %.4e s\n', dt);
        fprintf('  总时间步数: %d\n', actual_num_time_steps);

        T = ones(nx,ny,nz) * ambient_temp_C; 
        T_history = zeros(nx,ny,nz, num_time_steps_output); 
        output_time_indices = round(linspace(1, actual_num_time_steps, num_time_steps_output)); 
        output_counter = 1;

        for t_step = 1:actual_num_time_steps
            T_old = T;
            for i = 2:nx-1
                for j = 2:ny-1
                    for k = 2:nz-1
                        if VHC_map(i,j,k) <= 1e-9, continue; end 
                        term_x = (T_old(i+1,j,k) - 2*T_old(i,j,k) + T_old(i-1,j,k))/(dx^2);
                        term_y = (T_old(i,j+1,k) - 2*T_old(i,j,k) + T_old(i,j-1,k))/(dy^2);
                        term_z = (T_old(i,j,k+1) - 2*T_old(i,j,k) + T_old(i,j,k-1))/(dz^2);
                        
                        T(i,j,k) = T_old(i,j,k) + Alpha_map(i,j,k) * dt * (term_x + term_y + term_z) + ...
                                     (Q_gen_map_W_cm3(i,j,k) / VHC_map(i,j,k)) * dt;
                    end
                end
            end

            T(1,:,:)   = T(2,:,:); T(nx,:,:) = T(nx-1,:,:);
            T(:,1,:)   = T(:,2,:); T(:,ny,:) = T(:,ny-1,:);
            T(:,:,1)   = T(:,:,2); T(:,:,nz) = T(:,:,nz-1);
            
            if ismember(t_step, output_time_indices) && output_counter <= num_time_steps_output
                T_history(:,:,:,output_counter) = T;
                if mod(t_step, round(actual_num_time_steps/10)) == 0 || t_step == 1 
                    fprintf('  时间步: %d/%d, 当前模拟时间: %.2f s, 最大温度: %.2f °C\n', ...
                            t_step, actual_num_time_steps, t_step*dt, max(T(:)));
                end
                output_counter = output_counter + 1;
            end
        end
        disp('热传导模拟完成。');

        disp('--- 可视化结果 ---');
        if output_counter > 1 
            final_T_field = T_history(:,:,:,output_counter-1); 
        else
            final_T_field = T; 
        end

        slice_z_idx = ceil(nz/2);
        figure('Name', '光学MC与热模型: 温度分布切片');
        imagesc(x_coords_centers, y_coords_centers, squeeze(final_T_field(:,:,slice_z_idx))'); 
        axis image; colorbar;
        xlabel('X (cm)'); ylabel('Y (cm)');
        title(sprintf('温度分布 @ Z=%.2f cm, Time=%.1f s (Max T: %.2f °C)', ...
                      z_coords_centers(slice_z_idx), sim_duration_s, max(final_T_field(:)))); 
        set(gca, 'YDir', 'normal');
        disp('结果图已显示。');

    catch ME
        disp('错误发生');
        fprintf(2, '错误信息: %s\n', ME.message);
        if ~isempty(ME.stack)
            for k_stack = 1:length(ME.stack)
                fprintf(2, '  文件: %s, 函数: %s, 行: %d\n', ME.stack(k_stack).file, ME.stack(k_stack).name, ME.stack(k_stack).line);
            end
        end
        disp(ME); 
        
        disp('尝试在错误后恢复原始路径和目录...');
        if exist('original_dir_at_start', 'var') && ischar(original_dir_at_start) && ~isempty(original_dir_at_start) && isfolder(original_dir_at_start)
            try cd(original_dir_at_start); disp(['已返回原始目录: ', original_dir_at_start]); catch, end
        end
        if exist('original_path_at_start', 'var') && ischar(original_path_at_start) && ~isempty(original_path_at_start)
            try path(original_path_at_start); disp('已恢复原始MATLAB路径。'); catch, end
        end
        disp('错误处理结束');
        return; 
    end
    
    disp('--- 成功运行后清理 ---');
    if exist('original_dir_at_start', 'var') && ischar(original_dir_at_start) && ~isempty(original_dir_at_start) && isfolder(original_dir_at_start)
        try cd(original_dir_at_start); disp(['已返回原始目录: ', original_dir_at_start]); catch, end
    end
    if exist('original_path_at_start', 'var') && ischar(original_path_at_start) && ~isempty(original_path_at_start)
        try path(original_path_at_start); disp('已恢复原始MATLAB路径。'); catch, end
    end
    
    disp('光学与热分析流程结束。');
end 
