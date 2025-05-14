function mc_simulation_with_bvh_demo()


   
    vertices = [
        0 0 0; 1 0 0; 1 1 0; 0 1 0; 
        0 0 1; 1 0 1; 1 1 1; 0 1 1  
    ];
    faces = [ 
        1 2 3; 1 3 4; 
        5 6 7; 5 7 8; 
        1 2 6; 1 6 5; 
        2 3 7; 2 7 6; 
        3 4 8; 3 8 7;
        4 1 5; 4 5 8  
    ];
    
    scene_triangles_vertices = cell(size(faces,1), 1);
    for i = 1:size(faces,1)
        scene_triangles_vertices{i} = vertices(faces(i,:), :);
    end
  
    extra_tri_verts = [ -0.5 -0.5 0.5; 1.5 0.5 1.5; 0.5 1.5 0.2 ];
    scene_triangles_vertices{end+1} = extra_tri_verts;
    
    num_total_triangles = length(scene_triangles_vertices);
    fprintf('场景中定义的三角形数量: %d\n', num_total_triangles);

    % --- BVH 构建 ---
    disp('开始构建BVH...');
    tic;
    bvh_root = build_bvh_recursive(scene_triangles_vertices, 1:num_total_triangles, 0);
    toc;
    disp('BVH构建完成。');
    
   
    num_photons = 100;    
    max_steps = 100;          
    initial_weight = 1.0;
    weight_threshold = 1e-3;
    
    mua = 0.1; 
    mus = 10.0; 
    g_anisotropy = 0.9;
    n_medium = 1.33; 
    n_outside = 1.0; 
    
    mut = mua + mus;
    albedo = mus / mut; 

      source_pos = [0.5; 0.5; -0.1]; 
    source_dir = [0; 0; 1];    

    all_photon_paths = cell(num_photons, 1);
    total_absorbed_weight = 0;

 
    disp('开始蒙特卡洛光子追踪...');
    tic;
    for i_photon = 1:num_photons
        if mod(i_photon, round(num_photons/10))==0, fprintf('  追踪光子: %d/%d\n', i_photon, num_photons); end


        photon.pos = source_pos;
        photon.dir = source_dir;
        photon.weight = initial_weight;
        photon.alive = true;
        photon.path = photon.pos';

        for i_step = 1:max_steps
            if ~photon.alive, break; end

   
            step_length = -log(rand()) / mut;

           
            [hit_triangle_idx, dist_to_triangle, hit_barycentric_coords, hit_normal_geom] = ...
                intersect_scene_with_bvh(photon.pos, photon.dir, bvh_root, scene_triangles_vertices);

            
            if ~isempty(hit_triangle_idx) && dist_to_triangle < step_length
                
                photon.pos = photon.pos + photon.dir * dist_to_triangle; 
                photon.path = [photon.path; photon.pos'];
             
                if rand() < 0.5 
                    total_absorbed_weight = total_absorbed_weight + photon.weight;
                    photon.alive = false;
                else 
                    photon.dir = photon.dir - 2 * dot(photon.dir, hit_normal_geom) * hit_normal_geom;
                    photon.pos = photon.pos + photon.dir * 1e-4; 
                end
            else
              
                photon.pos = photon.pos + photon.dir * step_length;
                photon.path = [photon.path; photon.pos'];

  
                if any(photon.pos < [-1 -1 -0.2]') || any(photon.pos > [2 2 2]') 
                    photon.alive = false; %
                    continue;
                end

               
                delta_weight = photon.weight * (mua / mut);
                photon.weight = photon.weight - delta_weight;
                total_absorbed_weight = total_absorbed_weight + delta_weight;

                if photon.weight < weight_threshold 
                    m_rr = 10;
                    if rand() > (1/m_rr)
                        photon.alive = false;
                    else
                        photon.weight = photon.weight * m_rr;
                    end
                end
                if ~photon.alive, continue; end
                
                            if mus > 1e-6
                    if abs(g_anisotropy) < 1e-3, costheta_scatter = 2*rand() - 1;
                    else, costheta_scatter = (1/(2*g_anisotropy)) * (1 + g_anisotropy^2 - ((1-g_anisotropy^2)/(1-g_anisotropy+2*g_anisotropy*rand()))^2);
                    end
                    sintheta_scatter = sqrt(max(0,1 - costheta_scatter^2));
                    phi_scatter = 2*pi*rand();

                    ux_scatter = sintheta_scatter * cos(phi_scatter);
                    uy_scatter = sintheta_scatter * sin(phi_scatter);
                    uz_scatter = costheta_scatter;
                    
                   
                    w_vec = photon.dir;
                    if abs(w_vec(3)) > 0.99999
                        u_vec = [1;0;0]; 
                        if abs(dot(w_vec,u_vec)) > 0.99999, u_vec = [0;1;0]; end
                        v_vec = cross(w_vec,u_vec); v_vec = v_vec/norm(v_vec);
                        u_vec = cross(v_vec,w_vec); u_vec = u_vec/norm(u_vec);
                    else
                        u_vec = cross([0;0;1], w_vec); u_vec = u_vec/norm(u_vec);
                        v_vec = cross(w_vec, u_vec);   v_vec = v_vec/norm(v_vec);
                    end
                    
                    photon.dir = ux_scatter * u_vec + uy_scatter * v_vec + uz_scatter * w_vec;
                    photon.dir = photon.dir / norm(photon.dir);
                end
            end
        end 
        all_photon_paths{i_photon} = photon.path;
    end 
    toc;
    disp('蒙特卡洛追踪完成。');
    fprintf('总吸收权重: %.4f\n', total_absorbed_weight);

        figure('Name', 'BVH加速的MC光子追踪演示');
    ax_3d = subplot(1,2,1);
    hold(ax_3d, 'on');
    title(ax_3d, '3D场景与BVH包围盒');
    
   
    for i = 1:num_total_triangles
        verts = scene_triangles_vertices{i};
        patch(ax_3d, verts(:,1), verts(:,2), verts(:,3), 'c', 'FaceAlpha', 0.3, 'EdgeColor', 'k');
    end
    
   
    plot_bvh_recursive(ax_3d, bvh_root, 0);
    
    xlabel(ax_3d,'X (mm)'); ylabel(ax_3d,'Y (mm)'); zlabel(ax_3d,'Z (mm)');
    axis equal; view(3); grid on;
    
     ax_paths = subplot(1,2,2);
    hold(ax_paths, 'on');
    title(ax_paths, '部分光子路径 (3D)');
    num_paths_to_plot = min(num_photons, 20); 
    colors_path = lines(num_paths_to_plot);
    for i=1:num_paths_to_plot
        if ~isempty(all_photon_paths{i})
            path_data = all_photon_paths{i};
            plot3(ax_paths, path_data(:,1), path_data(:,2), path_data(:,3), '-', 'Color', colors_path(i,:), 'LineWidth', 0.8);
            plot3(ax_paths, path_data(1,1), path_data(1,2), path_data(1,3), 'go', 'MarkerFaceColor', 'g'); % 起点
            plot3(ax_paths, path_data(end,1), path_data(end,2), path_data(end,3), 'rx', 'MarkerSize', 8); % 终点
        end
    end
    xlabel(ax_paths,'X (mm)'); ylabel(ax_paths,'Y (mm)'); zlabel(ax_paths,'Z (mm)');
    axis equal; view(3); grid on;
    
   
    linkprop([ax_3d, ax_paths], {'CameraPosition', 'CameraUpVector', 'CameraTarget', 'CameraViewAngle'});
    
    disp('可视化完成。');
end


function node = build_bvh_recursive(all_triangles_verts, triangle_indices, depth)


    node = struct('aabb_min', [], 'aabb_max', [], 'is_leaf', false, ...
                  'left_child', [], 'right_child', [], 'tri_indices', []);

    num_node_triangles = length(triangle_indices);

    if num_node_triangles == 0
        node.is_leaf = true; 
        return;
    end

    node_aabb_min = [Inf; Inf; Inf];
    node_aabb_max = [-Inf; -Inf; -Inf];
    for i = 1:num_node_triangles
        idx = triangle_indices(i);
        tri_verts = all_triangles_verts{idx}; 
        node_aabb_min = min(node_aabb_min, min(tri_verts,[],1)');
        node_aabb_max = max(node_aabb_max, max(tri_verts,[],1)');
    end
    node.aabb_min = node_aabb_min;
    node.aabb_max = node_aabb_max;


    max_tris_in_leaf = 4; 
    max_bvh_depth = 15;
    if num_node_triangles <= max_tris_in_leaf || depth >= max_bvh_depth
        node.is_leaf = true;
        node.tri_indices = triangle_indices;
        return;
    end


    aabb_dims = node.aabb_max - node.aabb_min;
    [~, split_axis] = max(aabb_dims); 
    split_coord = node.aabb_min(split_axis) + aabb_dims(split_axis) / 2;

    left_indices = [];
    right_indices = [];
    for i = 1:num_node_triangles
        idx = triangle_indices(i);
        tri_verts = all_triangles_verts{idx};
      
        tri_center = mean(tri_verts, 1)';
        if tri_center(split_axis) < split_coord
            left_indices = [left_indices, idx];
        else
            right_indices = [right_indices, idx];
        end
    end
    
    if isempty(left_indices) || isempty(right_indices)
     
        node.is_leaf = true;
        node.tri_indices = triangle_indices;
        return;
    end


    node.left_child = build_bvh_recursive(all_triangles_verts, left_indices, depth + 1);
    node.right_child = build_bvh_recursive(all_triangles_verts, right_indices, depth + 1);
end


function [hit, t_near, t_far] = ray_aabb_intersect(ray_origin, ray_dir_inv, aabb_min, aabb_max)

    t1 = (aabb_min - ray_origin) .* ray_dir_inv;
    t2 = (aabb_max - ray_origin) .* ray_dir_inv;

    t_min_comp = min(t1, t2);
    t_max_comp = max(t1, t2);

    t_near = max(max(t_min_comp(1), t_min_comp(2)), t_min_comp(3));
    t_far  = min(min(t_max_comp(1), t_max_comp(2)), t_max_comp(3));

    hit = (t_near < t_far) && (t_far > 0); 
end

function [hit, t, u, v, normal_geom] = ray_triangle_intersect_moller_trumbore(ray_origin, ray_dir, v0, v1, v2)

    
    hit = false; t = Inf; u = 0; v = 0; normal_geom = [0;0;0];
    epsilon = 1e-7;
    edge1 = v1 - v0;
    edge2 = v2 - v0;
    
    pvec = cross(ray_dir, edge2);
    det_val = dot(edge1, pvec);

    if abs(det_val) < epsilon
        return;
    end
    
    inv_det = 1.0 / det_val;
    
    tvec = ray_origin - v0;
    u = dot(tvec, pvec) * inv_det;
    if u < 0 || u > 1
        return;
    end
    
    qvec = cross(tvec, edge1);
    v = dot(ray_dir, qvec) * inv_det;
    if v < 0 || u + v > 1
        return;
    end
    
    t = dot(edge2, qvec) * inv_det;
    
    if t > epsilon %
        hit = true;
        normal_geom = cross(edge1, edge2);
        normal_geom = normal_geom / norm(normal_geom);
     
        if dot(normal_geom, ray_dir) > 0
            normal_geom = -normal_geom;
        end
    else
        t = Inf; 
    end
end


function [closest_hit_tri_idx, min_dist, closest_bary_coords, closest_normal] = ...
    intersect_scene_with_bvh(ray_origin, ray_dir, bvh_root, all_triangles_verts)
    
    closest_hit_tri_idx = [];
    min_dist = Inf;
    closest_bary_coords = [];
    closest_normal = [];
    
    if isempty(bvh_root), return; end

    ray_dir_inv = 1 ./ ray_dir;
    
    
    stack = {}; 
    if ~isempty(bvh_root)
        stack{1} = bvh_root; 
    end
    
    while ~isempty(stack)
        node = stack{1}; 
        stack(1) = [];  
        
      
        [hit_aabb, t_near_aabb, ~] = ray_aabb_intersect(ray_origin, ray_dir_inv, node.aabb_min, node.aabb_max);
        
        if hit_aabb && t_near_aabb < min_dist 
            if node.is_leaf

                for i = 1:length(node.tri_indices)
                    tri_global_idx = node.tri_indices(i);
                    verts = all_triangles_verts{tri_global_idx};
                    v0 = verts(1,:)'; v1 = verts(2,:)'; v2 = verts(3,:)';
                    
                    [hit_tri, t_tri, u_tri, v_tri, normal_tri] = ...
                        ray_triangle_intersect_moller_trumbore(ray_origin, ray_dir, v0, v1, v2);
                    
                    if hit_tri && t_tri < min_dist
                        min_dist = t_tri;
                        closest_hit_tri_idx = tri_global_idx;
                        closest_bary_coords = [u_tri, v_tri];
                        closest_normal = normal_tri;
                    end
                end
            else
                
                if ~isempty(node.right_child)
                    stack = [{node.right_child}, stack]; 
                end
                if ~isempty(node.left_child)
                    stack = [{node.left_child}, stack];   
                end
            end
        end
    end
end


function plot_bvh_recursive(ax, node, depth)
    if isempty(node), return; end
    

    colors = lines(10); 
    color_idx = mod(depth, size(colors,1)) + 1;
        min_p = node.aabb_min;
    max_p = node.aabb_max;
    

    v = [min_p(1) min_p(2) min_p(3); max_p(1) min_p(2) min_p(3); max_p(1) max_p(2) min_p(3); min_p(1) max_p(2) min_p(3);
         min_p(1) min_p(2) max_p(3); max_p(1) min_p(2) max_p(3); max_p(1) max_p(2) max_p(3); min_p(1) max_p(2) max_p(3)];

    edges = [1 2; 2 3; 3 4; 4 1; 5 6; 6 7; 7 8; 8 5; 1 5; 2 6; 3 7; 4 8];
    
    for i = 1:size(edges,1)
        plot3(ax, [v(edges(i,1),1), v(edges(i,2),1)], ...
                  [v(edges(i,1),2), v(edges(i,2),2)], ...
                  [v(edges(i,1),3), v(edges(i,2),3)], ...
                  '-', 'Color', [colors(color_idx,:), 0.3], 'LineWidth', 0.5); % 带透明度
    end

    if ~node.is_leaf
        plot_bvh_recursive(ax, node.left_child, depth + 1);
        plot_bvh_recursive(ax, node.right_child, depth + 1);
    end
end
