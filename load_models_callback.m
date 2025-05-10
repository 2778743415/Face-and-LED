function load_models_callback(~, ~)
    % 加载STL模型（眼镜、人脸）
    try
        [g_skin_vertices, g_skin_triangle_indices] = import_stl('skin.stl');
        [g_glasses_vertices, g_glasses_faces] = import_stl('glasses.stl');
        
        % 更新图形界面
        update_axes_plot(g_skin_vertices, g_skin_triangle_indices);
    catch ME
        errordlg(['加载模型失败: ', ME.message], '错误');
    end
end

function [vertices, faces] = import_stl(filename)
    % 使用stlread等函数导入STL文件
    [vertices, faces] = stlread(filename);
end

function update_axes_plot(vertices, faces)
    % 在GUI中绘制模型
    axes(ax);
    trisurf(faces, vertices(:,1), vertices(:,2), vertices(:,3));
end