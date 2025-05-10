function stl_face_gui_mcmatlab_main
% 主函数：初始化GUI并注册回调
    % 创建主界面
    fig = create_main_figure();
    
    % 初始化全局变量
    initialize_appdata(fig);
    
    % 创建UI控件并绑定回调
    create_ui_controls(fig);
    
    % 设置关闭时的清理操作
    set(fig, 'CloseRequestFcn', @close_request_callback);
end

% ---------------------
% 子函数调用（定义在同文件或独立文件）
% ---------------------
function fig = create_main_figure()
    fig = figure('Name', 'STL眼镜模型与人脸曲面可视化 (蒙特卡洛 - 批量LED组合与多角度旋转)', ...
                 'NumberTitle', 'off', ...
                 'Position', [100 100 850 550]);
end

function initialize_appdata(fig)
    % 初始化全局变量（如历史记录、路径设置）
    setappdata(fig, 'all_input_points_history', {});
    setappdata(fig, 'mcmatlab_path', 'E:\toolbox\MCmatlab\MCmatlab'); % 示例路径
end

function create_ui_controls(fig)
    % 创建UI控件并绑定回调
    uicontrol('Style', 'pushbutton', 'String', '加载眼镜与人脸', ...
              'FontSize', 12, 'Position', [20 490 150 30], ...
              'Callback', @(~,~) load_models_callback(fig));
    
    uicontrol('Style', 'pushbutton', 'String', '改变皮肤颜色', ...
              'FontSize', 12, 'Position', [20 450 150 30], ...
              'Callback', @(~,~) change_surface_color_callback(fig));
    
    % 其他控件类似...
end

function close_request_callback(~, ~)
    % 清理资源（如删除临时文件、关闭进度条）
    close all;
end