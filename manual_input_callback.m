function manual_input_callback(~, ~)
    % 手动设置光源模块中心
    num_cuboids = str2double(inputdlg('输入光源模块数量'));
    if isempty(num_cuboids) || num_cuboids <= 0, return; end
    
    % 获取用户输入并保存
    centers = zeros(num_cuboids, 3);
    for i = 1:num_cuboids
        centers(i,:) = get_user_input_position();
    end
    
    % 保存到appdata
    setappdata(fig, 'manual_light_centers', centers);
end

function center = get_user_input_position()
    % 弹出对话框获取坐标
    coords = inputdlg('输入光源坐标 (X,Y,Z):');
    center = str2num(coords{:});
end