function run_monte_carlo_simulation(~, ~)
    % 配置MCmatlab参数
    mc_input = configure_mc_input();
    
    % 运行仿真
    try
        mc_output = mcmatlab(mc_input);
        visualize_results(mc_output);
    catch ME
        errordlg(['仿真出错: ', ME.message], 'MCmatlab错误');
    end
end

function mc_input = configure_mc_input()
    % 设置几何、光源、介质属性
    mc_input.G.Lx = 10; % 示例参数
    mc_input.MC.lightSource.sourceType = 5;
    % ...其他参数配置
end

function visualize_results(mc_output)
    % 可视化温度分布
    figure;
    imagesc(mc_output.Thermal.x_coords, mc_output.Thermal.y_coords, mc_output.Thermal.T(:,:,end));
    colorbar;
end