function run_thermal_simulation(~, ~)
    % 配置热力学参数
    mc_input.HS.active = true;
    mc_input.HS.durationOn = 60;
    
    % 运行MCmatlab仿真
    mc_output = mcmatlab(mc_input);
    
    % 可视化温度分布
    plot_temperature_distribution(mc_output);
end

function plot_temperature_distribution(mc_output)
    % 绘制温度切片
    final_T_slice_z_idx = ceil(size(mc_output.Thermal.T,3)/2);
    imagesc(mc_output.Thermal.x_coords, mc_output.Thermal.y_coords, ...
            squeeze(mc_output.Thermal.T(:,:,final_T_slice_z_idx)));
    title('温度分布');
end