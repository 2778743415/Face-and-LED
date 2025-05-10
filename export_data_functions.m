function export_light_intensity_callback(~, ~)
    % 导出光强数据到文件
    data = getappdata(fig, 'AccumulatedEnergy_SkinVertex');
    writematrix(data, 'light_intensity_data.csv');
    helpdlg('数据已保存到 light_intensity_data.csv');
end

function export_all_history_callback(~, ~)
    % 导出历史光源位置
    history = getappdata(fig, 'all_input_points_history');
    writematrix(cell2mat(history), 'history_light_positions.csv');
end