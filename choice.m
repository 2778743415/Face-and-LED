% MATLAB 脚本：读取光强数据文件，进行统计分析并三维可视化

% --- 用户配置参数 ---
basePath = 'C:\Users\WC\Desktop\test_output_matlab'; % 指定包含角度文件夹的根目录
angleDegs = [30, 40, 50, 60, 70]; % 要处理的角度值
% --------------------

% 初始化存储结果的容器
allStdDevs = [];           % 存储所有文件的标准差
allCVs = [];               % 存储所有文件的变异系数
allPercentagesInRange = []; % 存储所有文件在均值±10%范围内的数据占比
allAnglesForPlot = [];     % 存储每个数据点对应的角度，用于绘图区分颜色

% 定义角度文件夹名称格式
angleFolderNames = arrayfun(@(x) sprintf('angle_%d_deg', x), angleDegs, 'UniformOutput', false);
numAngleFolders = length(angleFolderNames);

% 为不同角度定义颜色 (MATLAB默认颜色序列)
plotColors = lines(numAngleFolders);

fprintf('开始处理数据...\n');

% 遍历每个角度文件夹
for i = 1:numAngleFolders
    currentAngle = angleDegs(i);
    currentFolderName = angleFolderNames{i};
    folderPath = fullfile(basePath, currentFolderName);

    fprintf('\n正在处理文件夹: %s (角度: %d°)\n', currentFolderName, currentAngle);

    % 检查文件夹是否存在
    if ~isfolder(folderPath)
        warning('文件夹 %s 不存在，跳过。\n', folderPath);
        continue;
    end

    % 构建文件搜索模式，例如 combo_*_angle30.txt
    filePattern = fullfile(folderPath, sprintf('combo_*_angle%d.txt', currentAngle));
    txtFiles = dir(filePattern);

    if isempty(txtFiles)
        fprintf('在 %s 中没有找到匹配的文件 (例如 combo_*_angle%d.txt)。\n', folderPath, currentAngle);
        continue;
    end

    fprintf('找到 %d 个文件进行处理。\n', length(txtFiles));

    % 遍历文件夹中的每个txt文件
    for j = 1:length(txtFiles)
        fileName = txtFiles(j).name;
        filePath = fullfile(folderPath, fileName);
        fprintf('  正在读取文件: %s ... ', fileName);

        try
            % 读取数据，假设数据是数值型，并且由空格、制表符或逗号分隔
            % readmatrix 会尝试自动检测分隔符
            data = readmatrix(filePath);

            if isempty(data)
                fprintf('文件为空或无法读取内容，跳过。\n');
                continue;
            end

            % 假设光强度数据在最后一列
            if size(data, 2) < 1
                fprintf('文件不包含任何数据列，跳过。\n');
                continue;
            end
            intensityData = data(:, end);

            if isempty(intensityData)
                fprintf('未能提取光强度数据，跳过。\n');
                continue;
            end
            
            % 确保intensityData是列向量
            if size(intensityData, 2) > 1
                intensityData = intensityData'; % 转置为列向量
            end
            
            % 移除NaN值，以防影响计算
            intensityData = intensityData(~isnan(intensityData));
            if isempty(intensityData)
                fprintf('移除NaN后无有效光强度数据，跳过。\n');
                continue;
            end


            % 1. 计算标准差 (Standard Deviation)
            stdVal = std(intensityData);

            % 2. 计算平均值和变异系数 (Coefficient of Variation)
            meanVal = mean(intensityData);
            if meanVal == 0
                cvVal = NaN; % 或者根据需求设为0或Inf
                fprintf('平均值为0，变异系数设为NaN。');
            else
                cvVal = stdVal / meanVal;
            end

            % 3. 计算在平均值±10.0%范围内的数据占比
            lowerBound = meanVal * (1 - 0.10);
            upperBound = meanVal * (1 + 0.10);
            
            numInRange = sum(intensityData >= lowerBound & intensityData <= upperBound);
            totalDataPoints = length(intensityData);
            
            if totalDataPoints == 0
                percentageInRange = NaN; % 如果没有数据点，则为NaN
                 fprintf('无有效数据点计算占比，设为NaN。');
            else
                percentageInRange = (numInRange / totalDataPoints) * 100;
            end

            % 存储计算结果
            allStdDevs = [allStdDevs; stdVal];
            allCVs = [allCVs; cvVal];
            allPercentagesInRange = [allPercentagesInRange; percentageInRange];
            allAnglesForPlot = [allAnglesForPlot; currentAngle];
            
            fprintf('处理完成 (Std=%.2f, CV=%.2f, PctInRange=%.2f%%)。\n', stdVal, cvVal, percentageInRange);

        catch ME
            fprintf('读取或处理文件 %s 时发生错误: %s。跳过此文件。\n', fileName, ME.message);
        end
    end
end

fprintf('\n所有文件处理完毕。\n');

% --- 三维可视化 ---
if isempty(allStdDevs)
    fprintf('没有成功处理任何数据，无法生成图像。\n');
else
    fprintf('正在生成三维散点图...\n');
    figure('Name', '光强度数据三维分析图', 'NumberTitle', 'off');
    hold on;
    
    legendEntries = {}; % 用于存储图例条目

    for i = 1:numAngleFolders
        currentAngleToPlot = angleDegs(i);
        % 找到属于当前角度的数据点的索引
        indices = (allAnglesForPlot == currentAngleToPlot);

        if any(indices) % 如果当前角度有数据点
            scatter3(allStdDevs(indices), allCVs(indices), allPercentagesInRange(indices), ...
                      50, plotColors(i,:), 'filled', 'MarkerFaceAlpha', 0.7); % 使用指定颜色和大小绘制点
            legendEntries{end+1} = sprintf('%d° 角', currentAngleToPlot);
        end
    end

    hold off;

    % 设置坐标轴标签和标题
    xlabel('标准差 (\sigma)', 'FontSize', 12);
    ylabel('变异系数 (CV = \sigma / \mu)', 'FontSize', 12);
    zlabel('均值\pm10%范围内数据占比 (%)', 'FontSize', 12);
    title('光强度数据统计特性三维可视化', 'FontSize', 14);

    % 添加图例
    if ~isempty(legendEntries)
        legend(legendEntries, 'Location', 'northeastoutside', 'FontSize', 10);
    end
    
    grid on;        % 显示网格
    axis tight;     % 紧缩坐标轴
    view(3);        % 设置为标准三维视角 (azimuth = -37.5, elevation = 30)
    rotate3d on;    % 允许用鼠标交互式旋转图像
    
    fprintf('三维散点图已生成。\n');
end

