function [R2, RMSE, MARE] = evaluate_metrics(y_true, y_pred)
    % 计算 R^2 (决定系数)
    SS_res = sum((y_true - y_pred).^2);
    SS_tot = sum((y_true - mean(y_true)).^2);
    R2 = 1 - (SS_res / SS_tot);

    % 计算 RMSE (均方根误差)
    RMSE = sqrt(mean((y_true - y_pred).^2));

    % 计算 MARE (平均绝对相对误差)
    MARE = mean(abs((y_true - y_pred) ./ y_true)) * 100; % 以百分比表示

end
