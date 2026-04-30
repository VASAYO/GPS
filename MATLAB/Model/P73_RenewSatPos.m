function OutSatPos = P73_RenewSatPos(SatPos, TProp, Params) %#ok<INUSD>
%
% Пересчёт координат спутника по новому значению времени распространения
% сигнала
%
% Входные переменные
%   SatPos - массив координат спутника (8х1) из P72_GetSatPos;
%   TProp - время распространения сигнала.
%
% Выходные переменные
%   OutSatPos - массив (3x1) скорректированных координат спутника.

%% ПЕРЕСЧЁТ КООРДИНАТ СПУТНИКА
% Извлечение параметров
    xs_k = SatPos(4);
    ys_k = SatPos(5);
    i_k  = SatPos(6);
    Omega_k  = SatPos(7);
    dOmega_e = SatPos(8);

% Поправка долготы восходящего узла
    Omega_k_new = Omega_k - dOmega_e * TProp;

% Вычисление координат
    x_k = xs_k*cos(Omega_k_new) - ys_k*cos(i_k)*sin(Omega_k_new);
    y_k = xs_k*sin(Omega_k_new) + ys_k*cos(i_k)*cos(Omega_k_new);
    z_k = ys_k*sin(i_k);

% Результат
    OutSatPos = [x_k; y_k; z_k];
end