function [SatPos, GPSTime, TProp] = P72_GetSatPos(Data, inGPSTime, ...
    inTProp, Params) %#ok<INUSD>
%
% Фунция производит вычисление координат спутника в момент времени
% inGPSTime и учитывает время распространения сигнала inTProp при переводе
% координат в систему ECEF
%
% Входные переменные
%   Data - структура, содержащая, как минимум, параметры подкадров 1, 2 и
%     3;
%   inGPSTime - время испускания сигнала;
%   inTProp - время распространения сигнала.
%
% Выходные переменные
%   SatPos - массив (8х1) координат и параметров спутника для пересчёта
%       координат:
%         [x; y; z; ... % координаты в прямоугольной системе координат
%         xs_k; ys_k; i_k; ... % исходные координаты спутника
%         Omega_k; % исходное значение Omega_k
%         ZaZa]; % параметр для пересчёта Omega_k в Omega_k_TProp
%   GPSTime - скорректированное время испускания сигнала;
%   TProp - скорректированное время распространения сигнала.

%% УСТАНОВКА ВЕЛИЧИН
sqrtA   = Data.sqrtA;
t_oe    = Data.t_oe;
Delta_n = Data.Delta_n * pi;
M_0     = Data.M_0 * pi;
e       = Data.e;
omega   = Data.omega * pi;
i_0     = Data.i_0 * pi;
Omega_0 = Data.Omega_0 * pi;
DOmega  = Data.DOmega * pi;

IDOT    = Data.IDOT * pi;

%% УСТАНОВКА КОНСТАНТ
% Гравитационная постоянная, м^3/с^2
    nu = 3.986005 * 1e14;

% Угловая скорость вращения Земли, рад/с
    dOmega_e = 7.2921151467 * 1e-5;

% Константа, необходимая для поправки часов, с / sqrt(м)
    F = -4.442807633 * 1e-10;

%% ВЫЧИСЛЕНИЕ КООРДИНАТ
% Большая полуось
    A = sqrtA^2;
% Среднее движение
    n_0 = sqrt(nu / A^3);
% Время, прошедшее с момента, для которого были рассчитаны эфемериды
    t_k = inGPSTime - t_oe;
    % Корректировка получившегося значения
        if t_k > 302400
            t_k = t_k - 604800;

        elseif t_k < -302400
            t_k = t_k + 604800;
        end
% Скорректированное среднее движение
    n = n_0 + Delta_n;
% Средняя аномалия
    M_k = M_0 + n * t_k;
% Вычисление эксцентрической аномалии
    % Функция
        fun = @(Ek) Ek - e * sin(Ek) - M_k;
    % Начальное приближение
        EkInit = 10e-5;
    % Вычисление
        E_k = fzero(fun, EkInit);
% Истинная аномалия
    % Числитель и знаменатель аргумента арктангенса
        Num = sqrt(1 - e^2)*sin(E_k) / (1 - e * cos(E_k) );
        Den = (cos(E_k)-e) / (1 - e * cos(E_k) );

    mu_k = atan2(Num, Den);
% Аргумент широты
    F_k = mu_k + omega;
% Вторые гармоническое возмущения
    delta_u_k = Data.C_us*sin(2*F_k) + Data.C_uc*cos(2*F_k);
    delta_r_k = Data.C_rs*sin(2*F_k) + Data.C_rc*cos(2*F_k);
    delta_i_k = Data.C_is*sin(2*F_k) + Data.C_ic*cos(2*F_k);
% Скорректированный аргумент широты
    u_k = F_k + delta_u_k;
% Скорректированный радиус
    r_k = A*(1 - e*cos(E_k) ) + delta_r_k;
% Скорректированное наклонение
    i_k = i_0 + delta_i_k + IDOT * t_k;
% Координаты в орбитальной плоскости
    xs_k = r_k * cos(u_k);
    ys_k = r_k * sin(u_k);
% Скорректированная долгота восходящего угла
    Omega_k = Omega_0 + (DOmega - dOmega_e)*t_k - dOmega_e*t_oe;
% Координаты относительно Земли
    x_k = xs_k*cos(Omega_k) - ys_k*cos(i_k)*sin(Omega_k);
    y_k = xs_k*sin(Omega_k) + ys_k*cos(i_k)*cos(Omega_k);
    z_k = ys_k*sin(i_k);

% Присвоение результата
    SatPos = zeros(8, 1);
    SatPos(1:8) = [x_k y_k z_k xs_k ys_k i_k Omega_k dOmega_e].';

%% КОРРЕКТИРОВКА ЧАСОВ СПУТНИКА И ВРЕМЕНИ РАСПРОСТРАНЕНИЯ
% Релятивистская поправка
    Delta_t_r = F*Data.e*Data.sqrtA*sin(E_k);
% Поправка времени испускания
    Delta_t_sv = Data.a_f0 + Data.a_f1*(inGPSTime-Data.t_oc) + ...
        Data.a_f2*(inGPSTime-Data.t_oc)^2 + Delta_t_r - Data.T_GD;

% Корректировки
    GPSTime = inGPSTime - Delta_t_sv;
    TProp   = inTProp + Delta_t_sv;
end