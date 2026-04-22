function UPos = P71_GetOneRXPos(Es, inGPSTimes, inTimeShifts, ...
    SampleNums, Params)
%
% Функция расчёта одного набора координат приёмника
%
% Входные переменные
%   Es - cell-массив с эфемеридами спутников;
%   inGPSTimes - моменты времени, в которые были испущены сигналы со
%       спутников;
%   inTimeShifts - отличия значений задержки распространения сигналов
%       спутников от общей постоянной составляющей задержки
%       распространения;
%   SampleNums - номера отсчётов, в которые пришли сигналы спутников с
%       метками времени inGPSTimes.
%
% Выходные переменные
%   UPos - результат-структура с полями:
%       x, y, z - координаты в прямоугольной ДСК
%       T0 - общая составляющая сдвига по времени
%       tGPSs, SampleNums - значения времени GPS для заданных отсчётов 
%       Lat, Lon, Alt - широта, долгота, высота
%       SatsPoses - координаты спутников, таблица со столбцами
%           x, y, z - координаты в прямоугольной ДСК;
%           xs_k, ys_k, i_k - координаты перед преобразованиями СК;
%           Lat, Lon, Alt - широта, долгота, высота;
%           El, Az - угол склонения и азимут;
%       NumIters, MaxNumIters - выполненное и максимальное число итераций;
%       Delta, MaxDelta - достигнутое и максимальное значение оценки
%           изменения координат пользователя между соседними итерациями
%           (м);
%       inGPSTimes, GPSTimes, inTimeShifts, TimeShifts - сохранение
%           параметров и скорректированных параметров.

%% УСТАНОВКА ПАРАМЕТРОВ
% Максимальное число итераций
    MaxNumIters = Params.P71_GetOneRXPos.MaxNumIters;
% Максимальное изменение координат пользователя между соседними
% итерациями (м). Если фактическое изменение меньше, то цикл
% останавливается
    MaxDelta = Params.P71_GetOneRXPos.MaxDelta;

%% РАСЧЁТ ПАРАМЕТРОВ
% Число спутников, используемых для вычисления координат
    NumSats = length(Es);

%% УСТАНОВКА КОНСТАНТ
% Скорость света, м/с
    c = 299792458;
% Радиус Земли, м
    R = 6356863;
% Начальное значение общей составляющей времени распространения сигналов
    T0init = 68e-3;

%% ОСНОВНАЯ ЧАСТЬ
% Инициализация результата
    UPos = struct( ...
        'x', [], ...
        'y', [], ...
        'z', [], ...
        'T0',    [], ...
        'tGPSs', [], ...
        'SampleNums', SampleNums, ...
        'Lat',   [], ...
        'Lon',   [], ...
        'Alt',   [], ...
        'SatsPoses',    [], ...
        'NumIters',     [], ...
        'MaxNumIters',  MaxNumIters, ...
        'Delta',        [], ...
        'MaxDelta',     MaxDelta, ...
        'inGPSTimes',   inGPSTimes, ...
        'GPSTimes',     [], ...
        'inTimeShifts', inTimeShifts, ...
        'TimeShifts',   [] ...
        );

% Инициализация исходных значений перед итерационным вычислением координат
% приёмника
    % Координаты спутников и скорректированные параметры
        GPSTimes   = zeros(size(inGPSTimes  ) );
        TimeShifts = zeros(size(inTimeShifts) );
        SatPoses   = zeros(8, NumSats);

        % Цикл по спутникам
        for k = 1 : NumSats
            % Вычисление координат спутника и получение скорректированных
            % параметрв
                [SatPoses(:, k), GPSTimes(k), TProp] = P72_GetSatPos( ...
                    Es{k}, ...
                    inGPSTimes(k), ...
                    T0init + inTimeShifts(k), ...
                    Params);
            % Скорректированные задержки распространения сигналов
                TimeShifts(k) = TProp - T0init;
        end
    % Постоянная составляющая времени распространения
        T0 = T0init;
    % Координаты приёмника
        RXPos = [0, 0, 0]';

% Итерационное вычисление координат
    % Флаг остановки вычислений
        isStop = false;
    % Счётчик количества итераций
        IterCount = 0;

    while ~isStop
        % Текущие оценки значений задержек распространения сигналов
            TimeShiftsIter = zeros(1, NumSats);

            % Цикл по спутникам
            for sat = 1 : NumSats
                TimeShiftsIter(sat) = sqrt(sum( (RXPos-SatPoses(1:3, sat) ).^2) ) / c - T0;
            end

        % Значения матрицы B
            B = zeros(NumSats, 1);

            % Цикл по спутникам
            for sat = 1 : NumSats
                B(sat) = c * (TimeShifts(sat) - TimeShiftsIter(sat) );
            end

        % Значения матрицы А
            A = zeros(NumSats, 4);
            A(:, 4) = -c;

            % Циклы по строкам и столбцам матрицы
            for row = 1 : NumSats
            for col = 1 : 3
                A(row, col) = -(SatPoses(col, row) - RXPos(col) ) / c / ...
                    (TimeShiftsIter(sat) + T0);
            end
            end
            clear row col sat;

        % Обратная или псевдообратная матрица к A
            if NumSats == 4
                invA = inv(A);

            elseif NumSats > 4
                invA = pinv(A);

            else
                error(['%s Невозможно рассчитать координаты т.к. ' ...
                    'используется менее 4 спутников'], datetime('now') );
            end

        % Матрица X
            X = invA * B;

        % Оценка отклонения от истинного значения
            buf = [X(1), X(2), X(3), X(4)*c];
            Delta = sqrt(sum(buf.^2) );

        % Увеличим счётчик итераций
            IterCount = IterCount + 1;

        % Определим, нужно ли выполнять следующую итерацию
            isStop = (Delta < MaxDelta) || (IterCount >= MaxNumIters);

        if isStop
        % Если полученное значение нас устраивает, принимаем в качестве
        % решения координаты RXPos

        else
        % Иначе обновим значения переменных, необходимых для расчёта, после
        % чего выполним еще одну итерацию

            RXPos = RXPos + X(1:3);
            T0    = T0 + X(4);

            % Обновим координаты спутников
            for sat = 1 : NumSats
                SatPoses(1:3, sat) = P73_RenewSatPos(SatPoses(:, sat), ...
                    T0 + TimeShifts(sat), Params);
            end
        end
    end

% Преобразуем результат в сферическую систему координат
    [Lat, Lon, Alt] = P74_Cartesian2Spherical(RXPos, Params);

% Сохранение результатов
    UPos.x = RXPos(1);
    UPos.y = RXPos(2);
    UPos.z = RXPos(3);
    UPos.T0 = T0;

    UPos.Lat = Lat;
    UPos.Lon = Lon;
    UPos.Alt = Alt;

    UPos.NumIters = IterCount;
    UPos.MaxNumIters = MaxNumIters;
    UPos.Delta = Delta;
    UPos.MaxDelta = MaxDelta;

    UPos.GPSTimes   = GPSTimes;
    UPos.TimeShifts = TimeShifts;

% P75_CalculateSatElAz
end