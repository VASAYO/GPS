function Res = P30_NonCohDemodSatsSigs(inRes, Params)
%
% Функция некогерентной демодуляции сигналов спутников
%
% Входные переменные
%   inRes - структура с результатами модели, объявленная в Main;
%
% Выходные переменные
%   Res - структура, которая отличается от inRes добавлением нового поля,
%       описание которого дано ниже в коде.

% Пересохранение результатов
    Res = inRes;

%% ИНИЦИАЛИЗАЦИЯ РЕЗУЛЬТАТА
    Demod = struct( ...
        'Bits', {cell(Res.Search.NumSats, 1)} ...
    );
    % Каждый элемент cell-массива Bits - массив 4хN значений 0/1, где
    %   N - количество демодулированных бит.

%% УСТАНОВКА ПАРАМЕТРОВ
% Необходимость прорисовки результатов
    isDraw = Params.Main.isDraw;

%% РАСЧЁТ ПАРАМЕТРОВ
    % Количество периодов CA-кода, приходящихся на один бит
        CAPerBit = 20;

%% ОСНОВНАЯ ЧАСТЬ ФУНКЦИИ - ЦИКЛ ПО НАЙДЕННЫМ СПУТНИКАМ
for k = 1 : Res.Search.NumSats
    % Строка состояния
        fprintf('%s Некогерентная демодуляция бит спутника №%02d (%d из %d) ...\n', ...
            datetime("now"), Res.Search.SatNums(k), k,  Res.Search.NumSats);

    % Массив корреляций CA-кодов
        CorVals = Res.Track.CorVals{k};

    % Синхронизация с началом бита и отбрасывание оконечных корреляций
        CorVals = CorVals(Res.BitSync.CAShifts(k)+1 : end);
        CorVals = CorVals(1 : end-mod(length(CorVals), CAPerBit) );

    % Разность фаз корреляций, отстоящих на 20 позиций друг относительно
    % друга
        dCorVals = CorVals(1+CAPerBit : end) .* conj(CorVals(1 : end-CAPerBit) );

    % Накопление корреляций на длительности бита
        dCorVals20 = reshape(dCorVals, CAPerBit, [] );
        dCorVals20 = sum(dCorVals20, 1);

    % Вторая разность
        ddCorVals = dCorVals20(2:end) .* conj(dCorVals20(1:end-1) );

    % Вторая разность битовой последовательности
        ddBits = pskdemod(ddCorVals, 2, 0, "gray", "OutputType", "bit").';

    % Варианты первой разности битовой последовательности
        % Заготовка с предположением о первом бите
            dBits = [ [0 1]; repmat(ddBits, 1, 2) ];

        dBits = cumsum(dBits, 1);
        dBits = mod(dBits, 2);

    % Варианты битовой последовательности
        % Заготовка с предположением о первом бите
            Bits = [[0 1 1 0]; repmat(dBits, 1, 2) ];

        Bits = cumsum(Bits, 1);
        Bits = mod(Bits, 2).';

    % Отбрасывание последовательностей, которые являются копиями
    % существующих с точностью до инверсии
        Bits( [2 3], :) = [];

    % Сохранение результата в структуру
        Demod.Bits{k} = Bits;

    % Строка состояния
        fprintf('%s   Завершено.\n', datetime("now") );

    % Прорисовка результатов и сохранение рисунков
        if isDraw > 0
            figure( ...
                Name=['P30_NonCohDemod_SatNum', num2str(Res.Search.SatNums(k) ) ], ...
                WindowStyle="docked" ...
                );

            plot(ddCorVals, '.'); grid on; axis equal;

            title( ['Вторая разность сигнального созвездия для спутника № ', num2str(Res.Search.SatNums(k) ) ] );
            xlabel('I');
            ylabel('Q');
        end
        if isDraw > 1
            saveas(gcf, ...
                cat(2, ...
                    Params.Main.SaveDirName, '/', ...
                    'P30_NonCohDemod_SatNum_', num2str(Res.Search.SatNums(k) ) ...
                ) ...
            );
        end
        if isDraw > 2
            close(gcf);
        end
end

% Добавление поля в структуру результатов
    Res.Demod = Demod;
