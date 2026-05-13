function Res = P10_CohSearchSats(inRes, Params)
%
% Функция когерентного поиска спутников в файле-записи
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
Search = struct( ...
    'NumSats',       0,  ... % Скаляр, количество найденных спутников
    'SatNums',       [], ... % массив 1хNumSats с номерами найденных
        ... % спутников
    'SamplesShifts', [], ... % массив 1хNumSats, каждый элемент -
        ... % количество отсчётов, которые нужно пропустить в файле-
        ... % записи до начала первого периода CA-кода соответствующего
        ... % спутника
    'FreqShifts',    [], ... % массив 1хNumSats со значениями частотных
        ... % сдвигов найденных спутников в Гц
    'CorVals',       [], ... % массив 1хNumSats вещественных значений
        ... % пиков корреляционных функций нормированных на среднее
        ... % значение, по которым были найдены спутники
    'AllCorVals',    zeros(1, 32) ... % массив максимальных значений
        ... % всех корреляционных функций
);

%% УСТАНОВКА ПАРАМЕТРОВ
% Количество периодов, учитываемых при обнаружении.
% Для когерентного обнаружения 1 <= NumCA2Search <= 10
    NumCA2Search = Params.P10_CohSearchSats.NumCA2Search;

% Массив центральных частот анализируемых диапазонов, Гц
    CentralFreqs = Params.P10_CohSearchSats.CentralFreqs;

% Порог обнаружения
    SearchThreshold = Params.P10_CohSearchSats.SearchThreshold;

% Флаг необходимости прорисовки результатов
    isDraw = Params.Main.isDraw;

%% СОХРАНЕНИЕ ПАРАМЕТРОВ
Search.NumCA2Search    = NumCA2Search;
Search.CentralFreqs    = CentralFreqs;
Search.SearchThreshold = SearchThreshold;

%% РАСЧЁТ ПАРАМЕТРОВ
% Количество рассматрвиаемых частотных диапазонов
    NumCFreqs = length(CentralFreqs);

% Длина CA-кода с учётом частоты дискретизации
    CALen = 1023 * Res.File.R;

% Число отсчётов записи, используемых для обнаружения
    NumSamples2Read = (2 * NumCA2Search + 1) * CALen - 1;

%% ОСНОВНАЯ ЧАСТЬ ФУНКЦИИ
% Лог
    fprintf('%s Когерентное обнаружение спутников...\n', datetime);

% Считаем сигнал из файла и разделим его на две части, чтобы построить два
% тела неопределённости
    Signal  = ReadSignalFromFile(Res.File, 0, NumSamples2Read);
    PartLen = (NumCA2Search + 1) * CALen - 1;
    Parts{1} = Signal(1:PartLen);
    Parts{2} = Signal(end-PartLen+1:end);

% Цикл по спутникам
for sat = 1 : 32
    % Лог
        fprintf('%s     Поиск спутника №%d: ', datetime, sat);

    % Банк опорных последовательностей для обнаружения
        refSeq = 1 - 2 * repelem(GenCACode(sat, NumCA2Search), Res.File.R);
        refLen = length(refSeq);

        refSeqs = repmat(refSeq, NumCFreqs, 1);
        F = repmat(CentralFreqs', 1, refLen);
        T = repmat( (0:refLen-1) / Res.File.Fs, NumCFreqs, 1);
        refSeqs = refSeqs .* exp(1j*2*pi * F .* T);

    % Построение двух тел неопределённости: для первой и второй пачки из 
    % NumCA2Search периодов CA-кода в сигнале
        CorBodies = {zeros(NumCFreqs, CALen), zeros(NumCFreqs, CALen)};

        % Цикл по частотным отстройкам
        for dfIdx = 1 : NumCFreqs
            % Цикл по двум частям сигнала
            for part = 1 : 2
                % Корреляция
                    buf = conv(Parts{part}, ...
                        fliplr(conj(refSeqs(dfIdx, :) ) ), "valid");
                % Когерентное накопление
                    bufAcc = reshape(buf, CALen, []);
                    bufAcc = abs(sum(bufAcc, 2) ).';

                CorBodies{part}(dfIdx, :) = bufAcc;
            end
        end

    % Решающие метрики об отсутствии или наличии сигнала
        Metrics = [0 0];
        for part = 1 : 2
            Metrics(part) = max(abs(CorBodies{part}), [], "all") / ...
                mean(abs(CorBodies{part}), "all");
        end
        [~, MaxMetPos] = max(Metrics);

    % Прорисовка результатов и сохранение рисунков
        if isDraw > 0
            figure(Name=['SatNum', num2str(sat) ] );
            surf(1:CALen, CentralFreqs, CorBodies{MaxMetPos});
            title('Тело неопределённости при некогерентном ', ...
                ['обнаружении спутника № ', num2str(sat) ] ...
            );
            xlabel('Частота');
            ylabel('Отсчёты');
            zlabel('Значение КФ');

            set(gcf, 'WindowStyle', 'docked');
        end
        if isDraw > 1
            saveas(gcf, ...
                cat(2, ...
                    Params.Main.SaveDirName, '/', ...
                    'P10_NonCohSearchSats_SatNum_', num2str(sat) ...
                ) ...
            );
        end
        if isDraw > 2
            close(gcf);
        end

    % Обновление структуры-результата
        Search.AllCorVals(sat) = Metrics(MaxMetPos);

    % Вынесение решения о наличии или отсутствии
        if any(Metrics >= SearchThreshold)
            fprintf('найден.\n');

        else
            fprintf('не найден.\n');
            continue;
        end

    % Определим временный и частотный сдвиг для тела неопределённости с
    % большим значением метрики
        [buf, MaxPosesByF] = max(CorBodies{MaxMetPos});
        [~,   MaxPosT    ] = max(buf);
        MaxPosF = MaxPosesByF(MaxPosT);

        SampleShift = MaxPosT - 1;
        FreqShift   = CentralFreqs(MaxPosF);

    % Обновление структуры-результата
        Search.NumSats              = Search.NumSats + 1;
        Search.SatNums      (end+1) = sat;
        Search.SamplesShifts(end+1) = SampleShift;
        Search.FreqShifts   (end+1) = FreqShift;
        Search.CorVals      (end+1) = Metrics(MaxMetPos);
end

% Ранжирование найденных спутников по убыванию мощности
    [Search.CorVals, indsSatsSorted] = sort(Search.CorVals, "descend");

    Search.SatNums       = Search.SatNums(indsSatsSorted);
    Search.SamplesShifts = Search.SamplesShifts(indsSatsSorted);
    Search.FreqShifts    = Search.FreqShifts(indsSatsSorted);

% Перезаписываем поле переменной Res
    Res.Search = Search;

% Лог
    fprintf('%s     Завершено.\n', datetime);
