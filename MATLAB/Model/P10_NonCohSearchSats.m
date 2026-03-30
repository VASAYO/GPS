function Res = P10_NonCohSearchSats(inRes, Params)
%
% Функция некогерентного поиска спутников в файле-записи
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
        'NumSats',        0, ... % Скаляр, количество найденных спутников
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
        NumCA2Search = Params.P10_NonCohSearchSats.NumCA2Search;

    % Массив центральных частот анализируемых диапазонов, Гц
        CentralFreqs = Params.P10_NonCohSearchSats.CentralFreqs;

    % Порог обнаружения
        SearchThreshold = Params.P10_NonCohSearchSats.SearchThreshold;

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

    % Число первых отсчётов записи, используемых для обнаружения
        NumFirstSigSamples = (NumCA2Search+1) * CALen - 1;

%% ОСНОВНАЯ ЧАСТЬ ФУНКЦИИ
    % Лог
        fprintf('%s Некогерентное обнаружение спутников...\n', datetime("now") );

    % Считываем сигнал из файла
        Signal = ReadSignalFromFile(Res.File, 0, NumFirstSigSamples).';

    % Цикл по C/A кодам
        for caIdx = 1 : 32
            % Лог
                fprintf('%s     Поиск спутника №%d: ', ...
                    datetime("now"), caIdx);

            % Генерация набора опорных последовательностей для обнаружения
                CACode1023 = 1 - 2*GenCACode(caIdx, 1).';
                CACodeR    = repelem(CACode1023, Res.File.R);

                refSeqs    = repmat(CACodeR, 1, NumCFreqs);

                expF = repmat(CentralFreqs, CALen, 1);
                expT = (0 : CALen-1)' / Res.File.Fs;
                expT = repmat(expT, 1, NumCFreqs);

                refSeqs = refSeqs .* ...
                    exp(1j*2*pi*expF .* expT);

            % Построение тела неопределённости
                CorrBody = zeros(CALen, NumCFreqs);

                for fIdx = 1 : NumCFreqs
                    % Корреляция сигнала с опорной последовательностью
                        buf = conv(Signal, flipud(conj(refSeqs(:, fIdx) ) ), ...
                            "valid");

                    % Некогерентное накопление
                        buf = reshape(buf, CALen, []);
                        CorrBody(:, fIdx) = sum(abs(buf), 2);
                end

            % Вынесение решения о наличии или отсутствии спутника
                % Решающая метрика
                    Peak2Aver = max(abs(CorrBody), [], "all") / ...
                        mean(abs(CorrBody), "all");

                % Сравнение значения метрики с пороговым
                    isSatFound = (Peak2Aver >= SearchThreshold);

            % Лог
                if isSatFound
                    fprintf('найден.\n');
                else
                    fprintf('не найден.\n');
                end

            % Обновляем структуру Search
                if isSatFound
                    Search.NumSats = Search.NumSats + 1;
                    Search.SatNums(end+1) = caIdx;

                    % Определение временного и частотного сдвигов
                        [buf, indsMaxDimT] = max(abs(CorrBody) );
                        [~,   indMaxDimF ] = max(buf);
                        indMaxDimT = indsMaxDimT(indMaxDimF);

                        Search.SamplesShifts(end+1) = indMaxDimT - 1;
                        Search.FreqShifts   (end+1) = ...
                            CentralFreqs(indMaxDimF);

                        Search.CorVals(end+1) = Peak2Aver;
                end

                Search.AllCorVals(caIdx) = Peak2Aver;

            % Прорисовка результатов и сохранение рисунков
                if isDraw > 0
                    figure(Name=['SatNum', num2str(caIdx) ] );
                    surf(CentralFreqs,  (1:CALen)', CorrBody);
                    title('Тело неопределённости при некогерентном ', ...
                        ['обнаружении спутника № ', num2str(caIdx) ] ...
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
                            'P10_NonCohSearchSats_SatNum_', num2str(caIdx) ...
                        ) ...
                    );
                end
                if isDraw > 2
                    close(gcf);
                end
        end

        % Ранжирование найденных спутников по убыванию их мощности в записи
            [Search.CorVals, indsSatsSorted] = sort(Search.CorVals, "descend");

            Search.SatNums       = Search.SatNums(indsSatsSorted);
            Search.SamplesShifts = Search.SamplesShifts(indsSatsSorted);
            Search.FreqShifts    = Search.FreqShifts(indsSatsSorted);

        % Перезаписываем поле переменной Res
            Res.Search = Search;

        % Лог
            fprintf('%s     Завершено.\n', datetime("now") );
end
