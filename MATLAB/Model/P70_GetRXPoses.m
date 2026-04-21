function Res = P70_GetRXPoses(inRes, Params)
%
% Функция сбора навигационной информации для спутников, у которых было
% найдено хотя бы одно значение TOW_Count_Message
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
% Positioning = struct(  ...
%   'RXPoses', cell(N, M), ...
%   'CAStep', CAStep, ...
%   'isCommonRxTime', isCommonRxTime ...
% );

% Количество строк N cell-массива RXPoses совпадает с количеством
%   подкадров спутников с одинаковыми значениями TOW. Количество
%   столбцов определяется количеством вычислений координат на
%   длительности одного подкадра и зависит от параметра CAStep,
%   определяемого ниже.
% CAStep - шаг в периодах CA-кода между соседними вычислениями
%   координат.
% isCommonRxTime - параметр, определяющий как вычислять параметры
%   спутников: в одно время приёмника или в разное.

%% УСТАНОВКА ПАРАМЕТРОВ
% Шаг в периодах CA-кода между соседними вычислениями координат. Всего
% в подкадре 6000 периодов CA-кода, поэтому, например, CAStep = 1000
% приведёт к вычислению 6 координат за один подкадр.
    CAStep = Params.P70_GetRXPoses.CAStep;
    
% Вариант вычисления координат.
% isCommonRxTime = 1 - координаты спутников вычисляются в одинаковый
%   момент  времени приёмника, соответствующий разным меткам
%   времени GPS
% isCommonRxTime = 0 - координаты спутников вычисляются в разные
%   моменты времени приёмника, соответствующие одинаковой метке
%   времени GPS
    isCommonRxTime = Params.P70_GetRXPoses.isCommonRxTime;

% Работа функции при isCommonRxTime = 0 пока не реализована
    if isCommonRxTime == 0
        error(['%s Работа функции при isCommonRxTime = 0 пока не ' ...
            'реализована'], datetime("now") );
    end
    
% Порядковые номера спутников, учитываемых при вычислении координат:
% 'all' - все спутники;
% 'firstX' - первые Х спутников, например 'first5';
% [1, 2, 5, 7] - конкретные номера.
    SatNums2Pos = Params.P70_GetRXPoses.SatNums2Pos;

% Поля структуры-результата
    Track = Res.Track;
    Ephemeris = Res.Ephemeris;

% Длительность одного периода CA-кода, с
    Tca = 1e-3;

%% РАСЧЁТ ПАРАМЕТРОВ
% Интервал дискретизации сигнала
    dt = 1/Res.File.Fs;

% Определим конкретные номера спутников
    if ischar(SatNums2Pos)
        if strcmp(SatNums2Pos, 'all')
            CurSatNums2Pos = 1:Res.Search.NumSats;
        else
            Buf = str2double(SatNums2Pos(6:end));
            CurSatNums2Pos = 1:Buf;
        end
    else
        CurSatNums2Pos = SatNums2Pos;
    end

% Выполним проверку спутников, выбранных для вычисления координат. Если для
% какого-либо спутника isSat2Use = 0, исключим его из списка
    isSat2Remove = false(size(CurSatNums2Pos) );
    for k = 1 : length(CurSatNums2Pos)
        if Res.SatsData.isSat2Use(CurSatNums2Pos(k) ) == 0
            isSat2Remove(k) = true;
        end
    end
    CurSatNums2Pos(isSat2Remove) = [];
    clear isSat2Remove;

% Число вычислений координат за длительность подкадра
    NumCalcsPerSF = floor(6000 / CAStep);

% Количество подкадров спутников с одинаковыми значениями TOW
    NumCommonSFs = size(Ephemeris, 1);

% Длина CA-кода с учётом частоты дискретизации
    CALen = 1023 * Res.File.R;
    
%% РАСЧЁТ КООРДИНАТ
% Строка состояния
    fprintf('%s Вычисление координат ...\n', datetime('now') );

% Инициализация результата
    Positioning = struct(  ...
      'RXPoses', {cell(NumCommonSFs, NumCalcsPerSF) }, ...
      'CAStep', CAStep, ...
      'isCommonRxTime', isCommonRxTime ...
    );

% Определим номера отсчётов сигнала, для которых будем вычислять координаты
    Poses2Calc = zeros(NumCommonSFs, NumCalcsPerSF);
    % Цикл по подкадрам
    for sfIdx = 1 : NumCommonSFs
        % Позиция первого отсчёта подкадра, принятого со спутника, первого
        % из списка спутников, используемых для вычисления координат
            SamplesShifts = Track.SamplesShifts{CurSatNums2Pos(1) };
            CANum = Ephemeris{sfIdx, CurSatNums2Pos(1) }.CANum;
            FirstPos = SamplesShifts(CANum) +1;

        Poses2Calc(sfIdx, :) = ...
            FirstPos + (0 : NumCalcsPerSF-1) * CAStep*CALen;
    end
    clear FirstPos SamplesShifts CANum;

% Вычислим координаты для каждого момента времени
for sfIdx = 1 : NumCommonSFs % Цикл по подкадрам с общим значением TOW для 
                             % всех спутников

    % Эфемериды спутников, использующиеся для вычисления координат
        Es = Ephemeris(sfIdx, CurSatNums2Pos);

    for tIdx = 1 : NumCalcsPerSF % Цикл по моментам времени, для которых 
                                 % вычисляюся координаты

        % Значение RefTOW
            RefTOW = (Es{1}.TOW - 1) * 6 + (tIdx-1) * CAStep * Tca;

        % Номер отсчёта, для которого вычисляются координаты
            Pos2Calc = Poses2Calc(sfIdx, tIdx);

        % Определим массив значений tGPS
            tGPSVals = zeros(1, length(CurSatNums2Pos) );

            % Цикл по спутникам
            for satIdx = 1 : length(CurSatNums2Pos)
                % Номера первых отсчётов CA-кодов спутника
                    SatCAPoses = Track.SamplesShifts{CurSatNums2Pos(satIdx) } +1;

                % Номер CA-кода, соответствующего значению RefTOW
                    RefTOWCANum = Es{CurSatNums2Pos(satIdx) }.CANum + ...
                                  (tIdx-1)*CAStep;

                % Вычислим tGPS спутника
                    tGPSVals(satIdx) = GettGPS(Pos2Calc, SatCAPoses, ...
                        RefTOWCANum, RefTOW, dt);
            end

        % Найдём отличия значений задержки распространения сигналов 
        % спутников от общей постоянной составляющей задержки 
        % распространения
            TimeShifts = tGPSVals - min(tGPSVals);

        % Номера отсчётов, в которые пришли сигналы спутников с метками 
        % времени inGPSTimes
            PosesRX = repelem(Pos2Calc, length(CurSatNums2Pos) );

        % Рассчитаем один набор координат
            Positioning.RXPoses{sfIdx, tIdx} = ...
                P71_GetOneRXPos(Es, tGPSVals, TimeShifts, PosesRX, Params);
    end
end

% Добавление нового поля в структуру-результат
    Res.Positioning = Positioning;

% Экспорт результатов в .kml файл
    P76_ExportResults(Positioning.RXPoses, Params);

% Строка состояния
    fprintf('%s     Завершено.\n', datetime('now') );
end

function tGPS = GettGPS(SampleNum, SamplesNums, RefCANum, RefTOW, dt)
%
% Функция определяет tGPS для отсчёта сигнала SampleNum
%
% Входные переменные
%   SampleNum - номер отсчёта записи, для которого надо расчитать время
%       GPS;
%   SamplesNums - номера первых отсчётов CA-кодов текущего спутника;
%   RefCANum - номер CA-кода, который является первым в подкадре, в котором
%       передаётся значение RefTOW;
%   RefTOW - значение RefTOW;
%   dt - интервал дискретизации записи.
%
% Выходные переменные
%   tGPS - время GPS в отсчёт SampleNum.

    % Номер первого отсчёта подкадра, для которого валидно значение RefTOW
        RefSampleNum = SamplesNums(RefCANum);
    
    % Разница между SampleNum и RefSampleNum
        DeltaSamples = SampleNum - RefSampleNum;
    
    tGPS = RefTOW - DeltaSamples * dt;
end