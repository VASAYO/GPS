function Res = P60_GatherSatsEphemeris(inRes, Params) %#ok<INUSD>
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
% Ephemeris = cell(N, Res.Search.NumSats);

% Количество строк cell-массива Ephemeris совпадает с количеством
% подкадров спутников с одинаковым значением TOW (естественно, 
% учитываются только те спутники, у которых SatsData.isSat2Use = 1).
% Элементами Ephemeris являются структуры, содержащие значения всех
% параметров первого, второго и третьего подкадров, а также порядковый
% номер подкадра, номер первого CA подкадра, передаваемое значение TOW,
% для которого верна эта информация. Если полную информацию собрать не
% удалось, то элемент cell-массива должен быть пустым.

%% УСТАНОВКА ПАРАМЕТРОВ

%% РАСЧЁТ ПАРАМЕТРОВ
% Имена всех полей структур, являющихся элементами Ephemeris
    ENames = { ...
        ... % Эти поля не относятся к навигационной информации
        'SFNum', ... %   Порядковый номер подкадра спутника,
        ...          % соответствующего текущей строке (подкадру) Ephemeris
        'CANum', ... %   Номер CA-кода спутника, с которого начинается
        ...          % подкадр с порядковым номером SFNum
        'TOW', ...   %   Значение TOW, передаваемое в подкадре с порядковым
        ...          % номером SFNum. Это значение одинаковое для всех
        ...          % элементов одной строки Ephemeris
        ...
        ... % Поля с навигационной информацией
        'WeekNumber', ...
        'CodesOnL2', ...
        'URA', ...
        'URA_in_meters', ...
        'SV_Health_Summary', ...
        'SV_Health', ...
        'IODC', ...
        'L2_P_Data_Flag', ...
        'T_GD', ...
        't_oc', ...
        'a_f2', ...
        'a_f1', ...
        'a_f0', ...
        'IODE', ...
        'C_rs', ...
        'Delta_n', ...
        'M_0', ...
        'C_uc', ...
        'e', ...
        'C_us', ...
        'sqrtA', ...
        't_oe', ...
        'Fit_Interval_Flag', ...
        'AODO', ...
        'C_ic', ...
        'Omega_0', ...
        'C_is', ...
        'i_0', ...
        'C_rc', ...
        'omega', ...
        'DOmega', ...
        'IDOT', ...
        };

%% ОСНОВНАЯ ЧАСТЬ ФУНКЦИИ - ЦИКЛ ПО НАЙДЕННЫМ СПУТНИКАМ
% Строка состояния
    fprintf('%s Сбор эфемерид ...\n', datetime('now') );

% Определим порядковые номера спутников, для котрых мы будем пытаться
% собирать эфемериды
    SatsInds = [];
    for k = 1 : Res.Search.NumSats
        if Res.SatsData.isSat2Use(k), SatsInds(end+1) = k; end %#ok<AGROW>
    end

% Определим значения TOW, общие для всех спутников
    CommonTOWVals = [];

    % Цикл по спутникам
    for k = SatsInds
        % Вытащим все значения TOW подкадров спутника
            TOWVals = GetTOWVals(Res.SatsData.HOW, k);
        % Если среди значений TOW есть ошибочно декодированные, 
        % восстановим их
            TOWVals = RecoverTOWVals(TOWVals);

        % Определим общие значения для всех спутников
            if isempty(CommonTOWVals)
                CommonTOWVals = TOWVals;
            else
                CommonTOWVals = intersect(CommonTOWVals, TOWVals);
            end
    end

% Число значений TOW, общих для всех спутников
    NumCommonTOWVals = length(CommonTOWVals);

% Для каждого спутника определим порядковый номер подкадра, в котором
% встречается первое значение TOW общее с остальными спутниками
    FirstSFNums = zeros(1, Res.Search.NumSats);

    for k = SatsInds
        % Вытащим все значения TOW подкадров спутника
            TOWVals = GetTOWVals(Res.SatsData.HOW, k);
        % Если среди значений TOW есть ошибочно декодированные, 
        % восстановим их
            TOWVals = RecoverTOWVals(TOWVals);

        FirstSFNums(k) = find(TOWVals == CommonTOWVals(1) );
    end

% Заготовим результат
    Ephemeris = cell(NumCommonTOWVals, Res.Search.NumSats);
        
% Теперь попробуем для каждого подкадра каждого спутника определить
% эфемериды
for satIdx = SatsInds % Цикл по спутникам
    % Создадим заготовку, в которую мы будем собирать эфемериды
        E = MakeEmptyE(ENames);

    % Выделим результаты парсинга текущего спутника
        SFData = struct( ...
            'HOW', Res.SatsData.HOW{satIdx}, ...
            'SF1', Res.SatsData.SF1{satIdx}, ...
            'SF2', Res.SatsData.SF2{satIdx}, ...
            'SF3', Res.SatsData.SF3{satIdx} ...
            );

    for sfIdx = 1 : NumCommonTOWVals % Цикл по подкадрам
        % Попробуем собрать эфемериды, валидные для значения TOW текущего
        % подкадра
            E = CheckAndAddE(E, FirstSFNums(satIdx) + sfIdx -1, ...
                SFData, ENames);

        % Заполним поля, не относящиеся к навигационным данным
            E.SFNum = FirstSFNums(satIdx) + sfIdx -1;
            E.CANum = Res.BitSync.CAShifts(satIdx) + ...
                      Res.SubFrames.BitShift(satIdx) * 20 + ...
                      (E.SFNum-1) * 300 * 20 + 1;
            E.TOW = CommonTOWVals(sfIdx);
            
        % Если удалось полностью собрать структуру, сохраняем результат
            if CheckE(E, ENames)
                Ephemeris{sfIdx, satIdx} = E;
            end
    end
end

% Добавим новое поле с результатами в Res
    Res.Ephemeris = Ephemeris;

% Строка состояния
    fprintf('%s     Завершено.\n', datetime('now') );
end

%% Подфункции
function E = MakeEmptyE(ENames)
% Создадим все поля

    E = struct();
    for k = 1 : length(ENames)
        E.(ENames{k} ) = [];
    end

    % Установим в поля, не относящиеся к навигационным данным, произвольные
    % параметры, чтобы тест isGathered проходил успешно (см. CheckAndAddE)
        E.SFNum = -1;
        E.CANum = -1;
        E.TOW   = -1;
end

function [outE, isNew] = CheckAndAddE(inE, SFNum, SFData, ENames)
% Функция выполняет сбор эфемерид спутника для подкадра с порядковым
% номером в записи SFNum
% 
% inE   - структура с полями, хранящими эфемериды. Имена полей определяются
%         массивом ENames;
% SFNum - порядковый номер подкадра в записи;
% SFData- структура с результатами парсинга первых трёх подкадров и HOW
%         спутника, для которого выполняется сбор эфемерид (содержит поля 
%         HOW, SF1, SF2, SF3).
% ENames - cell-массив названий всех эфемерид, которые необходимо собрать.

    % SFID подкадра, для которого собираем эфемериды
        SFID = SFData.HOW(SFNum).Subframe_ID;

    % Сравниваем/обновляем значения IODC, IODE
        switch SFID 
            case 1
                % Варианты различных ситуаций:
                %   1. Поле inE.IODC НЕПУСТОЕ, значение
                %   SFData.SF1(SFNum).IODC НЕ РАВНО NaN. В таком случае
                %   сравниваем 8 младших бит этих значений. Далее в 
                %   зависимости от результата:
                %     - Значения РАВНЫ: эфемериды актуальны, тогда
                %     приравниваем IODC = SFData.SF1(SFNum).IODC,
                %     переприсваиваем inE на выход, isNew = 0;
                %     - Значения НЕ РАВНЫ: создаём пустую заготовку, в
                %     ней выставляем значение IODC равным 
                %     SFData.SF1(SFNum).IODC, значение IODE равным 
                %     mod(IODC, 256) и пробуем собрать эфемериды для этих 
                %     значений, isNew = 1;
                %   
                %   2. Поле inE.IODC НЕПУСТОЕ, значение
                %   SFData.SF1(SFNum).IODC РАВНО NaN. В таком случае,
                %   присваиваем inE на выход т.к. не можем определить 
                %   актуальность эфемерид и просто воспользуемся 
                %   предыдущими, isNew = 0;
                % 
                %   3. Поле inE.IODC ПУСТОЕ, значение SFData.SF1(SFNum).IODC 
                %   НЕ РАВНО NaN. В таком случае, создаём пустую заготовку,
                %   в ней выставляем значение IODC равным 
                %   SFData.SF1(SFNum).IODC, значение IODE равным 
                %   mod(IODC, 256) и пробуем собрать эфемериды. 
                %   isNew = 1;
                %   
                %   4. Поле inE.IODC ПУСТОЕ, значение SFData.SF1(SFNum).IODC 
                %   РАВНО NaN. В таком случае, пробуем найти ближайшее
                %   значение IODC либо IODE (в зависимости от того, что 
                %   встретится раньше) в предыдущих подкадрах (или в 
                %   следующих, если значения не были найдены в предыдущих).
                %   В зависимости от результатов поиска: 
                %     - Ни IODC, ни IODE не было найдено: мы никак не можем
                %     собрать валидные эфемериды, присваиваем на выход
                %     пустую заготовку. isNew неважно;
                %     - Нашли значение IODC: в пустой заготовке
                %     инициализируем значение IODC, значение IODE как 8
                %     младших бит IODC и пробуем собрать эфемериды;
                %     - Нашли значение IODE: в пустой заготовке
                %     инициализируем значение IODE, значение IODC как IODE 
                %     и пробуем собрать эфемериды;

                if ~isempty(inE.IODC) && ~isnan(SFData.SF1(SFNum).IODC)
                % inE.IODC непустое, SFData.SF1(SFNum).IODC - не NaN

                    % Сравнение значений
                        if isequal(mod(inE.IODC,               256), ...
                                   mod(SFData.SF1(SFNum).IODC, 256) )

                            inE.IODC = SFData.SF1(SFNum).IODC;
                            outE  = inE;
                            isNew = 0;
                            return;

                        else
                            outE = MakeEmptyE(ENames);
                            outE.IODC = SFData.SF1(SFNum).IODC;
                            outE.IODE = mod(SFData.SF1(SFNum).IODC, 256);
                            isNew = 1;
                        end

                elseif ~isempty(inE.IODC) && isnan(SFData.SF1(SFNum).IODC)
                % inE.IODC непустое, SFData.SF1(SFNum).IODC - NaN

                    outE  = inE;
                    isNew = 0;
                    return;

                elseif isempty(inE.IODC) && ~isnan(SFData.SF1(SFNum).IODC)
                % inE.IODC пустое, SFData.SF1(SFNum).IODC - не NaN

                    outE = MakeEmptyE(ENames);
                    outE.IODC = SFData.SF1(SFNum).IODC;
                    outE.IODE = mod(SFData.SF1(SFNum).IODC, 256);
                    isNew = 1;

                elseif isempty(inE.IODC) && isnan(SFData.SF1(SFNum).IODC)
                % inE.IODC пустое, SFData.SF1(SFNum).IODC - NaN

                    outE = inE;
                    isNew = 0;
                    return;
                    % todo: НЕОБХОДИМО РЕАЛИЗОВАТЬ В БУДУЩЕМ
                end

            case {2, 3}
                % Логика работы аналогична работе алгоритма для SFID = 1,
                % только теперь выполняется сравнение не младших, а всех
                % бит IODE

                % Взятие IODE
                    if SFID == 2
                        IODEParsed = SFData.SF2(SFNum).IODE;
                    else
                        IODEParsed = SFData.SF3(SFNum).IODE;
                    end

                if ~isempty(inE.IODE) && ~isnan(IODEParsed)
                % inE.IODE непустое, IODEParsed - не NaN

                    % Сравнение значений
                        if isequal(inE.IODE, IODEParsed)

                            outE  = inE;
                            isNew = 0;
                            return;

                        else
                            outE = MakeEmptyE(ENames);
                            outE.IODE = IODEParsed;
                            outE.IODC = IODEParsed;
                            isNew = 1;
                        end

                elseif ~isempty(inE.IODE) && isnan(IODEParsed)
                % inE.IODE непустое, IODEParsed - NaN

                    outE  = inE;
                    isNew = 0;
                    return;

                elseif isempty(inE.IODE) && ~isnan(IODEParsed)
                % inE.IODE пустое, IODEParsed - не NaN

                    outE = MakeEmptyE(ENames);
                    outE.IODE = IODEParsed;
                    outE.IODC = IODEParsed;
                    isNew = 1;

                elseif isempty(inE.IODE) && isnan(IODEParsed)
                % inE.IODE пустое, IODEParsed - NaN

                    outE = inE;
                    isNew = 0;
                    return;
                    % todo: НЕОБХОДИМО РЕАЛИЗОВАТЬ В БУДУЩЕМ
                end
                
            otherwise
                % Варианты различных ситуаций: 
                %   1. Значения IODC, IODE в inE НЕПУСТЫЕ: эфемериды
                %   актуальны. outE = inE, isNew = 0;
                %   2. Значения IODC, IODE в inE ПУСТЫЕ: действуем
                %   аналогично варианту 4 при SFID = 1;

                if ~isempty(inE.IODC) && ~isempty(inE.IODE)
                % Значения IODC, IODE в inE НЕПУСТЫЕ

                    outE  = inE;
                    isNew = 0;
                    return;

                else
                % Значения IODC, IODE в inE ПУСТЫЕ
                    
                    % Найдем IODC или IODE, для которых будем собирать
                    % эфемериды
                        sfIdxs = fliplr(1 : SFNum-1);
                        sfIdxs = [sfIdxs, SFNum+1 : length(SFData.HOW) ];

                        for k = sfIdxs
                            if SFData.HOW(k).Subframe_ID == 1 && ...
                               ~isnan(SFData.SF1(k).IODC)

                                outE = MakeEmptyE(ENames);
                                outE.IODC = SFData.SF1(k).IODC;
                                outE.IODE = mod(SFData.SF1(k).IODC, 256);
                                break;

                            elseif SFData.HOW(k).Subframe_ID == 2 && ...
                                   ~isnan(SFData.SF2(k).IODE)

                                outE = MakeEmptyE(ENames);
                                outE.IODE = SFData.SF2(k).IODE;
                                outE.IODC = outE.IODE;
                                break;

                            elseif SFData.HOW(k).Subframe_ID == 3 && ...
                                   ~isnan(SFData.SF3(k).IODE)

                                outE = MakeEmptyE(ENames);
                                outE.IODE = SFData.SF3(k).IODE;
                                outE.IODC = outE.IODE;
                                break;
                            end
                        end

                end
        end

        % Сбор эфемерид
            % Цикл по подкадрам
            for k = 1 : length(SFData.SF1)
                switch SFData.HOW(k).Subframe_ID
                    case 1
                        % Названия эфемерид, содержащихся в подкадре
                            SFFieldNames = fieldnames(SFData.SF1(1));

                        % Добавление эфемерид в результат
                        for valIdx = 1 : length(SFFieldNames)
                            if SFData.SF1(k).IODC == outE.IODC && ...
                               any(~isnan(SFData.SF1(k).(SFFieldNames{valIdx} ) ) )

                                outE.(SFFieldNames{valIdx} ) = ...
                                    SFData.SF1(k).(SFFieldNames{valIdx} );
                            end
                        end

                    case 2
                        % Названия эфемерид, содержащихся в подкадре
                            SFFieldNames = fieldnames(SFData.SF2(1));

                        % Добавление эфемерид в результат
                        for valIdx = 1 : length(SFFieldNames)
                            if SFData.SF2(k).IODE == outE.IODE && ...
                               any(~isnan(SFData.SF2(k).(SFFieldNames{valIdx} ) ) )

                                outE.(SFFieldNames{valIdx} ) = ...
                                    SFData.SF2(k).(SFFieldNames{valIdx} );
                            end
                        end

                    case 3
                        % Названия эфемерид, содержащихся в подкадре
                            SFFieldNames = fieldnames(SFData.SF3(1));

                        % Добавление эфемерид в результат
                        for valIdx = 1 : length(SFFieldNames)
                            if SFData.SF3(k).IODE == outE.IODE && ...
                               any(~isnan(SFData.SF3(k).(SFFieldNames{valIdx} ) ) )

                                outE.(SFFieldNames{valIdx} ) = ...
                                    SFData.SF3(k).(SFFieldNames{valIdx} );
                            end
                        end

                    otherwise
                        continue;
                end
            end
end

function isGathered = CheckE(E, ENames)
% Проверим, остались ли пустые поля

    % Инициализация результата
        isGathered = true;

    % Проверка каждого поля
    for k = 1 : length(ENames)
        if isempty(E.(ENames{k} ) )
            isGathered = false;
            return;
        end
    end
end

function TOWVals = GetTOWVals(HOW, k)
% Вытащим все значения TOW подкадров спутника
% 
%   HOW - поле Res.SatsData, содержащее результат парсинга TOW подкадров
%         всех спутников;
%   k   - порядковый номер спутника, значения которого необходимо получить.

    % Инициализация результата
        TOWVals = [];

    % Проверка наличия данных парсинга
        if isempty(HOW{k} ), return; end

    TOWVals = zeros(1, length(HOW{k} ) );
    for idx = 1 : length(TOWVals)
        TOWVals(idx) = HOW{k}(idx).TOW_Count_Message;
    end
end

function TOWVals = RecoverTOWVals(inTOWVals)
% Восстановим ошибочно декодированные значения TOW_Count_Message
%
%   inTOWVals - массив-строка значений TOW_Count_Message последовательно
%               следующих подкадров. Если значение не было успешно 
%               декодированно, то оно должно быть равно NaN.

    % Пересохраним результат
        TOWVals = inTOWVals;

    % Проверка, на наличие ошибочно декодированных значений
        if ~any(isnan(TOWVals) ), return; end

    % Первое успешно декодированное значение и его позиция
        Poses  = find(~isnan(TOWVals) );
        Pos    = Poses(1);
        RefTow = TOWVals(Pos);

    % Число элементов в массиве до и после RefTow
        N1 = Pos - 1;
        N2 = length(TOWVals) - Pos;
        
    % Восстановление значений
        Data1 = RefTow-N1 : RefTow-1;
        Data2 = RefTow+1 : RefTow+N2;

        TOWVals = [Data1, RefTow, Data2];
end