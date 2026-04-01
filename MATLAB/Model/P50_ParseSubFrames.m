function Res = P50_ParseSubFrames(inRes, Params) %#ok<INUSD>
%
% Функция демодуляции сигналов спутников
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
    SatsData = struct( ...
        'isSat2Use', zeros(1, Res.Search.NumSats), ...
        'TLM', {cell(Res.Search.NumSats, 1)}, ...
        'HOW', {cell(Res.Search.NumSats, 1)}, ...
        'SF1', {cell(Res.Search.NumSats, 1)}, ...
        'SF2', {cell(Res.Search.NumSats, 1)}, ...
        'SF3', {cell(Res.Search.NumSats, 1)}, ...
        'SF4', {cell(Res.Search.NumSats, 1)}, ...
        'SF5', {cell(Res.Search.NumSats, 1)} ...
    );
    % Элементами всех cell-массивов (TLM, HOW, SF1, SF2, SF3, SF4, SF5)
    % являются структуры-массивы (1хN) с результатами парсинга, где N -
    % количество обработанных для спутника подкадров. Если какое то поле не
    % расшифровано из-за того, что не сошлось CRC, то его значение должно
    % быть установлено в nan. isSat2Use - массив флагов, указывающих,
    % было ли расшифровано хотя бы одно поле HOW.TOW_Count_Message, т.е.
    % имеет ли смысл в дальнейшем изучать содержимое подкадров (конечно,
    % isSat2Use = 0, если у этого спутника isSubFrameSync = 0).

%% УСТАНОВКА ПАРАМЕТРОВ

%% РАСЧЁТ ПАРАМЕТРОВ

%% ОСНОВНАЯ ЧАСТЬ ФУНКЦИИ - ЦИКЛ ПО НАЙДЕННЫМ СПУТНИКАМ С УСПЕШНОЙ
% ПОДКАДРОВОЙ СИНХРОНИЗАЦИЕЙ

% Строка состояния
    fprintf('%s Парсинг подкадров ...\n', datetime("now") );

for k = 1 : Res.Search.NumSats
    % Если для спутника не была выполнена подкадровая синхронизация,
    % пропустим его
        if ~Res.SubFrames.isSubFrameSync(k)
            continue;
        end

    % Выделим декодированные слова спутника
        Words = Res.SubFrames.Words{k};

    % Число подкадров
        NumSFs = size(Words, 1);

    % Парсинг слов TLM, HOW
        for sfIdx = 1 : NumSFs % Цикл по подкадрам
            % Парсинг TLM
                TLM = ParseTLM(Words{sfIdx, 1} );

            % Парсинг HOW
                HOW = ParseHOW(Words{sfIdx, 2} );

            % Если было успешно расшифровано HOW.TOW_Count_Message, то со 
            % спутником есть смысл работать дальше
                if ~isnan(HOW.TOW_Count_Message)
                    SatsData.isSat2Use(k) = 1;
                end
    
            % Присвоение результата
                SatsData.TLM{k}(end+1) = TLM;
                SatsData.HOW{k}(end+1) = HOW;
        end

    % Если ни разу не было успешно расшифровано поле HOW.TOW_Count_Message,
    % пропускаем спутник
        if ~SatsData.isSat2Use(k), continue; end

    % Массив номеров подкадров
        for sfIdx = 1 : NumSFs
            SFID = SatsData.HOW{k}(sfIdx).SubFrameID;

            % Если находим успешно декодированный номер подкадра, на его
            % основании определяем остальные номера
            if ~isnan(SFID)
                % Сколько номеров до и после найденного нужно восстановить
                    N1 = sfIdx-1;
                    N2 = NumSFs - sfIdx;

                % Восстановление остальных номеров
                    buf       = SFID - 1;
                    DataBe4   = buf-N1 : buf-1;
                    DataAfter = buf+1  : buf+N2;
                    SFIDs     = mod( [DataBe4, buf, DataAfter], 5) + 1;
                    break;
            end
        end

    % Парсинг остального содержимого подкадров
    for sfIdx = 1 : NumSFs % Цикл по подкадрам
        % Номер текущего подкадра
            SFID = SFIDs(sfIdx);

        % Объединим 3-10 слова подкадра в единую битовую последовательность
            [Bits, isCRC] = Words2BitFrame(Words(sfIdx, 3:10) );

        % По очереди выполним функции парсинга подкадров 1-5
        for funIdx = 1:5
            % Указатель на функцию парсинга
                Fun = str2func(['ParseSF', num2str(funIdx) ] );
            % Вызов функции
                ParseRes = Fun(Bits, isCRC, SFID);

            % Сохранение результата парсинга в структуру
                FieldName = ['SF', num2str(funIdx) ];
                SatsData.(FieldName){k}(end+1) = ParseRes;
        end
    end
end

% Присваивание структуре с результатами нового поля
    Res.SatsData = SatsData;

% Строка состояния
    fprintf('%s     Завершено.\n', datetime("now") );
end

function [Bits, isCRC] = Words2BitFrame(Words)
% Из (1х8) cell-массива Words составим кадр, т.е. добавим нулевые биты CRC
% и нулевые первые два слова. Это удобно для анализа кода по спецификации.
% Также составим массив флагов, указывающих на то, сошлось CRC в конкретном
% слове или нет

    % Инициализация резульата
        Bits  = zeros(1, 60);
        isCRC = zeros(1, 8);

    % Цикл по словам
        for k = 1 : 8
            % Проверка CRC
                if isempty(Words{k} ) % Не сошлось
                    buf = zeros(1, 30);

                else % Сошлось
                    buf = [Words{k}, zeros(1, 6) ];
                    isCRC(k) = 1;
                end

            Bits = [Bits, buf]; %#ok<AGROW>
        end
end

function SF = ParseSF1(Words, isCRC, SFID)
% Парсинг подкадра №1
%
% Входные аргументы:
%   Words - 1x300 массив бит подкадра с номером SFID. В массиве биты
%           чётности, а также биты первых двух слов равны нулю. Если для
%           слова не сошлось CRC, то позиции, соответствующие слову, 
%           заполнены нулями;
%   isCRC - 1х8 массив флагов, показывающих, сошлось ли CRC для слов №3-10 
%           подкадра;
%   SFID  - номер подкадра от 1 до 5.
%
% Выходные аргументы:
%   SF - структура, описание которой дано ниже. Если расчёт какого либо
%         поля невозможен потому что не сошлось CRC, необходимо установить
%         соответствующее значение в NaN.

    % Инициализация результата
        SF = struct( ...
            'WeekNumber',        [], ...
            'CodesOnL2',         [], ...
            'URA',               [], ...
            'URA_in_meters',     [], ...
            'SV_Health_Summary', [], ...
            'SV_Health',         [], ...
            'IODC',              [], ...
            'L2_P_Data_Flag',    [], ...
            'T_GD', [], ...
            't_oc', [], ...
            'a_f2', [], ...
            'a_f1', [], ...
            'a_f0', []  ...
            );

    % Если Words содержит данные другого подкадра, завершаем выполнение
        if SFID ~= 1, return; end

    % Парсинг
        % WeekNumber
            if isCRC(3 -2)
                buf = Words(61 + (0:10-1) );
                SF.WeekNumber = comp2de(buf, false) * 1;
            else
                SF.WeekNumber = NaN;
            end
        % CodesOnL2
            if isCRC(3 -2)
                buf = comp2de(Words(71 + (0:2-1) ), false);
                switch buf
                    case 0
                        SF.CodesOnL2 = 'Reserved';
                    case 1
                        SF.CodesOnL2 = 'P code ON';
                    case 2
                        SF.CodesOnL2 = 'C/A code ON';
                    otherwise
                        SF.CodesOnL2 = 3;
                end
            else
                SF.CodesOnL2 = NaN;
            end
        % URA
            if isCRC(3 -2)
                buf = Words(73 + (0:4-1) );
                SF.URA = comp2de(buf, false) * 1;
            else
                SF.URA = NaN;
            end
        % URA_in_meters
            if ~isnan(SF.URA)
                URAmetVals = [2 2.9 4.13 5.85 8.25 11.65 18.83 36 72 ...
                    144 288 576 1152 2300 4600 6144];

                SF.URA_in_meters = URAmetVals(SF.URA +1);
            else
                SF.URA_in_meters = NaN;
            end
        % SV_Health
            if isCRC(3 -2)
                buf = Words(77 + (0:6-1) );
                buf = comp2de(buf, false);
                if buf == 0
                    SF.SV_Health = 'All Signals OK';
                else
                    SF.SV_Health = buf;
                end
            else
                SF.SV_Health = NaN;
            end
        % SV_Health_Summary
            if ~isnan(SF.SV_Health)
                if isequal(SF.SV_Health, 'All Signals OK')
                    SF.SV_Health_Summary = 'All NAV data are OK';
                else
                    SF.SV_Health_Summary = 'Some or all NAV data are bad';
                end
            else
                SF.SV_Health_Summary = NaN;
            end
        % IODC
            if isCRC(3 -2) && isCRC(8 -2)
                buf = [Words(83 + (0:2-1) ), Words(211 + (0:8-1) )];

                SF.IODC = comp2de(buf, false);
            else
                SF.IODC = NaN;
            end
        % L2_P_Data_Flag
            if isCRC(4 -2)
                buf = Words(91);
                if buf == 1
                    SF.L2_P_Data_Flag = 'Data OFF on the L2 P-code';
                else
                    SF.L2_P_Data_Flag = 'Data ON on the L2 P-code';
                end
            else
                SF.L2_P_Data_Flag = NaN;
            end
        % T_GD
            if isCRC(7 -2)
                buf = Words(197 + (0:8-1) );
                SF.T_GD = comp2de(buf, true) * 2^-31;
            else
                SF.T_GD = NaN;
            end
        % t_oc
            if isCRC(8 -2)
                buf = Words(219 + (0:16-1) );
                SF.t_oc = comp2de(buf, false) * 2^4;
            else
                SF.t_oc = NaN;
            end
        % a_f2
            if isCRC(9 -2)
                buf = Words(241 + (0:8-1) );
                SF.a_f2 = comp2de(buf, true) * 2^-55;
            else
                SF.a_f2 = NaN;
            end
        % a_f1
            if isCRC(9 -2)
                buf = Words(249 + (0:16-1) );
                SF.a_f1 = comp2de(buf, true) * 2^-43;
            else
                SF.a_f1 = NaN;
            end
        % a_f0
            if isCRC(10 -2)
                buf = Words(271 + (0:22-1) );
                SF.a_f0 = comp2de(buf, true) * 2^-31;
            else
                SF.a_f0 = NaN;
            end
end

function SF = ParseSF2(Words, isCRC, SFID)
% Парсинг подкадра №2
%
% Входные аргументы:
%   Words - 1x300 массив бит подкадра с номером SFID. В массиве биты
%           чётности, а также биты первых двух слов равны нулю. Если для
%           слова не сошлось CRC, то позиции, соответствующие слову, 
%           заполнены нулями;
%   isCRC - 1х8 массив флагов, показывающих, сошлось ли CRC для слов №3-10 
%           подкадра;
%   SFID  - номер подкадра от 1 до 5.
%
% Выходные аргументы:
%   SF - структура, описание которой дано ниже. Если расчёт какого либо
%         поля невозможен потому что не сошлось CRC, необходимо установить
%         соответствующее значение в NaN.

    % Инициализация результата
        SF = struct( ...
            'IODE',    [], ...
            'C_rs',    [], ...
            'Delta_n', [], ...
            'M_0',     [], ...
            'C_uc',    [], ...
            'e',       [], ...
            'C_us',    [], ...
            'sqrtA',   [], ...
            't_oe',    [], ...
            'Fit_Interval_Flag', [], ...
            'AODO',    [] ...
            );

    % Если Words содержит данные другого подкадра, завершаем выполнение
        if SFID ~= 2, return; end

    % Парсинг
        % IODE
            if isCRC(3 -2)
                SF.IODE = comp2de(Words(61 : 61+8 -1), false);
            else
                SF.IODE = NaN;
            end
        % C_rs
            if isCRC(3 -2)
                buf = comp2de(Words(69 : 69+16 -1), true);
                SF.C_rs = buf * 2^-5;
            else
                SF.C_rs = NaN;
            end
        % Delta_n
            if isCRC(4 -2)
                buf = comp2de(Words(91 : 91+16 -1), true);
                SF.Delta_n = buf * 2^-43;
            else
                SF.Delta_n = NaN;
            end
        % M_0
            if isCRC(4 -2) && isCRC(5 -2)
                buf = [Words(107:107+8 -1), Words(121:121+24 -1) ];
                SF.M_0 = comp2de(buf, true) * 2^-31;
            else
                SF.M_0 = NaN;
            end
        % C_uc
            if isCRC(6 -2)
                buf = Words(151:151+16 -1);
                SF.C_uc = comp2de(buf, true) * 2^-29;
            else
                SF.C_uc = NaN;
            end
        % e
            if isCRC(6 -2) && isCRC(7 -2)
                buf = [Words(167:167+8 -1), Words(181:181+24 -1) ];
                SF.e = comp2de(buf, false) * 2^-33;
            else
                SF.e = NaN;
            end
        % C_us
            if isCRC(8 -2)
                buf = Words(211:211+16 -1);
                SF.C_us = comp2de(buf, true) * 2^-29;
            else
                SF.C_us = NaN;
            end
        % sqrtA
            if isCRC(8 -2) && isCRC(9 -2)
                buf = [Words(227:227+8 -1), Words(241:241+24 -1) ];
                SF.sqrtA = comp2de(buf) * 2^-19;
            else
                SF.sqrtA = NaN;
            end
        % t_oe
            if isCRC(10 -2)
                buf = Words(271:271+16 -1);
                SF.t_oe = comp2de(buf, false) * 2^4;
            else
                SF.t_oe = NaN;
            end
        % Fit_Interval_Flag
            if isCRC(10 -2)
                SF.Fit_Interval_Flag = Words(287);
            else
                SF.Fit_Interval_Flag = NaN;
            end
        % AODO
            if isCRC(10 -2)
                SF.AODO = comp2de(Words(288:288+5-1), false) * 900;
            else
                SF.AODO = NaN;
            end
end

function SF = ParseSF3(Words, isCRC, SFID)
% Парсинг подкадра №3
%
% Входные аргументы:
%   Words - 1x300 массив бит подкадра с номером SFID. В массиве биты
%           чётности, а также биты первых двух слов равны нулю. Если для
%           слова не сошлось CRC, то позиции, соответствующие слову, 
%           заполнены нулями;
%   isCRC - 1х8 массив флагов, показывающих, сошлось ли CRC для слов №3-10 
%           подкадра;
%   SFID  - номер подкадра от 1 до 5.
%
% Выходные аргументы:
%   SF - структура, описание которой дано ниже. Если расчёт какого либо
%         поля невозможен потому что не сошлось CRC, необходимо установить
%         соответствующее значение в NaN.

    % Инициализация результата
        SF = struct( ...
            'C_ic',    [], ...
            'Omega_0', [], ...
            'C_is',    [], ...
            'i_0',     [], ...
            'C_rc',    [], ...
            'omega',   [], ...
            'DOmega',  [], ...
            'IODE',    [], ...
            'IDOT',    [] ...
            );

    % Если Words содержит данные другого подкадра, завершаем выполнение
        if SFID ~= 3, return; end

    % Парсинг
        % C_ic
            if isCRC(3 -2)
                buf = Words(61 + (0:16 -1) );
                SF.C_ic = comp2de(buf, true) * 2^-29;
            else
                SF.C_ic = NaN;
            end
        % Omega_0
            if isCRC(3 -2) && isCRC(4 -2)
                buf = [Words(77 + (0:8 -1) ), Words(91 + (0:24 -1) ) ];
                SF.Omega_0 = comp2de(buf, true) * 2^-31;
            else
                SF.Omega_0 = NaN;
            end
        % C_is
            if isCRC(5 -2)
                buf = Words(121 + (0:16 -1) );
                SF.C_is = comp2de(buf, true) * 2^-29;
            else
                SF.C_is = NaN;
            end
        % i_0
            if isCRC(5 -2) && isCRC(6 -2)
                buf = [Words(137 + (0:8 -1) ), Words(151 + (0:24 -1) ) ];
                SF.i_0 = comp2de(buf, true) * 2^-31;
            else
                SF.i_0 = NaN;
            end
        % C_rc
            if isCRC(7 -2)
                buf = Words(181 + (0:16 -1) );
                SF.C_rc = comp2de(buf, true) * 2^-5;
            else
                SF.C_rc = NaN;
            end
        % omega
            if isCRC(7 -2) && isCRC(8 -2)
                buf = [Words(197 + (0:8 -1) ), Words(211 + (0:24 -1) ) ];
                SF.omega = comp2de(buf, true) * 2^-31;
            else
                SF.omega = NaN;
            end
        % DOmega
            if isCRC(9 -2)
                buf = Words(241 + (0:24 -1) );
                SF.DOmega = comp2de(buf, true) * 2^-43;
            else
                SF.DOmega = NaN;
            end
        % IODE
            if isCRC(10 -2)
                buf = Words(271 + (0:8 -1) );
                SF.IODE = comp2de(buf, false) * 1;
            else
                SF.IODE = NaN;
            end
        % IDOT
            if isCRC(10 -2)
                buf = Words(279 + (0:14 -1) );
                SF.IDOT = comp2de(buf, true) * 2^-43;
            else
                SF.IDOT = NaN;
            end
end

function SF = ParseSF4(Words, isCRC, SFID) %#ok<INUSD>
% Парсинг подкадра №4 - реализован только для (SV_Page_ID = 56)
%
% Входные аргументы:
%   Words - 1x300 массив бит подкадра с номером SFID. В массиве биты
%           чётности, а также биты первых двух слов равны нулю. Если для
%           слова не сошлось CRC, то позиции, соответствующие слову, 
%           заполнены нулями;
%   isCRC - 1х8 массив флагов, показывающих, сошлось ли CRC для слов №3-10 
%           подкадра;
%   SFID  - номер подкадра от 1 до 5.
%
% Выходные аргументы:
%   SF - структура, описание которой дано ниже. Если расчёт какого либо
%         поля невозможен потому что не сошлось CRC, необходимо установить
%         соответствующее значение в NaN.

    % Инициализация результата
        SF = struct( ...
            'Data_ID',    [], ...
            'SV_Page_ID', [], ...
            'alpha_0',    [], ...
            'alpha_1',    [], ...
            'alpha_2',    [], ...
            'alpha_3',    [], ...
            'beta_0',  [], ...
            'beta_1',  [], ...
            'beta_2',  [], ...
            'beta_3',  [], ...
            'A_1',     [], ...
            'A_0',     [], ...
            't_ot',    [], ...
            'WN_t',    [], ...
            'Delta_t_LS',  [], ...
            'WN_LSF',      [], ...
            'DN',          [], ...
            'Delta_t_LSF', []  ...
            );

    % Если Words содержит данные другого подкадра, завершаем выполнение
        if SFID ~= 4, return; end

    % Парсинг
        ...
end

function SF = ParseSF5(Words, isCRC, SFID) %#ok<INUSD>
% Парсинг подкадра №5
%
% Входные аргументы:
%   Words - 1x300 массив бит подкадра с номером SFID. В массиве биты
%           чётности, а также биты первых двух слов равны нулю. Если для
%           слова не сошлось CRC, то позиции, соответствующие слову, 
%           заполнены нулями;
%   isCRC - 1х8 массив флагов, показывающих, сошлось ли CRC для слов №3-10 
%           подкадра;
%   SFID  - номер подкадра от 1 до 5.
%
% Выходные аргументы:
%   SF - структура, описание которой дано ниже. Если расчёт какого либо
%         поля невозможен потому что не сошлось CRC, необходимо установить
%         соответствующее значение в NaN.

    % Инициализация результата
        SF = struct('SF5_parsing_to_be_done', [] );

    % Если Words содержит данные другого подкадра, завершаем выполнение
        if SFID ~= 5, return; end

    % Тут должен находиться парсинг
        SF.SF5_parsing_to_be_done = 1;
end

function TLM = ParseTLM(Word)
% Парсинг слова TLM
% 
% Word - 1х24 массив бит первого слова подкадра.

    % Инициализация результата
        TLM = struct( ...
            'Preamble',                  NaN, ...
            'TLM_Message',               NaN, ...
            'TLM_Integrity_Status_Flag', NaN  ...
            );

    % Если для данного слова не сошелся CRC, прекращаем парсинг
        if isempty(Word), return; end

    % Парсинг Preamble
        TLM.Preamble = comp2de(Word(1:8), false);

    % Парсинг TLM_Message
        TLM.TLM_Message = comp2de(Word(9:22), false);

    % Парсинг TLM_Integrity_Status_Flag
        TLM.TLM_Integrity_Status_Flag = Word(23);
end

function HOW = ParseHOW(Word)
% Парсинг слова HOW
% 
% Word - 1х24 массив бит второго слова подкадра.

    % Инициализация результата
        HOW = struct( ...
            'TOW_Count_Message', NaN, ...
            'Alert_Flag',        NaN, ...
            'Anti_Spoof_Flag',   NaN, ...
            'SubFrameID',        NaN  ...
            );

    % Если для данного слова не сошелся CRC, прекращаем парсинг
        if isempty(Word), return; end

    % Парсинг TOW_Count_Message
        HOW.TOW_Count_Message = comp2de(Word(1:17), false);

    % Парсинг Alert_Flag
        HOW.Alert_Flag = Word(18);

    % Парсинг Anti_Spoof_Flag
        HOW.Anti_Spoof_Flag = Word(19);

    % Парсинг SubFrameID
        HOW.SubFrameID = comp2de(Word(20:22), false);
end
                
function Out = comp2de(In, isSigned)
% Функция перевода двоичного дополнительного кода в десятичное число
%
% In       - массив-строка бит. In(1) соответствует MSB;
% isSigned - флаг, считается ли In знаковым числом;
% Out      - десятичное целое знаковое число.

    arguments
        In       double;
        isSigned logical = false;
    end

    % Число разрядов числа
        L = length(In);
    
    Out = sum(In .* 2.^fliplr(0 : L-1) );
    % Проверка знака числа
        if In(1) == 1 && isSigned
            Out = -2^L + Out;
        end
end
