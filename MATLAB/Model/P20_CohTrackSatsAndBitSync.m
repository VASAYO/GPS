  function Res = P20_CohTrackSatsAndBitSync(inRes, Params)
%
% Функция когерентного трекинга спутников и битовой синхронизации
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
    Track = struct( ...
        'SamplesShifts',     {cell(Res.Search.NumSats, 1)}, ... 
        'CorVals',           {cell(Res.Search.NumSats, 1)}, ...
        'HardSamplesShifts', {cell(Res.Search.NumSats, 1)}, ... 
        'FineSamplesShifts', {cell(Res.Search.NumSats, 1)}, ... 
        'EPLCorVals',        {cell(Res.Search.NumSats, 1)}, ...
        'DLL',               {cell(Res.Search.NumSats, 1)}, ...
        'FPLL',              {cell(Res.Search.NumSats, 1)} ...
    );
    % Каждая ячейка cell-массивов SamplesShifts, CorVals, HardSamplesShifts
    %   FineSamplesShifts является массивом 1xN, где N - количество
    %   периодов CA-кода соответствующего спутника, найденных в файле-
    %   записи (N может быть разным для разных спутников).
    % Каждый элемент массива SamplesShifts{k} - дробное количество
    %   отсчётов, которые надо пропустить в файле-записи до начала
    %   соответствующего периода CA-кода.
    % Каждый элемент массива CorVals{k} - комплексное значение корреляции
    %   части сигнала, содержащей соответствующий период CA-кода, с опорным
    %   сигналом.
    % Каждый элемент массивов HardSamplesShifts{k}, FineSamplesShifts{k} -
    %   соответственно дробная и целая части значений SamplesShifts{k}.
    % Каждая ячейка cell-массива EPLCorVals является массивом 3xN значений
    %   Early, Promt и Late корреляций. При этом: SamplesShifts{k} =
    %   EPLCorVals{k}(2, :).
    % DLL, FPLL - лог сопровождения фазы кода и частоты-фазы сигнала.

    BitSync = struct( ...
        'CAShifts', zeros(Res.Search.NumSats, 1), ... 
        'Cors', zeros(Res.Search.NumSats, 20) ...
    );
    % Каждый элемент массива CAShifts - количество периодов CA-кода,
    %   которые надо пропустить до начала бита.
    % Каждая строка массива Cors - корреляции, по позиции минимума которых
    %   определяется битовая синхронизация.

%% УСТАНОВКА ПАРАМЕТРОВ
    % Порядок фильтров
        DLL.FilterOrder = Params.P20_CohTrackSatsAndBitSync.DLL.FilterOrder;
        FPLL.FilterOrder = Params.P20_CohTrackSatsAndBitSync.FPLL.FilterOrder;
        
    % И DLL и FPLL имеют несколько режимов работы для каждого из них нужно
    % определить
        % Полосы фильтров
            DLL.FilterBands  = Params.P20_CohTrackSatsAndBitSync.DLL.FilterBands;
            FPLL.FilterBands = Params.P20_CohTrackSatsAndBitSync.FPLL.FilterBands;
            
        % Количество периодов накопления для фильтрации
            DLL.NumsIntCA  = Params.P20_CohTrackSatsAndBitSync.DLL.NumsIntCA;
            FPLL.NumsIntCA = Params.P20_CohTrackSatsAndBitSync.FPLL.NumsIntCA;

	% Определим количество периодов CA-кода, учитываемых для проверки
	% необходимости перехода между состояниями DLL и FPLL. Проверка
	% работает по принципу integrate and dump
        DLL.NumsCA2CheckState  = Params.P20_CohTrackSatsAndBitSync.DLL.NumsCA2CheckState;
        FPLL.NumsCA2CheckState = Params.P20_CohTrackSatsAndBitSync.FPLL.NumsCA2CheckState;
        
    % Граничные значения для перехода между состояниями
    % Если значение > HiTr, то переходим в следующее (более робастное)
    %   состояние
    % Если значение < LoTr, то переходим в предыдущее (более
    %   чувствительное)состояние
        DLL.HiTr = Params.P20_CohTrackSatsAndBitSync.DLL.HiTr;
        DLL.LoTr = Params.P20_CohTrackSatsAndBitSync.DLL.LoTr;
        
        FPLL.HiTr = Params.P20_CohTrackSatsAndBitSync.FPLL.HiTr;
        FPLL.LoTr = Params.P20_CohTrackSatsAndBitSync.FPLL.LoTr;

    % Период, с которым производится отображение числа обработанных
    % CA-кодов
        NumCA2Disp = Params.P20_CohTrackSatsAndBitSync.NumCA2Disp;

    % Максимальное число обрабатываемых CA-кодов (inf - до конца файла!)
        MaxNumCA2Process = Params.P20_CohTrackSatsAndBitSync.MaxNumCA2Process;

    % Количество бит, используемых для битовой синхронизации
        NBits4Sync = Params.P20_CohTrackSatsAndBitSync.NBits4Sync;

    % Частота дискретизации
        Fs = Res.File.Fs;

    % Необходимость прорисовки результатов
        isDraw = Params.Main.isDraw;

%% СОХРАНЕНИЕ ПАРАМЕТРОВ
    % Track.FPLL = FPLL; % не нужно, так как всё равно будет сделано в
    % Track.DLL = DLL;   % конце
    Track.MaxNumCA2Process = MaxNumCA2Process;

    BitSync.NBits4Sync     = NBits4Sync;

%% РАСЧЁТ ПАРАМЕТРОВ
    % Длина CA-кода с учётом частоты дискретизации
        CALen = 1023 * Res.File.R;

    % Количество периодов CA-кода, приходящихся на один бит
        CAPerBit = 20;

    % Длительность CA-кода, мс
        TCA = 10^-3;

%% ОСНОВНАЯ ЧАСТЬ ФУНКЦИИ - ТРЕКИНГ

% Пояснения от 13.05.2026
%
% 1. Дискриминатор: D = 0.5 * (E - L) / (E + L), 
%   E, L - early и late корреляции;
% 2. Петлевой фильтр: фильтрует сигнал ошибки;
% 3. NCO (Numerically Controlled Oscillator): устройство, выносящее решение
%   о сдвиге на +1/-1 отсчёт;
% 4. Решение с выхода NCO применяется к сигналу, который затем снова
%   поступит на вход дискриминатора;
% 5. Процедура выполняется каждые 20 мс (или можно задать параметром);
% 6. Перед таким трекингом нужно избавиться от доплера при помощи
%   использования FPLL;
% 7. Реализация петлевого фильтра: Understanding GPS, figure. 5.20.
%   Выполняется в ClassFilter;
% 8. Пока что - одно состояние: работает только DLL;
% 9. Битовую синхронизацию пока что захардкодить и считать известной;
%
% todo: сделано фото доски (папка Other).

% Строка состояния
    fprintf('%s Трекинг спутников...\n', datetime("now") );

% Число отсчётов сигнала, которые необходимо считать из файла для
% обработки
    if MaxNumCA2Process == inf
        NumSamples2Read = Res.File.SamplesLen;
    else
        NumSamples2Read = CALen * (1 + MaxNumCA2Process);
        NumSyncs = length(1 : 4 : MaxNumCA2Process); % todo убрать магические числа
        NumSamples2Read = NumSamples2Read + NumSyncs * 1;
        clear NumSyncs;

        % Проверка, чтобы рассчитанное значение не превышало количество
        % отсчётов в записи
            if NumSamples2Read > Res.File.SamplesLen
                NumSamples2Read = Res.File.SamplesLen;
            end
    end

% Считывание сигнала из файла
    Signal = ReadSignalFromFile(Res.File, 0, NumSamples2Read);

% EPL корряляции для каждого спутника
    EPLCorsSats = cell(1, Res.Search.NumSats);
% Значения на выходе дискриминатора для каждого спутника
    DoutSats = cell(1, Res.Search.NumSats);

% Цикл по спутникам
for k = 1:Res.Search.NumSats
    % Строка состояния
        fprintf('%s     Трекинг спутника №%02d (%d из %d) ...\n', ...
            datetime("now"), Res.Search.SatNums(k), k, ...
            Res.Search.NumSats);
        
    % Компенсация частотной отстройки
        df = Res.Search.FreqShifts(k);
        Signaldf = Signal .* exp(-1j*2*pi*df * (0:length(Signal)-1) / Fs);

    % Подготовка к трекингу
        % Указатель на очередной период CA-кода
            Ptr = Res.Search.SamplesShifts(k);
        % Счётчик обработанных CA-кодов
            CACount = 0;
        % Массив значений корреляций периодов CA-кода с опорной
        % последовательностью
            CorVals = [];
        % Массив числа отсчётов, которые необходимо пропустить в записи для
        % синхронизации с соответствующим периодом СА-кода
            SamplesShifts = [];
        % Опорная последовательность
            refSeq = GenCACode(Res.Search.SatNums(k), 1);
            refSeq = 1 - 2 * repelem(refSeq, Res.File.R);

        % С каким шагом выбираем и копим СА-коды
            NumCA2Integrate = DLL.NumsIntCA(1);

        % Значения дискриминатора
            Dout = [];

        % Выход петлевого фильтра
            LFOut = [];

        % Инициализируем и настроим петлевой фильтр
            LoopFilterDLL = ClassFilter();
            LoopFilterDLL.PrepareFilter(DLL.FilterOrder, DLL.FilterBands(1), TCA*NumCA2Integrate, 0, 0);

        % Значения NCO на каждом шаге
            NCOVal = 0;

    % Трекинг
    while (CACount < MaxNumCA2Process) && ( (Ptr + NumCA2Integrate*CALen) <= NumSamples2Read)

        % Выбираем отсчёты NumCA2Integrate подряд идущих СА-кодов + 1
        % отсчёт слева и справа
            buf = Signaldf(Ptr : Ptr + CALen*NumCA2Integrate + 1);

        % EPL корреляции для каждого CA-кода
            EPLCors = zeros(NumCA2Integrate, 3);
            for caIdx = 1 : NumCA2Integrate
                EPLCors(caIdx, :) = ( ...
                    conv(buf( (1:CALen+2) + (caIdx-1)*CALen), fliplr(conj(refSeq) ), "valid") ...
                    );
            end
            EPLCors = sum(abs(EPLCors), 1);
            EPLCorsSats{k}(end+1, :) = EPLCors;

        % Значение на выходе дискриминатора
            Dout(end+1) = 1/2 * (EPLCors(1)-EPLCors(3) ) / (EPLCors(1)+EPLCors(3) );

        % Выход петлевого фильтра
            LFOut(end+1) = LoopFilterDLL.Step(Dout(end) );

        % Накопление в NCO
            NCOVal(end+1) = NCOVal(end) + LFOut(end) * TCA * NumCA2Integrate;

        % Решение о подстройке синхронизации
            if NCOVal(end) > 1 / Res.File.R
                NCOVal(end) = NCOVal(end) - 2  / Res.File.R;
                Ptr = Ptr - 1;

            elseif NCOVal(end) < - 1 / Res.File.R
                NCOVal(end) = NCOVal(end) + 2  / Res.File.R;
                Ptr = Ptr + 1;
            end

        % Вычислим корреляции СА-кодов с опорной последовательностью и
        % дробные сдвиги в отсчётах до каждого из них
            buf = Signaldf(Ptr+1 : Ptr + NumCA2Integrate*CALen);
            for ca = 1 : NumCA2Integrate
                CorVals(end+1) = buf( (1:CALen) + (ca-1)*CALen ) * refSeq';
                SamplesShifts(end+1) = Ptr + (ca-1) * CALen - NCOVal(end);
            end

        % Инкремент счётчика
            CACount = CACount + NumCA2Integrate;
        % Обновление указателя на следующую группу CA-кодов
            Ptr = Ptr + NumCA2Integrate * CALen;

        % Строка состояния
            if mod(CACount, NumCA2Disp) == 0
                fprintf('%s         Обработано %d периодов CA-кода.\n', ...
                    datetime("now"), CACount);
            end
    end

    % Сохранение результатов трекинга спутника
        Track.CorVals{k}       = CorVals;
        Track.SamplesShifts{k} = SamplesShifts;
            DoutSats{k} = Dout;

    % Прорисовка результатов и сохранение рисунков
        if isDraw > 0
            figure( ...
                Name=['P20_NonCohTrack_SatNum', num2str(Res.Search.SatNums(k) ) ], ...
                WindowStyle="docked" ...
                );

            plot(abs(CorVals) ); grid on;

            title( ['Значения модуля корреляции периодов CA-кода с опорной последовательностью для спутника № ', num2str(Res.Search.SatNums(k) ) ] );
            xlabel('Период CA-кода в записи');
            ylabel('Корреляция');

            % Установка единого масштаба для всех рисунков
                if k == 1
                    YLim = get(gca, "YLim");
                    YLim(1) = 0;
                    
                    set(gca, "YLim", YLim);
                else
                    set(gca, "YLim", YLim);
                end
        end
        if isDraw > 1
            saveas(gcf, ...
                cat(2, ...
                    Params.Main.SaveDirName, '/', ...
                    'P20_NonCohTrack_SatNum_', num2str(Res.Search.SatNums(k) ) ...
                ) ...
            );
        end
        if isDraw > 2
            close(gcf);
        end

    % Строка состояния
        fprintf('%s         Завершено.\n', datetime("now") );
end

% Добавим новое поле с результатами в Res
    Res.Track = Track;

% Очищаем рабочее пространство
    clear Signal Signaldf CorVals SamplesShifts CACount EPLCors ...
        Ptr refSeq YLim  ...
        df buf;

% Строка состояния
    fprintf('%s     Завершено.\n', datetime("now") );

%% ОСНОВНАЯ ЧАСТЬ ФУНКЦИИ - БИТОВАЯ СИНХРОНИЗАЦИЯ
% Строка состояния
    fprintf('%s Битовая синхронизация спутников ...\n', datetime("now") );
    
% Цикл по спутникам
for k = 1:Res.Search.NumSats
    % Массив корреляций, полученных при трекинге
        CorVals = Track.CorVals{k};
        CorVals = CorVals(1 : CAPerBit * NBits4Sync + 1);

    % Дифференциальное созвездие
        dCorVals = CorVals(2:end) .* conj(CorVals(1:end-1) );
        dCorVals = reshape(dCorVals, CAPerBit, []);

    % Вычисление метрик для битовой синхронизации
        Metrics = abs(sum(dCorVals, 2) );

    % Число CA-кодов, которые необходимо пропустить до начала бита
        [~, ShiftCACodes] = min(abs(Metrics) );
        if ShiftCACodes == CAPerBit, ShiftCACodes = 0; end

    % Сохранение результатов в структуру
        BitSync.CAShifts(k) = ShiftCACodes;
        BitSync.Cors(k, :) = Metrics';

    % Прорисовка результатов и сохранение рисунков
        if isDraw > 0
            figure( ...
                Name=['P20_NonCohBitSync_SatNum', num2str(Res.Search.SatNums(k) ) ], ...
                WindowStyle="docked" ...
                );

            stem(Metrics); grid on;

            title( ['Метрики, использующиеся для битовой синхронизации спутника № ', num2str(Res.Search.SatNums(k) ) ] );
            xlabel('Сколько нужно пропустить CA-кодов до начала бита');
            ylabel('Значение метрики');
        end
        if isDraw > 1
            saveas(gcf, ...
                cat(2, ...
                    Params.Main.SaveDirName, '/', ...
                    'P20_NonCohBitSync_SatNum_', num2str(Res.Search.SatNums(k) ) ...
                ) ...
            );
        end
        if isDraw > 2
            close(gcf);
        end
end

% Добавим новое поле с результатами в Res
    Res.BitSync = BitSync;

% Строка состояния
    fprintf('%s     Завершено.\n', datetime("now") );
