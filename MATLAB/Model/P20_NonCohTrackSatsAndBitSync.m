function Res = P20_NonCohTrackSatsAndBitSync(inRes, Params)
%
% Функция некогерентного трекинга спутников и битовой синхронизации
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
        'SamplesShifts', {cell(Res.Search.NumSats, 1)}, ... 
        'CorVals',       {cell(Res.Search.NumSats, 1)} ...
    );
    % Каждая ячейка cell-массивов SamplesShifts и CorVals является массивом
    %   1xN, где N - количество периодов CA-кода соответствующего спутника,
    %   найденных в файле-записи (N может быть разным для разных
    %   спутников).
    % Каждый элемент массива SamplesShifts{k} - количество отсчётов,
    %   которые надо пропустить в файле-записи до начала соответствующего
    %   периода CA-кода.
    % Каждый элемент массива CorVals{k} - комплексное значение корреляции
    %   части сигнала, содержащей соответствующий период CA-кода, с опорным
    %   сигналом.

    BitSync = struct( ...
        'CAShifts', zeros(Res.Search.NumSats, 1), ... 
        'Cors', zeros(Res.Search.NumSats, 20) ...
    );
    % Каждый элемент массива CAShifts - количество периодов CA-кода,
    %   которые надо пропустить до начала бита.
    % Каждая строка массива Cors - корреляции, по позиции минимума которых
    %   определяется битовая синхронизация.

%% УСТАНОВКА ПАРАМЕТРОВ
    % Количество периодов CA-кода между соседними синхронизациями по
    % времени (NumCA2NextSync >= 1, NumCA2NextSync = 1 - синхронизация для
    % каждого CA-кода)
        NumCA2NextSync = Params.P20_NonCohTrackSatsAndBitSync.NumCA2NextSync;

    % Половина количества дополнительных периодов CA-кода, используемых для
    % синхронизации по времени
        HalfNumCA4Sync = Params.P20_NonCohTrackSatsAndBitSync.HalfNumCA4Sync;

    % Количество учитываемых значений задержки/набега синхронизации по
    % времени
        HalfCorLen = Params.P20_NonCohTrackSatsAndBitSync.HalfCorLen;

    % Период, с которым производится отображение числа обработанных
    % CA-кодов
        NumCA2Disp = Params.P20_NonCohTrackSatsAndBitSync.NumCA2Disp;

    % Максимальное число обрабатываемых CA-кодов (inf - до конца файла!)
        MaxNumCA2Process = Params.P20_NonCohTrackSatsAndBitSync.MaxNumCA2Process;

    % Количество бит, используемых для битовой синхронизации
        NBits4Sync = Params.P20_NonCohTrackSatsAndBitSync.NBits4Sync;

    % Частота дискретизации
        Fs = Res.File.Fs;

    % Необходимость прорисовки результатов
        isDraw = Params.Main.isDraw;

%% СОХРАНЕНИЕ ПАРАМЕТРОВ
    Track.NumCA2NextSync   = NumCA2NextSync;
    Track.HalfNumCA4Sync   = HalfNumCA4Sync;
    Track.HalfCorLen       = HalfCorLen;
    Track.MaxNumCA2Process = MaxNumCA2Process;

    BitSync.NBits4Sync     = NBits4Sync;

%% РАСЧЁТ ПАРАМЕТРОВ
    % Длина CA-кода с учётом частоты дискретизации
        CALen = 1023 * Res.File.R;

    % Количество периодов CA-кода, приходящихся на один бит
        CAPerBit = 20;

%% ОСНОВНАЯ ЧАСТЬ ФУНКЦИИ - ТРЕКИНГ
% Строка состояния
    fprintf('%s Трекинг спутников\n', datetime("now") );

% Число отсчётов сигнала, которые необходимо считать из файла для
% обработки
    if MaxNumCA2Process == inf
        NumSamples2Read = Res.File.SamplesLen;
    else
        NumSamples2Read = CALen * (1 + MaxNumCA2Process + HalfNumCA4Sync);
        NumSyncs = length(1 : NumCA2NextSync : MaxNumCA2Process);
        NumSamples2Read = NumSamples2Read + NumSyncs * HalfCorLen;
        clear NumSyncs;

        % Проверка, чтобы рассчитанное значение не превышало количество
        % отсчётов в записи
            if NumSamples2Read > Res.File.SamplesLen
                NumSamples2Read = Res.File.SamplesLen;
            end
    end

% Считывание сигнала из файла
    Signal = ReadSignalFromFile(Res.File, 0, NumSamples2Read);

% Цикл по спутникам
for k = 1:Res.Search.NumSats
    % Строка состояния
        fprintf('%s     Трекинг спутника №%02d (%d из %d) ...\n', ...
            datetime("now"), Res.Search.SatNums(k), k, ...
            Res.Search.NumSats);
        
    % Компенсация частотной отстройки
        df = Res.Search.FreqShifts(k);
        Signaldf = Signal .* exp(-1j*2*pi*df * (0:length(Signal)-1) / Fs );

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
        % Опорная последовательность для корреляции
            refSeqCor = GenCACode(Res.Search.SatNums(k), 1).';
            refSeqCor = 1 - 2 * repelem(refSeqCor, Res.File.R);
        % Опорная последовательность для синхронизации
            refSeqSync = GenCACode( ...
                Res.Search.SatNums(k), 1 + 2 * HalfNumCA4Sync);
            refSeqSync = 1 - 2 * repelem(refSeqSync, Res.File.R);

    % Трекинг
    while (CACount < MaxNumCA2Process) && ( (Ptr + CALen) <= NumSamples2Read)

        % Синхронизация при выполнении условия
        if mod(CACount, NumCA2NextSync) == 0 && CACount ~= 0

            % Позиции начала и конца отрезка сигнала, который необходимо
            % использовать для синхронизации
                P1 = Ptr + 1 - CALen * HalfNumCA4Sync - HalfCorLen;
                P2 = Ptr + CALen * (1+HalfNumCA4Sync) + HalfCorLen;

            % Если P1 или P2 выходят за диапазон [1; NumSamples2Read],
            % нужно предусмотреть добавление нулей для соответствия
            % размерностей
                ZerosBefore = 0;
                ZerosAfter  = 0;

                if P1 < 1
                    ZerosBefore = 1 - P1;
                    P1 = 1;
                end
                if P2 > NumSamples2Read
                    ZerosAfter = P2 - NumSamples2Read;
                    P2 = NumSamples2Read;
                end

            % Выбор отрезка сигнала для временной синхронизации
                buf = Signaldf(P1 : P2);
                buf = [zeros(1, ZerosBefore), buf, zeros(1, ZerosAfter) ]; %#ok<AGROW>

            % Корреляция с опорной последовательностью
                EPLCors = abs( ...
                    conv(buf, fliplr(conj(refSeqSync) ), "valid") ...
                    );

            % Определение ухода синхронизации по времени и её подстройка
                [~, PosMax] = max(EPLCors);
                drift = PosMax - median(1:length(EPLCors) );

                Ptr = Ptr + drift;
        end

        % Сохранение сдвига до CA-кода в массив
            SamplesShifts(end+1) = Ptr; %#ok<AGROW>

        % Корреляция периода СА-кода с опорной последовательностью
            CorVals(end+1) = Signaldf(Ptr + (1:CALen) ) * conj(refSeqCor); %#ok<AGROW>

        % Инкремент счётчика
            CACount = CACount + 1;
        % Обновление указателя на следующий период CA-кода
            Ptr = Ptr + CALen;

        % Строка состояния
            if mod(CACount, NumCA2Disp) == 0
                fprintf('%s       Обработано %d периодов CA-кода.\n', ...
                    datetime("now"), CACount);
            end
    end

    % Сохранение результатов трекинга спутника
        Track.CorVals{k}       = CorVals;
        Track.SamplesShifts{k} = SamplesShifts;

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
        fprintf('%s     Завершено.\n', datetime("now") );
end

% Добавим новое поле с результатами в Res
    Res.Track = Track;

% Очищаем рабочее пространство
    clear Signal Signaldf CorVals SamplesShifts CACount EPLCors ...
        P1 P2 Ptr refSeqCor refSeqSync PosMax YLim ZerosAfter  ...
        ZerosBefore drift df buf;

% Строка состояния
    fprintf('%s Завершено.\n', datetime("now") );

%% ОСНОВНАЯ ЧАСТЬ ФУНКЦИИ - БИТОВАЯ СИНХРОНИЗАЦИЯ
% Цикл по спутникам
for k = 1:Res.Search.NumSats

    % Строка состояния
        fprintf('%s Битовая синхронизация спутника №%02d (%d из %d) ...\n', ...
            datetime("now"), Res.Search.SatNums(k), k, ...
            Res.Search.NumSats);

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

    % Строка состояния
        fprintf('%s   Завершено.\n', datetime("now") );
end

% Добавим новое поле с результатами в Res
    Res.BitSync = BitSync;
