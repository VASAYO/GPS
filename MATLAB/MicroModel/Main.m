% Микромодель обработки сигналов GPS
    clc; clear;
    close all;

% Параметры
    % Опорная частота дискретизации
        Fbase = 1.023e6;
    % Коэффициент передискретизации
        sps = 2;
    % Длина C/A кода в чипах
        CACodeLen = 1023;
    % Число периодов C/A кода, исп-мых при обнаружении
        NumCACodePers = 20;

    % Число периодов C/A кодов, после которых необходимо подстраивать
    % синхронизацию по времени
        NumCA2Sync = 100;

    % Прорисовка результатов
        isDrawRes = false;

% Вычисляемые параметры

% Считывание сигнала из файла
    File.Name = '../../Signals/28_01_2019__17_02_51_x02_1ch_16b_15pos_200000ms.dat';
    File.HeadLenInBytes = 0;
    File.NumOfChannels  = 1;
    File.ChanNum  = 0;
    File.DataType = 'int16';
    File.Fs0      = Fbase * 2;
    File.dF       = 0;
    File.FsDown   = 1;
    File.FsUp     = 1;
    NumOfShiftedSamples = 0;
    NumOfNeededSamples  = 4 * File.Fs0 * File.FsUp;

    [Signal, File] = ReadSignalFromFile(File, NumOfShiftedSamples, ...
        NumOfNeededSamples);

% Обнаружение спутника
    CACodeNum = 1;

    % Отрезок сигнала в начале записи
        SignalShort = Signal(1:(NumCACodePers+1) * CACodeLen * sps - 1);

    % Массив значений сдвигов частоты
        FVals = 0 : 50 : 5200;
        FVals = [-fliplr(FVals(2:end) ), FVals];

    % Эталонный C/A код
        ethCACode = 1 - 2 * GenCACode(CACodeNum, 1);

    % Построение тела неопределённости
        CorrVals = zeros(length(FVals), CACodeLen * sps);

        for k = 1 : length(FVals)
            % Опорная последовательность
                refSeq = repelem(ethCACode, sps);
                refSeq = refSeq .* ...
                    exp(1j*2*pi*FVals(k) * (0:length(refSeq)-1) / File.Fs);

            % Корреляция
                buf  = conv(SignalShort, fliplr(conj(refSeq) ), "valid");

            % Некогерентное накопление результата
                buf1 = reshape(buf, CACodeLen * sps, [] ).';
                CorrVals(k, :) = sum(abs(buf1) );
        end

        if isDrawRes
            figure(CACodeNum)
            surf(CorrVals)
        end

    % Определение грубого сдвига частоты и сдвига до начала первого C/A
    % кода
        [buf, IndsMaxY] = max(CorrVals);
        [~, IndMaxX]    = max(buf);
        IndMaxY = IndsMaxY(IndMaxX);

        df = FVals(IndMaxY);
        OffsetSamples = IndMaxX - 1;

% Грубая подстройка частоты
    Signal = Signal .* exp(-1j*2*pi * df * (0:length(Signal)-1)/File.Fs);

% Вычисление корреляций C/A кодов с опорной последовательностью
    % Число полных C/A кодов в записи
        NumFullCACodes = floor( (length(Signal) - OffsetSamples) / length(refSeq) );

    % Указатель, на начало очередного периода C/A кода
        Ptr = OffsetSamples +1;

    % Значения корреляций
        PCorrs = zeros(NumFullCACodes, 1);

    % Опорная последовательность
        refSeq = repelem(ethCACode, sps).';

    % Счётчик периодов C/A кодов, обработанных после последней подстройки
    % синхронизации по времени
        CntrCA = 0;

    for k = 1 : NumFullCACodes
        % Обнуление переменной, отвечающей за подстройку символьной
        % синхронизации
            Ofst = 0;

        % Вычисление корреляции. При необходимости выполнение
        % Early-Late-Prompt корреляций и подстройка синхронизации по
        % времени
            if CntrCA == NumCA2Sync
                % Обнуление счётчика
                    CntrCA = 0;

                % Early-Prompt-Late корреляции
                    CorrValsEPL = zeros(1, 3);

                    CorrValsEPL(1) = Signal( (0:length(refSeq)-1) + Ptr - 1) * refSeq;
                    CorrValsEPL(2) = Signal( (0:length(refSeq)-1) + Ptr)     * refSeq;
                    CorrValsEPL(3) = Signal( (0:length(refSeq)-1) + Ptr + 1) * refSeq;

                % Выбор наибольшего значения
                    [PCorrs(k), Ofst] = max(abs(CorrValsEPL) );

                Ofst = Ofst - 2;

            else
                PCorrs(k) = Signal( (0:length(refSeq)-1) + Ptr) * refSeq;
            end

        % Инкремент счётчика
            CntrCA = CntrCA + 1;

        % Обновление значение указателя
            Ptr = Ptr + length(refSeq) + Ofst;
    end

% Синхронизация с началом бита
    % Дифференциальное созвездие P-корреляций C/A кодов
        PCorrsDiff = PCorrs(2 : end) .* conj(PCorrs(1 : end-1) );
    
    % Накопление
        PCorrsDiff = reshape(PCorrsDiff(1 : end - mod(length(PCorrsDiff), 20) ), ...
            20, [] ...
        );
        Metrics = abs(sum(PCorrsDiff, 2) );

    % Позиция корреляции, соответствующая началу первого бита в записи
        [~, OffsetCACodes] = min(Metrics);
        OffsetCACodes = OffsetCACodes - 1;

% Синхронизация с началом бита
    % Число бит в записи
        NumBits = floor(length(PCorrs(OffsetCACodes+1 : end) ) / 20);

    PCorrsSync = PCorrs(OffsetCACodes + (1:20*NumBits) );
    PCorrsSync = reshape(PCorrsSync, 20, []);

    % Когерентное накопление
        BitCorrVals = sum(PCorrsSync, 1);