clc; clear;
close all;

% Параметры
    % Опорная частота дискретизации
        Fbase = 1.023e6;
    % Коэффициент передискретизации
        sps = 6;
    % Длина C/A кода в чипах
        CACodeLen = 1023;
    % Число периодов C/A кода, исп-мых при обнаружении
        NumCACodePers = 20;

% Вычисляемые параметры

% Считывание сигнала из файла
    File.Name = '../Signals/28_01_2019__17_02_51_x02_1ch_16b_15pos_200000ms.dat';
    File.HeadLenInBytes = 0;
    File.NumOfChannels  = 1;
    File.ChanNum  = 0;
    File.DataType = 'int16';
    File.Fs0      = Fbase * 2;
    File.dF       = 0;
    File.FsDown   = 1;
    File.FsUp     = 3;
    NumOfShiftedSamples = 0;
    NumOfNeededSamples  = (NumCACodePers+1) * CACodeLen * sps - 1;

    [Signal, File] = ReadSignalFromFile(File, NumOfShiftedSamples, ...
        NumOfNeededSamples);

% Обнаружение спутника
    CACodeNum = 1;

    % Массив значений сдвигов частоты
        FVals = 0 : 400 : 5200;
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
                buf  = conv(Signal, fliplr(conj(refSeq) ), "valid");

            % Некогерентное накопление результата
                buf1 = reshape(buf, CACodeLen * sps, [] ).';
                CorrVals(k, :) = sum(abs(buf1) );
        end

        figure(CACodeNum)
        surf(CorrVals)
