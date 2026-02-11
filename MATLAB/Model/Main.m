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
    inFile.Name = '../Signals/28_01_2019__17_02_51_x02_1ch_16b_15pos_200000ms.dat';
    inFile.HeadLenInBytes = 0;
    inFile.NumOfChannels  = 1;
    inFile.ChanNum  = 0;
    inFile.DataType = 'int16';
    inFile.Fs0      = Fbase * 2;
    inFile.dF       = 0;
    inFile.FsDown   = 1;
    inFile.FsUp     = 3;
    NumOfShiftedSamples = 0;
    NumOfNeededSamples  = (NumCACodePers+1) * CACodeLen * sps - 1;

    [Signal, File] = ReadSignalFromFile(inFile, NumOfShiftedSamples, ...
        NumOfNeededSamples);

% Обнаружение спутников
    CACodeNum = 31;

    % Массив значений сдвигов частоты
        FVals = 0 : 720 : 720 * 7;
        FVals = [-fliplr(FVals(2:end) ), FVals];

    for i = 1:63
        % Эталонный C/A код
            ethCACode = 1 - 2 * GenCACode(i, 1);
    
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
    
            figure(i)
        surf(CorrVals)
    end