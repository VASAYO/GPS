% Исследование корреляционных свойств C/A кода

clc; clear;
close all;

addpath('..\Interfaces\');

% Отсчёты C/A кода
    CACode = GenCACode(1, 1);
    CACode = 1 - 2 * CACode;

% Корреляционная функция
    arg1 = repmat(CACode, 1, 30);
    arg1 = arg1(1:end-1) .* exp(1j * 2 * pi * 0 * (0:length(arg1)-1 -1) / 1.023e6 );
    CorrFun1Dim = conv(arg1, fliplr(conj(CACode) ), "valid");

% Накопление результата
    buf = reshape(CorrFun1Dim, 1023, []);
    AccRes = sum(abs(buf), 2);

% Прорисовка
    figure(1);
    stem(abs(CorrFun1Dim) ); grid on;

    figure(2);
    stem(AccRes)

%% Двумерная корреляционная функция
clc; clear;
close all;

addpath('..\Interfaces\');

% Чиповая скорость
    Rc = 1.023e6;

% C/A код
    CACode = GenCACode(1, 1);
    CACode = CACode * 2 - 1;

% Массив значений частот
    FVals = -5e3 : 1 : 5e3;

% Значения двумерного тела неопределенности
    CorrFun2Dim = zeros(length(FVals), length(CACode) );

% Моделирование тела неопределенности
    arg1 = [CACode, CACode];
    arg1 = arg1(1 : end - 1);

    for k = 1 : length(FVals)
        buf = arg1 .* exp(1j * 2 * pi * FVals(k) * (0:length(arg1)-1)/Rc);
        CorrFun2Dim(k, :) = conv(buf, fliplr(conj(CACode)), "valid");
    end

% Прорисовка
    figure(1);
    surf(abs(CorrFun2Dim) );

    figure(2)
    plot(FVals, abs(CorrFun2Dim(:, 1) ) );
    yline( max(abs(CorrFun2Dim(:, 1) ) ) / sqrt(2) );