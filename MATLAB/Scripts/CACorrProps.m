% Исследование корреляционных свойств CA-кодов
    clc; clear; close all;

% Параметры и исходные данные
    % Частота дискретизации, Гц
        Fs = 1.023e6;
    % Тип корреляционной функции:
    %   'aperiodic' - апериодическая;
    %   'periodic'  - периодическая;
        Type = 'aperiodic';
    % Номера СА-кодов, для которых исследуются взаимные корреляционные
    % свойства.
        CANum1 = 2;
        CANum2 = 2;
    % Массив частот для построения тела неопределённости
        dfVals = -5e3:10:5e3;
    % Число отсчётов одного периода CA-кода
        CALen = 1023;

% Генерация одного периода СА-кодов в двухполярной форме
    CA1 = 1 - 2*GenCACode(CANum1, 1);
    CA2 = 1 - 2*GenCACode(CANum2, 1);

% В зависимости от того, исследуем мы апериодическую или периодическую КФ,
% добавим к одному из СА-кодов по бокам нули или доп. периоды СА-кода
    if strcmp(Type, 'aperiodic')
        S1 = [zeros(1, 1022), CA1, zeros(1, 1022)]; 
    else
        S1 = [CA1(2:end), CA1, CA1(1:end-1)];
    end

%% Построение тела неопределённости
% Память под результат
    CorRes = zeros(length(dfVals), 1023*2-1);

% Цикл по частотным сдвигам
t = (0 : length(CA2)-1) / Fs;
for idf = 1 : length(dfVals)
    S2 = CA2 .* exp(1j*2*pi*dfVals(idf) * t);
    buf = conv(S1, fliplr(conj(S2) ), "valid");
    CorRes(idf, :) = buf;
end

%% Прорисовка результатов
dtVals = (-CALen+1 : CALen-1) / Fs;

figure(WindowStyle="docked");
surf(dtVals, dfVals, abs(CorRes) );
grid on;
title('Тело неопределённости')
xlabel("Задержка, с");
ylabel("Частота, Гц");
zlabel("Корреляция");

figure(WindowStyle="docked");
plot(dtVals, abs(CorRes(dfVals == 0, :) ) );
grid on;
title('Сечение тела неопределённости при нулевом сдвиге по частоте')
xlabel("Задержка, с");
ylabel("Корреляция");

figure(WindowStyle="docked");
plot(dfVals, abs(CorRes(:, median(1:size(CorRes, 2) ) ).' ) );
grid on;
title('Сечение тела неопределённости при нулевом сдвиге по времени')
xlabel("Частота, Гц");
ylabel("Корреляция");
