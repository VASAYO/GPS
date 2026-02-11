function CACode = GenCACode(SigNum, NumCycles)
%
% Функция генерирует заданное число периодов (циклов) C/A кода системы GPS,
% т.е. массив-строку длиной кратной 1023.
%  
%   SigNum    - номер кода, 1...63;
%   NumCycles - количество периодов C/A кода по 1023 отсчёта (1...1000), по
%               умолчанию NumCycles = 1.
%
%   CACode    - массив строка 1хNumCycles*1023, содержащая NumCycles 
%               периодов C/A кода.

arguments
    SigNum    (1, 1) double {mustBeInteger, ...
                             mustBeInRange(SigNum, 1, 63)} = 1
    NumCycles (1, 1) double {mustBeInteger, ...
                             mustBeInRange(NumCycles, 1, 1000)} = 1
end

% Номера отводов регистра G2 для генерации G2_i при SigNum <= 37
    Taps_Table = [ ...
        2, 6;
        3, 7;
        4, 8;
        5, 9;
        1, 9;
        2, 10;
        1, 8;
        2, 9;
        3, 10;
        2, 3;
        3, 4;
        5, 6;
        6, 7;
        7, 8;
        8, 9;
        9, 10;
        1, 4;
        2, 5;
        3, 6;
        4, 7;
        5, 8;
        6, 9;
        1, 3;
        4, 6;
        5, 7;
        6, 8;
        7, 9;
        8, 10;
        1, 6;
        2, 7;
        3, 8;
        4, 9;
        5, 10;
        4, 10;
        1, 7;
        2, 8;
        4, 10 ...
    ];

% Задержка последовательности G2 относительно G1 для генерации G2_i при
% SigNum > 37
    G2_Delays = [ ...
        67
        103
        91
        19
        679
        225
        625
        946
        638
        161
        1001
        554
        280
        710
        709
        775
        864
        558
        220
        397
        55
        898
        759
        367
        299
        1018 ...
    ];

% Генерация одного периода C/A кода
    RegG1 = ones(1, 10);
    RegG2 = ones(1, 10);
    G1   = zeros(1, 1023);
    G2_i = zeros(1, 1023);

    % Генерация G1
        for k = 1 : 1023
            G1(k) = RegG1(end);
            fb = mod(sum(RegG1( [3 10] ) ), 2);
            RegG1(2:end) = RegG1(1:end-1);
            RegG1(1) = fb;
        end

    % Генерация G2_i
        if SigNum <= 37 % Реализация при помощи двухотводного кодера
            for k = 1 : 1023
                fb = mod( sum(RegG2( [2 3 6 8 9 10] ) ), 2);
                G2_i(k) = mod(sum(RegG2(Taps_Table(SigNum, :) ) ), 2);
                RegG2(2:end) = RegG2(1:end-1);
                RegG2(1) = fb;
            end

        else % Реализация задержки напрямую
            G2Delay = G2_Delays(SigNum - 37);
            G2 = zeros(1, 1023);

            for k = 1 : 1023
                G2(k) = RegG2(end);
                fb = mod( sum(RegG2( [2 3 6 8 9 10] ) ), 2);
                RegG2(2:end) = RegG2(1:end-1);
                RegG2(1) = fb;
            end

            G2_i = circshift(G2, G2Delay);
        end

    % Поэлементное сложение последовательностей
        CACode = mod(G1 + G2_i , 2);

% Генерация нескольких периодов кода
    CACode = repmat(CACode, 1, NumCycles);
