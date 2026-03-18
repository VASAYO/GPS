function [Word, Success] = DecodeWord(CodeWord, D29Prev, D30Prev)
% Функция выполняет декодирование поданного на вход кодового слова
%
% CodeWord         - 30х1 массив бит кодового слова D1, D2, ..., D30;
% D29Prev, D30Prev - биты D29, D30 предыдущего кодового слова;
%
% Word    - 24х1 массив систематических бит после декодирования;
% Success - результат проверки на четность. Success = true означает 
%           отсутствие обнаруженных ошибок при декодировании.
%
% See IS-GPS-200H, p. 136-138

% Систематические и проверочные биты кодового слова
    SysBitsD = CodeWord(1:24);
    rxParityBits = CodeWord(25:30);

% Получение d1, d2, ..., d30, т.е. устранение возможной инверсии
% систематических бит
    SysBitsd = mod(SysBitsD + D30Prev, 2);
    Word = SysBitsd;

% Вычисление проверочных бит
    refParityBits = zeros(6, 1);

    refParityBits(1) = D29Prev + sum(SysBitsd( [1 2 3 5 6 10 11 12 13 14 17 18 20 23] ) );
    refParityBits(2) = D30Prev + sum(SysBitsd( [2 3 4 6 7 11 12 13 14 15 18 19 21 24] ) );
    refParityBits(3) = D29Prev + sum(SysBitsd( [1 3 4 5 7 8  12 13 14 15 16 19 20 22] ) );
    refParityBits(4) = D30Prev + sum(SysBitsd( [2 4 5 6 8 9  13 14 15 16 17 20 21 23] ) );
    refParityBits(5) = D30Prev + sum(SysBitsd( [1 3 5 6 7 9  10 14 15 16 17 18 21 22 24] ) );
    refParityBits(6) = D29Prev + sum(SysBitsd( [3 5 6 8 9 10 11 13 15 19 22 23 24] ) );

    refParityBits = mod(refParityBits, 2);

% Сравнение вычисленных и принятых проверочных бит
    Success = isequal(rxParityBits, refParityBits);
