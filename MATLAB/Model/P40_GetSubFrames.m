function Res = P40_GetSubFrames(inRes, Params)
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
    SubFrames = struct( ...
        'isSubFrameSync', zeros(Res.Search.NumSats, 1), ... 
        'BitSeqNum',      zeros(Res.Search.NumSats, 1), ...
        'BitShift',       zeros(Res.Search.NumSats, 1), ...
        'Words',          {cell(Res.Search.NumSats, 1)} ...
    );
    % Каждый элемент массива isSubFrameSync - флаг успешности подкадровой
    %   синхронизации.
    % Каждый элемент массива BitSeqNum - номер битового потока, в котором
    %   удалось выполнить синхронизацию с началом подкадра.
    % Каждый элемент массива BitShift - количество бит, которые надо
    %   пропустить от начала битового потока до начала первого подкадра.
    % Каждая ячейка cell-массива Words - cell-массив (Nx10), где N -
    %   количество обработанных подкадров, каждая ячейка - массив 1х24 бит
    %   декодированного слова, если CRC сошлось, и пустой массив, если CRC
    %   не сошлось.

%% УСТАНОВКА ПАРАМЕТРОВ
% Необходимость прорисовки результатов
    isDraw = Params.Main.isDraw;

% Папка для сохранения результатов
    SaveDirName = Params.Main.SaveDirName;

%% РАСЧЁТ ПАРАМЕТРОВ

%% ОСНОВНАЯ ЧАСТЬ ФУНКЦИИ - ЦИКЛ ПО НАЙДЕННЫМ СПУТНИКАМ
for k = 1 : Res.Search.NumSats
    % Выделение демодулированных битовых последовательностей спутника
        Bits = Res.Demod.Bits{k};

    % Определение номера валидной битовой последовательности и подкадровая
    % синхронизация
        [isSFSyncOk, BitSeqNum, BitShift] = ...
            SubFrameSync(Bits, isDraw, SaveDirName, Res.Search.SatNums(k) );

    % Присвоение результата
        SubFrames.isSubFrameSync(k) = isSFSyncOk;
        SubFrames.BitSeqNum(k)      = BitSeqNum;
        SubFrames.BitShift(k)       = BitShift;

    % Если подкадровая синхронизация не была выполнена, работать с данным
    % спутником более нет смысла
        if ~isSFSyncOk, continue; end

    % Выделение валидной последовательности
        Bits = Bits(BitSeqNum, :);

    % Синхронизация с началом подкадра с оставлением двух бит предыдщего
    % слова
        P1 = BitShift - 2 +1;
        if P1 < 1
            Bits = [zeros(1 - P1, 1), Bits]; %#ok<AGROW>
            P1   = 1;
        end
        BitsSync = Bits(P1 : end);

    % Выделение и декодирование слов
        Words = CheckFrames(BitsSync);

    % Присвоение результата
        SubFrames.Words{k} = Words;
end

% Присвоение нового поля структуре с результатом
    Res.SubFrames = SubFrames;
end

function Words = CheckFrames(Bits)
%
% Из битового потока выделяются все возможные кадры, в каждом кадре
% проверяется CRC каждого слова, если CRC сошлось, то сохраняется
% декодированное слово, в противном случае сохраняется пустой массив
%
% Bits - битовый поток. Третий бит потока совпадает с началом подкадра,
%        первые два бита принадлежат предыдущему слову.

    % Число бит в одном подкадре, слове
        SFLen   = 300;
        WordLen = 30;

    % Определение числа подкадров
        NumSFs = floor(length(Bits(3:end) ) / SFLen );

    % Инициализация результата
        Words = cell(NumSFs, SFLen / WordLen);

    % Выделение целого числа подкадров + 2 бита предыдущего слова
        Bits = Bits(1 : NumSFs * SFLen + 2);

    % Декодирование всех слов в битовом потоке
    for sfIdx = 1 : NumSFs % Цикл по подкадрам

        % Выделение бит подкадра
            SF = Bits( (1 : SFLen+2) + (sfIdx-1) * SFLen);

        for wIdx = 1 : 10 % Цикл по словам
            % Выделение слова
                EWord = SF( (1:WordLen+2) + (wIdx-1) * WordLen);

            % Декодирование
                [isOk, DWord] = CheckCRC(EWord);

            % Если CRC сошлось, присваиваем результат
                if isOk
                    Words{sfIdx, wIdx} = DWord;
                end
        end
    end
end

function [isOk, BitSeqNum, BitShift] = SubFrameSync(Bits, isDraw, ...
    SaveDirName, SatNum)
%
% Функция подкадровой синхронизации
%
% isOk      - флаг, указывающий, найдена синхронизация или нет, причём она
%   должна быть найдена только один раз!
% BitSeqNum - номер битовой последовательности, для которой найдена
%   синхронизация. т.е. последовательности, с которой надо дальше работать.
% BitShift  - количество бит, которые нужно пропустить в битовой
%   последовательности до начала подкадра.

    % Инициализация результата
        BitSeqNum = 1;
        BitShift = 0;

    % Количество и длина битовых последовательностей
        NumBitsSeqs = size(Bits, 1);
        BitSeqLen   = size(Bits, 2);

    % Длина подкадра, бит
        SFLen = 300;

    % Синхрослово начала подкадра
        refSeq    = [1 0 0 0 1 0 1 1];
        refSeq    = 1 - 2*refSeq;
        refSeqLen = length(refSeq);

    % Число периодов при накоплении
        NumAccPers = floor(BitSeqLen / SFLen);

    % Определим порог
        Threshold = refSeqLen * NumAccPers * 0.96;

    % Память под результат корреляции битовой последовательности и
    % синхрослова
        CorRes = zeros(NumBitsSeqs, BitSeqLen - refSeqLen + 1);
    % Массив метрик, сравниваемых с порогом, и их позиций в массиве
        Metrics     = zeros(NumBitsSeqs, 1);
        MetricsInds = zeros(NumBitsSeqs, 1);

    % Цикл по битовым последовательностям
    for k = 1 : NumBitsSeqs
        % Корреляция с синхрословом
            CorRes(k, :) = conv(1 - 2*Bits(k, :), fliplr(refSeq), "valid");

        % Накопление результата с периодом, равным длине подкадра
            CorResAcc = CorRes(k, 1 : end-mod(length(CorRes), SFLen) );
            CorResAcc = reshape(CorResAcc, SFLen, []);
            CorResAcc = sum(abs(CorResAcc), 2);

        % Максимальное значение результата накопления
            [Metrics(k), MetricsInds(k) ] = max(CorResAcc);
    end

    % Сравним рассчитанные метрики с порогом
        isThresholdExceeded = Metrics >= Threshold;

    % Определим, был ли превышен порог хотя бы для одной последовательности
        isOk = any(isThresholdExceeded);
    % Если порог не был превышен, завершаем выполнение функции
        if ~isOk, return; end

    % Присвоение результата
        if sum(isThresholdExceeded) > 1 % Если превышений порога больше 
                                        % одного, выберем результат с 
                                        % наибольшей метрикой
            [~, BitSeqNum] = max(Metrics);
        else
            BitSeqNum = find(isThresholdExceeded == 1);
        end
        BitShift  = MetricsInds(BitSeqNum) - 1;
end

function [isOk, DWord] = CheckCRC(EWord)
% Функиця осуществляет проверку CRC для одного слова навигационного
% сообщения
%
% See IS-GPS-200H, p. 136-138

% На входе:
%   EWord - слово (строка) с двумя битами предыдущего слова в начале, т.е.
%     всего 32 бита.

% На выходе: 
%   isOk - 1, если CRC сходится, 0 в противном случае.
%   DWord - декодированное слово (строка), т.е. всего 24 бита.

    % Выделение систематических, проыерочных бит, а также бит предедущего
    % слова
        D29Prev = EWord(1);
        D30Prev = EWord(2);
        SysBitsD     = EWord(2 + (1:24) );
        rxParityBits = EWord(end - 6 +1 : end);

    % Получение d1, d2, ..., d30, т.е. устранение возможной инверсии
    % систематических бит
        SysBitsd = mod(SysBitsD + D30Prev, 2);
        DWord = SysBitsd;

    % Вычисление проверочных бит
        refParityBits = zeros(1, 6);

        refParityBits(1) = D29Prev + sum(SysBitsd( [1 2 3 5 6 10 11 12 13 14 17 18 20 23]    ) );
        refParityBits(2) = D30Prev + sum(SysBitsd( [2 3 4 6 7 11 12 13 14 15 18 19 21 24]    ) );
        refParityBits(3) = D29Prev + sum(SysBitsd( [1 3 4 5 7 8  12 13 14 15 16 19 20 22]    ) );
        refParityBits(4) = D30Prev + sum(SysBitsd( [2 4 5 6 8 9  13 14 15 16 17 20 21 23]    ) );
        refParityBits(5) = D30Prev + sum(SysBitsd( [1 3 5 6 7 9  10 14 15 16 17 18 21 22 24] ) );
        refParityBits(6) = D29Prev + sum(SysBitsd( [3 5 6 8 9 10 11 13 15 19 22 23 24]       ) );

        refParityBits = mod(refParityBits, 2);

    % Сравнение вычисленных и принятых проверочных бит
        isOk = isequal(rxParityBits, refParityBits);
end
