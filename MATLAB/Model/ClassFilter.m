classdef ClassFilter < handle
% Описание класса для loop filter типа DLL или FLL-assisted PLL filter    
    properties
        % Вариант фильтра (DLL) или (FLL, PLL)
            isDLL;
        % Порядок фильтра [1], [2], [3] (DLL) или [1, 2], [2, 3] (FLL, PLL)
            Order;
        % Полоса фильтров DLL, FLL, PLL
            Bnd, Bnf, Bnp;
        % Шаг поступления отсчётов, сек
            T;
        % Аккумуляторы хранения скорости и ускорения
            VelocAcc, AccelAcc;
        % Коэффициенты фильтров DLL, FLL, PLL
            CoefsDLL, CoefsFLL, CoefsPLL;
    end
    methods
        function PrepareFilter(obj, Order, Bn, T, VelocAcc, AccelAcc)
        % Функция инициализации фильтра
        %
        % Входные параметры
        %   Order - 1 | 2 | 3 | [1, 2] | [2, 3];
        %   Bn    - 1 или 2 элемента в массиве, значения >0 в Гц.
        %   T     - значение >0 в секундах.
        %   VelocAcc, AccelAcc - необязательные переменные для
        %       инициализации аккумуляторов

            % Заглушка
                if ~isequal(Order, 1)
                    error("Реализован только фильтр первого порядка для DLL");
                end
        
            % Определим тип фильтра
                obj.isDLL = isscalar(Order);
            % Сохраним параметры фильтра внутри объекта
                obj.Order = Order;
                if obj.isDLL
                    obj.Bnd = Bn;
                else
                    obj.Bnf = Bn(1);
                    obj.Bnp = Bn(2);
                end
                obj.T = T;

            % Расчитаем коэффициенты фильтров DLL, FLL и PLL
                if obj.isDLL
                    obj.CoefsDLL = CalcCoefs(Order, Bn);
                else

                end
            % Инициализируем значения аккумуляторов скорости и ускорения
                obj.VelocAcc = VelocAcc;
                obj.AccelAcc = AccelAcc;
        end
        function ChangeParams(obj, Bn, T, VelocAcc, AccelAcc)
            % Сохраним параметры фильтра внутри объекта
                if obj.isDLL
                    obj.Bnd = Bn;
                else
                    obj.Bnf = Bn(1);
                    obj.Bnp = Bn(2);
                end
                obj.T = T;

            % Расчитаем коэффициенты фильтров DLL, FLL и PLL
                if obj.isDLL
                    obj.CoefsDLL = CalcCoefs(obj.Order, Bn);
                else

                end
            
            % Инициализируем значения аккумуляторов скорости и ускорения
                obj.VelocAcc = VelocAcc;
                obj.AccelAcc = AccelAcc;
        end
        function [Output, VelocAcc, AccelAcc] = Step(obj, Inp1, Inp2) %#ok<INUSD>
        % Функция выполнения действий одного шага фильтра
        % Формулы соответствуют Kaplan: page 181, fig 5.20
        %
        % Входные параметры
        %   Inp1, Inp2 - для isDLL = true используется только Inp1 и оно
        %       равно значению с выхода дискриминатора DLL, для isDLL =
        %       false соответственно Inp1 - значение с выхода
        %       дискриминатора FLL и Inp2 - значение с выхода
        %       дискриминатора PLL
        % Выходные переменные
        %   Output - значение, которое нужно подать на NCO
        %   VelocAcc, AccelAcc - текущие значения аккумуляторов

            % Шаг фильтра
                switch obj.Order
                    case 1
                        Output = Inp1 * obj.CoefsDLL;

                    otherwise
                end
            
            % Значения аккумуляторов после выполнения шага фильтра
                AccelAcc = obj.AccelAcc;
                VelocAcc = obj.VelocAcc;
        end
    end
end

function Coefs = CalcCoefs(Order, Bn)
% Функция расчёта коэффициентов фильтра по заданным порядку фильтра
% и его полосе. Формулы взяты из Kaplan: page 180, Table 5.6.
%
% Входные параметры
%   Order - 1 | 2 | 3;
%   Bn    - значение >0 в Гц.
% Выходные переменные
%   Coefs - массив [1xOrder] коэффициентов фильтра.

    switch Order
        case 1
            Coefs = Bn * 4;

        otherwise
    end
end
