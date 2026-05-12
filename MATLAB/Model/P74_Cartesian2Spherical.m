function [Lat, Lon, Alt] = P74_Cartesian2Spherical(xyz, Params)
% Функция выполняет преобразование координат из прямоугольной системы
% координат в сферическую систему координат
%
% Входные параметры
%   xyz - массив значений x, y, z - координат в прямоугольной системе 
%     координат в метрах.
%
% Выходные параметры
%   Latitude, Longitude, Altitude - широта и долгота в радианах, высота в
%     метрах.

    AlgType = Params.P74_Cartesian2Spherical.AlgType;
        % 0 - по стандарту РФ
        % 1 - по книге
        
    EllipseType = Params.P74_Cartesian2Spherical.EllipseType;
        % 0 - WGS84
        % 1 - ПЗ-90
        % 2 - Красовский - 1942
        
    % ellipticity
        switch EllipseType
            case 0
                ell = 1/298.257223563; % WGS84
            case 1
                ell = 1/298.257839303; % ПЗ-90
            case 2
                ell = 1/298.3; % Красовский - 1942
        end
        
    % semi-major axe of the earth
        switch EllipseType
            case 0
                a = 6378137; % WGS84
            case 1
                a = 6378136; % ПЗ-90
            case 2
                a = 6378245; % Красовский - 1942
        end

    % Координаты в прямоугольной системе
        x = xyz(1);
        y = xyz(2);
        z = xyz(3);

 % Вычисление долготы
    Lon = atan2(y, x);
    
% Вычисление широты и высоты в зависимости от выбранного алгоритма
    n = size(x, 1);
    Lat = zeros(n, 1);
    Alt = zeros(n, 1);
    
% Квадрат эксцентриситета
    e2 = 2*ell - ell^2;

if AlgType == 0
% Алгоритм по стандарту РФ (итерационный)
    for i = 1:n
        % Начальное приближение широты
            p = sqrt(x(i)^2 + y(i)^2);
            Lat0 = atan2(z(i), p * (1 - e2));
        % Итерационный процесс
        for iter = 1:10
            N = a / sqrt(1 - e2 * sin(Lat0)^2);
            h = p / cos(Lat0) - N;
            Lat_new = atan2(z(i), p * (1 - e2 * N / (N + h)));
            
            if abs(Lat_new - Lat0) < 1e-12
                Lat0 = Lat_new;
                break;
            end
            Lat0 = Lat_new;
        end
        
        Lat(i) = Lat0;
        N = a / sqrt(1 - e2 * sin(Lat(i))^2);
        Alt(i) = p / cos(Lat(i)) - N;
    end
    
else
% Алгоритм по книге (метод Боуринга - прямой без итераций)
    for i = 1:n
        p = sqrt(x(i)^2 + y(i)^2);
        
        % Вычисление геодезической широты
            r = sqrt(p^2 + z(i)^2);
            u = atan2(z(i) * (1 - ell), p);
        
        % Вспомогательные величины
            sin_u = sin(u);
            cos_u = cos(u);
            
            Lat0 = atan2(z(i) + e2 * a * sin_u^3, ...
                        p - e2 * a * cos_u^3);
        
        % Уточнение для лучшей точности (один-два шага)
        for iter = 1:3
            N = a / sqrt(1 - e2 * sin(Lat0)^2);
            h = p / cos(Lat0) - N;
            Lat_new = atan2(z(i), p * (1 - e2 * N / (N + h)));
            
            if abs(Lat_new - Lat0) < 1e-12
                Lat0 = Lat_new;
                break;
            end
            Lat0 = Lat_new;
        end
        
        Lat(i) = Lat0;
        N = a / sqrt(1 - e2 * sin(Lat(i))^2);
        Alt(i) = p / cos(Lat(i)) - N;
    end
end
