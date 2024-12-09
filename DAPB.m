% WAVELETS & ECG
clc;close all;

% Se generan todas las wavelets--------------------------------------------------------
wName = ["haar" "db4" "sym4" "coif4" "bior3.1" "rbio4.4" "meyr" "dmey" "gaus3" "mexh" ...
         "morl" "cgau4" "shan1-1.3" "fbsp2-0.5-1" "cmor1.5-1" "fk14"];
itr = 10; %Iteraciones de las wavelets

% Cargar ECG--------------------------------------------------------------------------
data = readmatrix(['NRH_2023A.txt']);
ECG = data(:); % Datos a vector columna

% Normalización del ECG
ECG = ECG / max(abs(ECG));

% Parametros de ventana del complejo PQRST--------------------------------------------
fs = 1000; % Frecuencia de muestreo 
Duracion = 0.6; % Duración segmento del complejo (seg)
Iniciopqrs = 1150;  % Comienzo del complejo PQRST (muestras)
finpqrs = Iniciopqrs + Duracion * fs - 1;
Ventanapqrs = ECG(Iniciopqrs:finpqrs);

figure('Name', 'WAVELETS VS ECG', 'NumberTitle', 'off');
filas = 4; columnas = 4; % Subplots (4x4 para 16 wavelets)

% Procesamiento de las wavelets------------------------------------------------------
for i = 1:length(wName) %Desde la primera hasta el tamaño de las wavelets definidas
    wname = wName(i);

    % Calculo función wavelet madre
    try
% phi: valor de la wavelet/ psi: Aproximación wavelet (vector)/ xval: Aproximaciones de la función escala y wavelets (vector) 
        [phi, psi, xval] = wavefun(wname, itr); % Tres salidas para wavelets
    catch
        try
            [psi, xval] = wavefun(char(wname), itr); % Intentento calculo la wavelet madre
        catch % Si no muestra una aviso de alerta para el usuario
            warning(['La wavelet ' char(wname) 'No es compatible o requiere otro toolbox :(']);
            continue;
        end
    end

    % Revisión si la wavelet es compleja
    if ~isreal(psi)
        psi = abs(psi); % Obtiene su magnitud
    end

% Caracteritcas de la wavelet & complejo pqrs
% Amplitud---------------------------------------------------------------------------------------------------
Wevelet_adaptado = psi * (max(Ventanapqrs) / max(abs(psi))); % Escalar la wavelet para tener misma amplitud del complejo

%Longitud
% Adaptar misma longitud que el complejo PQRS
Wevelet_adaptado = interp1(1:length(Wevelet_adaptado), Wevelet_adaptado, linspace(1, length(Wevelet_adaptado), length(Ventanapqrs)), 'linear', 0);

% Correlación cruzada-------------------------------------------------------------------------------
correlacion = xcorr(Ventanapqrs, Wevelet_adaptado, 'coeff');
[maxcorr, lag_idx] = max(abs(correlacion)); % Encontrar la mejor alineación
best_shift = lag_idx - length(Ventanapqrs);

% Alinear la wavelet rescalada con el complejo PQRS-----------------------------------------------
CentrarW = circshift(Wevelet_adaptado, best_shift);

% Graficas (16 totales)
subplot(filas, columnas, i);
plot(Ventanapqrs, 'k', 'DisplayName', 'Complejo Pqrs'); hold on;
plot(CentrarW, 'color', '#9be175','LineWidth', 1.5, 'DisplayName', [char(wname), ' Wavelet']);
ylim([-1.2 * max(abs(Ventanapqrs)), 1.2 * max(abs(Ventanapqrs))]); xlim([1, length(Ventanapqrs)]);
title([char(wname), ' / Corr Max: ', num2str(maxcorr, '%.2f')]);
xlabel('Muestras'); ylabel('Amplitud'); legend('show'); hold off;
end

% Título general de la figura
sgtitle('Comparativa Wavelets vs ECG');

