% DAYANA ANALY PACHECO BAÑUELOS
% TAREA 4

% MENU DE OPCIONES
clc; close all; clear all;
opcionp=0;
while opcionp~=6,
 opcionp=menu('TAREA 4: ECG/FOURIER',...
     '1. Cargar & filtrar ECG',...
     '2. Fourier armonicas manualmente',...
     '3. Findchangepts',...
     '4. Fourier armonicas automaticamente',...
     '5. Detección de onda PQRS',...
     '6. Salir')

% Opcion 1. Cargar & Filtrar ECG--------------------------------------------------------------------------
if opcionp==1,
file = ['DAPB.txt'];  %Se carga el ECG propio
ECG = load(file); 

%Parametros de la señal de ECG / Variables        
F_muestreo = 1000; % Frecuencia de muestreo de la señal va ser igual a 1000Hz
Freclinea = 60; % Define la frecuencia de línea de 60 Hz
Numuestras = size(ECG,1); % Variable que menciona que las muestras son del tamaño de la primera muetra del ECG

%Calculo de la transformada de fourier para el Filtro de 60Hz
Resultadofft= fft(ECG); % Se aplica transformada de fourier rapida a ECG
Resultadofft(1:round(length(Resultadofft)*Freclinea/Numuestras))=0;  %Elimina las frecuencias bajas (cercanas a 60 Hz) 
Resultadofft(end - round(length(Resultadofft)*Freclinea/Numuestras):end) = 0;
EcgF = real(ifft(Resultadofft));%Variable con la Señal filtrada en el dominio de tiempo

%Grafica las señales ecg vs filtrado 60 hz
figure('Name', 'Cargar ECG & Filtrar', 'NumberTitle', 'on');
subplot (1,2,1); plot(ECG(1:1000), 'color', 'k', 'LineWidth', 2); %ECG ORIGINAL
title('ECG'); xlabel('ms'); ylabel('mV'); hold on;
plot(EcgF(1:1000), 'color', 'g', 'LineWidth', 2); %ECG FILTRO 60 HZ
legend('ECG ORIGINAL', 'ECG FILTRADO 60 HZ'); hold off;
subplot (1,2,2); bar(abs(Resultadofft)/(Numuestras/2)); title('FFT');

% Opcion 2. Fourier armonicas manualmente-------------------------------------------------------------------
elseif opcionp==2,
%Parametros de la señal     
F_muestreo = 1000; % Frecuencia de muestreo de la señal va ser igual a 1000Hz
Freclinea = 60; % Define la frecuencia de línea de 60 Hz
Numuestras = size(ECG,1); % Variable que menciona que las muestras son del tamaño de la primera muetra del ECG

%Calculo de la transformada de fourier en el ECG filtrado
Resultadofft= fft(ECG); % Se aplica transformada de fourier rapida a ECG
Resultadofft(1:round(length(Resultadofft)*Freclinea/Numuestras))=0;  %Elimina las frecuencias bajas
Resultadofft(end - round(length(Resultadofft)*Freclinea/Numuestras):end) = 0;
correctedL = real(ifft(Resultadofft));%Variable con la Señal filtrada en el dominio de tiempo

% Obtención de la señal & Filtro manual
Sigkill=Resultadofft
Sigkill(1:1) = [0];% Para eliminar la componente de frecuencia más baja
Sigkill(2000:59000)=zeros(57001,1); % Establece en cero un rango de frecuencias del espectro de Fourier
EcgF=real((ifft(Sigkill))); %Se aplica filtro pasa bajas/ transformada inversa de fourier
%Graficas de las armonicas fourier manualmente
subplot (1,2,1); plot(correctedL(1:1000),'color','k','LineWidth',2); 
title('ECG'); xlabel('ms'); ylabel('mV'); hold on;
plot(EcgF(1:1000),'color','c','LineWidth',2); %ECG FILTRADO
legend('ECG FILTRADO','ECG FILTRADO MANUAL'); hold off;
subplot (1,2,2); bar(abs(Sigkill)/(Numuestras/2)); title('FFT');

% Opcion 3. Findchangepts-------------------------------------------------------------------
elseif opcionp==3,
%Parte de Findpeaks
%Parametros definidos
Distancia = 500; %Distancia entre picos es de 500 muestras
T = Numuestras / F_muestreo; % Se define el periodo de muestreo
t = 1/F_muestreo:1/F_muestreo:T; % Vector de tiempo
Arranque = 1000; Final=3850; %Modifica el tamaño del numero de muestras ploteados

%ECG original señal
figure('Name','ECG & FindPeaks', 'NumberTitle', 'on');
subplot (2,1,1); plot(correctedL(Arranque:Final),'color','#7E2F8E', 'LineWidth',1.5);
legend('ECG ORIGINAL'); 

%Se buscan los picos respectivos de las ondas R
subplot(2,1,2);[pksl,locsl]=findpeaks(correctedL(Arranque:Final),'MinPeakDistance',Distancia);
findpeaks(correctedL(Arranque:Final),'MinPeakDistance',Distancia);legend('FindPeaks');
plot(Arranque:Final, correctedL(Arranque:Final),'color', '#4DBEEE', 'LineWidth',1.5);hold on; % Dibuja los datos de la señal
plot(locsl + Arranque - 1, pksl, 'v', 'MarkerFaceColor','b', 'MarkerSize', 10); % Marcadores de picos 
legend('Señal ECG', 'FindPeaks'); hold off;

%Parte de Findchangepts
figure('Name','Findchangepts','NumberTitle', 'on');
ipt=findchangepts(correctedL(Arranque:Final),MaxNumChanges= size(pksl,1)*2);
findchangepts(correctedL(Arranque:Final), MaxNumChanges= size(pksl,1)*2);
legend('findchangepts');
figure ('Name','Spectrograma & Findchangepts');
[s,f,t,pxx]= spectrogram(correctedL(Arranque:Final),128,120,128,F_muestreo);
findchangepts(pxx,MaxNumChanges=(size(pksl,1)*2))

% Opcion 4. Fourier armonicas Autimatico-------------------------------------------------------------------
elseif opcionp==4,
%Para fourier automatico
%Definicion de parametros
Muestrasmitad=round(Numuestras/2);
cambios=20;
%Calculo armonicas de fourier
y1=abs(fft(correctedL));
Muescambitad=findchangepts(y1(1:Muestrasmitad),'MaxNumChanges',cambios); %Función findchangepts para detectar cambios en el dominio de Fourier
figure('Name','Detección de los cambios'); sgtitle ('Fourier');
%Transformada de fourier & filtro
y=(fft(correctedL));
y(1:min(Muescambitad))=zeros(min(Muescambitad),1);
y(max(Muescambitad):Muestrasmitad)=zeros(Muestrasmitad-max(Muescambitad)+1,1);
y(Muestrasmitad:Muestrasmitad+max(Muescambitad))=zeros(max(Muescambitad)+1,1);
y(Numuestras-min(Muescambitad):Numuestras)=zeros(min(Muescambitad)+1,1);
z1=real(ifft(y));
%Creacion de figuras de los cambios
subplot(2,3,[1 3]); hold on; plot(correctedL(1:500),'Color','k','LineWidth',2);
plot(z1(1:500),'color','#0072BD','LineWidth',2); 
plot(z1(1:500),'color','#77AC30','LineWidth',2);
legend('Original','Filtrada con "findchancepts"','Filtrada Manual');
sgtitle('Detección de los cambios');hold off;
subplot(2,3,4); y=abs(fft(correctedL))/(Numuestras/2);ymax=max(y);
bar(y);xlabel('Frecuencia');title('Transformada de Fourier original');ylabel('Amplitude');
subplot(2,3,5);y2=abs(fft(z1))/(Numuestras/2);ymax=max(y);
bar(y2);xlabel('Frecuencia');title('Transformada de Fourier finchangepts');xlabel('Amlitude'); 

%FFT Manual faltante en cambios
Sig=y; %Se define la variable para el filtro con fft
Sig(1:1) = [0]; Sig(2000:59000)=zeros(57001,1); %Se eliminan frecuencias manualmente
FiltroECG=real((ifft(Sig))); %Se aplica filtro pasa bajas/ transformada inversa de fourier
subplot (2,3,6); bar(abs(Sig)/(Numuestras/2)); title('Transformada de Fourier manual');
 
% Opcion 5. Detección de onda PQRS-------------------------------------------------------------------
elseif opcionp==5,
 % Definir parámetros
 ECG_Filtrado = EcgF;  
 Distancia_RR = 500; % Distancia mínima entre picos R
 Umbral_R = 0.5 * max(ECG_Filtrado);  % Umbral relativo para medir duración de la onda R
 [peaks_R, locs_R] = findpeaks(ECG_Filtrado, 'MinPeakDistance', Distancia_RR);% Detectar los picos R
 duracion_completa = []; % Para almacenar la duración de cada onda PQRS

    % Búsqueda de la duración de la onda completa (PQRS)
    for i = 1:length(locs_R)
        % Encontrar el inicio de la onda Q
        inicio_Q = find(ECG_Filtrado(1:locs_R(i)) < Umbral_R, 1, 'last');  % Último cruce del umbral antes del pico R
        
        % Encontrar el fin de la onda S
        if i < length(locs_R) % Para evitar errores en el último pico
            siguiente_R = locs_R(i + 1);
            fin_S = siguiente_R + find(ECG_Filtrado(siguiente_R:end) < Umbral_R, 1, 'first');  % Primer cruce del umbral después del siguiente pico R
            duracion_onda_completa = fin_S - inicio_Q; % Duración de la onda completa
            duracion_completa = [duracion_completa; duracion_onda_completa]; % Almacenar duración
        end
    end

 % Calcular la distancia mínima y máxima entre las ondas completas
 min_distancia = min(duracion_completa);  % Distancia mínima entre ondas
 max_distancia = max(duracion_completa);  % Distancia máxima entre ondas
 num_total_complejos = length(duracion_completa); % Número total de ondas completas

 % Graficar solo una onda PQRS completa (la primera)
 inicio_onda = max(1, locs_R(1) - 210);  % Aproximadamente 210 muestras antes del primer pico R
 fin_onda = min(length(ECG_Filtrado), locs_R(1) + 200);  % Aproximadamente 200 muestras después del primer pico R

 figure('Name', 'Detección de Onda Completa PQRS', 'NumberTitle', 'on');
 plot(ECG_Filtrado(inicio_onda:fin_onda), 'color', '#7E2F8E','LineWidth', 1.5); hold on;
 plot(locs_R(1) - inicio_onda + 1, peaks_R(1), 'v', 'MarkerFaceColor', '#7E2F8E');
 legend('Complejo pqrs', 'Primer Pico R');

 % Añadir texto con la información solicitada en la gráfica
 texto = sprintf('Onda distancia: %d - %d,  Ondas Totales: %d', min_distancia, max_distancia, num_total_complejos);
 text(100, peaks_R(1) - 100, texto, 'FontSize', 10, 'Color', 'k', 'HorizontalAlignment', 'center');
 xlabel('Muestras'); ylabel('Amplitud'); title('Detección de una Onda Completa PQRS');hold off;

 % Mostrar la distancia mínima y máxima, y el número de ondas en la consola
 disp('Distancia mínima entre ondas completas (en muestras):');
 disp(min_distancia);
 disp('Distancia máxima entre ondas completas (en muestras):');
 disp(max_distancia);
 disp('Número total de ondas completas detectadas:');
 disp(num_total_complejos);
end
end
