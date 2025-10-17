clearvars; 
close all; 
clc;
load m5d.mat   

% --- Extraer señales y tiempo ---
x = Acceleration.X;
y = Acceleration.Y;
z = Acceleration.Z;
timestamp = Acceleration.Timestamp;

% Tiempo en segundos desde el inicio
t = seconds(timestamp - timestamp(1));

% Magnitud combinada (magnitud de la aceleración)
g_mag = sqrt(x.^2 + y.^2 + z.^2);

% --- 1) Graficar componentes y magnitud ---
figure('Name','Componentes y Magnitud','NumberTitle','off');
subplot(3,1,1)
plot(t, x, 'm'); xlabel('Tiempo (s)'); ylabel('Aceleración X'); grid on; title('Componente X');

subplot(3,1,2)
plot(t, y, 'b'); xlabel('Tiempo (s)'); ylabel('Aceleración Y'); grid on; title('Componente Y');

subplot(3,1,3)
plot(t, z, 'g'); xlabel('Tiempo (s)'); ylabel('Aceleración Z'); grid on; title('Componente Z');

figure('Name','Magnitud de la aceleración','NumberTitle','off');
plot(t, g_mag);
xlabel('Tiempo (s)'); ylabel('Magnitud (m/s^2)'); title('Magnitud de la aceleración'); grid on;

figure(5)
plot(t, x, 'm')
hold on
plot(t, y, 'b')
plot(t, z, 'g')
hold off
xlabel('Tiempo (s)')
ylabel('Posiciónes')
legend('X','Y','Z')
grid on
title('Componentes de la aceleración')
% --- 2) Promedio por segmentos de tamaño fijo (tu enfoque original) ---
tamano_segmento = 270;   % <-- ajustable
num_segmentos = floor(length(g_mag) / tamano_segmento);

medias_segmento = zeros(num_segmentos,1);
varianzas_segmento = zeros(num_segmentos,1);
inicio_t = zeros(num_segmentos,1);
fin_t = zeros(num_segmentos,1);

for k = 1:num_segmentos
    idx_inicio = (k-1)*tamano_segmento + 1;
    idx_fin = k * tamano_segmento;
    segmento = g_mag(idx_inicio:idx_fin);

    medias_segmento(k) = mean(segmento);
    varianzas_segmento(k) = var(segmento);
    inicio_t(k) = t(idx_inicio);
    fin_t(k) = t(idx_fin);
end

% Figura: medias por segmento (índice de segmento)
figure('Name','Promedio por segmento (fijo)','NumberTitle','off');
plot(1:num_segmentos, medias_segmento, 'o-','LineWidth',1.2);
hold on; 
yline(9.81,'r--','LineWidth',1.0);
xlabel('Segmento'); ylabel('Aceleración promedio (m/s^2)');
title(sprintf('Promedio por segmento (tamaño = %d muestras)', tamano_segmento));
legend('Medias por segmento','Referencia g=9.81','Location','best');
grid on; hold off;

% --- 3) Promedio por "cara" del icosaedro (división por tiempo) ---
numCaras = 20;          % <-- ajustable: número de caras del icosaedro
tiempoTotal = t(end);
tiempoCara = tiempoTotal / numCaras;

g_medidas_porCara = nan(1, numCaras);
for i = 1:numCaras
    idx = (t >= (i-1)*tiempoCara) & (t < i*tiempoCara);
    if any(idx)
        g_medidas_porCara(i) = mean(g_mag(idx));
    else
        g_medidas_porCara(i) = NaN; % si no hay datos en ese intervalo
    end
end

% Figura: medias por cara (índice de cara)
figure('Name','Promedio por cara (tiempo)','NumberTitle','off');
hold on;
yline(9.81, 'g--', 'LineWidth', 1.2);
plot(1:numCaras, g_medidas_porCara, 'b-o','LineWidth',1.2);
xlabel('Cara del icosaedro');
ylabel('Aceleración promedio (m/s^2)');
title('Magnitud promedio de la aceleración por cara del icosaedro');
legend('Gravedad de referencia (9.81 m/s^2)', 'Mediciones obtenidas','Location','best');
grid on; hold off;

% --- 4) Tabla con resultados y guardado ---
% Para las filas por segmento fijo
T_segmentos = table((1:num_segmentos)', medias_segmento, varianzas_segmento, inicio_t, fin_t, ...
    'VariableNames', {'Segmento', 'Media_m_s2', 'Varianza', 'Inicio_t_s', 'Fin_t_s'});

% Para las filas por cara (tiempo)
Cara = (1:numCaras)';
Media_porCara = g_medidas_porCara';
Inicio_t_cara = ((0:numCaras-1)') * tiempoCara;
Fin_t_cara = ((1:numCaras)') * tiempoCara;
T_caras = table(Cara, Media_porCara, Inicio_t_cara, Fin_t_cara, ...
    'VariableNames', {'Cara','Media_m_s2','Inicio_t_s','Fin_t_s'});

% Mostrar tablas en consola
disp('--- Estadísticas por segmento (fijo) ---');
disp(T_segmentos);
disp('--- Estadísticas por cara (tiempo) ---');
disp(T_caras);

% Guardar ambas tablas en un archivo Excel (dos hojas)
writetable(T_segmentos, 'estadisticas_aceleracion.xlsx', 'Sheet', 'Segmentos_fijos');
writetable(T_caras, 'estadisticas_aceleracion.xlsx', 'Sheet', 'Por_Cara_Tiempo');

disp('--- Gráficas y tablas generadas y guardadas en estadisticas_aceleracion.xlsx ---');
