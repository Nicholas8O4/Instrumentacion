clear
clear all
load m5d.mat

x = Acceleration.X;
y = Acceleration.Y;
z = Acceleration.Z;
timestamp = Acceleration.Timestamp;

combinedAccel = sqrt(x.^2 + y.^2 + z.^2);

t = zeros(size(timestamp));
for n = 1 : length(timestamp)
    t(n) = seconds(timestamp(n) - timestamp(1));
end

% --- Graficar las componentes ---
plot(t, x, 'm')
hold on
plot(t, y, 'b')
plot(t, z, 'c')
hold off
xlabel('Tiempo (s)')
ylabel('Posición')
legend('X','Y','Z')
grid on
title('Componentes de la aceleración')

% --- Graficar cada componente por separado ---
figure(1)
plot(t, x)
xlabel('Tiempo (s)')
ylabel('Posición X')

figure(2)
plot(t, y)
xlabel('Tiempo (s)')
ylabel('Posición Y')

figure(3)
plot(t, z)
xlabel('Tiempo (s)')
ylabel('Posición Z')

figure(4)
plot(t, combinedAccel)
xlabel('Tiempo (s)')
ylabel('Aceleración (m/s^2)')
title('Magnitud de aceleración')

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

% --- Calcular media y varianza cada 520 datos ---
tamano_segmento = 270;
num_segmentos = floor(length(combinedAccel) / tamano_segmento);

medias = zeros(num_segmentos, 1);
varianzas = zeros(num_segmentos, 1);
inicio_t = zeros(num_segmentos, 1);
fin_t = zeros(num_segmentos, 1);

for k = 1:num_segmentos
    idx_inicio = (k-1)*tamano_segmento + 1;
    idx_fin = k * tamano_segmento;
    segmento = combinedAccel(idx_inicio:idx_fin);

    medias(k) = mean(segmento);
    varianzas(k) = var(segmento);
    inicio_t(k) = t(idx_inicio);
    fin_t(k) = t(idx_fin);
end

% --- Crear tabla con resultados ---
estadisticas = table((1:num_segmentos)', medias, varianzas, inicio_t, fin_t, ...
    'VariableNames', {'Numero De Cara', 'Media De La Gravedad (m/s^2)', 'Varianza', 'Inicio_tiempo', 'Fin_tiempo'});

% --- Mostrar tabla ---
disp(estadisticas)

% --- (Opcional) Guardar en un archivo Excel ---
writetable(estadisticas, 'estadisticas_aceleracion.xlsx')
