clearvars; close all; clc;
load 1.mat

gravedad = 9.81;
duracionMinima = 1.0;
ventanaSeg = 0.5;
umbralMovstd = 0.08;
umbralOutlier = 3;

ax = Acceleration.X;
ay = Acceleration.Y;
az = Acceleration.Z;
tiempoMarca = Acceleration.Timestamp;
tiempo = seconds(tiempoMarca - tiempoMarca(1));
A = [ax(:), ay(:), az(:)];
N = size(A,1);

dt = median(diff(tiempo));
if dt <= 0, dt = 0.01; end
fs = 1/dt;

aceleracionCombinada = sqrt(ax.^2 + ay.^2 + az.^2);
ventanaMuestras = max(1, round(ventanaSeg * fs));
movstd_acc = movstd(aceleracionCombinada, ventanaMuestras);
mascaraEstable = movstd_acc < umbralMovstd;

d = diff([0; mascaraEstable(:); 0]);
inicios = find(d == 1);
fines   = find(d == -1) - 1;

if isempty(inicios) || isempty(fines)
    duraciones = [];
else
    duraciones = tiempo(fines) - tiempo(inicios);
end

mantener = duraciones >= duracionMinima;
inicios = inicios(mantener);
fines   = fines(mantener);
duraciones = duraciones(mantener);

numTramos = length(inicios);
G = zeros(numTramos,3);
for k = 1:numTramos
    idx = inicios(k):fines(k);
    G(k,:) = [mean(ax(idx)), mean(ay(idx)), mean(az(idx))];
end

figure;
plot(tiempo, aceleracionCombinada, 'b-'); hold on;
yl = ylim;
for k = 1:numTramos
    patch([tiempo(inicios(k)) tiempo(fines(k)) tiempo(fines(k)) tiempo(inicios(k))], ...
          [yl(1) yl(1) yl(2) yl(2)], [0.9 0.9 0.9], 'EdgeColor','none', 'FaceAlpha', 0.3);
end
plot(tiempo, aceleracionCombinada, 'b-');
xlabel('Tiempo (s)'); ylabel('Aceleración (m/s^2)');
title('Magnitud de aceleración');
grid on; hold off;

if numTramos > 0
    X0 = zeros(numTramos,1);
    Y0 = zeros(numTramos,1);
    Z0 = zeros(numTramos,1);
    U  = G(:,1);
    V  = G(:,2);
    W  = G(:,3);

    figure;
    quiver3(X0, Y0, Z0, U, V, W, 0, 'LineWidth', 2, 'MaxHeadSize', 0.6);
    xlabel('X'); ylabel('Y'); zlabel('Z');
    axis equal; grid on;
    title('Vectores promedio');

    hold on;
    [sx, sy, sz] = sphere(40);
    r = gravedad;
    surf(r*sx, r*sy, r*sz, 'FaceAlpha', 0.05, 'EdgeColor','none');
    hold off;
end

intervalos = zeros(numTramos,2);
for k=1:numTramos
    intervalos(k,:) = [inicios(k), fines(k)];
end

if isempty(intervalos)
    numCaras = 20;
    tiempoTotal = tiempo(end);
    tiempoCara = tiempoTotal / numCaras;
    intervalos = zeros(numCaras,2);
    for i=1:numCaras
        s = find(tiempo >= (i-1)*tiempoCara, 1, 'first');
        e = find(tiempo < i*tiempoCara, 1, 'last');
        if isempty(s), s = 1; end
        if isempty(e), e = N; end
        intervalos(i,:) = [s,e];
    end
end

numIntervalos = size(intervalos,1);
A_prom = nan(numIntervalos,3);
conteo_intervalo = zeros(numIntervalos,1);
for i=1:numIntervalos
    s = intervalos(i,1); e = intervalos(i,2);
    idx = s:e;
    A_prom(i,:) = mean(A(idx,:),1);
    conteo_intervalo(i) = numel(idx);
end

magnitud_prom_inicial = sqrt(sum(A_prom.^2,2));

validos = ~any(isnan(A_prom),2);
X = A_prom(validos,:);

if size(X,1) < 10
    error('No hay suficientes intervalos.');
end

xv = X(:,1); yv = X(:,2); zv = X(:,3);
D = [xv.^2, yv.^2, zv.^2, xv.*yv, xv.*zv, yv.*zv, xv, yv, zv, ones(size(xv))];

[~, ~, V] = svd(D,0);
coef = V(:,end);

A11 = coef(1); A22 = coef(2); A33 = coef(3);
A12 = coef(4)/2; A13 = coef(5)/2; A23 = coef(6)/2;
B1  = coef(7)/2; B2 = coef(8)/2; B3 = coef(9)/2;
C   = coef(10);

A_cuad = [A11 A12 A13; A12 A22 A23; A13 A23 A33];
b_cuad = [B1; B2; B3];

centro = - (A_cuad \ b_cuad);
R = centro' * A_cuad * centro - C;

if R <= 0 || any(eig(A_cuad) <= 0)
    M_inicial = eye(3);
    b_inicial = zeros(3,1);
else
    L = sqrtm(A_cuad);
    M_inicial = (gravedad / sqrt(R)) * L;
    b_inicial = - M_inicial * centro;
end

p0 = [reshape(M_inicial,9,1); b_inicial];

indices_validos = find(validos);
A_prom_validos = A_prom(validos,:);

resfun = @(p) residuos_aceleracion(p, A_prom_validos, gravedad);

opciones = optimoptions('lsqnonlin','Display','iter','TolFun',1e-12,'TolX',1e-12,'MaxIter',1000,'Algorithm','levenberg-marquardt');
try
    p_estimado = lsqnonlin(resfun, p0, [], [], opciones);
catch
    p_estimado = p0;
end

M_estimado = reshape(p_estimado(1:9),3,3);
b_estimado = p_estimado(10:12);

A_cal = (M_estimado * A_prom_validos.' )' + repmat(b_estimado(:)', size(A_prom_validos,1), 1);
residuos = sqrt(sum(A_cal.^2,2)) - gravedad;
res_std = std(residuos);
res_med = mean(residuos);

idx_out = abs(residuos - res_med) > umbralOutlier * res_std;
if any(idx_out)
    mantener = ~idx_out;
    A_prom_refinado = A_prom_validos(mantener,:);
    resfun2 = @(p) residuos_aceleracion(p, A_prom_refinado, gravedad);
    try
        p_ref = lsqnonlin(resfun2, p_estimado, [], [], opciones);
        p_estimado = p_ref;
        M_estimado = reshape(p_estimado(1:9),3,3);
        b_estimado = p_estimado(10:12);
        A_cal_final = (M_estimado * A_prom_validos.' )' + repmat(b_estimado(:)', size(A_prom_validos,1), 1);
        residuos = sqrt(sum(A_cal_final.^2,2)) - gravedad;
    catch
        A_cal_final = A_cal;
    end
else
    A_cal_final = A_cal;
end

magnitud_prom_final = nan(numIntervalos,1);
magnitud_prom_final(validos) = sqrt(sum(((M_estimado * A_prom(validos,:).')' + repmat(b_estimado(:)', sum(validos), 1)).^2,2));

A_cal_todo = (M_estimado * A.' )' + repmat(b_estimado(:)', N, 1);
magnitud_todo_inicial = sqrt(sum(A.^2,2));
magnitud_todo_final  = sqrt(sum(A_cal_todo.^2,2));

rms_inicial = sqrt(mean((magnitud_prom_inicial(validos) - gravedad).^2));
rms_final   = sqrt(mean((magnitud_prom_final(validos)  - gravedad).^2));

parametros_cal.M = M_estimado;
parametros_cal.b = b_estimado;
save('calibracion_acelerometro_mejorada.mat','parametros_cal');

Cara = (1:numIntervalos)';
Media_inicial = magnitud_prom_inicial;
Media_final = magnitud_prom_final;
Inicio_idx = intervalos(:,1);
Fin_idx = intervalos(:,2);

T = table(Cara, Media_inicial, Media_final, Inicio_idx, Fin_idx, conteo_intervalo, ...
          'VariableNames', {'Cara','Media_inicial','Media_final','Inicio_idx','Fin_idx','Conteo'});
disp(T);
writetable(T,'estadisticas_aceleracion_calibrada_mejorada.xlsx','Sheet','Por_Intervalo');

figure;
hold on;
yline(gravedad,'k--','LineWidth',1.2);
plot(1:numIntervalos, Media_inicial, 'o-','LineWidth',1.2);
plot(1:numIntervalos, Media_final, 's-','LineWidth',1.2);
xlabel('Intervalo'); ylabel('Magnitud (m/s^2)');
legend('g','Inicial','Final'); grid on; hold off;

figure;
histogram(magnitud_todo_inicial, 60, 'Normalization','pdf'); hold on;
histogram(magnitud_todo_final,  60, 'Normalization','pdf'); hold off;
xlabel('Magnitud (m/s^2)'); ylabel('Densidad');
legend('Inicial','Final'); grid on;

sub = 1:round(max(1,floor(N/2000))):N;
figure;
plot3(A(sub,1), A(sub,2), A(sub,3), '.', 'MarkerSize', 6); hold on;
plot3(A_cal_todo(sub,1), A_cal_todo(sub,2), A_cal_todo(sub,3), '.', 'MarkerSize', 6);
[sx,sy,sz] = sphere(60);
surf(gravedad*sx, gravedad*sy, gravedad*sz, 'FaceAlpha', 0.10, 'EdgeColor','none');
xlabel('X'); ylabel('Y'); zlabel('Z'); axis equal; grid on;
legend('Inicial','Final','Esfera g'); hold off;

err_inicial = Media_inicial - gravedad;
err_final  = Media_final - gravedad;
figure;
plot(1:numIntervalos, err_inicial, 'o-'); hold on;
plot(1:numIntervalos, err_final, 's-'); hold off;
xlabel('Intervalo'); ylabel('Error (m/s^2)');
legend('Inicial','Final'); grid on;

figure;
bar(1:sum(validos), residuos); hold on;
yline(3*std(residuos),'r--','LineWidth',1.2);
xlabel('Intervalos validos'); ylabel('Residual');
grid on; hold off;

function r = residuos_aceleracion(p, A_prom_validos, gravedad)
    M = reshape(p(1:9),3,3);
    b = p(10:12);
    Acal = (M * A_prom_validos.' )' + repmat(b(:)', size(A_prom_validos,1), 1);
    normas = sqrt(sum(Acal.^2,2));
    r = normas - gravedad;
end
