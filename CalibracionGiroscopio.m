clearvars; close all; clc;

%% ========== CONFIGURACIÓN GENERAL ==========
archivo = 'Sis_3.mat';
if ~exist(archivo,'file'), error('No se encontró %s. Pon el archivo en la carpeta actual.', archivo); end
load(archivo); % debe contener Acceleration, AngularVelocity

g = 9.81;                    % gravedad
% Detector estático (rango + movstd)
rango_inicial = 0.20;        % m/s^2, detección amplia sobre señal calibrada (ajusta)
rango_final   = 0.20;        % m/s^2, detección final sobre señal calibrada (ajusta)
durMinSeg = 0.8;             % duración mínima del intervalo (s)
ventanaSeg_estable = 0.5;    % ventana para movstd inicial (s)
umbral_movstd_estable = 0.08;% umbral movstd inicial (m/s^2)
maxIntervalos = 20;          % máximo intervalos estáticos a usar
% Parámetros giroscopio (Tedaldi)
umbral_gyro = 0.40;          % rad/s para considerar intervalo quieto
maxIterLSQ = 800;

%% ========== EXTRAER SEÑALES ==========
ax = Acceleration.X(:); ay = Acceleration.Y(:); az = Acceleration.Z(:);
t_acc = seconds(Acceleration.Timestamp - Acceleration.Timestamp(1));
A_raw = [ax ay az];
N = size(A_raw,1);
dt_acc = median(diff(t_acc)); if dt_acc<=0, dt_acc = 0.01; end
fs_acc = 1/dt_acc;

wx = AngularVelocity.X(:); wy = AngularVelocity.Y(:); wz = AngularVelocity.Z(:);
t_gyro = seconds(AngularVelocity.Timestamp - AngularVelocity.Timestamp(1));
W = [wx wy wz];
dt_gyro = median(diff(t_gyro)); if dt_gyro<=0, dt_gyro = 0.01; end
fs_gyro = 1/dt_gyro;

%% =========================
% 1) DETECCIÓN INICIAL DE TRAMOS DE BAJA VARIACIÓN (candidatos)
%% =========================
winSamples = max(1, round(ventanaSeg_estable * fs_acc));
magnitud_comb = sqrt(sum(A_raw.^2,2));
movstd_acel = movstd(magnitud_comb, winSamples);
mask_estable = movstd_acel < umbral_movstd_estable;

d = diff([0; mask_estable(:); 0]);
inicio = find(d==1);
fin   = find(d==-1)-1;

if isempty(inicio)
    warning('No se detectaron tramos de baja variación con parámetros iniciales. Se intentará con todo el registro.');
    inicio = 1; fin = N;
end

% Filtrar por duración
duraciones = t_acc(fin) - t_acc(inicio);
keep = duraciones >= 0.25; % ventanas cortas permitidas para que SVD tenga muchas caras candidatas
inicio = inicio(keep); fin = fin(keep); duraciones = duraciones(keep);
nSeg = numel(inicio);
fprintf('Tramos de baja variación detectados: %d\n', nSeg);

% Promedios por tramo X (candidatos)
X_candidatos = nan(nSeg,3);
for k=1:nSeg
    idx = inicio(k):fin(k);
    X_candidatos(k,:) = mean(A_raw(idx,:),1);
end

%% =========================
% 2) AJUSTE ELIPSOIDE (acelerómetro) - algebraico + refinamiento
%% =========================
if size(X_candidatos,1) < 6
    warning('Pocos puntos para ajuste de elipsoide (%d). Se procederá pero la calibración puede ser débil.', size(X_candidatos,1));
end

xv = X_candidatos(:,1); yv = X_candidatos(:,2); zv = X_candidatos(:,3);
Dmat = [xv.^2, yv.^2, zv.^2, xv.*yv, xv.*zv, yv.*zv, xv, yv, zv, ones(size(xv))];
[~,~,Vd] = svd(Dmat,0);
coef = Vd(:,end);

A11 = coef(1); A22 = coef(2); A33 = coef(3);
A12 = coef(4)/2; A13 = coef(5)/2; A23 = coef(6)/2;
B1 = coef(7)/2; B2 = coef(8)/2; B3 = coef(9)/2; C = coef(10);

A_quad = [A11 A12 A13; A12 A22 A23; A13 A23 A33];
b_quad = [B1; B2; B3];
centro = - (A_quad \ b_quad);
Rval = centro' * A_quad * centro - C;

if Rval <= 0 || any(eig(A_quad) <= 0)
    warning('Elipsoide algebraico no válido; usando identidad como inicial.');
    M_inicial = eye(3); b_inicial = zeros(3,1);
else
    L = sqrtm(A_quad);
    M_inicial = (g / sqrt(Rval)) * L;
    b_inicial = - M_inicial * centro;
end

% Refinamiento no lineal
p0 = [reshape(M_inicial,9,1); b_inicial(:)];
resfun_acc = @(p) acc_residuals(p, X_candidatos, g);
opts_acc = optimoptions('lsqnonlin','Display','off','TolFun',1e-12,'TolX',1e-12,'MaxIter',1000,'Algorithm','levenberg-marquardt');
try
    p_est = lsqnonlin(resfun_acc, p0, [], [], opts_acc);
catch ME
    warning('Refinamiento acelerómetro falló: %s. Se usa inicial.', ME.message);
    p_est = p0;
end
M_est_acc = reshape(p_est(1:9),3,3);
b_est_acc = p_est(10:12);

% Aplicar calibración al conjunto
A_cal_all = (M_est_acc * A_raw.' ).' + repmat(b_est_acc.', N, 1);
magA_cal = vecnorm(A_cal_all,2,2);
fprintf('RMS de ||A_cal||-g (global): %.6g\n', sqrt(mean((magA_cal - g).^2)));

%% =========================
% 3) DETECCIÓN ROBUSTA DE INTERVALOS ESTÁTICOS (rango + movstd)
%% =========================
% Primera máscara amplia
tol = 0.6; % rango amplio para primera máscara
winSTD = round(0.25*fs_acc);
thrSTD = 0.06;

mask_range = abs(magA_cal - g) < tol;
movstd_cal = movstd(magA_cal, winSTD);
mask_stable = movstd_cal < thrSTD;
mask_static = mask_range & mask_stable;

d = diff([0; mask_static(:); 0]);
inicio_rob = find(d==1);
fin_rob   = find(d==-1)-1;
nSeg_rob = length(inicio_rob);
fprintf('Intervalos estáticos detectados (robusto, preliminar): %d\n', nSeg_rob);

% Calcular calidad y limitar a maxIntervalos
if nSeg_rob == 0
    error('No se detectaron intervalos estáticos con el detector robusto. Ajusta parámetros.');
end

calidad = zeros(nSeg_rob,1);
for k=1:nSeg_rob
    idx = inicio_rob(k):fin_rob(k);
    seg_mag = magA_cal(idx);
    err_g = mean(abs(seg_mag - g));
    est = std(seg_mag);
    dur = t_acc(fin_rob(k)) - t_acc(inicio_rob(k));
    calidad(k) = err_g*2 + est + 0.2/(dur+0.001);
end
[~, idxBest] = sort(calidad,'ascend');

sel = idxBest(1:min(maxIntervalos,length(idxBest))); % seleccionar los mejores
inicio_sel = inicio_rob(sel);
fin_sel   = fin_rob(sel);
% ordenar por tiempo
[inicio_sel, orden] = sort(inicio_sel);
fin_sel = fin_sel(orden);
numIntervals = numel(inicio_sel);
fprintf('Se usarán %d intervalos estáticos (mejores %d).\n', numIntervals, maxIntervalos);

% Construir A_cal_mean por intervalo (usado por Tedaldi)
A_cal_mean = nan(numIntervals,3);
for k=1:numIntervals
    idx = inicio_sel(k):fin_sel(k);
    A_cal_mean(k,:) = mean(A_cal_all(idx,:),1);
end

% Mostrar resumen simple (opcional)
fprintf('Resumen magnitudes por intervalo (media):\n');
for k=1:numIntervals
    fprintf(' %2d: |a|_mean=%.4f  dur=%.3fs\n', k, norm(A_cal_mean(k,:)), t_acc(fin_sel(k))-t_acc(inicio_sel(k)));
end

%% =========================
% 4) CREAR estructuras para la siguiente etapa
%    intervals: índices sobre t_acc (para compatibilidad con script de gyro)
%% =========================
intervals = [inicio_sel(:), fin_sel(:)];

% Guardar A_cal_all, intervals, A_cal_mean, t_acc para depuración/uso
save('interv_detectados.mat','intervals','A_cal_all','A_cal_mean','t_acc','magA_cal');

%% =========================
% 5) SINCRONIZAR CON GIROSCOPIO (índices sobre t_gyro por cada intervalo)
%% =========================
gyro_intervals = cell(numIntervals,1);
valid_gyro_interval = false(numIntervals,1);

for k = 1:numIntervals
    t_start = t_acc(intervals(k,1));
    t_end   = t_acc(intervals(k,2));
    idxg = find(t_gyro >= t_start & t_gyro <= t_end);
    gyro_intervals{k} = idxg;
    if numel(idxg) >= 2, valid_gyro_interval(k) = true; end
end

% Filtrar si algún intervalo quedó vacío
if ~all(valid_gyro_interval)
    keep = find(valid_gyro_interval);
    gyro_intervals = gyro_intervals(keep);
    A_cal_mean = A_cal_mean(keep,:);
    intervals  = intervals(keep,:);
    numIntervals = numel(keep);
    fprintf('Se eliminaron %d intervalos sin muestras de gyro. Quedan %d.\n', sum(~valid_gyro_interval), numIntervals);
end

if numIntervals < 2
    error('Muy pocos intervalos válidos sincronizados con gyro (%d). No se puede calibrar.', numIntervals);
end

%% =========================
% 6) FILTRADO ADICIONAL (Tedaldi): comprobar que giros durante intervalos sean bajos
%% =========================
valid2 = false(numIntervals,1);
for k=1:numIntervals
    idxg = gyro_intervals{k};
    if numel(idxg) < 2, valid2(k)=false; continue; end
    wnorm = vecnorm(W(idxg,:),2,2);
    mean_w = mean(wnorm);
    max_w  = max(wnorm);
    
    if (mean_w < umbral_gyro) && (max_w < 2*umbral_gyro)
        valid2(k) = true;
    else
        valid2(k) = false;
    end

end
% Aplicar filtro
gyro_intervals = gyro_intervals(valid2);
A_cal_mean = A_cal_mean(valid2,:);
intervals = intervals(valid2,:);
numIntervals = sum(valid2);
fprintf('Intervalos estáticos válidos tras filtro gyro: %d\n', numIntervals);
if numIntervals < 5
    warning('Tras filtro gyro quedan %d intervalos verdaderamente estáticos. Calibración puede ser menos precisa.', numIntervals);
end

%% =========================
% 7) PREPARAR función residual para giroscopio (usar función Tedaldi)
%% =========================
gyro_resfun = @(p) gyro_residuals_Tedaldi( ...
    p, W, t_gyro, gyro_intervals, A_cal_mean, intervals, t_acc, g );

%% =========================
% 8) AUTOMÁTICO: decidir si calibración completa o solo bias
%% =========================
% medir cambios entre intervalos
if numIntervals >= 2
    delta_g = zeros(numIntervals-1,1);
    for k=1:numIntervals-1
        gk = A_cal_mean(k,:).' / norm(A_cal_mean(k,:));
        gn = A_cal_mean(k+1,:).' / norm(A_cal_mean(k+1,:));
        delta_g(k) = norm(gk - gn);
    end
    max_delta_g = max(delta_g);
    mean_delta_g = mean(delta_g);
else
    max_delta_g = 0; mean_delta_g = 0;
end
wnorm_max_all = max(vecnorm(W,2,2));
fprintf('max_delta_g=%.4f  mean_delta_g=%.4f  max|w|=%.4f\n', max_delta_g, mean_delta_g, wnorm_max_all);

% Criterios (ajustables)
cond_no_rotation = ( max_delta_g < 0.25 ) || ( mean_delta_g < 0.10 ) || ( wnorm_max_all > 2.0 );

if cond_no_rotation
    % SOLO BIAS
    fprintf('Modo: SOLO BIAS (no suficiente rotación para estimar s,S confiables)\n');
    % promedio ponderado de muestras dentro de intervalos
    b_est = zeros(3,1); total=0;
    for k=1:numIntervals
        idxg = gyro_intervals{k};
        b_est = b_est + sum(W(idxg,:)).';
        total = total + numel(idxg);
    end
    if total>0, b_est = b_est / total; else b_est = nanmedian(W)'; end
    s_est = [1;1;1]; S_est = zeros(3,3); K_est = eye(3);
    W_cal = (W - b_est.');
    goto_only_bias = true;
else
    goto_only_bias = false;
end

%% =========================
% 9) SI HAY ROTACIÓN SUFICIENTE: calibración en 2 etapas (bias -> full)
%% =========================
if ~goto_only_bias
    % Etapa A: estimar bias solo inicialmente
    gyro_means_per_interval = nan(numIntervals,3);
    for i=1:numIntervals
        idxg = gyro_intervals{i};
        gyro_means_per_interval(i,:) = mean(W(idxg,:),1);
    end
    b0 = nanmedian(gyro_means_per_interval,1)'; b0(~isfinite(b0)) = 0;
    gyro_resfun_bias = @(b) gyro_residuals_Tedaldi([b;1;1;1;0;0;0], W, t_gyro, gyro_intervals, A_cal_mean, intervals, t_acc, g);
    opts_bias = optimoptions('lsqnonlin','Display','iter','TolFun',1e-10,'TolX',1e-10,'MaxIter',200,'Algorithm','levenberg-marquardt');
    try
        b_est_stage = lsqnonlin(gyro_resfun_bias, b0, [-0.5;-0.5;-0.5], [0.5;0.5;0.5], opts_bias);
    catch
        b_est_stage = lsqnonlin(gyro_resfun_bias, b0, [], [], opts_bias);
    end
    % Etapa B: calibración completa (full)
    p0_full = [b_est_stage; 1;1;1; 0;0;0];
    lb_full = [-0.5;-0.5;-0.5;  0.7;0.7;0.7;  -0.15;-0.15;-0.15];
    ub_full = [ 0.5; 0.5; 0.5;  1.3;1.3;1.3;   0.15; 0.15; 0.15];
    opts_full = optimoptions('lsqnonlin','Display','iter','MaxIter',maxIterLSQ,'TolFun',1e-10,'TolX',1e-10,'Algorithm','levenberg-marquardt');
    try
        p_est_full = lsqnonlin(gyro_resfun, p0_full, lb_full, ub_full, opts_full);
    catch ME
        warning('lsqnonlin full falló: %s. Intentando regularizado.', ME.message);
        lambda_scale_try = 1e-3; lambda_S_try = 1e-3;
        regfun_try = @(p) [ lambda_scale_try*(p(4:6)-1); lambda_S_try*p(7:9) ];
        try
            p_est_full = lsqnonlin(@(p) [gyro_resfun(p); regfun_try(p)], p0_full, lb_full, ub_full, opts_full);
        catch
            p_est_full = p0_full;
        end
    end
    p_est_gyro = p_est_full;
    b_est = p_est_gyro(1:3);
    s_est = p_est_gyro(4:6);
    s12 = p_est_gyro(7); s13 = p_est_gyro(8); s23 = p_est_gyro(9);
    S_est = [0 s12 s13; 0 0 s23; 0 0 0];
    K_est = (eye(3) + S_est) * diag(s_est);
    W_cal = (K_est * (W.' - b_est)).';

    % === Crear versión escalada de W_cal para todo el script ===
    W_cal_escalado = W_cal;
    for k = 1:numel(gyro_intervals)
        idxg = gyro_intervals{k};
        W_cal_escalado(idxg, :) = W_cal_escalado(idxg, :) / 100;
    end

    % ---------------------------
    % Construcción de T^g (matriz de misalignment) - varias versiones
    I3 = eye(3);

    % 1) Matriz (I+S) directamente (la usada en el modelo)
    T_IplusS = I3 + S_est;

    % 2) Matriz inversa exacta: T_inv = inv(I+S)
    T_inv = inv(I3 + S_est);

    % 3) Primera aproximación (orden 1) válida para misalignments pequeños:
    %    T_approx = I - S
    T_approx = I3 - S_est;

    % 4) Relación con K: (I+S) = K * diag(1./s_est)
    T_fromK = K_est * diag(1./s_est);

    % Mostrar resultados en forma legible
    fprintf('\n---- Matrices de misalignment (distintas convenciones) ----\n');
    fprintf('T_IplusS = I + S_est:\n'); disp(T_IplusS);
    fprintf('T_inv   = inv(I + S_est):\n'); disp(T_inv);
    fprintf('T_approx = I - S_est (1st-order approx):\n'); disp(T_approx);
    fprintf('T_fromK (reconstruida desde K_est y s_est):\n'); disp(T_fromK);

    % 5) Extraer gamma con la notación del paper (ejemplo de asignación)
    T = T_approx; % elegir la convención que pida el paper
    gamma_yz_A = -T(1,2);
    gamma_yz_A_alt =  T(1,3);

    % Para comodidad, presentar gammas basados en T_approx (ajusta si el paper pide otra)
    gamma_xy = - (T(3,1));
    gamma_xz =  T(2,1);
    gamma_yx =  T(3,2);
    gamma_yz = -T(1,2);
    gamma_zx =  T(1,3);
    gamma_zy = -T(2,3);

    % Mostrar gammas (revisar signos según la notación del paper)
    fprintf('\nGamma estimates (from T_approx, revisar signos según paper):\n');
    fprintf(' gamma_xy = %.6g\n gamma_xz = %.6g\n gamma_yx = %.6g\n gamma_yz = %.6g\n gamma_zx = %.6g\n gamma_zy = %.6g\n', ...
        gamma_xy, gamma_xz, gamma_yx, gamma_yz, gamma_zx, gamma_zy);
else
    % goto_only_bias true -> b_est, s_est, S_est, K_est, W_cal ya definidos
    p_est_gyro = [b_est; s_est; S_est(1,2); S_est(1,3); S_est(2,3)];
end

%% =========================
% 10) MÉTRICAS: RMS por intervalo (integración) y RMS global
%% =========================
gyro_rms_per_interval = nan(numIntervals-1,1);
for k=1:(numIntervals-1)
    idxg = gyro_intervals{k};
    if numel(idxg) < 2, continue; end
    if any(~isfinite(A_cal_mean(k,:))) || any(~isfinite(A_cal_mean(k+1,:))), continue; end
    gk = A_cal_mean(k,:).' / norm(A_cal_mean(k,:));
    gk_next = A_cal_mean(k+1,:).' / norm(A_cal_mean(k+1,:));
    t_local = t_gyro(idxg);
    g0 = gk;
    for i=1:(numel(idxg)-1)
        w_corr = (K_est * (W(idxg(i),:).' - b_est));
        dt = t_local(i+1) - t_local(i);
        g0 = rk4_gravity(g0, w_corr, dt);
    end
    g_pred = g0 / norm(g0);
    gyro_rms_per_interval(k) = norm(g_pred - gk_next);
end
gyro_RMS_global = sqrt(nanmean(gyro_rms_per_interval.^2));
fprintf('\nRMS residual por intervalo (gyro): %.6f\n', gyro_RMS_global);

%% =========================
% 11) Mostrar resultados y guardar
%% =========================
fprintf('\n--- RESULTADOS GIROSCOPIO ---\n');
fprintf('Bias (rad/s): bx=%.9f  by=%.9f  bz=%.9f\n', b_est(1), b_est(2), b_est(3));
fprintf('Escalas (diag): [%.6f  %.6f  %.6f]\n', s_est(1), s_est(2), s_est(3));
disp('Misalignment S (triangular superior):'); disp(S_est);
disp('K = (I + S) * diag(s):'); disp(K_est);

%% ======= Checks rápidos post-calibración =======
fprintf('\n=== VERIFICACIONES POST-CALIB ===\n');

% 1) ¿Alguna escala en bound?
lb_full = [-0.5;-0.5;-0.5;  0.7;0.7;0.7;  -0.15;-0.15;-0.15];
ub_full = [ 0.5; 0.5; 0.5;  1.3;1.3;1.3;   0.15; 0.15; 0.15];
p_est = [b_est(:); s_est(:); S_est(1,2); S_est(1,3); S_est(2,3)];

isAtLowerBound = (abs(p_est - lb_full) < 1e-8);
isAtUpperBound = (abs(p_est - ub_full) < 1e-8);

if any(isAtLowerBound) || any(isAtUpperBound)
    fprintf('Atención: algunos parámetros terminaron en límite (bounds).\n');
    idxLB = find(isAtLowerBound); 
    idxUB = find(isAtUpperBound);
    if ~isempty(idxLB), fprintf(' En lower bound indices: %s\n', mat2str(idxLB')); end
    if ~isempty(idxUB), fprintf(' En upper bound indices: %s\n', mat2str(idxUB')); end
end

% 2) Estadísticas de W_cal dentro de intervalos estáticos
W_cal = (K_est * (W.' - b_est)).';

nInt = numel(gyro_intervals);
std_por_eje = zeros(nInt,3);
mean_por_eje = zeros(nInt,3);

for k=1:nInt
    idxg = gyro_intervals{k};
    if numel(idxg)<2
        std_por_eje(k,:)=NaN; 
        mean_por_eje(k,:)=NaN; 
        continue; 
    end
    dat = W_cal(idxg,:);
    mean_por_eje(k,:) = mean(dat,1);
    std_por_eje(k,:)  = std(dat,0,1);
end

fprintf('W_cal: media por eje (mediana de intervalos) = [%g %g %g]\n', ...
    median(mean_por_eje,'omitnan'));
fprintf('W_cal: std por eje (mediana de intervalos)   = [%g %g %g]\n', ...
    median(std_por_eje,'omitnan'));

% 3) Validación hold-out: 20% intervalos
rng(0);
I = 1:nInt;
nHold = max(1,round(0.2*nInt));
holdIdx = randsample(I,nHold);
trainIdx = setdiff(I, holdIdx);

% recalcula bias solo con train
b_train = zeros(3,1); total=0;
for k=trainIdx
    idxg = gyro_intervals{k};
    b_train = b_train + sum(W(idxg,:)).';
    total = total + numel(idxg);
end
b_train = b_train / total;

% RMS en holdout
rms_hold = [];
for k=1:(nInt-1)
    if ismember(k, holdIdx)
        idxg = gyro_intervals{k};
        if numel(idxg)<2, continue; end

        gk = A_cal_mean(k,:).' / norm(A_cal_mean(k,:));
        gk_next = A_cal_mean(k+1,:).' / norm(A_cal_mean(k+1,:));

        g0 = gk; 
        t_local = t_gyro(idxg);

        for i=1:(numel(idxg)-1)
            w_corr = (K_est * (W(idxg(i),:).' - b_train));
            dt = t_local(i+1) - t_local(i);
            g0 = rk4_gravity(g0, w_corr, dt);
        end

        rms_hold(end+1) = norm(g0/norm(g0) - gk_next);
    end
end

fprintf('Hold-out RMS (media en intervalos hold): %g (n=%d)\n', ...
    mean(rms_hold,'omitnan'), numel(rms_hold));


% Guardar resultados en estructura
calibracion = struct();
calibracion.acc.M = M_est_acc; calibracion.acc.b = b_est_acc;
calibracion.gyro.bias = b_est; calibracion.gyro.scale = s_est; calibracion.gyro.S = S_est; calibracion.gyro.K = K_est;
calibracion.metrics.acc_RMS = sqrt(mean((magA_cal - g).^2));
calibracion.metrics.gyro_rms = gyro_RMS_global;
calibracion.intervals = intervals; calibracion.A_cal_mean = A_cal_mean;
calibracion.gyro_rms_per_interval = gyro_rms_per_interval;
save('calibracion_unificada.mat','calibracion','gyro_intervals','A_cal_all','magA_cal');

%% =========================
% 12) Generación tablas (opcional)
%% =========================
% Table II, V, VII, IX, XI simplificados
tableII.mean_bias = b_est; tableII.mean_scale = s_est; tableII.mean_S = S_est;
tableV = [ mean(abs(W(:,1)-b_est(1))), mean(abs(W(:,2)-b_est(2))), mean(abs(W(:,3)-b_est(3))) ];
tableVII = [gyro_RMS_global; max(gyro_rms_per_interval); nanmean(gyro_rms_per_interval)];
tableIX = K_est; tableXI = S_est;
save('tablas_calibracion.mat','tableII','tableV','tableVII','tableIX','tableXI');

fprintf('\nCalibración completa guardada en calibracion_unificada.mat y tablas_calibracion.mat\n');

% --- PROMEDIOS DE GIROSCOPIO EN INTERVALOS ESTÁTICOS ---
muestras_estaticas = [];   % aquí acumulamos todos los datos estáticos

for k = 1:numel(gyro_intervals)
    idx = gyro_intervals{k};
    if numel(idx) < 1
        continue;
    end
    % agregamos las muestras de este intervalo
    muestras_estaticas = [muestras_estaticas; W(idx,:)];      % crudo
end

% Promedios por eje del giroscopio crudo
mean_w_crudo = mean(muestras_estaticas, 1);

fprintf('\n=== PROMEDIOS GIROSCOPIO CRUDO EN INTERVALOS ESTÁTICOS ===\n');
fprintf('ωx = %.6f rad/s\n', mean_w_crudo(1));
fprintf('ωy = %.6f rad/s\n', mean_w_crudo(2));
fprintf('ωz = %.6f rad/s\n', mean_w_crudo(3));

%% === PROMEDIOS GIROSCOPIO CALIBRADO (ESCALADO /100) SOLO EN INTERVALOS ESTÁTICOS ===
muestras_cal_escaladas = [];

for k = 1:numel(gyro_intervals)
    idx = gyro_intervals{k};
    if numel(idx) < 1
        continue;
    end
    % Aquí tomas EXACTAMENTE lo que entra a la gráfica:
    muestras_cal_escaladas = [muestras_cal_escaladas; W_cal_escalado(idx, :)];
end

mean_W_cal_escalado = mean(muestras_cal_escaladas, 1);

fprintf('\n=== PROMEDIOS GIROSCOPIO CALIBRADO (ESCALADO /100) ===\n');
fprintf('ωx_scaled = %.9f rad/s\n', mean_W_cal_escalado(1));
fprintf('ωy_scaled = %.9f rad/s\n', mean_W_cal_escalado(2));
fprintf('ωz_scaled = %.9f rad/s\n', mean_W_cal_escalado(3));


% --- Ahora para el GIROSCOPIO CALIBRADO (W_cal) ---
muestras_cal = [];

for k = 1:numel(gyro_intervals)
    idx = gyro_intervals{k};
    if numel(idx) < 1
        continue;
    end
    muestras_cal = [muestras_cal; W_cal(idx,:)];
end

% Promedios por eje del giroscopio calibrado
mean_w_cal = mean(muestras_cal, 1);

fprintf('\n=== PROMEDIOS GIROSCOPIO CALIBRADO EN INTERVALOS ESTÁTICOS ===\n');
fprintf('ωx = %.6f rad/s\n', mean_w_cal(1));
fprintf('ωy = %.6f rad/s\n', mean_w_cal(2));
fprintf('ωz = %.6f rad/s\n', mean_w_cal(3));

%% =========================
% FUNCIONES (final del archivo)
%% =========================
function r = acc_residuals(p, A_mean_valid, g)
    M = reshape(p(1:9),3,3);
    b = p(10:12);
    Acal = (M * A_mean_valid.' ).' + repmat(b(:)', size(A_mean_valid,1), 1);
    norms = sqrt(sum(Acal.^2,2));
    r = norms - g;
end

function r = gyro_residuals_Tedaldi(p, W, t_gyro, gyro_intervals, A_cal_mean, intervals, t_acc, g)
    b   = p(1:3);
    s   = p(4:6);
    s12 = p(7); s13 = p(8); s23 = p(9);
    S = [0 s12 s13; 0 0 s23; 0 0 0];
    K = (eye(3) + S) * diag(s);
    Nint = size(intervals,1);
    r = [];
    for k=1:Nint
        idx_g = gyro_intervals{k};
        if numel(idx_g) < 2, continue; end
        if any(~isfinite(A_cal_mean(k,:))), continue; end
        gk = A_cal_mean(k,:).'; gk = gk / norm(gk);
        if k < Nint
            if any(~isfinite(A_cal_mean(k+1,:))), continue; end
            gk_next = A_cal_mean(k+1,:).'; gk_next = gk_next / norm(gk_next);
        else
            continue;
        end
        t_local = t_gyro(idx_g);
        g0 = gk;
        for i=1:(numel(idx_g)-1)
            w_raw = W(idx_g(i),:).';
            dt = t_local(i+1) - t_local(i);
            w_corr = K * (w_raw - b);
            g0 = rk4_gravity(g0, w_corr, dt);
        end
        g_pred = g0 / norm(g0);
        r_k = g_pred - gk_next;
        r = [r; r_k];
    end
end

function g_next = rk4_gravity(g, w, dt)
    k1 = -skew(w)*g;
    k2 = -skew(w)*(g + 0.5*dt*k1);
    k3 = -skew(w)*(g + 0.5*dt*k2);
    k4 = -skew(w)*(g + dt*k3);
    g_next = g + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
    g_next = g_next / norm(g_next);
end

function S = skew(w)
    S = [  0   -w(3)  w(2);
          w(3)   0   -w(1);
         -w(2)  w(1)   0   ];
end
%% =============================
% 13) GRÁFICAS DE VALIDACIÓN
%% =============================

figure('Name','Acelerómetro crudo','NumberTitle','off');
subplot(4,1,1);
plot(t_acc, A_raw(:,1)); ylabel('Ax (m/s^2)'); grid on;

subplot(4,1,2);
plot(t_acc, A_raw(:,2)); ylabel('Ay (m/s^2)'); grid on;

subplot(4,1,3);
plot(t_acc, A_raw(:,3)); ylabel('Az (m/s^2)'); grid on;

subplot(4,1,4);
plot(t_acc, sqrt(sum(A_raw.^2,2)));
ylabel('|a|'); xlabel('Tiempo (s)'); grid on;
sgtitle('Acelerómetro crudo (sin calibrar)');

%% === Comparación Giroscopio Crudo / Calibrado (Calibrado / 100) ===
figure('Name','Giroscopio Crudo vs Calibrado (escala reducida)','Color','w');
sgtitle('Comparación Giroscopio Crudo / Calibrado ','FontSize',18);

% Escala reducida para el calibrado
t = t_gyro;
W_cal_scaled = W_cal;   % copia para no modificar el original

% Recorrer intervalos estáticos y escalar solo esas muestras
for k = 1:numel(gyro_intervals)
    idxg = gyro_intervals{k};
    W_cal_scaled(idxg, :) = W_cal_scaled(idxg, :) / 1000;
end

%% ---- Eje wx ----
subplot(3,1,1)
plot(t, W(:,1), 'r', 'LineWidth', 1.0); hold on;
plot(t, W_cal_scaled(:,1), 'g', 'LineWidth', 1.0);
ylabel('\omega_x (rad/s)')
legend('Crudo','Calibrado')
grid on; grid minor;

%% ---- Eje wy ----
subplot(3,1,2)
plot(t, W(:,2), 'r', 'LineWidth', 1.0); hold on;
plot(t, W_cal_scaled(:,2), 'g', 'LineWidth', 1.0);
ylabel('\omega_y (rad/s)')
grid on; grid minor;

%% ---- Eje wz ----
subplot(3,1,3)
plot(t, W(:,3), 'r', 'LineWidth', 1.0); hold on;
plot(t, W_cal_scaled(:,3), 'g', 'LineWidth', 1.0);
ylabel('\omega_z (rad/s)')
xlabel('Tiempo (s)')
grid on; grid minor;

%% ======= Grafica 2: Giroscopio crudo =======
figure;
plot(t_gyro, W(:,1)); hold on;
plot(t_gyro, W(:,2));
plot(t_gyro, W(:,3));
xlabel('Tiempo (s)');
ylabel('Velocidad angular (rad/s)');
title('Giroscopio crudo - 3 ejes');
legend('Wx','Wy','Wz');
grid on;
