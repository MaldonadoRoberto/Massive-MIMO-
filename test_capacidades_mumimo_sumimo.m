%%
%
% Script para comprobar que el sum-rate obtenido mediante (1+log2(SINR) es equivalente
% a obtenerlo mediante la formula de un sistema MIMO single user.

clear all
clc


%
% Voy a suponer un sistema mimo multi usuario de 3 antenas y 2 usuarios, cuyo
% equivalente en un sistema mimo single user consiste en 3 antenas transmisoras (ntx)
% y 3 antenas receptoras (nrx)
M = 80;
K = 10;

ntx = M;
nrx = K;

H = (randn(M,K)+1i*randn(M,K))/sqrt(2);  %(orden M x K)

% Como en este caso no se tienen en cuenta el latge scale fading pd es la
% snr recibida media. Supondremos que tiene un valor de 10db;

SNR_media_dB = 10;
SNR_media = 10^(10/10);

%% Primer metodo: Utilizaré en este primer método el proceso que se usa en
% la tesis por el cual se obtiene el sum-rate mediante la formula
% log2(1+SINR). Utilizaré la precodificación MRT:

W = conj(H);  % Matriz de precodifcaicón (orden M x K)

alpha =  1/trace(W*W'); % Factor de normalización

P_transmitida = 40; %40 W de potrencia transmitida;
P_ruido_dBm = -100;
P_ruido = 10^(-100/10)/1000; %dBm -> W

% En este bucle se recorren todos los usuarios del sistema para calcular la
% SINR de cada uno de ellos:
for t = 1 : K
    
    % La parte del denominador de la ecuación de SINR tiene una sumatoria
    % que incluye el vector de canal del usuario k y las columnas de la
    % matriz de precodificación de los demás usuarios diferentes de k.
    
    % Para eliminar la contribución del usuario k, creamos una matriz de
    % precodificación W_diff igual que W pero con una columna de 0 que se
    % corresponde a la precodificaicón del usuario k
    W_diff = W;
    W_diff(:,t) = zeros(1,M); % matriz de canal en la que se no tiene en cuenta el vector de coeficientes del usuario k
    
    v = 0;
    for m = 1:K % recorremos las filas de la matriz
        v(m) = norm(H(:,t).'*W_diff(:,m))^2;
    end
    
    sum_v = sum(v);  %Calculamos la sumatoria
    
    SINR = ((alpha)*(SNR_media)*norm(H(:,t).'*W(:,t))^2)/(sum_v*alpha*SNR_media+1);  % Ecuación (2.28)
    
    R(t) = log2(1+SINR);  % Rate para cada usuario
    
end

C_SINR = sum(R);  % Sum-rate total


%% Segundo metodo: En este caso se utilizará la ecuación de mimo single user
% para obtener la ¿¿capacidad?? con la ecuación (2) del tutorial de MIMO
% masivo que nos pasó Javier:

% C = log2(det(I + (1/No)*H*Kx*H^H)


% De este mismo documento, se conoce la matriz de covarianzas Kx está
% formada por las siguientes matrices:

% Kx = Q*P*Q^H (siendo Q la matriz de precodificación y P la matriz de asignación de potencias para cada antena transmisora)


% En el caso de MIMO punto a punto, la matriz de canal tiene dimensiones
% nrx X ntx por lo que primero se debe de trasponer la matriz H:

H = H.';  %orden (nrx X ntx)

% Ahora bien conociendo las dimensiones de H y por tando las de H^H (ntx X
% nrx). Quedan definidas las dimensiones de la matriz Kx que deben ser
% [ntx X ntx]

% Matriz de asignación de potencias (la asignación se hace con respecto a
% las antenas transmisoras):

P = eye(nrx).*P_transmitida/nrx;  % definida en la diapotiva 5 del massive mimo tutorial

% Sabiendo que Kx debe ser de dimension [ntx X ntx] y que P tiene dimension
% [ntx x ntx] la matriz de precodificación Q deberia ser de dimensiones
% [ntx X ntx]


% En el caso de utiliza como matriz de precodificación la mtriz V
% procedente de la descomposición de valores singulares de la matriz H, las
% dfimensiones coinciden:

[U,D,V] = svd(H);
fprintf('\nLas dimensiones de V son [%d,%d]',size(V));

% Si utilizamos como matriz de precodificación la matriz MRT (conj(H))
% tenemos que las dimensiones no concuerdan ya que aunque se pueda multiplicar
% el resultado de la matriz Kx no es de orden ntx X ntx sino de orden nrx X nrx
% y por tanto no puede ser multiplicado por H:
% W = conj(H);
% fprintf('\nLas dimensiones de W son [%d,%d]',size(W));
%
% Kx = W*P*W';
% fprintf('\nLas dimensiones de Kx son [%d,%d]',size(Kx));


% Si usamos W = H^H; la matriz Kx no puede ser obtenida. Pruebo ambas
% variatnes para MRT porque en un sistema MU-MIMO con H de orden M x K la
% matriz de precodificaicón con MRT es conj(H) pero como hemos traspuesto
% la matriz de canal y ahora es de orden K x M = nrx X ntx; he probado
% tambien a usar la matriz de precodficaición H^H.

W = H';
fprintf('\nLas dimensiones de W son [%d,%d]',size(W));
Kx = W*P*W';
% Aun asi las matrices tampoco concuerdan.


C_SUMIMO_con_matrices = log2(det(eye(nrx)+(1/P_ruido)*H*Kx*H'));

% A no ser que se me haya pasado algo, no se puede utilizar precodificacion
% lineal en MIMO sngle user (solamente SVD o matriz identidad)

%% Tercer metodo:

% Voy a utilizar la misma ecuación que el metodo anterior, pero con la
% diferencia de que ahora voy a calcular la matriz de covarianza
% directamente del vector de simbolos transmitidos generados según la
% ecuación (2.24) de la tesis:

% x = sqrt(alpha) * W * q; siendo q un vector de simbolos para todos los
% usuarios y de potencia unidad;

% Para realizar el calculo, lo primero que hago es generar simbolos de
% información (QAM) para cada usuario, calcular x y obtener su correlación
% ya que, al tener media nula, cov = corr

for n = 1:100
    q = (randn(1,K)+1i*randn(1,K))/sqrt(2);
    W = H';  % pongo la matriz hermitica, porque si es conjugada no hay forma de multiplicar dicha matriz
    % por el vector de información q (orden (1 X nrx) o (nrx x 1)
    alpha = 1/trace(W*W');
    x(:,n) = sqrt(alpha)*W*q.';
    correlation_iter(:,:,n) = x(:,n)*x(:,n)';
end

m_correlation = mean(correlation_iter,3);

Kx = m_correlation; % En este caso trace(Kx) = 1;

C_SUMIMO = log2(det(eye(K)+SNR_media*H*Kx*H'));


%% Cuarto metodo:

% En este metodo voy a utilizar la formula de la cpacidad para MU-MIMO
% (maximización del a información mutua) por lo que se necesita realizar
% una optimización de potencias asignadas a los usuarios

H = H.'; % Tomamos H como la H de un sistema MU MIMO (tamaño M x K)
Pmax = 1;
cvx_begin
cvx_quiet(true); % this suppresses screen output from the solver
variable p(K);
maximize (-det_inv(eye(M)+ SNR_media*conj(H)*diag(p)*H.'));
subject to
sum(p)<=Pmax
min(p)>=0
cvx_end
C_MUMIMO = log2(det(eye(M)+SNR_media*conj(H)*diag(p)*H.'));  %Formula (2.6) de la tesis

potencia_asignada = wf(H,SNR_media);
C_MUMIMO_WF = log2(det(eye(M)+SNR_media*H*diag(potencia_asignada)*H'));

% Si realizamos el mismo método pero sin optmización de potencias tenemos:
P_MUMIMO = eye(K)*1/K;
C_MUMIMO_equipotencia = log2(det(eye(M)+SNR_media*conj(H)*P_MUMIMO*H.'));

% en el caso de iid tiene sentido ya que los usuarios estan situados al a
% misma distancia y euipoitencia y optimización tienen resultados muy
% similares.

fprintf('\n\n\n\n--------RESUMEN DE LA SIMULACIÓN--------')
fprintf('\n\n Sum rate obtenido con (1+SINR) = %2f (bits/s/Hz)', C_SINR);
fprintf('\n Capacidad obtenida con single user MIMO con Kx = QPQ^H = %2f (bits/s/Hz)', C_SUMIMO_con_matrices);
fprintf('\n Capacidad obtenida con single user MIMO = %2f (bits/s/Hz)', C_SUMIMO);
fprintf('\n Capacidad obtenida con MultiUSER MIMO (optimación de potencias con CVX) = %2f (bits/s/Hz)', C_MUMIMO);
fprintf('\n Capacidad obtenida con MultiUSER MIMO (optimación de potencias con WF) = %2f (bits/s/Hz)', C_MUMIMO_WF);
fprintf('\n Capacidad obtenida con MultiUSER MIMO (equipotencias) = %2f (bits/s/Hz)', C_MUMIMO_equipotencia);