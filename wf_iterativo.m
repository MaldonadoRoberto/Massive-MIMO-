%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                  Funcion: wf_iterativo                            %%%%
%%%%                                                                   %%%%
%%%%       Funcion en la que se obtienen las potencias optimas en el   %%%%
%%%%       canal MIMO MAC dual con las que se maximizan la capacidad   %%%%
%%%%        de BC.                                                     %%%%
%%%%                                                                   %%%%
%%%%     El metodo utilizado esta descrito en el apartado 3.4.3        %%%%
%%%%                                                                   %%%%
%%%%      PARAMETROS DE ENTRADA:                                       %%%%
%%%%                                                                   %%%%
%%%%      H: matriz de canal                                           %%%%
%%%%                                                                   %%%%
%%%%      SNRlin : SNR en escala lineal                                %%%%
%%%%                                                                   %%%%
%%%%      PARAMETROS DE SALIDA:                                        %%%%
%%%%                                                                   %%%%
%%%%      sumsum: sum rate obtenido tras aplicar el algoritmo          %%%%
%%%%                                                                   %%%%
%%%%      Q_total: vector de potencias asignadas a cada usuario        %%%%
%%%%                                                                   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [sumsum,Q_total] = wf_iterativo(H,SNRlin)

M = size(H,2);
K = size(H,1);
Q_user = 0;

Q_total = repmat(Q_user,1,K); % iniciamos matrices de cov

H_eff = [];
sumsum = 0;
delta = 100;
while delta > 0.01
    
    for i=1:K % calculo de los canales efectivos para cada user;
        H_sin_user = H;
        H_sin_user(i,:) = zeros(1,M);
        summ = 0;
        for j =1:K
            summ = H_sin_user(j,:)'*Q_total(j)*H_sin_user(j,:) + summ;  % hermiticos al reves
        end
        
        H_eff(i,:) = H(i,:)*((eye(M) + summ)^(-1/2));
        
    end
    
    % calculo de los valores singulares de cada vector efectivo;
    
    for i=1:K
        [u,d,v] = svd(H_eff(i,:));
        valores_singularess(i) = d(1);
    end
    
    for i = 1:K
        eval(sprintf('efect_%d(1,:) = H_eff(i,:);',i));  % Sumatoria del rate de todos los usuarios para cada iteración
    end
    
    % matriz diagonal para 8 usuarios, en el caso de que se simule con otro
    % número de usuarios hay que modificar esta linea
    block = blkdiag(efect_1,efect_2,efect_3,efect_4,efect_5,efect_6,efect_7,efect_8);   
    
    [pot] = waterFiling(block,SNRlin);
    potencia = zeros(1,K);
    Q_nueva = zeros(1,K);
    
    
    % Este bucle se hace para relacionar la potencia obtenida al aplicar el
    % algoritmo de wf simple a los usuarios del sistema.
    sv_sort = sort(valores_singularess,'descend');
    
    for i=1:K
        [a,b] = find(valores_singularess == sv_sort(i));
        potencias(b) = pot(i);
        Q_nueva(b) = (1/K)*potencias(b) + ((K-1)/K)*Q_total(b);
    end
    A = diag(Q_nueva);
    
    sum_rate = log2(det(eye(M) + SNRlin*H'*A*H));
    delta = abs(sumsum - sum_rate);
    
    sumsum = sum_rate;
    Q_total = Q_nueva;
end


