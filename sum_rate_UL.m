%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                 Función: sum_rate_UL                              %%%%
%%%%                                                                   %%%%
%%%%      Función en las que se aplica las ecuaciones de la SINR       %%%%
%%%%      del capitulo de MIMO masivo para el caso del enlace UL.      %%%%
%%%%      Se realiza la simulación para un vector de SNRs              %%%%
%%%%                                                                   %%%%
%%%%    Este script utiliza la matriz de canal H (small scale fading)  %%%%
%%%%                                                                   %%%%
%%%%      PARAMETROS DE ENTRADA:                                       %%%%
%%%%                                                                   %%%%
%%%%      H: matriz de canal que solo contiene small scale fading      %%%%
%%%%                                                                   %%%%
%%%%      snrdB: vector de SNR en dB                                   %%%%
%%%%                                                                   %%%%
%%%%                                                                   %%%%
%%%%      PARAMETROS DE SALIDA:                                        %%%%
%%%%                                                                   %%%%
%%%%      R_total_mrc: vector que contiene el sum rate para cada SNR   %%%%
%%%%      para el detector MRC                                         %%%%
%%%%                                                                   %%%%
%%%%      R_total_zf: vector que contiene el sum rate para cada SNR    %%%%
%%%%      para el detector ZF                                          %%%%
%%%%                                                                   %%%%
%%%%      R_total_mmse: vector que contiene el sum rate para cada SNR  %%%%
%%%%      para el detector MMSE                                        %%%%
%%%%                                                                   %%%%
%%%%      sum_iterativo: vector que contiene el sum rate optimo        %%%%
%%%%      para cada SNR                                                %%%%
%%%%                                                                   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [R_total_mrc, R_total_zf, R_total_mmse, sumrate] = sum_rate_UL(H,SNRdB)

nrx = size(H,2);
ntx = size(H,1);

HH = conj(H);
precoders = {'mrc' 'zf' 'mmse'};
SNRlin = 10.^(SNRdB./10);


for r = 1:length(SNRdB)
    SNR = SNRlin(r);
    
    for x=1:3
        precoder = char(precoders(x));
        switch (precoder)
            case 'mrc'
                W = H;
                
            case 'zf'
                W = H*(H'*H)^-1;
                
            case'mmse'
                for g = 1:nrx
                    W(:,g) = sqrt(SNR)*(SNR*H*H'+eye(ntx))^-1*H(:,g);
                end                
        end
        
        for t = 1 : nrx
            H_diff = H;
            H_diff(:,t) = zeros(ntx,1);% matriz de canal en la
            % que se no tiene en cuenta el vector de coeficientes del usuario k
            
            v = 0;
            for m = 1:nrx % recorremos las columnas de la matriz
                v = norm(W(:,t)'*H_diff(:,m))^2 + v;
            end
            
            SINR = (SNR*(norm(W(:,t)'*H(:,t))^2))/(SNR*v+norm(W(:,t))^2);
            R(t) = log2(1+SINR);
        end
        
        eval(sprintf('R_total_%s(r) = sum(R);',precoder));  % Sumatoria del rate de todos los usuarios para cada iteración
        clear W
        
    end
    
    %%
    % Para cada iteración se va a obtener tambien el valor de la capacidad
    % usando la ecuación de la rewguión de capacidad en el caso de que
    % estemos ante un sistema en el que los usuarios tengan una sola antena
    % Segun se dice en los papers, en este caso, la capacidad obtenida con
    % mmse y la obtenida con la formula de la capacidad deben ser iguales
    sumrate(r) = real(log2(det(eye(ntx)+SNR*H*H')));
end

end
