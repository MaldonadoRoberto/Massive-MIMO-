%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                 Función: sum_rate_DL                              %%%%
%%%%                                                                   %%%%
%%%%      Función en las que se aplica las ecuaciones de la SINR       %%%%
%%%%      del capitulo de MIMO masivo para el caso del enlace DL.      %%%%
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
%%%%      R_total_mrt: vector que contiene el sum rate para cada SNR   %%%%
%%%%      para el precodificador MRT                                   %%%%
%%%%                                                                   %%%%
%%%%      R_total_zf: vector que contiene el sum rate para cada SNR    %%%%
%%%%      para el precodificador ZF                                    %%%%
%%%%                                                                   %%%%
%%%%      R_total_mmse: vector que contiene el sum rate para cada SNR  %%%%
%%%%      para el precodificador MMSE                                  %%%%
%%%%                                                                   %%%%
%%%%      sum_iterativo: vector que contiene el sum rate optimo        %%%%
%%%%      para cada SNR, en este caso, con water filing iterativo      %%%%
%%%%                                                                   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [R_total_mrt, R_total_zf, R_total_mmse,sum_iterativo] =   sumrate_DL(H,SNRdB)

nrx = size(H,2);
ntx = size(H,1);

precoders = {'mrt' 'zf' 'mmse'};
SNRlin = 10.^(SNRdB./10);

for r = 1:length(SNRdB)
    SNR = SNRlin(r);
    
    for x=1:3
        precoder = char(precoders(x));
        
        switch (precoder)
            case 'mrt'
                W = conj(H);
                
            case 'zf'
                W = conj(H)*(H.'*conj(H))^(-1);
                
            case'mmse'
                W = conj(H)*(H.'*conj(H)+(nrx/SNR).*eye(nrx))^(-1);
        end
        % Factor de normalización
        alpha =  1/trace(W*W');
        for t = 1 : nrx
            W_diff = W;
            W_diff(:,t) = zeros(1,ntx);% matriz de canal en la
            % que se no tiene en cuenta el vector de coeficientes del usuario k
            
            v = 0;
            for m = 1:nrx % recorremos las filas de la matriz
                v(m) = norm(H(:,t).'*W_diff(:,m))^2;
            end
            sum_v = sum(v);
            
            SINR = ((alpha)*(SNR)*norm(H(:,t).'*W(:,t))^2)/(sum_v*(SNR)*alpha+1);
            R(t) = real(log2(1+SINR));
            
        end
        
        eval(sprintf('R_total_%s(r) = sum(R);',precoder));  % Sumatoria del rate de todos los usuarios para cada iteración
        
    end
    
     [sum_iterativo(r), potencia] = wf_iterativo(H.',SNR);
end