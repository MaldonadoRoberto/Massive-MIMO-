%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                Funcion: getConditionNumber                        %%%% 
%%%%                                                                   %%%%
%%%%       Función en la que se calcula el número de condición para    %%%%
%%%%         una matriz de canal dada. Se descompone la matriz en      %%%%
%%%%        valores singulares, y se obtiene la relación entre el      %%%%
%%%%        máximo y el mínimo. Se realiza tanto para la matriz de     %%%%
%%%%              Quadriga como para la matriz i.i.d                   %%%%
%%%%                                                                   %%%%
%%%%      PARAMETROS DE ENTRADA:                                       %%%%
%%%%                                                                   %%%%
%%%%      H_struct: struct que contiene el canal de quadriga           %%%%
%%%%       normalizado                                                 %%%%
%%%%                                                                   %%%%
%%%%      H_struct_iid: struct que contiene el canal iid               %%%%
%%%%       normalizado                                                 %%%%
%%%%                                                                   %%%%                                                                                                                               %%%%
%%%%      PARAMETROS DE SALIDA:                                        %%%%
%%%%                                                                   %%%%
%%%%      condition_number: número de condicion para matriz qudriga    %%%%
%%%%                                                                   %%%%
%%%%      condition_number_iid: numero de condicion para iid           %%%%
%%%%                                                                   %%%%
%%%%            Resultados mostrados en sección: 7.4.2                 %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [condition_number,condition_number_iid] = get_condition_number(H_struct,H_struct_iid)


ntx = size(H_struct{1},1);
nrx = size(H_struct{1},2);

for index = 1:size(H_struct,2)
    % selecion de los datos
    H = H_struct{index};
    H_iid = H_struct_iid{index};

    %% Una vez se tienen normalizada la matriz del canal se procede
    % a obtener el spread de los valores signulares, para lo que primeramente
    % se debe de descomponer con SVD la matriz H:
    [U,D,V] = svd(H);
    D_diagonal = diag(D);
    
    % Se eliminan los posibles valores nulos
    D_diagonal(D_diagonal==0) = [];
    SV_spread(index) = 10*log10(D_diagonal(1,1)/min(D_diagonal));
    
    % SV spread idd  
    [U_iid,D_iid,V_iid] = svd(H_iid);
    D_iid_diagonal = diag(D_iid);
    
    max_iid_singularvalue = max(D_iid_diagonal);
    min_iid_singularvalue = min(D_iid_diagonal);
    
    SV_spreadiid(index) = 10*log10(max_iid_singularvalue/min_iid_singularvalue);
    
end

condition_number=mean(SV_spread);
condition_number_iid = mean(SV_spreadiid);
