%HF bulb, unlimited dividing potential, FGF signal

delete(gcp('nocreate'));
parpool

NBool_Vec= [50 100 200];
% this gives different gene regulation info. updating frequency
s_str_C_Vec = [0.2 0.5 1]; 

% n_NBool = 1;
% s_str = 1;
% s_str_C = 0.015;

parfor nsim = 1:9

    n_NBool = ceil(nsim/3);
    n_s_str_C = mod(nsim, 3);
    if n_s_str_C == 0
        n_s_str_C = 3;
    end

    s_str = s_str_C_Vec(n_s_str_C);
    NBool = NBool_Vec(n_NBool);
% s_str = s_str_C_Vec;
% NBool = NBool_Vec;


   TwoGenes_ZeroAngle_FGFHomoBoo(s_str, NBool)

end