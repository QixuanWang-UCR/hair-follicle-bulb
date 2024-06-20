%HF bulb, unlimited dividing potential, FGF signal

% delete(gcp('nocreate'));
% parpool

NBool_Vec= 150;
% this gives different gene regulation info. updating frequency
s_str_C_Vec = 0.8; 

% n_NBool = 1;
% s_str = 1;
% s_str_C = 0.015;

for nsim = 1:1

%     n_NBool = ceil(nsim/2);
%     n_s_str_C = mod(nsim, 3);
%     if n_s_str_C == 0
%         n_s_str_C = 3;
%     end
% 
%     s_str = s_str_C_Vec(n_s_str_C);
%     NBool = NBool_Vec(n_NBool);
s_str = s_str_C_Vec;
NBool = NBool_Vec;


   TwoGenesCL_Boo_WTBasDiv(s_str, NBool)

end