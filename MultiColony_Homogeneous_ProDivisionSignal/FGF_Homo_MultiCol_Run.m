%HF bulb, unlimited dividing potential, FGF signal

delete(gcp('nocreate'));
parpool

NSim = 20;
%rand(3,3);
parfor nsim = 1:NSim

    FGF_Homo_MultiCol(nsim)

end