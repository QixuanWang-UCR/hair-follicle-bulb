%HF bulb, unlimited dividing potential, FGF signal

delete(gcp('nocreate'));
parpool

NSim = 5;

parfor nsim = 1:1:NSim

    CL_ZeroAngle_FGF(nsim)

end