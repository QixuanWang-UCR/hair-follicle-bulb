%HF bulb, unlimited dividing potential, no signal

delete(gcp('nocreate'));
parpool

NSim = 30;

rand(1,2024);

parfor nsim = 1:1:NSim

    ZeroAngle_WT(nsim)

end
