%HF bulb, unlimited dividing potential, no signal

delete(gcp('nocreate'));
parpool

NSim = 20;

parfor nsim = 1:1:NSim

    HFBulb_WT_MultiCol(nsim)

end