example_params;

if (length(A)-rank(obsv(A,C))) == 0  % sprawdzenie obserwowalności
    sprintf("System obserwowalny")
else
    sprintf("System nie obserwowalny")
    return
end

time_step = 0.01;
startpoint = [1 0.2 0.5];
ploton = true;
searchiter = 50;

[G1, G2, time_samples, betas, T, x, fval, exitflag, output] = getoptiG1G2(A, B, C, time_step, startpoint, ploton, searchiter);
plotG1G2(G1, G2, time_samples);

