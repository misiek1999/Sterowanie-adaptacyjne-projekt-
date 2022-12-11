example_params;

if (length(A)-rank(obsv(A,C))) == 0  % sprawdzenie obserwowalności
    sprintf("System obserwowalny")
else
    sprintf("System nie obserwowalny")
    return
end

options = optimset('MaxIter',500,'PlotFcns',@optimplotfval);

params0 = [1, 0.2, 0.2];  % T0, beta1, beta2, ..., beta_n_B - długość zależy od wymiaru macierzy B
optifhand = @(p)optiG1G2(A, B, C, p);

[x,fval,exitflag,output] = fminsearch(optifhand, params0, options);
