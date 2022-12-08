example_params;

options = optimset('MaxIter',500,'PlotFcns',@optimplotfval);

params0 = [1, 0.5, 0.5];  % T0, beta1, beta2 - długość zależy od wymiaru macierzy A
optifhand = @(p)optiG1G2(A, B, C, p);

[x,fval,exitflag,output] = fminsearch(optifhand, params0, options);
