function []=LaboratoryWork5()
clear; clc;
t = 0:pi/50:2*pi;
q = exp(1i*t);
lambda = -14;

figure(1);
xlim([-1 0])
grid on
explicitB=[3/2, -1/2];
explicitAdams(explicitB,q,'green');
explicitB=[23/12, -16/12, 5/12];
explicitAdams(explicitB,q,'magenta');
explicitB=[55/24, -59/24, 37/24, -9/24];
explicitAdams(explicitB,q,'cyan');

function []= explicitAdams(explicitB,q,color)
    res1 = 0;
    for m = 1:length(explicitB)
        res1 = res1 + explicitB(m)*q.^(length(explicitB) - m);
    end
    res2 = (q.^length(explicitB) - q.^(length(explicitB) - 1)) ./ res1;
    patch(real(res2),imag(res2),color);
end
title("Explicit Adams method 2-4 order")
figure(2);
xlim([-7 0]);
grid on
implicitB=[5/12, 8/12, -1/12];
implicitAdams(implicitB,q,'green');
implicitB=[9/24, 19/24, -5/24, 1/24];
implicitAdams(implicitB,q,'magenta');
implicitB=[251/720, 646/720, -264/720, 106/720, -19/720];
implicitAdams(implicitB,q,'cyan');

function []= implicitAdams(implicitB,q,color)
    res1 = 0;
    for m = 1:length(implicitB)
        res1 = res1 + implicitB(m)*q.^(length(implicitB) - m);
    end
    res2 = (q.^(length(implicitB) - 1) - q.^(length(implicitB) - 2)) ./ res1;
    patch(real(res2),imag(res2),color);
end
title("Implicit Adams method 2-4 order")
disp('Critical tau step:')
disp(0.3/abs(lambda))
end