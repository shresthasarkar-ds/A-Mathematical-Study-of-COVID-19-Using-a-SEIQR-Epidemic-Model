// Euler's method against RK4 for the SEIQR model (Case 1)
// Part 1: largest error in I(t) for different step sizes.  Part 2: a very coarse step (h = 5 days).
// (Not tested by the author of this file.)
A = 100; mu = 0.02; beta = 4e-5;
alpha = 0.1; gamma = 0.05; sigma = 0.08; theta = 0.07;
k1 = alpha + mu + gamma; k2 = sigma + mu; k3 = theta + mu;

function dy = seiqr(t, y)
    S = y(1); E = y(2); I = y(3); Q = y(4); R = y(5);
    inf = beta * S * (E + I);
    dy = [A - mu*S - inf;
          inf - k1*E;
          alpha*E - k2*I;
          gamma*E + sigma*I - k3*Q;
          theta*Q - mu*R];
endfunction

y0 = [4990; 10; 0; 0; 0];

// ---- Part 1: error against step size (360 days) ----
T = 360;
hs = [5 2 1 0.5 0.2 0.1 0.05 0.01];
errE = []; errK = [];
for h = hs
    n = round(T/h);
    t = (0:n)*h;
    yref = [y0, ode(y0, 0, t(2:$), 1e-12, 1e-12, seiqr)];   // accurate reference solution
    ye = zeros(5, n+1); yk = zeros(5, n+1);
    ye(:,1) = y0; yk(:,1) = y0;
    for i = 1:n
        y = ye(:,i);
        ye(:,i+1) = y + h*seiqr(0, y);                       // Euler step
        y = yk(:,i);
        s1 = seiqr(0, y);
        s2 = seiqr(0, y + h/2*s1);
        s3 = seiqr(0, y + h/2*s2);
        s4 = seiqr(0, y + h*s3);
        yk(:,i+1) = y + h/6*(s1 + 2*s2 + 2*s3 + s4);         // RK4 step
    end
    errE = [errE, max(abs(ye(3,:) - yref(3,:)))];
    errK = [errK, max(abs(yk(3,:) - yref(3,:)))];
    mprintf("h = %g: Euler error %g, RK4 error %g\n", h, errE($), errK($));
end
scf(1); clf();
plot2d("ll", hs', [errE' errK'], [5 2]);
xlabel("step size h (days)"); ylabel("largest error in I(t) (people)");
legend("Euler", "RK4");
// expected at h = 0.1: Euler about 5.0, RK4 about 1e-7

// ---- Part 2: a very coarse step, h = 5 days ----
h = 5; T = 200; n = T/h; t = (0:n)*h;
ye = zeros(5, n+1); yk = zeros(5, n+1);
ye(:,1) = y0; yk(:,1) = y0;
for i = 1:n
    y = ye(:,i);
    ye(:,i+1) = y + h*seiqr(0, y);
    y = yk(:,i);
    s1 = seiqr(0, y);
    s2 = seiqr(0, y + h/2*s1);
    s3 = seiqr(0, y + h/2*s2);
    s4 = seiqr(0, y + h*s3);
    yk(:,i+1) = y + h/6*(s1 + 2*s2 + 2*s3 + s4);
end
tf = 0:0.5:T;
yf = [y0, ode(y0, 0, tf(2:$), 1e-12, 1e-12, seiqr)];
scf(2); clf();
plot(tf, yf(3,:), "k");
plot(t, ye(3,:), "r-o");
plot(t, yk(3,:), "b-s");
xlabel("time (days)"); ylabel("infected I(t)");
legend("accurate solution", "Euler, h = 5", "RK4, h = 5");
// expected peaks: Euler about 800, RK4 about 696, accurate solution about 698
