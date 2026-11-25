// SEIQR COVID-19 model: the five compartments (Case 1: beta = 4e-5; Case 2: beta = 1e-5)
// Run in Scilab: File > Open, then click Run.  (Not tested by the author of this file.)
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

t = 0:0.1:365;
y0 = [4990; 10; 0; 0; 0];
y = ode(y0, 0, t, seiqr);

names = ["S" "E" "I" "Q" "R"];
for i = 1:5
    scf(i);
    plot(t, y(i,:));
    xlabel("time (days)");
    ylabel(names(i));
end
[pk, ip] = max(y(3,:));
mprintf("Peak of I = %.1f on day %.1f\n", pk, t(ip));   // Case 1: about 698 on day 63.5
