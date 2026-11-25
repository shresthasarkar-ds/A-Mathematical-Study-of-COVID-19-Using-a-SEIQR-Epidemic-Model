// SEIQR phase plane: S against I from five different starting points
// casenr = 1: R0 = 2.35 (curves spiral into the endemic point); casenr = 2: R0 = 0.59 (curves go to the disease-free point)
// (Not tested by the author of this file.)
casenr = 1;
A = 100; mu = 0.02;
alpha = 0.1; gamma = 0.05; sigma = 0.08; theta = 0.07;
k1 = alpha + mu + gamma; k2 = sigma + mu; k3 = theta + mu;
if casenr == 1 then
    beta = 4e-5; T = 500; Sstar = 2125; Istar = 338.2;      // endemic point C*
else
    beta = 1e-5; T = 150; Sstar = 5000; Istar = 0;          // disease-free point C0
end

function dy = seiqr(t, y)
    S = y(1); E = y(2); I = y(3); Q = y(4); R = y(5);
    inf = beta * S * (E + I);
    dy = [A - mu*S - inf;
          inf - k1*E;
          alpha*E - k2*I;
          gamma*E + sigma*I - k3*Q;
          theta*Q - mu*R];
endfunction

t = 0:0.1:T;
starts = [4990 4900 3000 1000 4000;      // S(0)
            10  100  300   50   20;      // E(0)
             0    0    0    0    0;
             0    0    0    0    0;
             0    0    0    0    0];
clf();
for j = 1:5
    y = ode(starts(:,j), 0, t, seiqr);
    plot(y(1,:), y(3,:));
    plot(y(1,1), y(3,1), "ko");              // black dot = starting point
end
plot(Sstar, Istar, "r*");                    // red star = equilibrium point
xlabel("S (susceptible)");
ylabel("I (infected)");
