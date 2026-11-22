// SEIQR lockdown scenarios (Case 1 parameters): infected population I(t)
// beta is cut by the fraction c between day ts and day te.  (Not tested by the author of this file.)
global bb
A = 100; mu = 0.02; beta = 4e-5;
alpha = 0.1; gamma = 0.05; sigma = 0.08; theta = 0.07;
k1 = alpha + mu + gamma; k2 = sigma + mu; k3 = theta + mu;

function dy = seiqrL(t, y)
    global bb
    S = y(1); E = y(2); I = y(3); Q = y(4); R = y(5);
    inf = bb * S * (E + I);
    dy = [A - mu*S - inf;
          inf - k1*E;
          alpha*E - k2*I;
          gamma*E + sigma*I - k3*Q;
          theta*Q - mu*R];
endfunction

dt = 0.1;
y0 = [4990; 10; 0; 0; 0];
//          no lockdown | day 30 kept | 30-120 lifted | day 60 kept | weaker (40%)
c_  = [0      0.6    0.6    0.6    0.4];
ts_ = [30     30     30     60     30];
te_ = [120    365    120    365    365];
cols = ["k" "c" "b" "m" "r"];

clf();
for s = 1:5
    c = c_(s); ts = ts_(s); te = te_(s);
    bb = beta;
    tA = dt:dt:ts;
    yA = ode(y0, 0, tA, seiqrL);
    tt = [0 tA]; yy = [y0 yA];
    bb = beta*(1 - c);
    if te > ts then
        tB = ts+dt:dt:te;
        yB = ode(yy(:,$), ts, tB, seiqrL);
        tt = [tt tB]; yy = [yy yB];
    end
    bb = beta;
    if te < 365 then
        tC = te+dt:dt:365;
        yC = ode(yy(:,$), te, tC, seiqrL);
        tt = [tt tC]; yy = [yy yC];
    end
    plot(tt, yy(3,:), cols(s));
    [pk, ip] = max(yy(3,:));
    mprintf("scenario %d: peak I = %.1f on day %.1f\n", s, pk, tt(ip));
end
xlabel("time (days)"); ylabel("infected I(t)");
legend("No lockdown", "Lockdown from day 30, kept", "Lockdown days 30-120, then lifted", "Late lockdown from day 60, kept", "Weaker lockdown from day 30 (40%)");
// expected peaks: 698 (day 63.5), 146 (day 40), 636 (day 170), 690 (day 61), 269 (day 76)
