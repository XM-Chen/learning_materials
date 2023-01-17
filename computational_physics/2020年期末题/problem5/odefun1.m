function dydt = odefun1(t, v)
dydt = -(1 + t) * v + (4 * rand - 2);
end

