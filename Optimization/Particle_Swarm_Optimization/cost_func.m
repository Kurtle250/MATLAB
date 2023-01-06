function f = cost_func(x)
    h = x(1);
    L = x(2);
    t = x(3);
    b = x(4);
    f = 1.1047*L*h^2+0.04811*t*b*(14.0+L);
end