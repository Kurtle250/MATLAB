function T = Test_func(x0)
    x = x0(1);
    y = x0(2);
    T = (x+2*y-7)^2+(2*x+y-5)^2;
end