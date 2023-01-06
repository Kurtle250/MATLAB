
%% Initial cost function and initial design points
% syms X1 X2
% syms f(X1,X2)
% f(X1,X2) = (6*(X1)^2)+(2*(X2)^2)-(6*(X1)*(X2))-(X1)-(2*(X2));
% x0 = [1,1,1,1];
% DFP_func(f,x0,0.001,20);

function Newton_DFP(cost_func,x0,eps,max_iter)
    arguments
       cost_func
       x0        double
       eps       double
       max_iter  double
    end
    syms X1 X2 X3 X4
    syms f(X1,X2,X3,X4)
    f(X1,X2,X3,X4) = cost_func;
    A = {eye(2,2),0,0,0,0,0};B = {0,0,0,0,0,0};C = {0,0,0,0,0,0};
    cval = {0,0,0,0,0,0};y = {0,0,0,0,0,0};s = {0,0,0,0,0,0};
    z = {0,0,0,0,0,0};d ={0,0,0,0,0,0};
    fval = {0,0,0,0,0,0};c_norm = {0,0,0,0,0,0};x ={0,0,0,0,0,0};
    c_norm{1} = 5;
    iterater = 1;
    x{1} = x0;
    fval{1} = f(x0(1),x0(2),x0(3),x0(4));
    fprintf("cost function: %s \n",char(cost_func));
    fprintf("Starting x0 design points:(%f,%f) \n",x{iterater}(1),x{iterater}(2));
    %fprintf("cost function at x0(%d,%d): %f \n",char(x{1}(1),x{1}(2),f{1});
    while(iterater < max_iter)
        if(c_norm{iterater}(1) < eps)
            fprintf("-----------------------------------------------------\n");
            fprintf("Iteration number: %d\n",iterater);
            fprintf("Norm of c is: %f \n",c_norm{iterater});
            fprintf("x is:(%f,%f) \n",x{iterater}(1),x{iterater}(2));
            fprintf("Cost function at x* is:%f \n",f(x{iterater}(1),x{iterater}(2)));
            break;
        end
        fprintf("-----------------------------------------------------\n");
        fprintf("Iteration number: %d\n",iterater);
        syms a
        syms g(a)
        dc = [diff(f, X1), diff(f, X2)];
        cval{iterater} = dc(x{iterater}(1),x{iterater}(2));
        c_norm{iterater} = norm(cval{iterater});
        fprintf("Norm of c is: %f \n",c_norm{iterater});
        fprintf("x is:(%f,%f) \n",x{iterater}(1),x{iterater}(2));
        fprintf("Cost function at x is:%f \n",f(x{iterater}(1),x{iterater}(2)));
        if(iterater == 1)
            d{iterater}=(-1).*cval{iterater};
            x{iterater+1}=[(x{iterater}(1)+d{iterater}(1)*a),(x{iterater}(2)+d{iterater}(2)*a)];
            g(a) = subs(f,[X1,X2],[x{iterater+1}(1),x{iterater+1}(2)]);
            dg_da = diff(g,a);
            alpha = solve(dg_da,a);
            x{iterater+1}=[(x{iterater}(1)+d{iterater}(1)*alpha),(x{iterater}(2)+d{iterater}(2)*alpha)];
            s{iterater} = alpha.*d{iterater};
            cval{iterater+1} = dc(x{iterater+1}(1),x{iterater+1}(2));
            y{iterater} = cval{iterater+1} - cval{iterater};
            z{iterater} = y{iterater};
            B{iterater} = (s{iterater}.*s{iterater}.')/dot(s{iterater},y{iterater});
            C{iterater} = (z{iterater}.*z{iterater}.')/dot(z{iterater},y{iterater});
            A{iterater+1} = A{iterater} + B{iterater} +C{iterater};
            c_norm{iterater+1} = norm(cval{iterater+1});
        else
            d{iterater}=(-1).*cval{iterater}*A{2};
            x{iterater+1}=[(x{iterater}(1)+d{iterater}(1)*a),(x{iterater}(2)+d{iterater}(2)*a)];
            g(a) = subs(f,[X1,X2],[x{iterater+1}(1),x{iterater+1}(2)]);
            dg_da = diff(g,a);
            alpha = solve(dg_da,a);
            x{iterater+1}=[(x{iterater}(1)+d{iterater}(1)*alpha),(x{iterater}(2)+d{iterater}(2)*alpha)];
            s{iterater} = alpha*d{iterater};
            cval{iterater+1} = dc(x{iterater+1}(1),x{iterater+1}(2));
            y{iterater} = cval{iterater+1} - cval{iterater};
            z{iterater} = y{iterater};
            B{iterater} = (s{iterater}.*s{iterater}.')/dot(s{iterater},y{iterater});
            C{iterater} = (z{iterater}.*z{iterater}.')/dot(z{iterater},y{iterater});
            A{iterater+1} = A{iterater} + B{iterater} +C{iterater};
            c_norm{iterater+1} = norm(cval{iterater+1});
        end
        iterater=iterater+1;
    end 
end