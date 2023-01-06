function H = penalty_func(x)
    %% input varibles for X
    h = x(1);
    L = x(2);
    t = x(3);
    b = x(4);
    %% constant values
    P = 6000;
    W=14;
    E=30*exp(1)^6;
    G=12*exp(1)^6;
    Tau_max=13600;
    Sig_max=30000;
    Delta_max =0.25;
    %% constrain sub functions
    M = P*(W+(L/2));
    J = 2*(sqrt(2)*h*L*(L^2/12 +((h+t)/2)^2));
    R = sqrt((L^2/4)+((h+t)/2)^2);
    Tau_pp = (M*R)/J;
    Tau_p = P/(sqrt(2)*h*L);
    Tau = sqrt((Tau_p)^2+2*Tau_p*Tau_pp*(L/2*R)+(Tau_pp)^2);
    Sig = 6*P*W/b*t^2;
    Delta = (4*P*W^3)/(E*t^3*b);
    Pc = ((4.013*E*sqrt((t^2*b^6)/36))/W^2)*(1-(t/2*W)*sqrt(R/(4*G)));
    %% constraint functions in std form
    g(1) = Tau - Tau_max;                        % 𝜏 - 𝜏max ≤ 0
   % g(1) = Tau/Tau_max - 1;
    g(2) = Sig - Sig_max;                        % 𝜎 - 𝜎max ≤ 0 
    %g(2) = Sig/Sig_max -1;
    g(3) = h - b;                                % ℎ -b ≤ 0
    g(4) = 0.1047*h^2+0.04811*t*b*(14.0+L) -5;   % 0.10471ℎ^2 + 0.04811*𝑡*𝑏*(14.0 + 𝐿)-5 ≤ 0
    g(5) = 0.125 - h;                            % 0.125 -h ≤ 0
    g(6) = Delta - Delta_max;                            % 𝛿 - 𝛿max ≤ 0
    g(7) = P - Pc;                             % 𝑃 - 𝑃c ≤ 0 
    q = [max(g,0)];
    
    %% multi-stage assignment function.
    for i=1:length(q)
        if q(i) <= 0.001
            theta(i) = 10;
        elseif q(i) <= 0.01
            theta(i) = 20;
        elseif q(i) <= 1
            theta(i) = 100;
        else
            theta(i) = 300;
        end
    end
    %% power of the penalty function.
    for i=1:length(q)
        if q(i) < 1
            gamma(i) = 1;
        else
            gamma(i) = 2;
        end
        
    end
    for i=1:length(q)
        H_tmp(i) = theta(i)*q(i)^gamma(i);
    end
    H = sum(H_tmp);
end



