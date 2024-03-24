
%李相辰 202104060303

function[] = Lagrange(x,f,x0)
n = length(x) ;
m = length(x0);
    for i = 1:m
    D = x0(i);
    y = 0.0;
    for k = 1:n
    l = 1.0;
        for j = 1:n

            if j~=k
            l = l*(D-x(j))/(x(k)-x(j));
            end
        end
        %Pn(x)
        y = y + l*f(k);
    end
        xx = num2str(D,'%.4f');
        y =num2str(y,'%.4f');
        disp('f(x)的近似值点坐标为：');
        disp(['(',xx,',',y,')']);
    end
end