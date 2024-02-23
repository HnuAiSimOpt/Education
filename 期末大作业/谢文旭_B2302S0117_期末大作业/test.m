lc1 = [-1/sqrt(3) 1/sqrt(3)];
count = 1;
for i = 1 : 2
    for j = 1 : 2
        for k = 1 : 2
            A(count, 1) = 0.125 * (1 - lc1(i)) * (1 - lc1(j)) * (1 + lc1(k));
            A(count, 2) = 0.125 * (1 - lc1(i)) * (1 + lc1(j)) * (1 + lc1(k));
            A(count, 3) = 0.125 * (1 + lc1(i)) * (1 + lc1(j)) * (1 + lc1(k));
            A(count, 4) = 0.125 * (1 + lc1(i)) * (1 - lc1(j)) * (1 + lc1(k));
            A(count, 5) = 0.125 * (1 - lc1(i)) * (1 - lc1(j)) * (1 - lc1(k));
            A(count, 6) = 0.125 * (1 - lc1(i)) * (1 + lc1(j)) * (1 - lc1(k));
            A(count, 7) = 0.125 * (1 + lc1(i)) * (1 + lc1(j)) * (1 - lc1(k));
            A(count, 8) = 0.125 * (1 + lc1(i)) * (1 - lc1(j)) * (1 - lc1(k));
            count = count + 1;
        end
    end
end
k = inv(A);

a1 = [1/6 4/6 1/6];
a2 = [1/6 1/6 4/6];
b = [-1/sqrt(3) 1/sqrt(3)];
count = 1;
for i = 1 : 3
    for j = 1 : 2
        B(count, 1) = 0.5 * (1 - b(j)) * (1 - a1(i) - a2(i));
        B(count, 2) = 0.5 * (1 - b(j)) * a1(i);
        B(count, 3) = 0.5 * (1 - b(j)) * a2(i);
        B(count, 4) = 0.5 * (1 + b(j)) * (1 - a1(i) - a2(i));
        B(count, 5) = 0.5 * (1 + b(j)) * a1(i);
        B(count, 6) = 0.5 * (1 + b(j)) * a2(i);
        count = count + 1;
    end
end
k2 = inv(B);

count = 1;
for i = 1 : 2
    for j = 1 : 2
        C(count, 1) = 0.25 * (1 + lc1(i)) * (1 + lc1(j));
        C(count, 2) = 0.25 * (1 - lc1(i)) * (1 + lc1(j));
        C(count, 3) = 0.25 * (1 - lc1(i)) * (1 - lc1(j));
        C(count, 4) = 0.25 * (1 + lc1(i)) * (1 - lc1(j));
        count = count + 1;
    end
end
k3 = inv(C);



        
