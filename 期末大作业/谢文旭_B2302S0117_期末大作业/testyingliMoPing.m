a = (1 + 1 / sqrt(3)) * (1 + 1 / sqrt(3)) / 4;
b = (1 - 1 / sqrt(3)) * (1 - 1 / sqrt(3)) / 4;
c = (1 + 1 / sqrt(3)) * (1 - 1 / sqrt(3)) / 4;
A = [a c b c;
     c a c b;
     b c a c;
     c b c a;];
 inv(A)
