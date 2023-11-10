function df = Defoc4order(par0,z);

wo = par0(1);
c = par0(2);
d = par0(3);
A = par0(4);
B = par0(5);
df = wo.*sqrt(1 + ((z-c)./d).^2 + A.*((z-c)./d).^3 + B.*((z-c)./d).^4 );