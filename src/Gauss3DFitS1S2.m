function g = Gauss3DFitS1S2(par0,mp);

mu1 = par0(1);
mu2 = par0(2);
s1 = par0(3);
s2 = par0(4);
A = par0(5);
[X1,Y1] = meshgrid(1:1:mp);
g = A*exp(-((X1-mu1).^2)./(2*s1^2)-((Y1-mu2).^2)./(2*s2^2));

