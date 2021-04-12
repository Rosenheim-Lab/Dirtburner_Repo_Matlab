function out=Gauss2(x,p)
out= p(1)*exp(-((x-p(2))/p(3)).^2);