function out=Gauss2(x,p)
out= p(1)*exp(-((x-p(2))/p(3)).^2)+p(4)*exp(-((x-p(5))/p(6)).^2);