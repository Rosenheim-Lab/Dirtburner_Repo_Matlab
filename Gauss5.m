function out=Gauss5(x,p)
out= p(1)*exp(-((x-p(2))/p(3)).^2)+p(4)*exp(-((x-p(5))/p(6)).^2)+p(7)*exp(-((x-p(8))/p(9)).^2)+p(10)*exp(-((x-p(11))/p(12)).^2)+p(13)*exp(-((x-p(14))/p(15)).^2);