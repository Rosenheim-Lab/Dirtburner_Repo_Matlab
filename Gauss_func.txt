function out=Gauss3(x,p)
for j=1:length(p)/3
	terms(:,j) = p(3*j-2)*exp(-((x-p(3*j-1))/p(3*j)).^2);
end

out = sum(terms, 2)