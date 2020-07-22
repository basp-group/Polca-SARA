function [y,y_vec] = Phi(A,x,M)

xx = zeros(4*M,1);

for i =1:4
xg{i} = A{i}(x{i});
xx = xx + xg{i};
end

for i = 1:4
 y{i,1} = xx(i:4:end);  
 y_vec(:,i) = y{i,1};
end

end