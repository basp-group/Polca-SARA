function [b,s] = Phit(At,y,scale,L,Lt)

% L = [1,0,0,1;1,0,0,-1;0,1,1,0]; % Conversion matrix
% Lt = 0.5*(conj(L))'; %Adjoint conversion matrix

yy = y.';
yy = yy(:);

for i = 1:4
  b{i} = At{i}(yy);  
  b{i} = scale * b{i}(:);
end

 s = mat2cell([b{1:4}]*(Lt),size(b{1},1), [1 1 1])';

 for i =1:3
     s{i} = real(s{i});
 end

end