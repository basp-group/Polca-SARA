Ut1_ = U1{t};
Ut2_ = U2{t};
z1 = size(Ut1_(1,:,:));
z2 = size(Ut2_(1,:,:));
na_t = na(t);
id = true(na_t,1);
Yt_ = Y;


% To check for U1
H{1} = @(u1,u2,a) (H1t{1,a}+ H1t{3,a})*u1 + (H1t{2,a}+ H1t{4,a})*u2;
H{2} = @(u1,u2,a) (H1t{5,a}+ H1t{7,a})*u1 + (H1t{6,a}+ H1t{8,a})*u2;
H{3} = H{1};
H{4} = H{2};

a = 1;
id(a) = false;
for i=1:4
    Yt{i} = Yt_{i,t};
    Yt_alpha{i} = nonzeros(Yt{i}(a,id));  % transpose of the rows of Y /!\ nonzeros systematically returns a column vector
    Ut1{i} = reshape(Ut1_(i,:,:),[z1(2:end) 1]);
    Ut2{i} = reshape(Ut2_(i,:,:),[z2(2:end) 1]);
    u1t{i} = flipud(Ut1{i}(:,a));
    u2t{i} = flipud(Ut2{i}(:,a));
end

% 
a1 = ((H{1}(u1t{1},u1t{2},a)));
a2 = ((H{2}(u1t{1},u1t{2},a)));
a3 = ((H{3}(u1t{3},u1t{4},a)));
a4 = ((H{4}(u1t{3},u1t{4},a)));

a1 = H1t{1,a}*u1t{1};
a2 = H1t{2,a}*u1t{2};
a2 = H1t{3,a}*u1t{1};
a4 = H1t{4,a}*u1t{2};


% To check for U2
H{1} = @(u1,u2,a) (H2t{1,a}+ H2t{2,a})*u1 + (H2t{3,a}+ H2t{4,a})*u2;
H{3} = @(u1,u2,a) (H2t{5,a}+ H2t{6,a})*u1 + (H2t{7,a}+ H2t{8,a})*u2;
H{2} = H{1};
H{4} = H{3};

for i = 1:4
    Yt_alpha{i} = nonzeros(Yt{i}(id,a));
    Ut1{i} = reshape(Ut1_(i,:,:),[z1(2:end) 1]);
    Ut2{i} = reshape(Ut2_(i,:,:),[z2(2:end) 1]);
    u1t{i} = conj(Ut1{i}(:,a));
    u2t{i} = conj(Ut2{i}(:,a));
end
% 
b1 = ((H{1}(u2t{1},u2t{3},a)));
b2 = ((H{2}(u2t{2},u2t{4},a)));
b3 = ((H{3}(u2t{1},u2t{3},a)));
b4 = ((H{4}(u2t{3},u2t{4},a)));

b1 = H2t{1,a}*u2t{1};
b2 = H2t{2,a}*u2t{1};
b3 = H2t{3,a}*u2t{3};
b4 = H2t{4,a}*u2t{3};


