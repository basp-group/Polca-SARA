% To check the operators
na = 27;

% For U1
H1a = Ha;
% a = 2;
Ua1_ = reshape(U1(:, a,:, :), [S2, P, 4]);
Ua2_ = reshape(U2(:, a,:, :), [S2, P, 4]);
id = true(na, 1);
id(a) = false;
[Ya1] = create_Ya1(Y, id,a);
for i = 1:4
    Ua1{i} = squeeze(Ua1_(:,:,i)); %reshape(Ut1_(i,:,:),[z1(2:end) 1]);
    Ua2{i} = squeeze(Ua2_(:,:,i)); %reshape(Ut2_(i,:,:),[z2(2:end) 1]);
    u1a{i} = flipud(Ua1{i});
    u2a{i} = flipud(Ua2{i});
end
id_a = find(Ya{a}{1} == 0);

h11 = squeeze(H1a(1,:,:,:) + H1a(3,:,:,:));
h12 = squeeze(H1a(2,:,:,:) + H1a(4,:,:,:));
h13 = squeeze(H1a(5,:,:,:) + H1a(7,:,:,:));
h14 = squeeze(H1a(6,:,:,:) + H1a(8,:,:,:));

H1{1} = @(u1,u2) direct_operator2(u1,u2, h11, h12, F, T, Gt, scale_t, id_a);
H1{2} = @(u1,u2) direct_operator2(u1,u2, h13, h14, F, T, Gt, scale_t, id_a);
H1{3} = @(u1,u2) direct_operator2(u1,u2, h11, h12, F, T, Gt, scale_t, id_a);
H1{4} = @(u1,u2) direct_operator2(u1,u2, h13, h14, F, T, Gt, scale_t, id_a);



a1 = H1{1}(u1a{1},u1a{2});
a2 = H1{2}(u1a{1},u1a{2});
a3 = H1{3}(u1a{3},u1a{4});
a4 = H1{4}(u1a{3},u1a{4});


% For U2
H2a = Ha;

id = true(na, 1);
id(a) = false;
[Ya2] = create_Ya2(Y, id,a);
for i = 1:4
    Ua1{i} = squeeze(Ua1_(:,:,i)); %reshape(Ut1_(i,:,:),[z1(2:end) 1]);
    Ua2{i} = squeeze(Ua2_(:,:,i)); %reshape(Ut2_(i,:,:),[z2(2:end) 1]);
    u1a{i} = conj(Ua1{i}); % [S2, P]
    u2a{i} = conj(Ua2{i});
end
id_a = find(Ya{a}{1} == 0);


h11 = squeeze(H2a(1,:,:,:) + H2a(2,:,:,:));
h12 = squeeze(H2a(3,:,:,:) + H2a(4,:,:,:));
h13 = squeeze(H2a(5,:,:,:) + H2a(6,:,:,:));
h14 = squeeze(H2a(7,:,:,:) + H2a(8,:,:,:));

H2{1} = @(u1,u2) direct_operator2(fliplr(u1),fliplr(u2), h11, h12, F, T, Gt, scale_t, id_a); % do not forget the fliplr (property F(x*)[n] = {F(x)*}[-n]) 
H2{3} = @(u1,u2) direct_operator2(fliplr(u1),fliplr(u2), h13, h14, F, T, Gt, scale_t, id_a);
H2{2} = @(u1,u2) direct_operator2(fliplr(u1),fliplr(u2), h11, h12, F, T, Gt, scale_t, id_a);
H2{4} = @(u1,u2) direct_operator2(fliplr(u1),fliplr(u2), h13, h14, F, T, Gt, scale_t, id_a);

b1 = H2{1}(u2a{1},u2a{3});
b2 = H2{2}(u2a{2},u2a{4});
b3 = H2{3}(u2a{1},u2a{3});
b4 = H2{4}(u2a{3},u2a{4});



%%
% Test direct-adjoint operators
Ha = ones(na-1, S*S, T);
x0 = randn(S*S, P);
Hx = direct_operator(x0, Ha, F, T, [], ones(T, 1), []);
y0 = randn(na-1, T);
Hty = adjoint_operator(y0, Ha, F, T, Gt, ones(T, 1), P);
a = sum(sum(Hty.*x0));
b = sum(sum(Hx.*y0)); % a = conj(b)
norm(b - conj(a))

