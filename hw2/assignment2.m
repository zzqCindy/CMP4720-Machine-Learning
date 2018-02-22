load('data_class3.mat');

class1 = Data{1,1};
class2 = Data{1,2};
class3 = Data{1,3};
[u1,sigma1] = MeanCov(class1');
[u2,sigma2] = MeanCov(class2');
[u3,sigma3] = MeanCov(class3');

u = {u1,u2,u3};
sigma = {sigma1,sigma2,sigma3};

v = [1,3,2;4,6,1;7,-1,0;-2,6,5];
prior = [0.6,0.2,0.2];
g = zeros(3,3);
    
% discriminant
for i = 1:3
    for j = 1:4
        g(i,j) = discriFunction(v(j,:),sigma{1,i},u{1,i},3,prior(1,i));
    end
end

[max_g,index] = max(g,[],1);
disp(['(1,3,2) is classified to class ',num2str(index(1,1))]);
disp(['(4,6,1) is classified to class ',num2str(index(1,2))]);
disp(['(7,-1,0) is classified to class ',num2str(index(1,3))]);
disp(['(-2,6,5) is classified to class ',num2str(index(1,4))]);

% Question 2
% First
mu1 = [8,2];
mu2 = [2,8];
sigma = [4.1,0;0,2.8];
r1 = mvnrnd(mu1,sigma,1000);
r2 = mvnrnd(mu2,sigma,1000);
%scope
scop1_l = min(mu1(1,1),mu2(1,1))-sqrt(sigma(1,1))*4;
scop1_r = max(mu1(1,1),mu2(1,1))+sqrt(sigma(1,1))*4;
scop2_l = min(mu1(1,2),mu2(1,2))-sqrt(sigma(2,2))*4;
scop2_r = max(mu1(1,2),mu2(1,2))+sqrt(sigma(2,2))*4;
scope_x = scop1_l:0.2:scop1_r;
scope_y = scop2_l:0.2:scop2_r;
%boundary
syms x;
g_solve1 = bound(mu1,mu2,sigma,sigma,4/5,1/5);
boundary = double(subs(g_solve1,x,scope_x));
figure('Name','First Random Point');
plot(r1(:,1),r1(:,2),'.r',r2(:,1),r2(:,2),'.k',scope_x,boundary,'-b');
%3D 
[x1,y1] = meshgrid(scope_x,scope_y);
xy = [x1(:), y1(:)];
p1 = mvnpdf(xy,mu1,sigma);
p2 = mvnpdf(xy,mu2,sigma);
P1 = reshape(p1,size(x1));
P2 = reshape(p2,size(x1));
figure('Name','First 3D')
mesh(x1,y1,P1);
hold on;
mesh(x1,y1,P2);
%posterior probability
px = P1 * 4/5 + P2 * 1/5;
post1 = (P1*4/5) ./ px;
post2 = (P2*1/5) ./ px;
figure('Name','First Posterior1');
mesh(x1,y1,post1);
figure('Name','First Posterior2');
mesh(x1,y1,post2);

% Second
sigma = [4.1,0.4;0.4,2.8];
r1 = mvnrnd(mu1,sigma,1000);
r2 = mvnrnd(mu2,sigma,1000);
%boundary
g_solve2 = bound(mu1,mu2,sigma,sigma,4/5,1/5);
boundary = double(subs(g_solve2,x,scope_x));
figure('Name','Second Random Point');
plot(r1(:,1),r1(:,2),'.r',r2(:,1),r2(:,2),'.k',scope_x,boundary,'-b');
%3D
p1 = mvnpdf(xy,mu1,sigma);
p2 = mvnpdf(xy,mu2,sigma);
P1 = reshape(p1,size(x1));
P2 = reshape(p2,size(x1));
figure('Name','Second 3D')
mesh(x1,y1,P1);
hold on;
mesh(x1,y1,P2);
%posterior probability
px = P1 * 4/5 + P2 * 1/5;
post1 = (P1*4/5) ./ px;
post2 = (P2*1/5) ./ px;
figure('Name','Second Posterior1');
mesh(x1,y1,post1);
figure('Name','Second Posterior2');
mesh(x1,y1,post2);

% Third
sigma1 = [2.1,1.5;1.5,3.8];
sigma2 = [4.1,0.4;0.4,2.8];
r1 = mvnrnd(mu1,sigma1,1000);
r2 = mvnrnd(mu2,sigma2,1000);
%scope
scop1_l = min(mu1(1,1)-sqrt(sigma1(1,1))*4,mu2(1,1)-sqrt(sigma2(1,1))*4);
scop1_r = max(mu1(1,1)+sqrt(sigma1(1,1))*4,mu2(1,1)+sqrt(sigma2(1,1))*4);
scop2_l = min(mu1(1,2)-sqrt(sigma1(2,2))*4,mu2(1,2)-sqrt(sigma2(2,2))*4);
scop2_r = max(mu1(1,2)+sqrt(sigma1(2,2))*4,mu2(1,2)+sqrt(sigma2(2,2))*4);
scope_x = scop1_l:0.2:scop1_r;
scope_y = scop2_l:0.2:scop2_r;
%boundary
g_solve3 = bound(mu1,mu2,sigma1,sigma2,4/5,1/5);
boundary = double(subs(g_solve3,x,scope_x));
figure('Name','Third Random Point');
plot(r1(:,1),r1(:,2),'.r',r2(:,1),r2(:,2),'.k',scope_x,boundary,'-b');
ylim([scop2_l,scop2_r]);
%3D
[x1,y1] = meshgrid(scope_x,scope_y);
xy = [x1(:), y1(:)];
p1 = mvnpdf(xy,mu1,sigma1);
p2 = mvnpdf(xy,mu2,sigma2);
P1 = reshape(p1,size(x1));
P2 = reshape(p2,size(x1));
figure('Name','Third 3D')
mesh(x1,y1,P1);
hold on;
mesh(x1,y1,P2);
%posterior probability
px = P1 * 4/5 + P2 * 1/5;
post1 = (P1*4/5) ./ px;
post2 = (P2*1/5) ./ px;
figure('Name','Third Posterior1');
mesh(x1,y1,post1);
figure('Name','Third Posterior2');
mesh(x1,y1,post2);

% Fourth
%boundary
g_solve4 = bound(mu1,mu2,sigma1,sigma2,1/2,1/2);
boundary = double(subs(g_solve4,x,scope_x));
figure('Name','Fourth Random Point');
plot(r1(:,1),r1(:,2),'.r',r2(:,1),r2(:,2),'.k',scope_x,boundary,'-b');
ylim([scop2_l,scop2_r]);
%3D
figure('Name','Fourth 3D')
mesh(x1,y1,P1);
hold on;
mesh(x1,y1,P2);
%posterior probability
px = P1 * 1/2 + P2 * 1/2;
post1 = (P1*1/2) ./ px;
post2 = (P2*1/2) ./ px;
figure('Name','Fourth Posterior1');
mesh(x1,y1,post1);
figure('Name','Fourth Posterior2');
mesh(x1,y1,post2);


% Question 1 Function
% hw1 mean and covariance
function [ u, sigma ] = MeanCov(x)

[number,dimension] = size(x);
u = zeros(1,dimension);
sigma = zeros(dimension,dimension);

for i = 1:dimension
    for j = 1:number
        u(1,i) = u(1,i) + x(j,i)/number;
    end
end

for i=1:dimension
    for j=1:dimension
        sigma(i,j)=sum((x(:,i)-mean(x(:,i))).*(x(:,j)-mean(x(:,j))))/(number-1);
    end
end

end

% Mahalanobis distance
function [ distance ] = mahalanobis(x, sigma, u)

[row_x,~] = size(x);
distance = sum((x - repmat(u,row_x,1)) * inv(sigma).*...
    (x - repmat(u,row_x,1)),2);

end

% discriminant function
function [ discriValue ] = discriFunction(x, sigma, u, d, prior)

mahal = mahalanobis(x,sigma,u);
discriValue = -1/2*mahal - d/2*log(2*pi)- 1/2*log(det(sigma)) + log(prior);

end

% Question 2 Function
% solve the decision boundary
function [ g_solve ] = bound(mu1,mu2,sigma1,sigma2,prior1,prior2)
    syms x y;
    g_1 = (-0.5) * ([x-mu1(1,1),y-mu1(1,2)]) * inv(sigma1) * ([x-mu1(1,1);y-mu1(1,2)]) - 0.5 * log(det(sigma1)) + log(prior1); 
    g_2 = (-0.5) * ([x-mu2(1,1),y-mu2(1,2)]) * inv(sigma2) * ([x-mu2(1,1);y-mu2(1,2)]) - 0.5 * log(det(sigma2)) + log(prior2);
    g = g_1 - g_2;
    g_solve = solve(g,'y');
end