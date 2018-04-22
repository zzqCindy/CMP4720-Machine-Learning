clear all

load('test_train_data_class3');

TrainPts = Data.train;
TestPts = Data.test;

% load and transpose the matrix
train1 = TrainPts{1,1}';
train2 = TrainPts{1,2}';
train3 = TrainPts{1,3}';
test1 = TestPts{1,1}';
test2 = TestPts{1,2}';
test3 = TestPts{1,3}';

% a) maximum likelihood estimatesmax_p1
mu1 = mean(train1);
mu2 = mean(train2);
mu3 = mean(train3);
mu = repmat(mu1,size(train1,1),1);
var1 = (train1-mu)'*(train1-mu) ./ size(train1,1);
mu = repmat(mu2,size(train2,1),1);
var2 = (train2-mu)'*(train2-mu) ./ size(train2,1);
mu = repmat(mu3,size(train3,1),1);
var3 = (train3-mu)'*(train3-mu) ./ size(train3,1);

% c) classify
mu = [mu1;mu2;mu3];
var = [var1;var2;var3];
gx1 = classify(test1,var,mu);
gx2 = classify(test2,var,mu);
gx3 = classify(test3,var,mu);
[max_g1,index1] = max(gx1,[],2);
[max_g2,index2] = max(gx2,[],2);
[max_g3,index3] = max(gx3,[],2);

% plot discimination function
%syms x;
%scope_x = -4:0.1:4;
%g_solve1 = bound(mu1,mu2,var1,var2,1/3,1/3);
%boundary1 = double(subs(g_solve1,x,scope_x));
%g_solve2 = bound(mu1,mu3,var1,var3,1/3,1/3);
%boundary2 = double(subs(g_solve2,x,scope_x));
%g_solve3 = bound(mu2,mu3,var2,var3,1/3,1/3);
%boundary3 = double(subs(g_solve3,x,scope_x));

%figure(1);
%plot(train1(:,1),train1(:,2),'.r',train2(:,1),train2(:,2),'.b',train3(:,1),train3(:,2),'.k',...
%    scope_x,boundary1,'-r',scope_x,boundary2,'-b',scope_x,boundary3,'-k');
%axis equal,
%hold off;

% confusion matrix
confu = zeros(3,3);
for i = 1:size(index1,1)
    confu(index1(i,1),1) = confu(index1(i,1),1) + 1;
end

for i = 1:size(index2,1)
    confu(index2(i,1),2) = confu(index2(i,1),2) + 1;
end

for i = 1:size(index3,1)
    confu(index3(i,1),3) = confu(index3(i,1),3) + 1;
end

% error rate
error = zeros(3,1);
for i = 1:3
    error(i,1) = 1 - confu(i,i)/sum(confu(:,i));
end

% d) Bayesian estimates
[thetaTrain1, rhoTrain1] = cart2pol(train1(:,1),train1(:,2));
[thetaTrain2, rhoTrain2] = cart2pol(train2(:,1),train2(:,2));
[thetaTrain3, rhoTrain3] = cart2pol(train3(:,1),train3(:,2));
[thetaTest1, rhoTest1] = cart2pol(test1(:,1),test1(:,2));
[thetaTest2, rhoTest2] = cart2pol(test2(:,1),test2(:,2));
[thetaTest3, rhoTest3] = cart2pol(test3(:,1),test3(:,2));

% scatter figure after tranform
figure(2);
scatter(thetaTrain1,rhoTrain1);
hold on;
scatter(thetaTrain2,rhoTrain2);
hold on;
scatter(thetaTrain3,rhoTrain3);
axis equal,
hold off;

rhoMean1 = mean(rhoTrain1);
rhoMean2 = mean(rhoTrain2);
rhoMean3 = mean(rhoTrain3);

% estimate for mu
n = size(rhoTrain1,1);
rhoVar1 = (0.25 * 100) / (n*100 + 0.25);
rhoMu1 = (n * rhoMean1 / 0.25) * rhoVar1;
n = size(rhoTrain2,1);
rhoVar2 = (0.25 * 100) / (n*100 + 0.25);
rhoMu2 = (n * rhoMean2 / 0.25) * rhoVar2;
n = size(rhoTrain3,1);
rhoVar3 = (0.25 * 100) / (n*100 + 0.25);
rhoMu3 = (n * rhoMean3 / 0.25) * rhoVar3;

% posterier function
sym mu
post1 = 1 / sqrt(2*pi*rhoVar1) * exp(-1/2 * power((mu-rhoMu1),2) / rhoVar1);
post2 = 1 / sqrt(2*pi*rhoVar2) * exp(-1/2 * power((mu-rhoMu2),2) / rhoVar2);
post3 = 1 / sqrt(2*pi*rhoVar3) * exp(-1/2 * power((mu-rhoMu3),2) / rhoVar3);

% classify
mun = [rhoMu1; rhoMu2; rhoMu3];
varn = [rhoVar1; rhoVar2; rhoVar3];
posterior1 = bayesClassify(rhoTest1,mun,varn);
posterior2 = bayesClassify(rhoTest2,mun,varn);
posterior3 = bayesClassify(rhoTest3,mun,varn);
[max_p1,index1] = max(posterior1,[],2);
[max_p2,index2] = max(posterior2,[],2);
[max_p3,index3] = max(posterior3,[],2);

confu2 = zeros(3,3);
for i = 1:size(index1,1)
    confu2(index1(i,1),1) = confu2(index1(i,1),1) + 1;
end

for i = 1:size(index2,1)
    confu2(index2(i,1),2) = confu2(index2(i,1),2) + 1;
end

for i = 1:size(index3,1)
    confu2(index3(i,1),3) = confu2(index3(i,1),3) + 1;
end

% error rate
error2 = zeros(3,1);
for i = 1:3
    error2(i,1) = 1 - confu2(i,i)/sum(confu2(:,i));
end

% b) HW2 Q1
% Mahalanobis distance
function [ distance ] = mahalanobis(x, sigma, u)

[row_x,~] = size(x);
distance = sum((x - repmat(u,row_x,1)) * inv(sigma).*...
    (x - repmat(u,row_x,1)),2);

end

% discriminant function
function [ discriValue ] = discriFun(x, sigma, u, d, prior)

mahal = mahalanobis(x,sigma,u);
discriValue = -1/2*mahal - d/2*log(2*pi)- 1/2*log(det(sigma)) + log(prior);

end

% classify
function [ gx ] = classify(test,var,mu)

gx = zeros(size(test,1),3);
for i = 1:size(gx,1)
    for j = 1:3
        gx(i,j) = discriFun(test(i,:),var(2*j-1:2*j,:),mu(j,:),2,1/3);
    end
end

end

% bound
function [ g_solve ] = bound(mu1,mu2,sigma1,sigma2,prior1,prior2)
    syms x y;
    g_1 = (-0.5) * ([x-mu1(1,1),y-mu1(1,2)]) * inv(sigma1) * ([x-mu1(1,1);y-mu1(1,2)]) - 0.5 * log(det(sigma1)) + log(prior1); 
    g_2 = (-0.5) * ([x-mu2(1,1),y-mu2(1,2)]) * inv(sigma2) * ([x-mu2(1,1);y-mu2(1,2)]) - 0.5 * log(det(sigma2)) + log(prior2);
    g = g_1 - g_2;
    g_solve = solve(g,'y');
end

% density estimate
function [ postrior ] = bayesClassify(test,mu,var)

postrior = zeros(size(test,1),3);
for i = 1:size(test,1)
    for j = 1:3
        postrior(i,j) = 1 / (2*pi*sqrt(var(j,1)*0.25)) *...
            exp(-1/2 * power((test(i,1)-mu(j,1)),2) / (var(j,1)+0.25));
    end
end

end