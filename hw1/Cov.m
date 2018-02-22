% load Data Set
load('.\data_class4.mat');

% 4 classes
class1 = Data{1,1};
class2 = Data{1,2};
class3 = Data{1,3};
class4 = Data{1,4};
% transposed matrix
class1 = class1';
class2 = class2';
class3 = class3';
class4 = class4';

% Question 1.
% a)covariance, mean
mean1 = mean(class1);
mean2 = mean(class2);
mean3 = mean(class3);
mean4 = mean(class4);
cov1 = cov(class1);
cov2 = cov(class2);
cov3 = cov(class3);
cov4 = cov(class4);
% covariance 'by hand'
% initialize the matrix
c1 = zeros(2); co1 = zeros(2);
c2 = zeros(2); co2 = zeros(2);
c3 = zeros(2); co3 = zeros(2);
c4 = zeros(2); co4 = zeros(2);
% Cov(X,Y) = E((X-?)(Y-v))
for i=1:size(class1,2)
    for j=1:size(class1,2)
        c1(i,j)=sum((class1(:,i)-mean(class1(:,i))).*...
        (class1(:,j)-mean(class1(:,j))))/(size(class1,1)-1);
    end
end

for i=1:size(class2,2) 
    for j=1:size(class2,2) 
        c2(i,j)=sum((class2(:,i)-mean(class2(:,i))).*...
        (class2(:,j)-mean(class2(:,j))))/(size(class2,1)-1);
    end 
end

for i=1:size(class3,2) 
    for j=1:size(class3,2) 
        c3(i,j)=sum((class3(:,i)-mean(class3(:,i))).*...
        (class3(:,j)-mean(class3(:,j))))/(size(class3,1)-1);
    end 
end

for i=1:size(class4,2) 
    for j=1:size(class4,2) 
        c4(i,j)=sum((class4(:,i)-mean(class4(:,i))).*...
        (class4(:,j)-mean(class4(:,j))))/(size(class4,1)-1);
    end 
end
% E((X-?)(Y-v)) = E(XY)-?v
for i=1:size(class1,2)
    for j=1:size(class1,2)
        co1(i,j) = mean(class1(:,i).*class1(:,j))-mean(class1(:,i)).*mean(class1(:,j));
    end
end

for i=1:size(class2,2)
    for j=1:size(class2,2)
        co2(i,j) = mean(class2(:,i).*class2(:,j))-mean(class2(:,i)).*mean(class2(:,j));
    end
end

for i=1:size(class3,2)
    for j=1:size(class3,2)
        co3(i,j) = mean(class3(:,i).*class3(:,j))-mean(class3(:,i)).*mean(class3(:,j));
    end
end

for i=1:size(class4,2)
    for j=1:size(class4,2)
        co4(i,j) = mean(class4(:,i).*class4(:,j))-mean(class4(:,i)).*mean(class4(:,j));
    end
end

% b) eigenvectors and eigenvalues
[eigenValue1,vector1] = eigen(class1);
[eigenValue2,vector2] = eigen(class2);
[eigenValue3,vector3] = eigen(class3);
[eigenValue4,vector4] = eigen(class4);

% c) plot
% k is the coefficient to make the eigevector more visible.
k = 3;
x1 = [mean1(1,1),mean1(1,1);mean1(1,1)+k*vector1(1,1),mean1(1,1)+k*vector1(2,1)];
y1 = [mean1(1,2),mean1(1,2);mean1(1,2)+k*vector1(1,2),mean1(1,2)+k*vector1(2,2)];
x2 = [mean2(1,1),mean2(1,1);mean2(1,1)+k*vector2(1,1),mean2(1,1)+k*vector2(2,1)];
y2 = [mean2(1,2),mean2(1,2);mean2(1,2)+k*vector2(1,2),mean2(1,2)+k*vector2(2,2)];
x3 = [mean3(1,1),mean3(1,1);mean3(1,1)+k*vector3(1,1),mean3(1,1)+k*vector3(2,1)];
y3 = [mean3(1,2),mean3(1,2);mean3(1,2)+k*vector3(1,2),mean3(1,2)+k*vector3(2,2)];
x4 = [mean4(1,1),mean4(1,1);mean4(1,1)+k*vector4(1,1),mean4(1,1)+k*vector4(2,1)];
y4 = [mean4(1,2),mean4(1,2);mean4(1,2)+k*vector4(1,2),mean4(1,2)+k*vector4(2,2)];
plot(class1(:,1),class1(:,2),'.r',...
     class2(:,1),class2(:,2),'.g',...
     class3(:,1),class3(:,2),'.b',...
     class4(:,1),class4(:,2),'.c');
hold on
plot(x1,y1,'-m','LineWidth',2);
hold on
plot(x2,y2,'-m','LineWidth',2);
hold on
plot(x3,y3,'-m','LineWidth',2);
hold on
plot(x4,y4,'-m','LineWidth',2);
axis equal; 

% function to compute both eigenvector and eigenvalue.
function [value,vector] = eigen(class)
% covariance matrix
covMatrix = zeros(2);
% covMatrix = cov(class);
for i=1:size(class,2) 
    for j=1:size(class,2) 
        covMatrix(i,j)=sum((class(:,i)-mean(class(:,i))).*...
        (class(:,j)-mean(class(:,j))))/(size(class,1)-1);
    end 
end
% eigenvalue
b = covMatrix(1,1) + covMatrix(2,2);
% c = det(covMatrix);
c = covMatrix(1,1)*covMatrix(2,2)-covMatrix(1,2)*covMatrix(2,1);
value = [(b+sqrt(b^2-4*c))/2;(b-sqrt(b^2-4*c))/2];
% eigenvector
v1 = covMatrix(1,2)/(value(1,1)-covMatrix(1,1));
v2 = covMatrix(1,2)/(value(2,1)-covMatrix(1,1));
vector = [v1/sqrt(v1^2+1),1/sqrt(v1^2+1);...
          v2/sqrt(v2^2+1),1/sqrt(v2^2+1)];
end