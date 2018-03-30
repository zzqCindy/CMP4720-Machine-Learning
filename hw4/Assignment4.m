clear all

% data of the table in the textbook
pxw1 = [.42,-.087,.58;-.2,-3.3,-3.4;
        1.3,-.32,1.7;.39,.71,.23;
        -1.6,-5.3,-.15;-.029,.89,-4.7;
        -.23,1.9,2.2;.27,-.3,-.87;
        -1.9,.76,-2.1;.87,-1,-2.6];
pxw2 = [-.4,.58,.089;-.31,.27,-.04;
        .38,.055,-.035;-.15,.53,.011;
        -.35,.47,.034;.17,.69,.1;
        -.011,.55,-.18;-.27,.61,.12;
        -.065,.49,.0012;-.12,.054,-.063];
pxw3 = [.83,1.6,-.014;1.1,1.6,.48;
        -.44,-.41,.32;.047,-.45,1.4;
        .28,.35,3.1;-.39,-.48,.11;
        .34,-.079,.14;-.3,-.22,2.2;
        1.1,1.2,-.46;.18,-.11,-.49];

% a) three features xi of category w1
[mu1,sigma1] = ML_estimate(pxw1(:,1));
[mu2,sigma2] = ML_estimate(pxw1(:,2));
[mu3,sigma3] = ML_estimate(pxw1(:,3));
disp('a) estimation for mean and variance for w1')
disp(['mean for x1 = ', num2str(mu1)]);
disp(['variance for x1 = ',num2str(sigma1)]);
disp(['mean for x2 = ', num2str(mu2)]);
disp(['variance for x2 = ',num2str(sigma2)]);
disp(['mean for x3 = ', num2str(mu3)]);
disp(['variance for x3 = ',num2str(sigma3)]);

% b) two-dimension
[mu1,sigma1] = ML_estimate2(pxw1(:,1:2));
[mu2,sigma2] = ML_estimate2(pxw1(:,2:3));
[mu3,sigma3] = ML_estimate2(pxw1(:,[1,3]));
disp('b) 2-dimension estimation for mean and variance for w1')
disp('mean for x1,x2 = '); disp(mu1);
disp('variance for x1,x2 = '); disp(sigma1);
disp('mean for x2,x3 = '); disp(mu2);
disp('variance for x2,x3 = '); disp(sigma2);
disp('mean for x1,x3 = '); disp(mu3);
disp('variance for x1,x3 = '); disp(sigma3);

% c) three-dimension
[mu,sigma] = ML_estimate2(pxw1);
disp('c) 3-dimension estimation for mean and variance for w1')
disp('mean for x1,x2,x3 of w1 = '); disp(mu);
disp('variance for x1,x2,x3 of w1 = '); disp(sigma);

% d) diagnoal
[mu,sigma] = ML_diagonal(pxw2);
disp('d) mean and diagnoal variance for w2')
disp('mean for x1,x2,x3 of w2 = '); disp(mu);
disp('variance for x1,x2,x3 of w2 = '); disp(sigma);

% e/f) mean and variance of each feature
[mu1,sigma1] = ML_estimate2([pxw1(:,1),pxw2(:,1),pxw3(:,1)]);
[mu2,sigma2] = ML_estimate2([pxw1(:,2),pxw2(:,2),pxw3(:,2)]);
[mu3,sigma3] = ML_estimate2([pxw1(:,3),pxw2(:,3),pxw3(:,3)]);
disp('e/f) mean and variance of each feature x1,x2,x3 using function ML_estimate2')
disp('mean for x1 = '); disp(mu1);
disp('variance for x1 = '); disp(sigma1);
disp('mean for x2 = '); disp(mu2);
disp('variance for x2 = '); disp(sigma2);
disp('mean for x3 = '); disp(mu3);
disp('variance for x3 = '); disp(sigma3);
[mu1,sigma1] = ML_diagonal([pxw1(:,1),pxw2(:,1),pxw3(:,1)]);
[mu2,sigma2] = ML_diagonal([pxw1(:,2),pxw2(:,2),pxw3(:,2)]);
[mu3,sigma3] = ML_diagonal([pxw1(:,3),pxw2(:,3),pxw3(:,3)]);
disp('e/f) mean and variance of each feature x1,x2,x3 using function ML_diagonal')
disp('mean for x1 = '); disp(mu1);
disp('variance for x1 = '); disp(sigma1);
disp('mean for x2 = '); disp(mu2);
disp('variance for x2 = '); disp(sigma2);
disp('mean for x3 = '); disp(mu3);
disp('variance for x3 = '); disp(sigma3);

    
% a) maximum-likelihood
function [mu,sigma] = ML_estimate(x)
    len = length(x);
    mu = sum(x) / len;
    sigma = sum((x - mu).^2) / len;
end

% b/c) multi-dimensional Gaussian
function [mu,sigma] = ML_estimate2(x)
    len = size(x,1);
    mu = sum(x) / len;
    tmp = x - repmat(mu,len,1);
    sigma = (tmp'*tmp) / len;
end

% d) diagonal component of coveriance
function [mu,sigma] = ML_diagonal(x)
    len = size(x,1);
    mu = sum(x) / len;
    tmp = x - repmat(mu,len,1);
    sigma = (tmp'*tmp) / len;
    sigma = diag(sigma);
end