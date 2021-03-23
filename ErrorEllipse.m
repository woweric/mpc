 function ErrorEllipse(datamatrix,p,N,A)
% 自定义MATLAB函数ErrorEllipse：打印二维数据的置信椭圆
% 修改自原代码http://www.visiondummy.com/2014/04/draw-error-ellipse-representing-covariance-matrix/
% 输入值为一个n×2的数值矩阵和置信概率p
% 修改2019/10/09
% 修改2019/08/06
% T2控制限T2_limit=(N-1)A/(N-A)F(V,N-V,a)
% N为样本的数量
% A为提取的成分数量
% V为X的维度
% a为置信度,用finv函数表示时P取0.95
data = datamatrix;
% 计算协方差矩阵、特征向量、特征值
covariance = cov(data);
[eigenvec, eigenval ] = eig(covariance);
% 求取最大特征向量
[largest_eigenvec_ind_c, r] = find(eigenval == max(max(eigenval)));
largest_eigenvec = eigenvec(:, largest_eigenvec_ind_c);
% 求取最大特征值
largest_eigenval = max(max(eigenval));
% 计算最小特征向量和最小特征值
if(largest_eigenvec_ind_c == 1)
    smallest_eigenval = max(eigenval(:,2));
    smallest_eigenvec = eigenvec(:,2);
else
    smallest_eigenval = max(eigenval(:,1));
    smallest_eigenvec = eigenvec(1,:);
end
% 计算X轴和最大特征向量直接的夹角，值域为[-pi,pi]
angle = atan2(largest_eigenvec(2), largest_eigenvec(1));
% 当夹角为负时，加2pi求正值
if(angle < 0)
    angle = angle + 2*pi;
end
% 计算数据的两列均值，格式为2乘1的矩阵
avg = mean(data);
% 配置置信椭圆的参数，包括卡方值、旋转角度、均值、长短轴距
% chisquare_val = sqrt(chi2inv(p,2));
% f_val=finv(p,V,N-V);
f_val=finv(p,A,N-A);
% t_limit=sqrt((N-1)*A/(N-A)*f_val);
t_limit=sqrt(A*((N^2-1)/(N*(N-A)))*f_val);
theta_grid = linspace(0,2*pi);
phi = angle;
X0=avg(1);
Y0=avg(2);
a=t_limit*sqrt(largest_eigenval);
b=t_limit*sqrt(smallest_eigenval);
% 将椭圆投射到直角坐标轴中 
ellipse_x_r  = a*cos( theta_grid );
ellipse_y_r  = b*sin( theta_grid );
% 旋转矩阵
R = [ cos(phi) sin(phi); -sin(phi) cos(phi) ];
% 相乘，旋转椭圆
r_ellipse = [ellipse_x_r;ellipse_y_r]' * R;
% 打印置信椭圆
plot(r_ellipse(:,1) + X0,r_ellipse(:,2) + Y0,'--k')
hold on;
% 打印原始数据
scatter(data(:,1), data(:,2), 'b');
% xlabel('LV1 Score');
% ylabel('LV2 Score');
% legend('T^2 limit','observations');
hold on;