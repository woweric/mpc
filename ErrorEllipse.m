 function ErrorEllipse(datamatrix,p,N,A)
% �Զ���MATLAB����ErrorEllipse����ӡ��ά���ݵ�������Բ
% �޸���ԭ����http://www.visiondummy.com/2014/04/draw-error-ellipse-representing-covariance-matrix/
% ����ֵΪһ��n��2����ֵ��������Ÿ���p
% �޸�2019/10/09
% �޸�2019/08/06
% T2������T2_limit=(N-1)A/(N-A)F(V,N-V,a)
% NΪ����������
% AΪ��ȡ�ĳɷ�����
% VΪX��ά��
% aΪ���Ŷ�,��finv������ʾʱPȡ0.95
data = datamatrix;
% ����Э���������������������ֵ
covariance = cov(data);
[eigenvec, eigenval ] = eig(covariance);
% ��ȡ�����������
[largest_eigenvec_ind_c, r] = find(eigenval == max(max(eigenval)));
largest_eigenvec = eigenvec(:, largest_eigenvec_ind_c);
% ��ȡ�������ֵ
largest_eigenval = max(max(eigenval));
% ������С������������С����ֵ
if(largest_eigenvec_ind_c == 1)
    smallest_eigenval = max(eigenval(:,2));
    smallest_eigenvec = eigenvec(:,2);
else
    smallest_eigenval = max(eigenval(:,1));
    smallest_eigenvec = eigenvec(1,:);
end
% ����X��������������ֱ�ӵļнǣ�ֵ��Ϊ[-pi,pi]
angle = atan2(largest_eigenvec(2), largest_eigenvec(1));
% ���н�Ϊ��ʱ����2pi����ֵ
if(angle < 0)
    angle = angle + 2*pi;
end
% �������ݵ����о�ֵ����ʽΪ2��1�ľ���
avg = mean(data);
% ����������Բ�Ĳ�������������ֵ����ת�Ƕȡ���ֵ���������
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
% ����ԲͶ�䵽ֱ���������� 
ellipse_x_r  = a*cos( theta_grid );
ellipse_y_r  = b*sin( theta_grid );
% ��ת����
R = [ cos(phi) sin(phi); -sin(phi) cos(phi) ];
% ��ˣ���ת��Բ
r_ellipse = [ellipse_x_r;ellipse_y_r]' * R;
% ��ӡ������Բ
plot(r_ellipse(:,1) + X0,r_ellipse(:,2) + Y0,'--k')
hold on;
% ��ӡԭʼ����
scatter(data(:,1), data(:,2), 'b');
% xlabel('LV1 Score');
% ylabel('LV2 Score');
% legend('T^2 limit','observations');
hold on;