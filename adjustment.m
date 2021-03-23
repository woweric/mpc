load datamat.mat

%% Derive operated trajectories

X11 = X.Fs.y(k-1);
X12 = X.DO2.y(k-1);
X1 = zeros(1,2*(k-1));
for i = 1:k-1
    X1(:,2*i-1:2*i) = [X11(i) X12(i)];
end

D1k = D(:,1:k-1);D2k = D(:,k:end);
H1k = H(:,1:k-1);H2k = H(:,k:end);

GM = X1-D1k; 
H1e = [H1k GM];


%% adjuestment
% if rank(H1e) ~= rank(H1k)
% fun = @(l) norm(l'*H1k-GM);
% x0 = zeros(4,1);
% lambda = fmin(fun,x0,[],[]);
% elseif rank(H1e) == rank(H1k)
    lambda = GM*H1k'*(H1k*H1k')^(-1);
% elseif rank(H1e) == rank(H1k)
%         
% end

%% future trajectory
X2 = D2k+lambda'*H2k;
kk = length(X2);
Fs_list = 1:2:kk;
DO2_list = 2:2:kk;
Fs_full = X2(Fs_list);
DO2_full = X2(DO2_list);
Fs = Fs_full(1);
DO2 = DO2_full(1);
