function [S,W,w1,w2] = FP_negen(X,X_mean)

num_components = size(X,1);

%% G functions

%G is log(cosh(x));

g = @(x) tanh(x); %g = G'
g_der = @(x) 1 - g(x).*g(x); % g_der = g' = G''




%% Fast ICA algorithm
W = zeros(num_components);
max_iter = 100;
tol = 10^-6;

for i = 1:num_components
    %initialize random weights
    w = 2*rand(num_components, 1) - 1; %generate random vector with components between -1 and 1
    w = w/norm(w); %project onto unit circle 
    if i == 1
        w1 = w; %First row of weights
    elseif i == 2
        w2 = w; %Second row of weights
    end
    j = 1;
    distance = 0;
    while j < max_iter && distance < tol 
        %update weights 
        w_new = mean(X.*g(w'*X),2)-mean(g_der(w'*X),2)* w;
        if i>= 2
            %decorrelate weights
            w_new = w_new - sum((w_new'*W(1:i-1,:)').*W(1:i-1,:)',2);
        end
        %normalize weights
        w_new = w_new/norm(w_new); 
        
        if j > 1
            distance = abs(abs((w'*w_new))-1);
        end
        w = w_new;
        
        if i == 1
            w1 = [w1 w];
        elseif i == 2
            w2 = [w2 w];
        end
        j = j + 1;
    end
    W(i,:) = w;
end
%unmix signals
S = W*X;
%add back mean from signals
S = S + X_mean;

end

