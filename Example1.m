clear all
n = 200;
y = 57;
initial_lambda = [5;45];
alpha_star = y+1; beta_star = n-y+1; % true values
lambda_star = [alpha_star;beta_star];

%------------------------- VB with ordinary gradient -------------------------------%
alpha = initial_lambda(1); beta = initial_lambda(2);
lambda = [alpha;beta];
niter = 40000;iter = 1;
stop = false;
while ~stop
    lambda_old = lambda;
    
    stepsize = 1/(1+iter);
    alpha = alpha-stepsize*((alpha-y)*psi(1,alpha)-(alpha+beta-n)*psi(1,alpha+beta));
    beta = beta-stepsize*((beta-n+y)*psi(1,beta)-(alpha+beta-n)*psi(1,alpha+beta));
    
    if alpha<=0 alpha = .001; end
    if beta<=0 beta = .001; end    
    lambda = [alpha;beta];  
    
    iter = iter+1;
    if (mean(abs(lambda-lambda_old))<=0.00001)||iter>niter stop = true; end
end
lambda1 = [alpha,beta];

%------------------------- Exact NGVB -------------------------------%
lambda = initial_lambda;
niter = 40000;iter = 1;
stop = false;
MSE_NGVB = 0;
while ~stop    
    lambda_old = lambda;
    stepsize = 1/(1+iter);
    
    alpha = lambda(1); beta = lambda(2);
    aux = psi(1,alpha+beta);    
    I_F = [psi(1,alpha)-aux, -aux; -aux, psi(1,beta)-aux]; % exact Fisher matrix
    Grad_LB = [(alpha_star-alpha)*(psi(1,alpha)-aux)-(beta_star-beta)*aux;
               (beta_star-beta)*(psi(1,beta)-aux)-(alpha_star-alpha)*aux];
    lambda = lambda+stepsize*(I_F\Grad_LB);   
    lambda(lambda<=0) = 0.01;   
    
    MSE_NGVB(iter) = norm(lambda-lambda_star);
    
    iter = iter+1;
    if (mean(abs(lambda-lambda_old))<=0.0001)||iter>niter stop = true; end 
end
lambda_NGVB = lambda;

%------------------------- IFVB -------------------------------%
lambda = initial_lambda;
d_lambda = 2;
S = 10; % the number of Monte Carlo samples to estimate Fisher
alpha_stepsize = 0.8; % use in the step size
niter = 100000;iter = 1;
stop = false;
H_inverse = 0.01*eye(d_lambda);
H_inverse_scale = 0; % to count the number of samples in estimating inverse Fisher
MSE_IFVB = 0;
while ~stop    
    lambda_old = lambda;
    stepsize = 10/(1+iter)^alpha_stepsize;
    
    alpha = lambda(1); beta = lambda(2);
    aux = psi(1,alpha+beta);    
    % gradient of lower bound
    Grad_LB = [(alpha_star-alpha)*(psi(1,alpha)-aux)-(beta_star-beta)*aux;
               (beta_star-beta)*(psi(1,beta)-aux)-(alpha_star-alpha)*aux];
    % update inversion-free inverse Fisher       
    for s = 1:S
        theta = betarnd(alpha,beta);    
        grad_log_q_lambda = [psi(alpha+beta)-psi(alpha)+log(theta);
            psi(alpha+beta)-psi(beta)+log(1-theta)];
        aux = H_inverse*grad_log_q_lambda;
        H_inverse = H_inverse - 1/(1+grad_log_q_lambda'*aux)*aux*(aux');
    end
    H_inverse_scale = H_inverse_scale+S;
    InverseFreeFisher = H_inverse_scale*H_inverse;
    NatGradient = InverseFreeFisher*Grad_LB;
   
    lambda = lambda+stepsize*NatGradient;    
    lambda(lambda<=0) = 0.01;   
    
    MSE_IFVB(iter) = norm(lambda-lambda_star);
    
    iter = iter+1;
    if (mean(abs(lambda-lambda_old))<=0.0001)||iter>niter stop = true; end
end
lambda_IFVB = lambda;


%------------------------- AIFVB -------------------------------%
lambda = initial_lambda;
lambda_bar = initial_lambda;
alpha_stepsize = 0.55; %use in the step size
d_lambda = 2;
S = 10; % the number of Monte Carlo samples to estimate Fisher
w_AIFVB = 2; % the parameter w used in the weighted averaged scheme
sum_AIFVB = 0;

niter = 100000;iter = 1;
stop = false;
H_inverse = 0.01*eye(d_lambda);
H_inverse_scale = 0;
MSE_AIFVB = 0;
while ~stop    
    lambda_old = lambda_bar;
    stepsize = 10/(1+iter)^alpha_stepsize;
    
    alpha = lambda(1); beta = lambda(2);
    aux = psi(1,alpha+beta);    
    % gradient of lower bound
    Grad_LB = [(alpha_star-alpha)*(psi(1,alpha)-aux)-(beta_star-beta)*aux;
               (beta_star-beta)*(psi(1,beta)-aux)-(alpha_star-alpha)*aux];
    % update inversion-free inverse Fisher.       
    % Note: we set c_beta = 0 in this example, i.e. no Gaussian random variables
    % are added in the estimation of inverse Fisher 
    alpha_bar = lambda_bar(1); beta_bar = lambda_bar(2);
    for s = 1:S
        theta = betarnd(alpha_bar,beta_bar);    
        grad_log_q_lambda = [psi(alpha_bar+beta_bar)-psi(alpha_bar)+log(theta);
            psi(alpha_bar+beta_bar)-psi(beta_bar)+log(1-theta)];
        aux = H_inverse*grad_log_q_lambda;
        H_inverse = H_inverse - 1/(1+grad_log_q_lambda'*aux)*aux*(aux');
    end
    H_inverse_scale = H_inverse_scale+S;
    InverseFreeFisher = H_inverse_scale*H_inverse;
    NatGradient = InverseFreeFisher*Grad_LB;
   
    lambda = lambda+stepsize*NatGradient;    
    lambda(lambda<=0) = 0.01;   
    
    % updated weighted averaged estimate
    sum_AIFVB = sum_AIFVB + (log(iter+1))^w_AIFVB;
    lambda_bar = lambda_bar + ((log(iter+1))^w_AIFVB/sum_AIFVB)*(lambda-lambda_bar);

    MSE_AIFVB(iter) = norm(lambda_bar-lambda_star);
    
    iter = iter+1;
    if (mean(abs(lambda_bar-lambda_old))<=0.0001)||iter>niter stop = true; end
end
lambda_AIFVB = lambda_bar;


figure(1)
fontsize = 20;
plot(MSE_IFVB,'-','LineWidth',3);
hold on
plot(MSE_AIFVB,'--','LineWidth',3);
hold off
legend('IFVB','AIFVB')
xlabel('Iteration');
ylabel('MSE');


figure(2)
x = 0:.01:1;
yy_VB = betapdf(x,lambda1(1),lambda1(2));
yy_NGVB = betapdf(x,lambda_NGVB(1),lambda_NGVB(2));
yy_IFVB = betapdf(x,lambda_IFVB(1),lambda_IFVB(2));
yy_AIFVB = betapdf(x,lambda_AIFVB(1),lambda_AIFVB(2));
yy_exact = betapdf(x,y+1,n-y+1);
plot(x,yy_exact,'-',x,yy_VB,'-',x,yy_NGVB,'-',x,yy_IFVB,'-',x,yy_AIFVB,'-','LineWidth',3);
h_legend = legend('Exact','Ordinary gradient','NGVB','IFVB','AIFVB'); 
set(h_legend,'FontSize',16)


