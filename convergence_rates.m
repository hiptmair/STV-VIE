function [order, orderH1] = convergence_rates(errors, errorsH1, Niter, steps)


for count = 1:Niter-1


if count > 1
order(1,count) = -log(errors(1,count)/errors(1,count-1))/log(steps(1,count-1)/steps(1,count)); 

order(2,count) = -log(errors(2,count)/errors(2,count-1))/log(steps(1,count-1)/steps(1,count)); 

order(3,count) = -log(errors(3,count)/errors(3,count-1))/log(steps(1,count-1)/steps(1,count)); 
% 
% order(4,count) = -log(errors(4,count)/errors(4,count-1))/log(steps(1,count-1)/steps(1,count)); 
% 
% order(5,count) = -log(errors(5,count)/errors(5,count-1))/log(steps(1,count-1)/steps(1,count)); 

orderH1(1,count) = -log(errorsH1(1,count)/errorsH1(1,count-1))/log(steps(1,count-1)/steps(1,count)); 

orderH1(2,count) = -log(errorsH1(2,count)/errorsH1(2,count-1))/log(steps(1,count-1)/steps(1,count)); 

orderH1(3,count) = -log(errorsH1(3,count)/errorsH1(3,count-1))/log(steps(1,count-1)/steps(1,count)); 

% orderH1(4,count) = -log(errorsH1(4,count)/errorsH1(4,count-1))/log(steps(1,count-1)/steps(1,count)); 
% 
% orderH1(5,count) = -log(errorsH1(5,count)/errorsH1(5,count-1))/log(steps(1,count-1)/steps(1,count)); 

end

end