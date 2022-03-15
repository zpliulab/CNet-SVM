function logisError = my_error(Xm,y,theta_hat)
%% the function to get errors between y and y_hat
%% this function returns the L2 norm of the error

logisError=norm((transpose(y)-(1./(1+exp(-transpose(theta_hat)*Xm'))>.5)),2);

fprintf('The L2 norm of y-y^ is: %5.2f \n',logisError)

