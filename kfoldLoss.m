function rmse = kfoldLoss(x, y, cv, numHid, lr)
% Adopted from https://www.mathworks.com/matlabcentral/answers/428298-neural-network-hyperparameter-tuning

% Train net.
net = feedforwardnet(numHid, 'traingd');
net.trainParam.lr = lr;
net = train(net, x(:,:,:,cv.training), y(cv.training));
% Evaluate on validation set and compute rmse
ypred = net(x(:,:,:,cv.test));
rmse = sqrt(mean((ypred - y(cv.test)).^2));
end