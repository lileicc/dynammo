[X, model_true] = sample_clds(100, 5, 2);
Y = real(X);
figure;
plot(Y');

model_train = learn_clds(Y, 'Hidden', 2, 'MaxIter', 100);

Xhat = sample_clds(model_train, 100);
Yhat = real(Xhat);
figure;
subplot(2, 1, 1);
plot(Y');
subplot(2, 1, 2);
plot(Yhat');
