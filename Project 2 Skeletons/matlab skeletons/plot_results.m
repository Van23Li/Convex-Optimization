function [] = plot_results(y,prediction, d, titlename)
%PLOT_RESULTS plot the results
figure('Name',titlename);
subplot(1,3,1);
scatter(prediction, y, 5,'filled');
hold on;
plot([0,72],[0,72]);
xlabel('predicted');
ylabel('actual');
mae = mean(abs(y-prediction));
title(sprintf('MAE %f', mae));

subplot(1,3,2);
scatter(prediction(d==1,:), y(d==1,:), 5,'filled');
hold on;
plot([0,72],[0,72]);
xlabel('predicted');
ylabel('actual');
mae = mean(abs(y(d==1,:)-prediction(d==1,:)));
title(sprintf('MAE (Churned) %f', mae));

subplot(1,3,3);
scatter(prediction(d==0,:), y(d==0,:), 5,'filled');
hold on;
plot([0,72],[0,72]);
xlabel('predicted');
ylabel('actual');
mae = mean(abs(y(d==0,:)-prediction(d==0,:)));
title(sprintf('MAE (Not Churned) %f', mae));

figure('Name',titlename);
subplot(1,2,1)
histogram(y(d==1)-prediction(d==1));
title('Error if Churned');
grid on;
xlabel('y-f(x)');
subplot(1,2,2)
histogram(y(d==0)-prediction(d==0));
title('Error if Not Churned');
grid on;
xlabel('y-f(x)');

end

