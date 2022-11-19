function [] = compute_metrics(y, prediction,d)
%COMPUTE_METRICS compute the survival analysis metrics of interest
died_y = y(d==1,:);
died_p = prediction(d==1,:);
disp('C-Index')
disp(sum(sum((died_y<=transpose(y))&(died_p<=transpose(prediction))))/sum(sum(died_y<=transpose(y))));

disp('Average Underestimated Survival (Churned)')
disp(-mean(mean(prediction((prediction<y)&(d==1))-y((prediction<y)&(d==1)))))

disp('Average Underestimated Survival (Not Churned)')
disp(-mean(mean(prediction((prediction<y)&(d==0))-y((prediction<y)&(d==0)))))
end

