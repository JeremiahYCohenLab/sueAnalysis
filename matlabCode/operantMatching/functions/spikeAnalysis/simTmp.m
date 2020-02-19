clear all
y(1:5) = 0; 
for currT = 5:350
    y(currT+1) = 0.2*y(currT) + 0.2*y(currT-1) + 0.1*y(currT-2) + 0.1*y(currT-3) + 0.1*y(currT-4) + normrnd(0,1);
end
figure; 
subplot(1,3,1)
plot(y(6:end))

subplot(1,3,2)
histogram(y(6:end))

subplot(1,3,3); hold on;
aC = autocorr(y(6:end), 30);
plot(aC)


clear all
y(1:3) = 0; 
for currT = 3:350
    y(currT+1) = 0.2*y(currT) + 0.2*y(currT-1) + 0.1*y(currT-2) + normrnd(0,1);
end
figure; 
subplot(1,3,1)
plot(y)

subplot(1,3,2)
histogram(y)

subplot(1,3,3); hold on;
aC = autocorr(y, 30);
plot(aC)


%% 
clear all
for currS = 1:300
    y = [];
    y(1:5) = 0; 
    for currT = 5:350
        y(currT+1) = 0.2*y(currT) + 0.2*y(currT-1) + 0.1*y(currT-2) + 0.1*y(currT-3) + 0.1*y(currT-4) + normrnd(0,1);
    end
    aC(currS, :) = autocorr(y(6:end), 30);
end

x = [1:31];
figure; subplot(1,2,2); plotFilled(x, aC, 'b')



