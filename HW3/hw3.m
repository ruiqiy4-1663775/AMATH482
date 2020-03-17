%% Load the data (videos) for Test 1/Case 1
clear all; close all; clc

load('cam1_1.mat')
load('cam2_1.mat')
load('cam3_1.mat')
% implay(vidFrames3_1)

%%
filter = zeros(480 ,640);  
filter(190:440, 290:410) = 1; 

data1 = load_cropped_data(vidFrames1_1, filter, 240);

%%
filter = zeros(480 ,640);  
filter(90:400, 230:360) = 1;  

data2 = load_cropped_data(vidFrames2_1, filter, 240);

%%
filter = zeros(480 ,640);  
filter(220:340, 250:490) = 1;  

data3 = load_cropped_data(vidFrames3_1, filter, 240);

%%
collected_data = collect(data1, data2, data3);

%% Results for Case 1
[m,n]= size(collected_data);  

collected_data = collected_data - repmat(mean(collected_data, 2),1,n);  
[U,S,V]= svd(collected_data'/sqrt(n-1));  
lambda = diag(S).^2; 
Y = collected_data' * V;  

figure(1)
plot(1:6, lambda/sum(lambda), 'rx', 'Linewidth', 3);
title("Case 1: Variance and energy");
xlabel("Variances "); ylabel("Proportion of energy");

figure(2)
subplot(2,1,1)
plot(1:216, collected_data(2,:), 1:216, collected_data(1,:), 'Linewidth', 2)
ylabel("Displacement (pixels)"); xlabel("Time(frames)");
title("Case 1: Displacement along Z axis and XY - plane");
legend("Z", "XY")
subplot(2,1,2)
plot(1:216, Y(:,1),'g','Linewidth', 2)
ylabel("Displacement in pixels)"); xlabel("Time (frames)");
title("Case 1: Displacement along principal component directions");
legend("principle component 1")






%% Load the data for Test 2/Case 2
clear all; close all; clc;

load('cam1_2.mat')
load('cam2_2.mat')
load('cam3_2.mat')

%%
filter = zeros(480 ,640);  
filter(200:400, 300:430) = 1; 

data1 = load_cropped_data(vidFrames1_2, filter, 240);

%%
filter = zeros(480 ,640);  
filter(50:420, 180:440) = 1;  

data2 = load_cropped_data(vidFrames2_2, filter, 240);

%%
filter = zeros(480 ,640);  
filter(180:330, 280:470) = 1; 

data3 = load_cropped_data(vidFrames3_2, filter, 240);

%%
collected_data = collect(data1, data2, data3);

%% Results for Test 2
[m,n]= size(collected_data); 

collected_data = collected_data - repmat(mean(collected_data, 2),1,n);  
[U,S,V]= svd(collected_data'/sqrt(n-1));  
lambda = diag(S).^2;  
Y = collected_data' * V;  


figure (3)
plot (1:6,lambda/sum(lambda),'bx','Linewidth', 2) ;
title("Case 2: Energy of each Diagonal Variance ") ;
xlabel("Variances ") ; ylabel("Proportion of energy") ;
figure(4)
subplot(2 ,1 ,1)
plot(1:295, collected_data(2 ,:), 1:295, collected_data(1 ,:) ,'Linewidth', 2)
ylabel(" Displacement in pixels"); xlabel ("Time ( frames )") ;
legend("Z" , " XY ")
title("Case 2: Original displacement along Z axis and XY - plane ( cam 1)") ;
subplot(2,1,2)

plot(1:295 , Y (:,1) ,1:295 , Y(:,2) ,'g','Linewidth', 2)
ylabel("Displacement in pixels") ; xlabel("Time ( frames )") ;
title("Case 2: Displacement along principal component directions") ;
legend("principle component 1 " , " principle component 2")





%% Load the data for Test 3/Case 3
clear all; close all; clc;

load('cam1_3.mat')
load('cam2_3.mat')
load('cam3_3.mat')
%%
filter = zeros(480 ,640);  
filter(240:420, 280:380) = 1;

data1 = load_cropped_data(vidFrames1_3, filter, 240);

%%
filter = zeros(480 ,640);  
filter(180:380, 240:400) = 1; 

data2 = load_cropped_data(vidFrames2_3, filter, 240);


%%
filter = zeros(480 ,640);  
filter(180:330, 240:480) = 1;  

data3 = load_cropped_data(vidFrames3_3, filter, 240);

%%
collected_data = collect(data1, data2, data3);

%% Results for Test 3
[m,n]= size(collected_data);  

collected_data = collected_data - repmat(mean(collected_data, 2),1,n);  
[U,S,V]= svd(collected_data'/sqrt(n-1)); 
lambda = diag(S).^2;  
Y = collected_data' * V;

figure (5)
plot(1:6 , lambda / sum ( lambda ) , 'bx', 'Linewidth', 2) ;
title("Case 3: Energy of each Diagonal Variance") ;
xlabel("Variances ") ; ylabel ("Proportion of energy") ;

figure (6)
subplot (2 ,1 ,1)
plot (1:235 , collected_data(2 ,:), 1:235, collected_data(1,:),'Linewidth', 2)
ylabel (" Displacement in pixels") ; xlabel("Time (frames)") ;
legend("Z" , " XY ")
title (" Case 3: Original displacement along Z axis and XY - plane ") ;
subplot (2 ,1 ,2)
plot (1:235 , Y(: ,1) , 1:235 , Y(: ,2) , 1:235 , Y(: ,3) ,1:235 , Y(: ,4) ,'r','Linewidth', 2)
ylabel(" Displacement in pixels ) ") ; xlabel("Time( frames )") ;
title(" Case 3: Displacement along principal component directions ") ;
legend(" principle component 1 " , " principle component 2 " , " principle component 3 " , " principle component 4 ")




%% Load the data for Test 4/Case 4
clear all; close all; clc;

load('cam1_4.mat')
load('cam2_4.mat')
load('cam3_4.mat')
%%
filter = zeros(480 ,640);  
filter(230:440, 330:460) = 1;  

data1 = load_cropped_data(vidFrames1_4, filter, 240);

%%
filter = zeros(480 ,640);  
filter(100:360, 230:410) = 1;  

data2 = load_cropped_data(vidFrames2_4, filter, 240);


%%
filter = zeros(480 ,640);  
filter(140:280, 320:500) = 1;  

data3 = load_cropped_data(vidFrames3_4, filter, 230);

%%
[M,I] = min(data1(1:20,2));
data1 = data1(I:end,:);
[M,I] = min(data2(1:20,2));
data2 = data2(I:end,:);
[M,I] = min(data3(1:20,2));
data3 = data3(I:end,:);

data1 = data1(1:length(data3), :);
data2 = data2(1:length(data3), :);

collected_data = [data1'; data2'; data3'];

%% Results for Test 4
[m,n]= size(collected_data);  

collected_data = collected_data - repmat(mean(collected_data, 2),1,n);  
[U,S,V]= svd(collected_data'/sqrt(n-1)); 
lambda = diag(S).^2; 
Y = collected_data' * V;  

figure(7)
plot(1:6 , lambda / sum(lambda) , 'bx', 'Linewidth', 2) ;
title(" Case 4: Energy of each Diagonal Variance") ;
xlabel("Variances") ; ylabel ("Proportion of Energy") ;
figure(8)
subplot(2 ,1 ,1)
plot(1:376 , collected_data(2 ,:) , 1:376, collected_data(1,:),'Linewidth', 2)
ylabel("Displacement in pixels") ; xlabel("Time ( frames )") ;
title("Case 4: Original displacement along Z axis and XY - plane ") ;
subplot(2 ,1 ,2)
plot(1:376 , Y(: ,1) , 1:376 , Y(: ,2) , 1:376 , Y(: ,3) , 'Linewidth', 2)
ylabel("Displacement in pixels )") ; xlabel("Time ( frames )") ;
title("Case 4: Displacement along principal component directions") ;
legend("principle component 1" , "principle component 2" , "principle component 3")