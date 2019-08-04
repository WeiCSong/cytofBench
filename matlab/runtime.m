clear;

method={'Accense','ACDC','Depeche','Flowmean','Flowsom','kmeans','LDA','Phenograph','Xshift'};

%% Cell-cycle data
t=[3.1709 2.4123 2.3798 2.6068 0.7310 -0.0482 -0.3423 2.2105 2.5370];
t=10.^t;
%ac=[0.3529 0.8342 0.2808 0.3055 0.6605 0.3969 0.9095 0.2309 0.3622];
fm=[0.3500 0.8466 0.3551 0.3506 0.6897 0.4224 0.9110 0.2789 0.3970];
scatter(fm,t,60,'filled')
for i=1:9
    text(fm(i),t(i),method{i},'FontSize',12)
end
xlabel('Fmeasure','FontSize',14);
ylabel('runtime(second)','FontSize',14);
title('Runtime vs. Fmeasure: Cell-cycle data','FontSize',14)
figure;

%% Colon data
t=[3.1887 2.1829 2.4951 2.4804 0.6990 -0.5211 -0.3656 1.7898 2.1908];
t=10.^t;
fm=[0.3495 0.7874 0.6898 0.5901 0.5888 0.4951 0.8587 0.3994 0.3094];
scatter(fm,t,60,'filled')
for i=1:9
    text(fm(i),t(i),method{i},'FontSize',12)
end
xlabel('Fmeasure','FontSize',14);
ylabel('runtime(second)','FontSize',14);
title('Runtime vs. Fmeasure: Colon data','FontSize',14)
figure;

%% Levine32 data
t=[3.1768 2.2671 2.3003 2.8573 0.8518 -0.2924 -0.3768 1.8950 2.2342];
t=10.^t;
fm=[0.6008 0.9939 0.9231 0.9279 0.8979 0.6748 0.9807 0.7062 0.7706];
scatter(fm,t,60,'filled')
for i=1:9
    text(fm(i),t(i),method{i},'FontSize',12)
end
xlabel('Fmeasure','FontSize',14);
ylabel('runtime(second)','FontSize',14);
title('Runtime vs. Fmeasure: Levine32dim data','FontSize',14)
figure;

%% Muscle data
t=[3.1247 2.2287 2.2073 3.0570 0.6448 0.0338 -0.3989 1.9556 2.2611];
t=10.^t;
fm=[0.4072 0.8784 0.8850 0.7928 0.8286 0.6207 0.9238 0.4336 0.5133];
scatter(fm,t,60,'filled')
for i=1:9
    text(fm(i),t(i),method{i},'FontSize',12)
end
xlabel('Fmeasure','FontSize',14);
ylabel('runtime(second)','FontSize',14);
title('Runtime vs. Fmeasure: Muscle data','FontSize',14)
figure;

%% Samusik data
t=[3.1832 2.5164 2.3815 2.7195 1.0667 -0.2366 -0.3098 2.0074 2.2815];
t=10.^t;
fm=[0.6376 0.9731 0.8677 0.9099 0.8417 0.5535 0.9759 0.9249 0.9091];
scatter(fm,t,60,'filled')
for i=1:9
    text(fm(i),t(i),method{i},'FontSize',12)
end
xlabel('Fmeasure','FontSize',14);
ylabel('runtime(second)','FontSize',14);
title('Runtime vs. Fmeasure: Samusik01 data','FontSize',14)
figure;

%% Levine13 data
t=[3.1650 2.3380 2.3594 1.7178 0.8778 -0.0969 -0.6198 1.6708 2.0413];
t=10.^t;
fm=[0.6745 0.9332 0.8010 0.8842 0.8760 0.6293 0.9586 0.9393 0.7606];
scatter(fm,t,60,'filled')
for i=1:9
    text(fm(i),t(i),method{i},'FontSize',12)
end
xlabel('Fmeasure','FontSize',14);
ylabel('runtime(second)','FontSize',14);
title('Runtime vs. Fmeasure: Levine13dim data','FontSize',14)




