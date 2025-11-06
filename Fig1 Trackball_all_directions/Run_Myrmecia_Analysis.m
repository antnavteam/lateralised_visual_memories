load('Myrmecia_data.mat')

fps=30;
calib_0 = 1550; % Get it in ball rotation
calib_1 = 1105;
ball_diameter = 4.87; %in cm

%% get nb_frame per data

for i=1:length(data)
    nb_frame(i)=size(data{i},1);
end
hist(nb_frame)
f=find(nb_frame==min(nb_frame));
nb_frame(f)
filename{f}

%%
% the frame_window does not include the time before p
frame_window = frame_before_p + 360; % 450 =15 sec. (we made sure the last 15s were without cleaning)

turn_dyn ={};
speed_dyn ={};



for i=1:length(data)
    cur_data=data{i};
    
    X0=cur_data(1 : frame_window,1)/calib_0; % get the last Xsec (we made sure the last 15s were without cleaning)
    X1=cur_data(1:frame_window ,3)/calib_1;    
    Xs= (X0 + X1) / 2; % Make average of X0 and X1 (they should get same signal)
    
    Y0=cur_data(1 : frame_window,2) /calib_0;
    Y1=cur_data(1 : frame_window,4) /calib_1;
    
    % Path
    cum_angle{i} = pi2pi(cumsum(Xs) * 2* pi); % in radian (per frame)
    
    % Turns 

    Xs_radsec= Xs * 2*pi * fps;  % Put them in rad/sec
    turn_dyn{i} = Xs_radsec;
    
    speed = sqrt(Y0.^2 + Y1.^2) * (ball_diameter*pi) * fps; % put in cm/sec
    speed_dyn{i} = speed;

end


%% Plot Turn dynamic 
% figure()

x_axis = ((1:length(turn_dyn{1}))/fps) - 1/fps - frame_before_p/fps; % so that zero is when the ring is lifted 


cur_cond=strcmp(cond,'U'); %'U'for Unfamiliar or 'R' for familiar Route
rot_list = unique(abs(rot)); % abs, to pool +45 -45 / +90 -90 / +135 -135

     k=0;   
     for r = rot_list'; k=k+1; % go through all rotation
         
        dth=[];
        dr=[];
         j=0; 
         for a=unique(ant(cur_cond))'; j=j+1; % go through all individual
             f=find(cur_cond & ant == a & rot == r); %get the index
                dth(j,:) = turn_dyn{f}' *180/pi; % in deg/sec
                dr(j,:) = speed_dyn{f}' / fps; % so as to be back in cm/frame
              
              % now check for the opposite (-45 -90 and -135)
              if r>0 & r<3.14
              j=j+1;
              f=find(cur_cond & ant == a & rot == -r); %get the index
                dth(j,:) = - turn_dyn{f}' *180/pi; % in deg/sec
                dr(j,:) = speed_dyn{f}' / fps; % so as to be back in cm/frame
              end
              
         end
        
        subplot(2,5,k); 
        plot(x_axis, dth'); hold on
        plot(x_axis, fastsmooth(median(dth),3,3,1),'k','linewidth',2); hold on
        
        plot([min(x_axis) max(x_axis)], [0 0], 'k'); %mark the zero
        plot([0 0] , [min(min(dth)) max(max(dth))],'k'); % mark t0, when the ring is lifted
        ylim([-500 +500])
        title([num2str(r*180/pi) 'deg']);
        ylabel('tur nvel (deg/s)');
        xlabel('time (s)');
        
        subplot(2,5,k+5); %k+5 to make U second row
        a=dth(:,frame_before_p:end); % get read of frame before p
        bins=-600:6:600;
        hist(reshape(a,[size(a,1)*size(a,2), 1]),bins);
        text(-600,20,['mean=' num2str(mean(a(a<0)))]);
        text(0,40, ['mean=' num2str(mean(a(a>0)))]);
        text(-600,60, ['mean abs=' num2str(mean(mean((abs(a)))))]);
% To see each individuals clearly        
%          figure (100+k)
%          for i=1:size(dth,1)
%             subplot(4,4,i);
%             plot(dth(i,:)); hold on
%             plot(mean(dth),'k','linewidth',2)
%          end
        

     end
     
     
     
     %% Plot Turn dynamic neat example
 figure()

x_axis = ((1:length(turn_dyn{1}))/fps) - 1/fps - frame_before_p/fps; % so that zero is when the ring is lifted 
cur_cond=strcmp(cond,'R'); %'U' or 'R'
r = -45*pi/180
         
        dth=[];
        dr=[];
         j=0; 
         for a=unique(ant(cur_cond))'; j=j+1; % go through all individual
             f=find(cur_cond & ant == a & rot == r); %get the index
                dth(j,:) = turn_dyn{f}' *180/pi; % in deg/sec
                dr(j,:) = speed_dyn{f}' / fps; % so as to be back in cm/frame
          
         end
        
        plot(x_axis, dth'); hold on
        plot(x_axis, fastsmooth(median(dth),3,3,1),'k','linewidth',2); hold on
        
        plot([min(x_axis) max(x_axis)], [0 0], 'k'); %mark the zero
        plot([0 0] , [min(min(dth)) max(max(dth))],'k'); % mark t0, when the ring is lifted
        ylim([-500 +500])
        title([num2str(r*180/pi) 'deg']);
        ylabel('tur nvel (deg/s)');
        xlabel('time (s)');

    
     
       
         
%% --------------------------------------------------------ANALYSIS OF AVERAGES

%              first_frames_or_not=1
%   frame_window = 360 + frame_before_p; % 450 =15 sec. (we made sure the last 15s were without cleaning)
    
first_frames_or_not=0; % analyse the last frames because they warrant no much pauses
frame_window = 360; % 450 =15 sec. (we made sure the last 15s were without cleaning)

         
for i=1:length(data)
    cur_data=data{i};
    
    % Take last X frames or first X frames
    if first_frames_or_not == 1
        X0=cur_data(frame_before_p:frame_window,1)/calib_0; % get the last Xsec (we made sure the last 15s were without cleaning)
        X1=cur_data(frame_before_p:frame_window,3)/calib_1;    
    else
        X0=cur_data(end-frame_window+1 : end,1)/calib_0; % get the last Xsec (we made sure the last 15s were without cleaning)
        X1=cur_data(end-frame_window+1 : end,3)/calib_1;  
    end
    
    Xs= (X0 + X1) / 2; % Make average of X0 and X1 (they should get same signal)
    
    Y0=cur_data(end-frame_window+1 : end,2) /calib_0;
    Y1=cur_data(end-frame_window+1 : end,4) /calib_1;
    
    % Path
    cum_angle{i} = pi2pi(cumsum(Xs) * 2* pi); % in radian (per frame)
    
    % Turns 

    Xs_radsec= Xs * 2*pi * fps;  % Put them in rad/sec
    
    turn_l(i)= sum(Xs_radsec(Xs_radsec>0));
    turn_r(i)= sum(Xs_radsec(Xs_radsec<0));
    turn_ratio(i) = ( abs(turn_l(i)) - abs(turn_r(i)) ) / ( abs(turn_l(i)) + abs(turn_r(i)) ); % time spent turning L/R: between +1 and -1 
    turn_av(i)= mean(Xs_radsec);
    turn_dyn{i} = Xs_radsec;
    
   
    time_l(i)= sum(Xs>0) / fps; %time spent turning l or r , in sec
    time_r(i)= sum(Xs<0) / fps;
    time_ratio(i) = (time_l(i)- time_r(i))/ (time_l(i)+time_r(i)); % time spent turning L/R: between +1 and -1 

    turn_vel_l(i)= mean(Xs_radsec(Xs>0)); %how much it turn when turning in that direction
    turn_vel_r(i)= mean(Xs_radsec(Xs<0));
    
    %turn vel for prefered side and unprefered side
    if time_l(i) >= time_r(i)
        turn_vel_prefered(i) = abs(turn_vel_l(i));
        turn_vel_unprefered(i) = abs(turn_vel_r(i));
%         turn_prefered(i) = abs(turn_l(i));
%         turn_unprefered(i) = abs(turn_r(i));
%         time_prefered(i) = time_l(i);
%         time_unprefered(i) = time_r(i);
    else
        turn_vel_prefered(i) = abs(turn_vel_r(i));
        turn_vel_unprefered(i) = abs(turn_vel_l(i));
%         turn_prefered(i) = abs(turn_r(i));
%         turn_unprefered(i) = abs(turn_l(i));
%         time_prefered(i) = time_r(i);
%         time_unprefered(i) = time_l(i);
    end
    
    % Speed: fixed ant, inversion of sign means backward

    
%         % Assume no backward (so at least put them at 0)
%         if mean(Y0)>0 
%             Y0(Y0<0)=0;
%         else Y0(Y0>0)=0;
%         end
% 
%         if mean(Y1)>0 
%             Y1(Y1<0)=0;
%         else Y1(Y1>0)=0;
%         end
%          
% %     plot(Y0);hold on; plot(Y1,'r');    
%  
    speed = sqrt(Y0.^2 + Y1.^2) * (ball_diameter*pi) * fps; % put in cm/sec
    speed_av(i) = mean(speed);
    speed_dyn{i} = speed;

end

%% 

%-------------------
figure()
var=abs(time_ratio);

cur_cond=strcmp(cond,'R'); %'U' or 'R'
colour = 'r';
[t_time_ratio, rotlist] = rotation_plotter (rot,ant,cur_cond,var,colour,'time bias')

cur_cond=strcmp(cond,'U'); %'U' or 'R'
colour = 'b';
[t_time_ratio_U, rotlist] = rotation_plotter (rot,ant,cur_cond,var,colour,'time bias');

title('Time bias: Route vs Unfam')  
xlim([-180 180])
xticks([-180:45:180])


%------------------- 
figure()
var=turn_ratio;

cur_cond=strcmp(cond,'R'); %'U' or 'R'
colour = 'r';
[t_turnratio, rotlist] = rotation_plotter (rot,ant,cur_cond,var,colour,'turn bias')

cur_cond=strcmp(cond,'U'); %'U' or 'R'
colour = 'b';
[t_turnratio_U, rotlist] = rotation_plotter (rot,ant,cur_cond,var,colour,'turn bias');


title('Turn bias: Route')  
xlim([-180 180])
xticks([-180:45:180])




N=size(t_turnratio,1)


%-------------------
figure()
var=turn_av * 180/pi; %to get in degree /sec (instead of rad/sec)

cur_cond=strcmp(cond,'R'); %'U' or 'R'
colour = 'r';
rotation_plotter (rot,ant,cur_cond,var,colour,'turn deg/sec');

cur_cond=strcmp(cond,'U'); %'U' or 'R'
colour = 'b';
rotation_plotter (rot,ant,cur_cond,var,colour,'turn deg/sec');

title('Turn average: Route vs Unfam')  
xlim([-180 180])
xticks([-180:45:180])

%% Stats for time ratio


% do for all direction independantlyp=[];
P=[];
k=0;
for r= rotlist 
 k=k+1;
[P(k,1),H,STATS] = signrank(t_time_ratio(:,rotlist==r),  t_time_ratio_U(:,rotlist==r)   ,'method','exact','tail','right');
% [H,P(k,1),ci,STATS] = ttest(t_time_ratio(:,rotlist==r),  t_time_ratio_U(:,rotlist==r), 'tail','right');
% t(k,1)=STATS.tstat;
end
[rotlist',P]%,t]


%stats for turn ratio
P=[];
 [P(1,1),H,STATS] = signrank(t_turnratio(:,rotlist==+45),  t_turnratio(:,rotlist==-45)   ,'method','exact','tail','left')
 [P(2,1),H,STATS] = signrank(t_turnratio(:,rotlist==+90),  t_turnratio(:,rotlist==-90)   ,'method','exact','tail','left')
 [P(3,1),H,STATS] = signrank(t_turnratio(:,rotlist==+135),  t_turnratio(:,rotlist==-135)   ,'method','exact','tail','left')

 [P(1,1),H,STATS] = signrank(t_turnratio_U(:,rotlist==+45),  t_turnratio_U(:,rotlist==-45)   ,'method','exact','tail','left')
 [P(2,1),H,STATS] = signrank(t_turnratio_U(:,rotlist==+90),  t_turnratio_U(:,rotlist==-90)   ,'method','exact','tail','left')
 [P(3,1),H,STATS] = signrank(t_turnratio_U(:,rotlist==+135),  t_turnratio_U(:,rotlist==-135)   ,'method','exact','tail','left')

P=[];
 [H,P(1,1),STATS] = ttest(t_turnratio(:,rotlist==+45),  t_turnratio(:,rotlist==-45)  ,'tail','left')
 [H,P(2,1),STATS] = ttest(t_turnratio(:,rotlist==+90),  t_turnratio(:,rotlist==-90)  ,'tail','left')
 [H,P(3,1),STATS] = ttest(t_turnratio(:,rotlist==+135),  t_turnratio(:,rotlist==-135)   ,'tail','left')

 [H,P(1,1),STATS] = ttest(t_turnratio_U(:,rotlist==+45),  t_turnratio_U(:,rotlist==-45)  ,'tail','left')
 [H,P(2,1),STATS] = ttest(t_turnratio_U(:,rotlist==+90),  t_turnratio_U(:,rotlist==-90)   ,'tail','left')
 [H,P(3,1),STATS] = ttest(t_turnratio_U(:,rotlist==+135),  t_turnratio_U(:,rotlist==-135)  ,'tail','left')


%% Analysis turn velocity -------------------------------------
figure()

var=speed_av; % time spent turning L/R: between +1 and -1 

cur_cond=strcmp(cond,'R'); %'U' or 'R'
colour = 'r';
rotation_plotter (rot,ant,cur_cond,var,colour,'speed fwd');

cur_cond=[strcmp(cond,'U') & seq_nb<=8]; %'U' or 'R'
colour = 'b';
rotation_plotter (rot,ant,cur_cond,var,colour,'speed fwd');


title('speed fwd: Route vs Unfam')  
xlim([-180 180])
xticks([-180:45:180])


%% ---------------- turn vel prefered vs unprefered side
%---------------- turn vel prefered vs unprefered side
figure()


subplot(2,2,1);
cur_cond=strcmp(cond,'R'); %'U' or 'R'

colour = 'r';
var=turn_vel_prefered * 180/pi; % time spent turning L/R: between +1 and -1 
rotation_plotter (rot,ant,cur_cond,var,colour,'turn vel');
 
colour = 'y';
var=turn_vel_unprefered * 180/pi; % time spent turning L/R: between +1 and -1 
rotation_plotter (rot,ant,cur_cond,var,colour,'turn vel');

title('Route: turn vel pref vs. unprefered side (deg/sec)')
xlim([-180 180])
xticks([-180:45:180])
ylim([0 350])
 
   
subplot(2,2,2);
cur_cond=strcmp(cond,'U'); %'U' or 'R'

colour = 'b';
var=turn_vel_prefered * 180/pi; % time spent turning L/R: between +1 and -1 
rotation_plotter (rot,ant,cur_cond,var,colour,'turn vel');
 
colour = 'c';
var=turn_vel_unprefered * 180/pi; % time spent turning L/R: between +1 and -1 
rotation_plotter (rot,ant,cur_cond,var,colour,'turn vel');
     
title('Unfam: turn vel pref vs. unprefered side (deg/sec)')  
xlim([-180 180])
xticks([-180:45:180])
ylim([0 350])



subplot(2,2,3);
var=(turn_vel_prefered * 180/pi) - (turn_vel_unprefered * 180/pi); % time spent turning L/R: between +1 and -1 

cur_cond=strcmp(cond,'R'); %'U' or 'R'
colour = 'r';
rotation_plotter (rot,ant,cur_cond,var,colour,'turn vel diff');

cur_cond=strcmp(cond,'U'); %'U' or 'R'
colour = 'b';
rotation_plotter (rot,ant,cur_cond,var,colour,'turn vel diff');

title('turn vel pref-unpref: Route vs Unfam')  
xlim([-180 180])
xticks([-180:45:180]) 


subplot(2,2,4);
var=(turn_vel_prefered * 180/pi); % time spent turning L/R: between +1 and -1 

cur_cond=strcmp(cond,'R'); %'U' or 'R'
colour = 'r';
rotation_plotter (rot,ant,cur_cond,var,colour,'turn vel diff');

cur_cond=strcmp(cond,'U'); %'U' or 'R'
colour = 'b';
rotation_plotter (rot,ant,cur_cond,var,colour,'turn vel diff');

title('turn vel pref: Route vs Unfam')  
xlim([-180 180])
xticks([-180:45:180]) 
