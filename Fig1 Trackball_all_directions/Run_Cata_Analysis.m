load('Cata_Data.mat')
frame_before_p = 40;

cur_cond=strcmp(cond,'R');

ant=ant(cur_cond)
data=data(cur_cond)
filename = filename(cur_cond)
rot=rot(cur_cond)
cond=cond(cur_cond)

fps=30;
calib_0 = 1550; % Get it in ball rotation
calib_1 = 1105;
ball_diameter = 4.87; %in cm

%% get nb_frame per data

for i=1:length(data)
    nb_frame(i)=size(data{i},1);
end
hist(nb_frame,30)
f=find(nb_frame==min(nb_frame));
nb_frame(f)
filename{f}

%%
frame_window = 360 + frame_before_p; % 450 =15 sec. (we made sure the last 15s were without cleaning)

turn_dyn ={};
speed_dyn ={};



for i=1:length(data)
    cur_data=data{i};
    
    % Take first  X frames or first X frames

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
        plot([min(x_axis) max(x_axis)], [0 0], 'k'); hold on%mark the zero
        plot([0 0] , [-800 800],'k'); hold on% mark t0, when the ring is lifted
        ylim([-800 +800])
        title([num2str(r*180/pi) 'deg']);
        ylabel('tur nvel (deg/s)');
        xlabel('time (s)');
        
        subplot(2,5,k+5); %k+5 to make U second row
        a=dth(:,frame_before_p:end); % get read of frame before p
        bins=-1000:10:1000;
        hist(reshape(a,[size(a,1)*size(a,2), 1]),bins);
        text(-1000,20,['mean=' num2str(mean(a(a<0)))]);
        text(0,40, ['mean=' num2str(mean(a(a>0)))]);
        text(-1000,60, ['mean abs=' num2str(mean(mean((abs(a)))))]);
        
      

     end
%% Plot Turn dynamic keeping + and - separated
 figure()

x_axis = ((1:length(turn_dyn{1}))/fps) - 1/fps - frame_before_p/fps; % so that zero is when the ring is lifted 


rot_list = unique(abs(rot)); % abs, to pool +45 -45 / +90 -90 / +135 -135

     k=0;   
     for r = rot_list'; k=k+1; % go through all rotation
         
        dth=[];
        dr=[];
        dthn=[];
        drn=[];
         j=0; i=0;
         for a=unique(ant(cur_cond))'; j=j+1; % go through all individual
             f=find(cur_cond & ant == a & rot == r); %get the index
                dth(j,:) = turn_dyn{f}' *180/pi; % in deg/sec
                dr(j,:) = speed_dyn{f}' / fps; % so as to be back in cm/frame
              
              % now check for the opposite (-45 -90 and -135)
              if r>0 & r<3.14
              i=i+1;
              f=find(cur_cond & ant == a & rot == -r); %get the index
                dthn(i,:) = turn_dyn{f}' *180/pi; % in deg/sec
                drn(i,:) = speed_dyn{f}' / fps; % so as to be back in cm/frame
              end
              
         end
        
        subplot(2,5,k); 
        plot(x_axis, dth','b'); hold on
        if r>0 & r<3.14
        plot(x_axis, dthn','r'); hold on
        plot(x_axis, fastsmooth(median(dthn),3,3,1),'k','linewidth',2); hold on
        end
        plot(x_axis, fastsmooth(median(dth),3,3,1),'k','linewidth',2); hold on
        plot([min(x_axis) max(x_axis)], [0 0], 'k'); hold on%mark the zero
        plot([0 0] , [-800 800],'k'); hold on% mark t0, when the ring is lifted
        ylim([-800 +800])
        title([num2str(r*180/pi) 'deg']);
        ylabel('tur nvel (deg/s)');
        xlabel('time (s)');
           
      

     end

     
 %% Plot Turn dynamic single individual
for a = unique(ant(cur_cond))'
    figure(a)

x_axis = ((1:length(turn_dyn{1}))/fps) - 1/fps - frame_before_p/fps; % so that zero is when the ring is lifted 


rot_list = unique(abs(rot)); % abs, to pool +45 -45 / +90 -90 / +135 -135

     k=0;   
     for r = rot_list'; k=k+1; % go through all rotation
         
        dth=[];
        dr=[];
         j=0; 
         
         %for a=unique(ant(cur_cond))';
         j=j+1; % go through all individual
             f=find(cur_cond & ant == a & rot == r); %get the index
                dth(j,:) = turn_dyn{f}' *180/pi; % in deg/sec
                dr(j,:) = speed_dyn{f}' / fps; % so as to be back in cm/frame
              
              % now check for the opposite (-45 -90 and -135)
              if r>0 & r<3.14
              j=j+1;
              f=find(cur_cond & ant == a & rot == -r); %get the index
                dth(j,:) = turn_dyn{f}' *180/pi; % in deg/sec
                dr(j,:) = speed_dyn{f}' / fps; % so as to be back in cm/frame
              end
              
         %end
        
        subplot(1,5,k); 
        plot(x_axis, dth'); hold on
        plot([min(x_axis) max(x_axis)], [0 0], 'k'); hold on%mark the zero
        plot([0 0] , [-800 800],'k'); hold on% mark t0, when the ring is lifted
        ylim([-800 +800])
        title([num2str(r*180/pi) 'deg']);
        ylabel('tur nvel (deg/s)');
        xlabel('time (s)');
        
 
      

     end
end    
 %% Plot Turn dynamic neat example
 figure()

x_axis = ((1:length(turn_dyn{1}))/fps) - 1/fps - frame_before_p/fps; % so that zero is when the ring is lifted 
r = 45*pi/180
         
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
        ylim([-1000 +1000])
        title([num2str(r*180/pi) 'deg']);
        ylabel('tur nvel (deg/s)');
        xlabel('time (s)');
        
        
        % To see each individuals clearly        
         figure (100+1)
         for i=1:size(dth,1)
            subplot(4,5,i);
            plot(x_axis, dth(i,:)); hold on
%             plot(x_axis, fastsmooth(median(dth),3,3,1),'k','linewidth',2); hold on
            plot([min(x_axis) max(x_axis)], [0 0], 'k'); %mark the zero
            plot([0 0] , [min(min(dth)) max(max(dth))],'k'); % mark t0, when the ring is lifted
            ylim([-1000 +1000])
         end
         
         
         
         
%% --------------------------------------------------------ANALYSIS OF AVERAGES

%              first_frames_or_not=1
%   frame_window = 360 + frame_before_p; % 450 =15 sec. (we made sure the last 15s were without cleaning)
   
first_frames_or_not=0; % analyse the last frames like myrmecia
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

    else
        turn_vel_prefered(i) = abs(turn_vel_r(i));
        turn_vel_unprefered(i) = abs(turn_vel_l(i));

    end
   
    speed = sqrt(Y0.^2 + Y1.^2) * (ball_diameter*pi) * fps; % put in cm/sec
    speed_av(i) = mean(speed);
    speed_dyn{i} = speed;

end

%% 

%-------------------
figure()
var=abs(time_ratio);

colour = 'r';
[t_timeratio, rotlist] = rotation_plotter (rot,ant,cur_cond,var,colour,'time bias');

title('Time bias: Route')  
xlim([-180 180])
xticks([-180:45:180])

%------------------- 
figure()
var=turn_ratio;

colour = 'r';
[t_turnratio, rotlist] =rotation_plotter (rot,ant,cur_cond,var,colour,'turn bias')
title('Turn bias: Route')  
xlim([-180 180])
xticks([-180:45:180])


%-------------------
figure()
var=turn_av * 180/pi; %to get in degree /sec (instead of rad/sec)

colour = 'r';
[t_turnav, rotlist] = rotation_plotter (rot,ant,cur_cond,var,colour,'turn deg/sec');

title('Turn average (deg/s): Route')  
xlim([-180 180])
xticks([-180:45:180])



%% ----- stats

%stats for turn ratio
P=[];
 [P(1,1),H,STATS] = signrank(t_turnratio(:,rotlist==+45),  t_turnratio(:,rotlist==-45)   ,'method','exact','tail','left')
 [P(2,1),H,STATS] = signrank(t_turnratio(:,rotlist==+90),  t_turnratio(:,rotlist==-90)   ,'method','exact','tail','left')
 [P(3,1),H,STATS] = signrank(t_turnratio(:,rotlist==+135),  t_turnratio(:,rotlist==-135)   ,'method','exact','tail','left')

P=[];
 [H,P(1,1),STATS] = ttest(t_turnratio(:,rotlist==+45),  t_turnratio(:,rotlist==-45)  ,'tail','left')
 [H,P(2,1),STATS] = ttest(t_turnratio(:,rotlist==+90),  t_turnratio(:,rotlist==-90)  ,'tail','left')
 [H,P(3,1),STATS] = ttest(t_turnratio(:,rotlist==+135),  t_turnratio(:,rotlist==-135)   ,'tail','left')

%% --------------simple stats
t_var = t_turnratio;
% t_var = t_turnav
[N,P,Z]=deal([]);

k=0;
for i= rotlist
k=k+1;
        data=t_var(:,k);
   
    if i>0 & i<180; % one tail

        [P(k,1),H,STATS] = signrank (data,0,'tail','left');
        N(k,1)=length(data);

    elseif i<0 % other tail

        [P(k,1),H,STATS] = signrank (data,0,'tail','right');
        N(k,1)=length(data);

    else % means 0 or 180; then both tails

        [P(k,1),H,STATS] = signrank (data,0,'tail','both');
        N(k,1)=length(data);

    end
end

t_all_rot = table(rotlist',N,P)

%% -------------- stats by overlaying negatif and positve

[N,P,Z]=deal([]);
rotlist = [180 -135 -90 -45 0 45 90 135]
rot_list_deg_stat = unique(abs(rotlist));

k=0;
for i= rot_list_deg_stat
k=k+1;
    
    % overlay 
    if i>0 & i<180
        fpos=find(rotlist == i);
        fneg=find(rotlist == -i);
        data=[t_var(:,fpos) ; -t_var(:,fneg)];

        [P(k,1),H,STATS] = signrank (data,0,'tail','left');
        Z(k,1)=STATS.zval;
        N(k,1)=length(data);

    else % means it 0 or 180
        f=find(rotlist == i);
        f=f(1); %to make sure it doesn't double 180
        data=t_var(:,f) ;

        [P(k,1),H,STATS] = signrank (data,0,'tail','both');
        Z(k,1)=0;
        N(k,1)=length(data);

    end
end

t_all_rot = table(rot_list_deg_stat',N,P)


%% Forward speed -------------------------------------
figure()

var=speed_av; % time spent turning L/R: between +1 and -1 

colour = 'r';
rotation_plotter (rot,ant,cur_cond,var,colour,'speed fwd');

title('speed fwd (cm/s): Route ')  
xlim([-180 180])
xticks([-180:45:180])


%% ---------------- turn vel prefered vs unprefered side
figure()


subplot(2,1,1);
colour = 'r';
var=turn_vel_prefered * 180/pi; % time spent turning L/R: between +1 and -1 
rotation_plotter (rot,ant,cur_cond,var,colour,'turn vel');
 
colour = 'y';
var=turn_vel_unprefered * 180/pi; % time spent turning L/R: between +1 and -1 
rotation_plotter (rot,ant,cur_cond,var,colour,'turn vel');

title('Route: turn vel pref vs. unprefered side (deg/sec)')
xlim([-180 180])
xticks([-180:45:180])
ylim([0 550])
 
   

subplot(2,1,2);
var=(turn_vel_prefered * 180/pi) - (turn_vel_unprefered * 180/pi); % time spent turning L/R: between +1 and -1 

cur_cond=strcmp(cond,'R'); %'U' or 'R'
colour = 'r';
rotation_plotter (rot,ant,cur_cond,var,colour,'turn vel diff');


title('turn vel pref-unpref: Route ')  
xlim([-180 180])
xticks([-180:45:180]) 




