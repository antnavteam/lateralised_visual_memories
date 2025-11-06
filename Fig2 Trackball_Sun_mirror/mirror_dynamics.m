
load('miror_data_fordynamics.mat')

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
frame_window = 150 % 150 is good otherwise bucket or mirror off. This corresponds to the window before and the same is added after

turn_dyn ={};
speed_dyn ={};

mir_trial = 1

for i=1:length(data)
    
    % get the mirror index
    cur_mirror = mirror{i}; 
    mir_on=find( [cur_mirror(1:end-1)-cur_mirror(2:end)] < 0); 
    mir_off=find( [cur_mirror(1:end-1)-cur_mirror(2:end)] > 0); 

    if mir_trial ==1
        window = [mir_on(1)-frame_window : mir_on(1)+ frame_window]; % here look at the first mirroring event

        if max(window) > mir_off(1)
            warning([cond{i} , num2str(ant(i)), 'window to big, the end span over  mirror off'])
        end
        
    elseif mir_trial ==2
        
        window = [mir_on(2)-frame_window : mir_on(2)+ frame_window]; % here look at the first mirroring event

        if min(window) < mir_off(1)
            warning([cond{i} , num2str(ant(i)), 'window to big, the begining span over previous mirror on'])
        end
    end
        
    
    cur_data=data{i};
    
    X0=cur_data(window ,1)/calib_0; % get the last Xsec (we made sure the last 15s were without cleaning)
    X1=cur_data(window ,3)/calib_1;    
    Xs= (X0 + X1) / 2; % Make average of X0 and X1 (they should get same signal)
    
    Y0=cur_data(window,2) /calib_0;
    Y1=cur_data(window,4) /calib_1;
    
    % Path
    cum_angle{i} = pi2pi(cumsum(Xs) * 2* pi); % in radian (per frame)
    
    % Turns 

    Xs_radsec= Xs * 2*pi * fps;  % Put them in rad/sec
    turn_dyn{i} = Xs_radsec;
    
        

    speed = sqrt(Y0.^2 + Y1.^2) * (ball_diameter*pi) * fps; % put in cm/sec
    speed_dyn{i} = speed;

    
    % get averages before

    Xs_before = Xs_radsec (1:frame_window);
    turn_l_before(i)= sum(Xs_before(Xs_before>0));
    turn_r_before(i)= sum(Xs_before(Xs_before<0));
    turn_ratio_before(i) = ( abs(turn_l_before(i)) - abs(turn_r_before(i)) ) / ( abs(turn_l_before(i)) + abs(turn_r_before(i)) ); % time spent turning L/R: between +1 and -1 
    turn_av_before(i)= mean(Xs_before);
    turn_med_before(i)= median(Xs_before);

    
        % get averages after

    Xs_after = Xs_radsec (frame_window:end);
    turn_l_after(i)= sum(Xs_after(Xs_after>0));
    turn_r_after(i)= sum(Xs_after(Xs_after<0));
    turn_ratio_after(i) = ( abs(turn_l_before(i)) - abs(turn_r_after(i)) ) / ( abs(turn_l_after(i)) + abs(turn_r_after(i)) ); % time spent turning L/R: between +1 and -1 
    turn_av_after(i)= mean(Xs_after);
    turn_med_after(i)= median(Xs_after);

    
    
end


%% Plot Turn dynamic 
% figure()

nb_frames = length(turn_dyn{1});
x_axis = [[1:nb_frames] - ceil(nb_frames/2)]   /fps; % so that zero is when the ring is lifted 

fcond={};
fcond{1}='Trained_side_Realease_route';
fcond{2}='Trained_side_Realease_unfam';
fcond{3}='Trained_straight_Realease_route';

dth=[];
dr=[];
dth_return_up=[];

figure()

for c=1:length(fcond)
    
    cur_cond=strcmp(cond,fcond(c)); % get the condition

         j=0; 
         for a=unique(ant(cur_cond))'; j=j+1; % go through all individual
             f=find(ant == a); %get the index
                dth(j,:) = turn_dyn{f}' *180/pi; % in deg/sec
                dr(j,:) = speed_dyn{f}' / fps; % so as to be back in cm/frame
                
                % return the ant to align them base on their average turn before the mirror
             if turn_av_before(a) < 0  %return only if turnbefore_below zero
                dth_return_up(j,:) = -turn_dyn{f}' *180/pi; % in deg/sec
             else
                dth_return_up(j,:) = turn_dyn{f}' *180/pi; % in deg/sec
             end
                
                
         end
         

        
        % Plot raw data
        subplot(2,3,c); 
                plot(x_axis, dth'); hold on

%                error = fastsmooth (std(dth)/sqrt(size(dth,1)) , 3,3,1)  ;
%        patch_vector =  [mean(dth)+error, flip(mean(dth)-error)];

        error = fastsmooth (iqr(dth)/2, 3,3,1)  ; %/2 because we patch half above and half below
        patch_vector =  [median(dth)+error, flip(median(dth)-error)];

        patch([x_axis flip(x_axis)],patch_vector,[0.5 0.5 0.5]); hold on
        
        plot(x_axis, fastsmooth(median(dth),3,3,1),'k','linewidth',2); hold on
        
        plot([min(x_axis) max(x_axis)], [0 0], 'k'); %mark the zero
        plot([0 0] , [min(min(dth)) max(max(dth))],'k'); % mark t0, when the ring is lifted
        ylim([-500 +500])
        title(fcond(c));
        ylabel('turn vel (deg/s)');
        xlabel('time (s)');
        
        % Plot returned data
        subplot(2,3,c+3); 
        plot(x_axis, dth_return_up'); hold on

%        error = fastsmooth (std(dth_return_up)/sqrt(size(dth_return_up,1)) , 3,3,1)  ;
%        patch_vector =  [mean(dth_return_up)+error, flip(mean(dth_return_up)-error)];

        error = fastsmooth (iqr(dth_return_up)/2, 3,3,1)  ; %/2 because we patch half above and half below
        patch_vector =  [median(dth_return_up)+error, flip(median(dth_return_up)-error)];
        
        patch([x_axis flip(x_axis)],patch_vector,[0.5 0.5 0.5]); hold on
        
        plot(x_axis, fastsmooth(median(dth_return_up),3,3,1),'k','linewidth',2); hold on
        
        plot([min(x_axis) max(x_axis)], [0 0], 'k'); %mark the zero
        plot([0 0] , [min(min(dth_return_up)) max(max(dth_return_up))],'k'); % mark t0, when the ring is lifted
        ylim([-500 +500])
        title(['Returned up' , fcond(c)]);
        ylabel('turn vel (deg/s)');
        xlabel('time (s)');
end

     
%% Plot the means 
 
fcond={};
fcond{1}='Trained_side_Realease_route';
fcond{2}='Trained_side_Realease_unfam';
fcond{3}='Trained_straight_Realease_route';

 before = turn_av_before /pi*180 ; %put in deg/sec
 after = turn_av_after /pi*180 ;
%  before = turn_med_before /pi*180 ; %put in deg/sec
%  after = turn_med_after /pi*180 ;
  
ymin =-250;
ymax = 250;

figure()

for c=1:length(fcond)
    
    subplot(1,3,c)
    
    f=strcmp(cond,fcond(c));
    boxplot([before(f)',after(f)']); hold on;
    plotSpread ([before(f)',after(f)']); hold on;
    plot([before(f);after(f)],'color',[0.5 0.5 0.5]); hold on;
    plot([0.5 2.5],[0 0],'--k')
    ylim([ymin ymax])
    title(fcond{c})
    
    
end


% Stats 


%----- First test if group are oriented to the right before the mirror
    
    f=find(strcmp(cond,'Trained_side_Realease_route') ==1)
    [p,h,stats] = signrank(before(f),0,'tail','left')
%     signedrank: 21
%     p = 0.0156 

    f=find(strcmp(cond,'Trained_side_Realease_unfam') ==1)
    [p,h,stats] = signrank(before(f),0,'tail','left')
%     signedrank: 64
%     p = 0.9788


    f=find(strcmp(cond,'Trained_straight_Realease_route') ==1)
    [p,h,stats] = signrank(before(f),0,'tail','left')
%     signedrank: 66
%     p = 0.9866

%----------- Then test if different before after, makes sense only for
%group 1, because it is oriented

    f=find(strcmp(cond,'Trained_side_Realease_route') ==1)
    [p,h,stats] = signrank(before(f),after(f),'tail','left')
%     signedrank: 21
%     p = 0.0156 


%----- Then test if mirror reverses them
    
    f=find(strcmp(cond,'Trained_side_Realease_route') ==1)
    [p,h,stats] = signrank(after(f).*before(f),0,'tail','left')
%     signedrank: 2
%     p = 0.0469 

    f=find(strcmp(cond,'Trained_side_Realease_unfam') ==1)
    [p,h,stats] = signrank(after(f).*before(f),0,'tail','left')
%     signedrank: 35
%     p = 0.3955

    f=find(strcmp(cond,'Trained_straight_Realease_route') ==1)
        [p,h,stats] = signrank(after(f).*before(f),0,'tail','left')
%     signedrank: 15
%     p =  0.0320

%% make little csv table of the mean

t=table(ant, cond, before', after', 'VariableNames', {'ant','cond','before','mirror'});

writetable(t,'table_mirror_means.csv')
