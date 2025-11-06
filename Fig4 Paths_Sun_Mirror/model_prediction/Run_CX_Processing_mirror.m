%% The Ant Agent mirrored sun
% Note that this fonction will run agent simulations, which implies noise,
% Therefore, the average graphs resulting from these new paths may differ
% slightly across successive runs.


time=120;
pontine =1;

gain=1;
motor_noise=10;
decay_CPU4 = 0.2; %how much ratio CPU4 loses per step
mem_bias=45;
mem_noise=10;

sun_rotation_time=60;% first 90deg then 180deg
sun_rotation=120;


nb_ant=20;
[X, Y, th, turn] = deal([]);
for i=1:nb_ant
    [X(i,:),Y(i,:),th(i,:),turn(i,:),a(i),b(i)] =CX_Agent_Ant_mirrorsun(time, mem_bias, decay_CPU4, gain, motor_noise, mem_noise, pontine, sun_rotation_time, sun_rotation);
    X(i,:) =  X(i,:) - X(i,sun_rotation_time); %aligne them for mirror at 0
end

figure();
plot(X(:,1:sun_rotation_time)',Y(:,1:sun_rotation_time)','b'); hold on;
plot(X(:,sun_rotation_time:end)',Y(:,sun_rotation_time:end)','r'); hold on;

scatter(0,0);
axis equal
text((0),min(min(Y)),{['gain=' num2str(gain)];...
['motor noise=' num2str(motor_noise)];...
['decay CPU4=' num2str(decay_CPU4)] ;...
['mem bias=' num2str(mem_bias)];...
['mem noise=' num2str(mem_noise)]})

% figure(); plot(abs(turn)'*gain);
% figure(); plot(mean(abs(turn)*gain));




% % ------------make the same analysis as real ant

bearing_dist =1; % radial dist use to segment the path, on the real ant figure the ratio is 5 (bearing_dist=12; window=60)
window_cm = 5; %how much step before and after mirror event
xlim ([-window_cm*4 window_cm*4])

fillpath_step =0.1;


    
 for i=1:nb_ant
   x=X(i,:)';
   y=Y(i,:)';
   
   x_after=x(sun_rotation_time:end);
   y_after=y(sun_rotation_time:end);

   x_before=flipud(x(1:sun_rotation_time));
   y_before=flipud(y(1:sun_rotation_time));
     
   [xfilled, yfilled] = fillpath (x_before, y_before, fillpath_step);
   [~,~, bearings_before_flipped] = sinuosity (xfilled, yfilled, bearing_dist);

   [xfilled, yfilled] = fillpath (x_after, y_after, fillpath_step);
   [~,~, bearings_after] = sinuosity (xfilled, yfilled, bearing_dist);
   
   bearing_all= [flipud(pi2pi(bearings_before_flipped - pi)); bearings_after]' ;
   Bearings{i,1}= bearing_all;
   Bearings_idx{i,1}= [(-length(bearings_before_flipped): -1) , 1 : length(bearings_after)];
 
 
   %get turn angle for successive bearings
    turn=[];
    for j = 1 : length (bearing_all) - 1
        turn(j) = abs(pi2pi(bearing_all(j) - bearing_all (j+1)));
    end
   Turns{i,1} = turn;
   %Get the sequence, 0 correspond to the turn angle betwen the segment just before and
   % just after the mirror
   Turns_idx{i,1} = (-length(bearings_before_flipped)+1) : (length(bearings_after)-1);

 end




% window_cm = 30;
window = ceil(window_cm / bearing_dist); % how many segment before and after the mirror will be considered
table_idx = -window:window;

Turn_table=NaN(nb_ant,length(table_idx));
Bearing_table=NaN(nb_ant,length(table_idx));


for i=1:nb_ant

      a=Turns{i,1};
      b=Bearings{i,1};
    
    % fill turn table
    k=0;%nb of column (time along the table)
    for j=table_idx; k=k+1;
        f=find(Turns_idx{i,1} == j);
        if isempty(f)== 0 % if the path is long enough, otherwise it leaves it as a NaN
            Turn_table(i,k) = a(f);
        end
        
         f=find(Bearings_idx{i,1} == j);
        if isempty(f)== 0 % if the path is long enough, otherwise it leaves it as a NaN
            Bearing_table(i,k) = b(f);
        end                   

        
    end

end
        

% % -------------- plot group figures ------


figure()

subplot(2,1,1);
    errorbar(table_idx*bearing_dist, nanmean(Turn_table)*180/pi, nanstd(Turn_table)*180/pi / sqrt(nb_ant) ,'r','linewidth',2); hold on
    plot([0 0] , [0 30],'--k');
    xlabel('relative to onset (step)')
    ylabel('turns (deg)')
    title(['sgmt dist = ' , num2str(bearing_dist)]);
   
subplot(2,1,2);
    errorbar(table_idx*bearing_dist, nanmean(Bearing_table)*180/pi, nanstd(Bearing_table)*180/pi / sqrt(nb_ant) ,'r','linewidth',2); hold on
    plot([min(table_idx*bearing_dist) max(table_idx*bearing_dist)] , [0 0],'--k');
    plot([0 0] , [-20 20],'--k');
    xlabel('relative to onset (step)')
    ylabel('Bearing (deg)')
    title(['sgmt dist = ' , num2str(bearing_dist)]);    



%% Draw 

time=40;
pontine =1

gain=1;
motor_noise=10;
decay_CPU4 = 0.2; %how much ratio CPU4 loses per step
mem_bias=45;
mem_noise=10;

sun_rotation_time=20;% first 90deg then 180deg
sun_rotation=135;


nb_ant=1;
[X, Y, th, turn] = deal([]);

    [X,Y,th,turn,a,b,CPU4_L,CPU4_R,PB] =CX_Agent_Ant_mirrorsun(time, mem_bias, decay_CPU4, gain, motor_noise, mem_noise, pontine, sun_rotation_time, sun_rotation);


figure();
plot(X(1:sun_rotation_time)',Y(1:sun_rotation_time)','b'); hold on;
plot(X(sun_rotation_time:end)',Y(sun_rotation_time:end)','r'); hold on;

axis equal
text((0),min(min(Y)),{['gain=' num2str(gain)];...
['motor noise=' num2str(motor_noise)];...
['decay CPU4=' num2str(decay_CPU4)] ;...
['mem bias=' num2str(mem_bias)];...
['mem noise=' num2str(mem_noise)]})


%plot the CPU state

thref = 0:pi/4:2*pi-0.1;

fix_list = [5 10 19 20 21 22 23 25 30 35 40];
fix_nb = length(fix_list)
scatter(X(fix_list),Y(fix_list),'k');

figure()
k=0;
for i=fix_list
    k=k+1;

    subplot(3,fix_nb,k);
%     bar(PB(i,:),'k');
%     polarplot(thref,PB(i,:), 'k-o'); hold on;
      poloarplot2rose(thref,PB(i,:));
    title(num2str(i - sun_rotation_time));

    subplot(3,fix_nb,k+fix_nb);
%     bar(CPU4_L(i,:),'r');
%         polarplot(thref,CPU4_L(i,:), 'r-o'); hold on;
        poloarplot2rose(thref,CPU4_L(i,:));

    subplot(3,fix_nb,k+fix_nb*2);
%     bar(CPU4_R(i,:),'b');
%         polarplot(thref,CPU4_R(i,:), 'b-o'); hold on;
        poloarplot2rose(thref,CPU4_R(i,:));

end