
%% get bearings and speed from path

load('mirror_cata_paths.mat')

fillpath_step = 0.1;%fill the path to make bearings more accurate
bearing_dist = 12; %radius of the circle for cutting the path into sgmt
which_sector =12;

[Bearings, Bearings_idx, Turns, Turns_idx, Speed, Speed_idx] = deal([]);
for i=1:length(X)
   x=X{i,1};
   y=Y{i,1};
      
   
   if which_sector == 12
       % Remove third sector
       f=find(Sector{i,1}==3,1,'first');
       x=x(1:f);
       y=y(1:f);

       %Start from onset towards the end
       f=find(Sector{i,1}==2,1,'first');
       x_after=x(f:end);
       y_after=y(f:end);
   
   elseif which_sector == 23
       %Start from onset towards the end
       f=find(Sector{i,1}==3,1,'first');
       x_after=x(f:end);
       y_after=y(f:end);
       
       % Remove first sector
       f=find(Sector{i,1}==1,1,'last');
       x=x(f:end);
       y=y(f:end);

   else
       warning('which_sector must be 12(1-2) or 23(2-3)')
   end
   
      % Start from onset towards the begining
   % fliplr the path before, so that bearings start exactly at the path
   x_before= flipud(x(1:f));
   y_before=flipud(y(1:f));
    
%       % get nb_clicks per chunk (and thus speed)
%    [nb_clicks_after] = nb_click_per_chunk (x_after, y_after, bearing_dist, fillpath_step);
%    [nb_clicks_before_flipped] = nb_click_per_chunk (x_before, y_before, bearing_dist, fillpath_step);
%    nb_clicks_all = [flipud(nb_clicks_before_flipped)'; nb_clicks_after']' ;
%    
%    speed_all = 1./nb_clicks_all * fps * bearing_dist;
%    Speed{i,1}= speed_all;
%    Speed_idx{i,1}= [(-length(nb_clicks_before_flipped): -1) , 1 : length(nb_clicks_after)];

   
   % get bearings
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

 clear ('turn', 'nb_clicks_all', 'speed_all', 'nb_clicks_before_flipped', 'nb_clicks_after', 'bearing_all','bearings_after', 'bearings_before_flipped', 'f', 'i','j','x','y','x_after','y_after','xfilled','yfilled','x_before','y_before','fillpath_step')

%    save ('mirror_cata_bearings.mat')
 
 
%% ----------------------- plot turns and bearing history


%   load('mirror_cata_bearings.mat')

%how much before and after the change do we look at
window_cm = 60;
window = ceil(window_cm / bearing_dist); % how many segment before and after the mirror will be considered

nb_ant=length(X);
table_idx = -window:window;

Turn_table=NaN(nb_ant,length(table_idx));
Bearing_table=NaN(nb_ant,length(table_idx));


for i=1:nb_ant

      a=Turns{i,1};
      b=Bearings{i,1};
%     c=Speed{i,1};
    
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

    f_SM = strcmp(cond,'SMS');
    f_ctr = strcmp(cond,'CTR');

subplot(2,1,1);
    errorbar(table_idx*bearing_dist, nanmean(Turn_table(f_SM,:))*180/pi, nanstd(Turn_table(f_SM,:))*180/pi / sqrt(sum(f_SM)) ,'r','linewidth',2); hold on
    errorbar(table_idx*bearing_dist, nanmean(Turn_table(f_ctr,:))*180/pi, nanstd(Turn_table(f_ctr,:))*180/pi / sqrt(sum(f_ctr)) , 'K','linewidth',1); hold on
    plot([0 0] , [0 30],'--k');
    xlabel('relative to onset (cm)')
    ylabel('turns (deg)')
    title(['sgmt dist = ' , num2str(bearing_dist)]);
   
subplot(2,1,2);
    errorbar(table_idx*bearing_dist, nanmean(Bearing_table(f_SM,:))*180/pi, nanstd(Bearing_table(f_SM,:))*180/pi / sqrt(sum(f_SM)) ,'r','linewidth',2); hold on
    errorbar(table_idx*bearing_dist, nanmean(Bearing_table(f_ctr,:))*180/pi, nanstd(Bearing_table(f_ctr,:))*180/pi / sqrt(sum(f_ctr)) , 'K','linewidth',1); hold on
    plot([min(table_idx*bearing_dist) max(table_idx*bearing_dist)] , [0 0],'--k');
    plot([0 0] , [-20 20],'--k');
    xlabel('relative to onset (cm)')
    ylabel('Bearing (deg)')
    title(['sgmt dist = ' , num2str(bearing_dist)]);    
    

%% plot ant path

figure();

for i=1:nb_ant
%     n= nest(i); %to separate both nest
    n=1 ; % to pool the nest
    x=X{i} ;
    y=Y{i} ;
    sector=Sector{i} ;
    
    if strcmp(cond{i},'CTR')
        subplot(2,1,n); title('CTR');
        plot(x(sector==1),y(sector==1),'b'); hold on;
        plot(x(find(sector==1,1,'last'): find(sector==2,1,'last')),y(find(sector==1,1,'last'): find(sector==2,1,'last')),'r'); hold on;
        
    elseif strcmp(cond{i},'SMS')
        subplot(2,1,n+1); title('SMS');   
        plot(x(sector==1),y(sector==1),'b'); hold on;
        plot(x(find(sector==1,1,'last'): find(sector==2,1,'last')),y(find(sector==1,1,'last'): find(sector==2,1,'last')),'r'); hold on;

    end
    
end

