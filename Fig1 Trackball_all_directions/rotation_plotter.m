function [return_data, return_rotlist] = rotation_plotter (rot,ant,cur_cond,var,colour,label)



% rot_list = unique(rot);
rot_list = [max(rot) ; unique(rot)]; %dupplicate 180 by adding it at the begining
     j=0; 
     for a=unique(ant(cur_cond))'; j=j+1; % go through all individual
         k=0;   
         for r = rot_list'; k=k+1; % go through all rotation

                f=find(cur_cond & ant == a & rot == r); %get the index
%                 t_turn_av(a,k) = turn_av(f);
%                 t_turn_l(a,k) = turn_l(f);
%                 t_turn_r(a,k) = turn_r(f);
%                 
%                 t_time_ratio(a,k) = (time_l(f)- time_r(f))/ (time_l(f)+time_r(f)); % time spent turning L/R: between +1 and -1 
                  t_var(j,k) = var(f);

         end
     end
 
% --------- plot it
     
rot_list = [-max(rot) ; unique(rot)];
x = rot_list/pi*180;
ym = nanmean(t_var);
ybar = nanstd(t_var)/sqrt(size(t_var,1));

errorbar(x , ym , ybar , 'color', colour , 'Linewidth', 2 ); hold on       
plot(x , t_var', 'color', colour ); hold on  

%plot rep�re:
ymax=max(max(t_var));
ymin=min(min(t_var));
if ymin>0; ymin=0;end

plot([0 0]',[ymin ymax]','k'); hold on
plot([x(1) x(end)]',[0 0]','k'); hold on

ylabel(label)
xlabel('angle away from goal (deg)') 

return_data = t_var(:,1:end-1);
return_rotlist = x(1:end-1)';


