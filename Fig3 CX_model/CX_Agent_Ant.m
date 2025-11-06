
function [x,y,th,turn,err_dir,straightness] =CX_Agent_Ant(time, mem_bias, decay_CPU4, gain, motor_noise_deg, mem_noise_deg, pontine)

%(x,y) are the agent coordinate
%th is the current heading direction;

% Basic parameters
step_length = 1;
max_angle = 180/180*pi; %boundary condition for the oscillation Careful below of the sign if not 180 and zero!
min_angle = 0/180*pi;

PB_thref = 0:pi/4:2*pi-0.1;
bump_pwr = 4; % the higher the more the winner take all, the smaller the more linear
% PB=((pi-abs(pi2pi(PB_thref - th))).^pwr) / pi.^pwr; % example

mem_L= mem_bias/180*pi;
mem_R= -mem_bias/180*pi;

mem_noise=mem_noise_deg/180*pi;
motor_noise=motor_noise_deg/180*pi;

% Allocation variable cumulative
turn=nan([1,time]);
th=nan([1,time+1]);
x=nan([1,time+1]);
y=nan([1,time+1]);

% Initial conditions
x(1)=0; 
y(1)=0; 
th(1)=pi22pi(rand*2*pi); % th is the actual agent bearing
CPU4_L = zeros([1,8]);
CPU4_R = zeros([1,8]);

for i=1:time %i refer to time before movement. 
    
   % Get PB activation (0-1)
      PB_raw =((pi-abs(pi2pi(PB_thref -  th(i)))).^bump_pwr);
      PB = (PB_raw / max(PB_raw) );
    
 %--------------- Update CPU4 memory -------------

     % Add noise to the direction of the visual memory
    mem_Ln = mem_L + normrnd(0,mem_noise);
    mem_Rn = mem_R + normrnd(0,mem_noise);
     % Get View Familiarities (0-1)
    Fam_L = ((pi - abs(pi2pi(mem_Ln - th(i)))) /pi)^2;
    Fam_R = ((pi - abs(pi2pi(mem_Rn - th(i)))) /pi)^2;
%   fam_R =0; % Covered eye = always max unfamiliair (pi)
    

    
     % Get CPU4 input (0 to 1) = Fam inhibited by circhsift PB
    CPU4_L_input = Fam_L - circshift(PB,+1);
    CPU4_R_input = Fam_R - circshift(PB,-1);
    
    CPU4_L_input(CPU4_L_input < 0) = 0;
    CPU4_R_input(CPU4_R_input < 0) = 0;
    
    % CPU4 input added to previous CPU4 decayed.
    CPU4_L = CPU4_L_input + (CPU4_L - CPU4_L*decay_CPU4);
    CPU4_R = CPU4_R_input + (CPU4_R - CPU4_R*decay_CPU4);
    
%--------------- Steering -------------
 
    % Get CPU1 activation
    if pontine == 0
        CPU1_L = CPU4_L - PB ;
        CPU1_R = CPU4_R - PB ;
    elseif pontine == 1
        pontine_L = CPU4_L;
        pontine_R = CPU4_R;
        CPU1_L = CPU4_L - PB - circshift(pontine_L,4); 
        CPU1_R = CPU4_R - PB - circshift(pontine_R,4);
    else
        warning('pontine must be 1(on) or 0(off)')
    end
    
        CPU1_L(CPU1_L<0) = 0;
        CPU1_R(CPU1_R<0) = 0;


    
    % Get the turning rate (mean instead of sum to stay between 0-1)
    turn_r = sum(CPU1_L); % Note that a good left match (high CPU1_L) means turn right
    turn_l = sum(CPU1_R);
    
    % Will need to change when the oscillator comes into play
    turn(i) = turn_l - turn_r;
    
    angle_modul = turn(i)*gain ;
        % boundary
        if abs(angle_modul) < min_angle; angle_modul = min_angle * sign(angle_modul); end   
        if abs(angle_modul) > max_angle; angle_modul = max_angle * sign(angle_modul); end 
  
        
        
        
    th(i+1)= th(i)+ angle_modul; %set the new heading

    %add some motor noise (sensory noise could be also added)
    th(i+1) = th(i+1) + normrnd(0,motor_noise);
    
%     % calculate the stepsize
    stepsize = step_length;
%     stepsize = (step_length* (1-fwd_modul)) / nb_step;
%     if stepsize<0 ; stepsize=0;
%     end
    
    %--------move to time i+1------------   
     
        [X, Y]=pol2cart(th(i+1) , stepsize); %move along this direction
        x(i+1)= x(i)+X;
        y(i+1)= y(i)+Y;
    
end

[th_end, r_end]=cart2pol(x(end),y(end));
err_dir=abs(th_end);
straightness = r_end/time;
