function MPC_out = MPC_alfa_tracking( dri_trq_act, alfa_desire, tot_trq_desire, oldu1)

persistent controller N aa

% if t == 0
if isempty(aa)
    aa = 1;
     
    % optimization problem parameters
    N = 6;%prediction horizon
    T_s = 0.01;%sampling time
    nx = 1;% number of state varibles
    nu = 1;% number of inputs
        
    t_alfa = 1*(1/0.01)^2;
    t_u1 = 1*(1/20)^2;
    t_dtrq = 35*(1/1/T_s)^2*1.4;
%     t_dtrq = 10*(1/1/T_s)^2*1.4;
    t_dalfa = 0*(1/0.01)^2;
    tao_driver = .5;
    
    K = 7;
    
    
    
    yalmip('clear');
    % optimization varibles define
    uvar = sdpvar(repmat(nu,1,N+1),repmat(1,1,N+1));
    xvar = sdpvar(repmat(nx,1,N+1),repmat(1,1,N+1));
        
    alfa_tgt = sdpvar(N+1,1);
    
    tot_trq_tgt = sdpvar(N+1,1);
    
    pastu1 = sdpvar(1);
        
    s_alfa = sdpvar(repmat(1,1,N),repmat(1,1,N));
       
    % problem formulation
    constraints=[];
    objective=0;
    
    for k=1:N
               
%         s_alfa{k} = xvar{k}(1) - alfa_tgt(k)*tot_trq_tgt(1);
        
        s_alfa{k} = xvar{k}(1)/tot_trq_tgt(1) - alfa_tgt(k);
        
        objective = objective + t_alfa*s_alfa{k}^2  + t_u1*uvar{k}(1)^2 ; %  + t_u2*uvar{k}(2)^2
        
        constraints = [constraints, xvar{k+1}(1) == (1 - T_s/tao_driver)*xvar{k}(1) + K*T_s/tao_driver*uvar{k}(1)];
                
        constraints = [constraints,-30 <= uvar{k}(1) <= 30];
        
        constraints = [constraints,-30*T_s <= uvar{1}(1) - pastu1 <= 30*T_s];
        
%         constraints = [constraints,-1*T_s <= uvar{k+1}(1) - uvar{k}(1) <= 1*T_s];
                
    end
    
     for j = 2:N
        
        objective = objective + t_dtrq*(uvar{j} - uvar{j-1})^2; %
        objective = objective + t_dalfa*(s_alfa{k} - s_alfa{k-1})^2;
    
    end
    
    objective = objective + t_dtrq*(uvar{1} - pastu1)^2;
    
    parameters_in = {xvar{1},alfa_tgt,tot_trq_tgt,pastu1};%,pastu2
    
    solutions_out = uvar{1};
    
    %     assign(uvar{1},delta_ego);
    
    controller = optimizer(constraints,objective,sdpsettings('solver','gurobi','usex0',1),parameters_in,solutions_out);
    
    [MPC_out,diagnostic,~,~,controller] = controller{{ dri_trq_act,alfa_desire,tot_trq_desire,oldu1}};%,ax_ego
    
    %     MPC_out = controller{{[X_ego;Y_ego;vy_ego;fai_ego;dfai_ego],X_desire*ones(N+1,1),Y_desire*ones(N+1,1),fai_desire*ones(N+1,1),vx_desire*ones(N+1,1),X_obs*ones(N+1,1),Y_obs*ones(N+1,1),delta_ego,A_R,B_R,vy_desire*ones(N+1,1),delta_Y_act}};
else
    %     assign(uvar{1},delta_ego);
    [MPC_out,diagnostic,~,~,controller] = controller{{ dri_trq_act,alfa_desire,tot_trq_desire,oldu1}};%,ax_ego,vx_ego
    
    %     MPC_out = controller{{[X_ego;Y_ego;vy_ego;fai_ego;dfai_ego],X_desire*ones(N+1,1),Y_desire*ones(N+1,1),fai_desire*ones(N+1,1),vx_desire*ones(N+1,1),X_obs*ones(N+1,1),Y_obs*ones(N+1,1),delta_ego,A_R,B_R,vy_desire*ones(N+1,1),delta_Y_act}};
end
end