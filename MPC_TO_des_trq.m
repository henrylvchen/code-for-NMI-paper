function MPC_out = MPC_TO_des_trq( Y_ego, vy_ego, fai_ego, dfai_ego, Y_desire, fai_desire, oldu1, delta_ego, ddelta_ego, sin_fai, cos_fai)

persistent controller N aa

% if t == 0
if isempty(aa)
    aa = 1;
    
    % optimization problem parameters
    N = 6;%prediction horizon
    T_s = 0.01;%sampling time
    nx = 8;% number of state varibles
    nu = 1;% number of inputs
    
    vx = 13;
    m = 1830;
    lf = 1.152;
    lr = 1.693;
    Iz = 3477;
    K_f = 40703*2;
    K_r = 64495*2;
    
    K_sw = 2.29*6;
    B_sw = 1.56/2;
    J_sw = 0.172;
    is = 14.5;
    K_load = 30;
    K_sw_load = K_load*0.4;
    
    t_Y = 50*(1/0.1)^2;
    t_fai = 20*(1/0.01)^2;
    t_trq = .1*(1/20)^2;
    t_dtrq = .00001*(1/8/T_s)^2;
        
    yalmip('clear');
    % optimization varibles define
    uvar = sdpvar(repmat(nu,1,N+1),repmat(1,1,N+1));
    xvar = sdpvar(repmat(nx,1,N+1),repmat(1,1,N+1));
    
    
    Y_tgt = sdpvar(N+1,1);
    
    yaw_angle_tgt = sdpvar(N+1,1);
    
    pastu1 = sdpvar(1);
    
    s_Y = sdpvar(repmat(1,1,N),repmat(1,1,N));
    
    s_fai = sdpvar(repmat(1,1,N),repmat(1,1,N));
    
    alfa_f = sdpvar(repmat(1,1,N),repmat(1,1,N));
    
    alfa_r = sdpvar(repmat(1,1,N),repmat(1,1,N));
    
    T_load = sdpvar(repmat(1,1,N),repmat(1,1,N));
    
    % problem formulation
    constraints=[];
    objective=0;
    
    for k=1:N
        
        s_Y{k} = xvar{k}(3) - Y_tgt(k);
        
        s_fai{k} = xvar{k}(8) - yaw_angle_tgt(k);
        
        alfa_f{k} = (xvar{k}(1) + lf*xvar{k}(2))/vx - xvar{k}(6)/is;
        
        alfa_r{k} = (xvar{k}(1) - lr*xvar{k}(2))/vx;
        
        T_load{k} = -K_sw_load*alfa_f{k}/is*180/3.14;
        
        objective = objective + t_Y*s_Y{k}^2  + t_fai*s_fai{k}^2  + t_trq*uvar{k}(1)^2 ; % 
        
        constraints = [constraints, xvar{k+1}(1) == xvar{k}(1) + T_s*(-vx*xvar{k}(2) + 1/m*(-K_f*alfa_f{k} + (-K_r)*alfa_r{k}))];
        
        constraints = [constraints, xvar{k+1}(2) == xvar{k}(2) + T_s/Iz*(lf*(-K_f)*alfa_f{k} - lr*(-K_r)*alfa_r{k})];
        
        constraints = [constraints, xvar{k+1}(3) == xvar{k}(3) + T_s*(vx*xvar{k}(4) + xvar{k}(1)*xvar{k}(5))];
        
        constraints = [constraints, xvar{k+1}(4) == xvar{k}(4) + T_s*xvar{k}(5)*xvar{k}(2)];
        
        constraints = [constraints, xvar{k+1}(5) == xvar{k}(5) - T_s*xvar{k}(4)*xvar{k}(2)];
        
        constraints = [constraints, xvar{k+1}(6) == xvar{k}(6) + T_s*xvar{k}(7)];
        
        constraints = [constraints, xvar{k+1}(7) == xvar{k}(7) + T_s*(-B_sw/J_sw*xvar{k}(7) - K_sw/J_sw*xvar{k}(6) + 1/J_sw*(uvar{k} - T_load{k}))];
        
        constraints = [constraints, xvar{k+1}(8) == xvar{k}(8) + T_s*xvar{k}(2)];
        
        
        
        constraints = [constraints,-30 <= uvar{k} <= 30];
        
        constraints = [constraints,-700*T_s <= uvar{k+1}(1) - uvar{k}(1) <= 700*T_s];
        
        constraints = [constraints,-1500*T_s <= uvar{1}(1) - pastu1 <= 1500*T_s];
        
    end
    
    for j = 2:N
        
        objective = objective + t_dtrq*(uvar{j} - uvar{j-1})^2; %
    
    end
    
    objective = objective + t_dtrq*(uvar{1} - pastu1)^2;
    
    parameters_in = {xvar{1},Y_tgt,yaw_angle_tgt,pastu1};%,pastu2
    
    solutions_out = [uvar{1};uvar{2};uvar{3};uvar{4};uvar{5};uvar{6};uvar{7}];
%     solutions_out = uvar{1};
    %     assign(uvar{1},delta_ego);
    
    controller = optimizer(constraints,objective,sdpsettings('solver','ipopt','usex0',1),parameters_in,solutions_out);
    
    [MPC_out,diagnostic,~,~,controller] = controller{{[vy_ego; dfai_ego; Y_ego; sin_fai; cos_fai; delta_ego; ddelta_ego; fai_ego],Y_desire,fai_desire,oldu1}};%,ax_ego
    
    %     MPC_out = controller{{[X_ego;Y_ego;vy_ego;fai_ego;dfai_ego],X_desire*ones(N+1,1),Y_desire*ones(N+1,1),fai_desire*ones(N+1,1),vx_desire*ones(N+1,1),X_obs*ones(N+1,1),Y_obs*ones(N+1,1),delta_ego,A_R,B_R,vy_desire*ones(N+1,1),delta_Y_act}};
else
    %     assign(uvar{1},delta_ego);
    [MPC_out,diagnostic,~,~,controller] = controller{{[vy_ego; dfai_ego; Y_ego; sin_fai; cos_fai; delta_ego; ddelta_ego; fai_ego],Y_desire,fai_desire,oldu1}};%,ax_ego,vx_ego
    
    %     MPC_out = controller{{[X_ego;Y_ego;vy_ego;fai_ego;dfai_ego],X_desire*ones(N+1,1),Y_desire*ones(N+1,1),fai_desire*ones(N+1,1),vx_desire*ones(N+1,1),X_obs*ones(N+1,1),Y_obs*ones(N+1,1),delta_ego,A_R,B_R,vy_desire*ones(N+1,1),delta_Y_act}};
end
end