function y=slipFunction(x,u,tdPos)
%%slipFunction与之前版本做了相应改变（最后一行，求半个周期最后的状态）
 
    g=9.8;    m=80; L0=1;
    k=30000;
    
    N = 1000;      % Number of Samples
    tSim = 2;     % Max Sim time
        
    tdPosL = [0 0.15 0];
    tdPosR = [0 0 0];
    tdParams = [pi/6 0];
    heightThreshold = 0;
    A=eye(6,6);
    A(2,2) = -1;
    A(5,5) = -1;
    
%     u=[0.2796356 0.3995468 18474.621];
    function [value isterminal direction] = leftTouchDownEvents(t,state)
        value = state(3)-heightThreshold;
        isterminal = 1;	direction = -1;
    end

    function [value isterminal direction] = RightTouchDownEvents(t,state)
        value = state(3)-heightThreshold;
        isterminal = 1;	direction = -1;
    end

    function [value isterminal direction] = leftLiftOffEvents(t,state)
        relPos = state(1:3)' - tdPosL;
        len = norm(relPos);
        value = len-L0;
        isterminal = 1;	direction = +1;
    end

    function [value isterminal direction] = rightLiftOffEvents(t,state)
        relPos = state(1:3)' - tdPosR;
        len = norm(relPos);
        value = len-L0;
        isterminal = 1;	direction = +1;
    end

    function [value isterminal direction] = LowestHighEvent(t,state)
        value = state(6);
        isterminal = 1;	direction = +1;
    end

    function [value isterminal direction] = msStanceEvent(t,state)
        value = state(6);
        isterminal = 1;	direction = -1;
    end

    function dState  = leftStanceDynamics(t,state)
        dState = 0*state;
        dState(1:3) = state(4:6);
        relPos = state(1:3)' - tdPosL;
        len = norm(relPos);
        relPos = relPos/len;
        force = (L0-len)*k;
        dState(4:6) = relPos*force/m;
        dState(6) = dState(6) - g;
        %state
        %dState
    end
 
    function dState  = rightStanceDynamics(t,state)
        dState = 0*state;
        dState(1:3) = state(4:6);
        relPos = state(1:3)' - tdPosR;
        len = norm(relPos);
        relPos = relPos/len;
        force = (L0-len)*k;
        dState(4:6) = relPos*force/m;
        dState(6) = dState(6) - g;
        %state
        %dState
    end
 
    function dState  = twoStanceDynamics(t,state)
        dState = 0*state;
        dState(1:3) = state(4:6);
        relPosL = state(1:3)' - tdPosL;
        relPosR = state(1:3)' - tdPosR;
        lenL = norm(relPosL);
        lenR = norm(relPosR);
        relPosL = relPosL/lenL;
        relPosR = relPosR/lenR;
        forceL = (L0-lenL)*k;
        forceR = (L0-lenR)*k;
        dState(4:6) = (relPosL*forceL+relPosR*forceR)/m;
        dState(6) = dState(6) - g;
        %state
        %dState
    end
 
    function [T_out STATE_out FOOT_out tf statef EN_out lowestState tdPosRNow] = halfSimulatePeriod(t0,stateInit)
    tspan = linspace(t0,t0+tSim,N);
    state0 = stateInit;
    rightTouchDownOptions = odeset('Events',@RightTouchDownEvents,'AbsTol',1e-9,'RelTol',1e-9);
    [T STATE tf statef] = ode45(@leftStanceDynamics,tspan,state0,rightTouchDownOptions);
    
    foot=0;
    en=0;
     
    STATE_out = STATE;
    T_out = T;
    FOOT_out = foot;
    EN_out = en;
    lowestState = zeros(6,1);
    
    if length(tf) == 0
            'break after 1'
            return;
    end
    
    tdPosR = [STATE(end,1)+cos(tdParams(2))*sin(tdParams(1)) STATE(end,2)-sin(tdParams(2))*sin(tdParams(1)) 0];
    tdPosRNow=tdPosR;
    
    tspan = linspace(tf,tf+tSim,N);
    LowestHighOptions = odeset('Events',@LowestHighEvent,'AbsTol',1e-9,'RelTol',1e-9);
    [T STATE tf statef] = ode45(@twoStanceDynamics,tspan,statef,LowestHighOptions);
    
    foot=0;
    en=0;
    
    if length(tf) == 0
            'break after 2'
            return;
    end
    
    lowestState = statef;
    T_out = [T_out ; T];
    STATE_out = [STATE_out ; STATE];
    FOOT_out = [ FOOT_out ; foot];
    EN_out = [EN_out ; en];
       
    tspan = linspace(tf,tf+tSim,N);
    leftLiftOffOptions = odeset('Events',@leftLiftOffEvents,'AbsTol',1e-9,'RelTol',1e-9);
    [T STATE tf statef] = ode45(@twoStanceDynamics,tspan,statef,leftLiftOffOptions);
    
    foot=0;
    en=0;
    
    if length(tf) == 0
            'break after 2'
            return;
    end 
    
    T_out = [T_out ; T];
    STATE_out = [STATE_out ; STATE];
    FOOT_out = [ FOOT_out ; foot];
    EN_out = [EN_out ; en];
    
    tspan = linspace(tf,tf+tSim,N);
    msStanceOptions = odeset('Events',@msStanceEvent,'AbsTol',1e-9,'RelTol',1e-9);
    [T STATE tf statef] = ode45(@rightStanceDynamics,tspan,statef,msStanceOptions);

    foot=0;
    en=0;
    
    if length(tf) == 0
            'break after 3'
            return;
    end
    
    T_out = [T_out ; T];
    STATE_out = [STATE_out ; STATE];
    FOOT_out = [ FOOT_out ; foot];
    EN_out = [EN_out ; en];
    end
    

    stateInit = x + [tdPos 0 0 0]';
%     stateInit = zeros(6,1); 
%     stateInit(2) = 0.003955567;
%     stateInit(3) = 0.994369144;
%     stateInit(4) = x; 
    tdParams=[u(1) u(2)];
    k=u(3);
%     tdPosL = [0 0.1+stateInit(2) 0];
    tdPosL = tdPos;
    tdPosR = [0 0 0];
    heightThreshold = L0*cos(tdParams(1));
    t0=0;
    [T_out STATE_out FOOT_out tf statef EN_out lowestState tdPosRNow]=halfSimulatePeriod(t0,stateInit);
%     y=statef';
    y=A*(statef'-[tdPosRNow 0 0 0]');
    
end