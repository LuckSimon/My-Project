function run3dSLIPofWalking()
%%初始化参数
    g=9.8;    m=80; L0=1;
    k=30000;
    deltaK=0;
    
    N = 1000;      % Number of Samples
    tSim = 2;     % Max Sim time
        
    tdPosL = [0 0.15 0];
    tdPosR = [0 0 0];
    tdParams = [pi/6 0];
    heightThreshold = 0;
    stepLength = 0.35;
    stepWidth = 0.2;
    
    stepTime = 0;%%后期获取
    
    
    %%最优控制
    xbest = [];
    ubest = [];
    lqrCoefficient = [];
    Jx = [];
    Ju = [];
    
    
%     Ce=[1 4 2 1 2 2]; 
%     Ce=[2 5 4 5 1 1]; %%滑模面参数 
    Ce=[6 16 6 10 8 4]; %%滑模面参数 
%     Ce=[3 5 4 5 1 1];
    
    A=eye(6,6);
    A(2,2) = -1;
    A(5,5) = -1;
    B=eye(3,3);
    B(2,2) = -1;
    stepNumbers=0;
    stepNumbers=0;
    
    vx =1;
    
%% 事件切换函数定义
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
%% 动力学方程
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
%% 一个周期仿真
    function [T_out STATE_out FOOT_out tf statef EN_out TD_params] = simulatePeriod(t0,stateInit)
    tspan = linspace(t0,t0+tSim,N);
    state0 = stateInit;
    rightTouchDownOptions = odeset('Events',@RightTouchDownEvents,'AbsTol',1e-12,'RelTol',1e-12);
    [T STATE tf statef] = ode45(@leftStanceDynamics,tspan,state0,rightTouchDownOptions);
    
    foot=[ones(size(STATE,1),1)*tdPosL,STATE(:,1:3)];
    en=0;
    
    figure(1)
    subplot(411)
    hold on
%     plot3(STATE(:,1),STATE(:,2),STATE(:,3),'b-.');
    plot(T,STATE(:,3));
%     plot(T,STATE(:,6));
%     plot3(STATE(:,4),STATE(:,5),STATE(:,6),'r');
    subplot(412)
    hold on
    plot(T,STATE(:,2));
%     hold on
%     plot(tdPosL(1),tdPosL(2),'ro');
    subplot(413)
    hold on
    plot(T,STATE(:,4));
%     relPos = STATE(1:3) - tdPosL;
%     plot(T,sqrt(transpose(sum(transpose(relPos.^2))))','r')
    subplot(414)
    hold on
    plot(T,-k*(sum((STATE(:,1:3)-tdPosL).^2,2)-1)*cos(tdParams(1)));
    
    figure(3)
    hold on
    plot3(STATE(:,1),STATE(:,2),STATE(:,3));
    
    STATE_out = STATE;
    T_out = T;
    FOOT_out = foot;
    EN_out = en;
    TD_params=tdPosL;
    
%     
%     tdParams(2)=-1*tdParams(2);
    tdPosR = [STATE(end,1)+cos(tdParams(2))*sin(tdParams(1)) STATE(end,2)-sin(tdParams(2))*sin(tdParams(1)) 0];
    tspan = linspace(tf,tf+tSim,N);
    LowestHighOptions = odeset('Events',@LowestHighEvent,'AbsTol',1e-12,'RelTol',1e-12);
    [T STATE tf statef] = ode45(@twoStanceDynamics,tspan,statef,LowestHighOptions);
    
    foot=[ones(size(STATE,1),1)*tdPosL,ones(size(STATE,1),1)*tdPosR];
    en=0;
    
    figure(1)
    subplot(411)
    hold on
    plot(T,STATE(:,3));
%     plot(T,STATE(:,6));
%     plot3(STATE(:,1),STATE(:,2),STATE(:,3),'b-.');
%     plot3(STATE(:,4),STATE(:,5),STATE(:,6),'r');
    subplot(412)
    hold on
    plot(T,STATE(:,2));
%     relPos = STATE(1:3) - tdPosL;
%     plot(T,sqrt(transpose(sum(transpose(relPos.^2))))','r')
    subplot(413)
    hold on
    plot(T,STATE(:,4));
    
    subplot(414)
    hold on
    plot(T,-k*(sum((STATE(:,1:3)-tdPosL).^2,2)-1)*cos(tdParams(1)));
    
    figure(3)
    hold on
    plot3(STATE(:,1),STATE(:,2),STATE(:,3));
    
    T_out = [T_out ; T];
    STATE_out = [STATE_out ; STATE];
    FOOT_out = [ FOOT_out ; foot];
    EN_out = [EN_out ; en];
    TD_params=[TD_params;tdPosR];
    
    
    tspan = linspace(tf,tf+tSim,N);
    leftLiftOffOptions = odeset('Events',@leftLiftOffEvents,'AbsTol',1e-12,'RelTol',1e-12);
    [T STATE tf statef] = ode45(@twoStanceDynamics,tspan,statef,leftLiftOffOptions);
    
    foot=[ones(size(STATE,1),1)*tdPosL,ones(size(STATE,1),1)*tdPosR];
    en=0;
    
    
    figure(1)
    subplot(411)
    hold on
    plot(T,STATE(:,3));
    
    subplot(412)
    hold on
    plot(T,STATE(:,2));
%     relPos = STATE(1:3) - tdPosL;
%     plot(T,sqrt(transpose(sum(transpose(relPos.^2))))','r')
    subplot(413)
    hold on
    plot(T,STATE(:,4));
    
    subplot(414)
    hold on
    plot(T,-k*(sum((STATE(:,1:3)-tdPosL).^2,2)-1)*cos(tdParams(1)));
    
    figure(3)
    hold on
    plot3(STATE(:,1),STATE(:,2),STATE(:,3));
    
    T_out = [T_out ; T];
    STATE_out = [STATE_out ; STATE];
    FOOT_out = [ FOOT_out ; foot];
    EN_out = [EN_out ; en];
    
    
    
    tspan = linspace(tf,tf+tSim,N);
    msStanceOptions = odeset('Events',@msStanceEvent,'AbsTol',1e-12,'RelTol',1e-12);
    [T STATE tf statef] = ode45(@rightStanceDynamics,tspan,statef,msStanceOptions);
    
    foot=[STATE(:,1:3),ones(size(STATE,1),1)*tdPosR];
    en=0;
   
    
    figure(1)
    subplot(411)
    hold on
    plot(T,STATE(:,3));
    
    subplot(412)
    hold on
    plot(T,STATE(:,2));
%     relPos = STATE(1:3) - tdPosL;
%     plot(T,sqrt(transpose(sum(transpose(relPos.^2))))','r')
    subplot(413)
    hold on
    plot(T,STATE(:,4));
    
    subplot(414)
    hold on
    plot(T,0*k*(sum((STATE(:,1:3)-tdPosL).^2,2)-1));
    
    figure(3)
    hold on
    plot3(STATE(:,1),STATE(:,2),STATE(:,3));
    
    T_out = [T_out ; T];
    STATE_out = [STATE_out ; STATE];
    FOOT_out = [ FOOT_out ; foot];
    EN_out = [EN_out ; en];
    
    
%%滑模控制
    stepNumbers=stepNumbers+1;
    rightMsStanceState = statef-[tdPosR 0 0 0];
    stateDistance = rightMsStanceState'-(A^stepNumbers)*xbest;
    pinv(Ce*Ju);
    stateDistance=(A^stepNumbers)*stateDistance;
    Sn=Ce*stateDistance
    ueq=-pinv(Ce*Ju)*Ce*(Jx-eye(6,6))*stateDistance;
    [FDn,s1,s2,s3]=FD(stateDistance);
    
    if 2*Sn*Ce*Ju*FDn*stateDistance+(Ce*Ju*FDn*stateDistance)^2>0
        'break after 2'
    end
    
    Un=ubest+ueq+FDn*stateDistance;
    ueq+FDn*stateDistance;
    Snn=Ce*Ju*(ueq+FDn*stateDistance)+Ce*Jx*stateDistance;
    tdParams=Un(1:2);
    k=Un(3);
    
    tspan = linspace(tf,tf+tSim,N);
    leftTouchDownOptions = odeset('Events',@leftTouchDownEvents,'AbsTol',1e-12,'RelTol',1e-12);
    [T STATE tf statef] = ode45(@rightStanceDynamics,tspan,statef,leftTouchDownOptions);
    
    foot=[STATE(:,1:3),ones(size(STATE,1),1)*tdPosR];
    en=0;
   
    
    figure(1)
    subplot(411)
    hold on
    plot(T,STATE(:,3));
    
    subplot(412)
    hold on
    plot(T,STATE(:,2));
%     relPos = STATE(1:3) - tdPosL;
%     plot(T,sqrt(transpose(sum(transpose(relPos.^2))))','r')
    subplot(413)
    hold on
    plot(T,STATE(:,4));
    
    subplot(414)
    hold on
    plot(T,0*k*(sum((STATE(:,1:3)-tdPosL).^2,2)-1));
   
    figure(3)
    hold on
    plot3(STATE(:,1),STATE(:,2),STATE(:,3));
    
    T_out = [T_out ; T];
    STATE_out = [STATE_out ; STATE];
    FOOT_out = [ FOOT_out ; foot];
    EN_out = [EN_out ; en];
    
    
%     tdParams(2)=-1*tdParams(2);
    tdPosL = [STATE(end,1)+cos(tdParams(2))*sin(tdParams(1)) STATE(end,2)+sin(tdParams(2))*sin(tdParams(1)) 0];
%     subplot(312)
%     hold on
%     plot(tdPosL(1),tdPosL(2),'ro');
    
    tspan = linspace(tf,tf+tSim,N);
    rightLiftOffOptions = odeset('Events',@rightLiftOffEvents,'AbsTol',1e-12,'RelTol',1e-12);
    [T STATE tf statef] = ode45(@twoStanceDynamics,tspan,statef,rightLiftOffOptions);
    
    foot=[ones(size(STATE,1),1)*tdPosL,ones(size(STATE,1),1)*tdPosR];
    en=0;
    
    figure(1)
    subplot(411)
    hold on
    plot(T,STATE(:,3));
    
    subplot(412)
    hold on
    plot(T,STATE(:,2));
%     relPos = STATE(1:3) - tdPosL;
%     plot(T,sqrt(transpose(sum(transpose(relPos.^2))))','r')
    subplot(413)
    hold on
    plot(T,STATE(:,4));
    
    subplot(414)
    hold on
    plot(T,-k*(sum((STATE(:,1:3)-tdPosL).^2,2)-1)*cos(tdParams(1)));
    
    figure(3) 
    hold on
    plot3(STATE(:,1),STATE(:,2),STATE(:,3));
    
    T_out = [T_out ; T];
    STATE_out = [STATE_out ; STATE];
    FOOT_out = [ FOOT_out ; foot];
    EN_out = [EN_out ; en];
%     TD_params=[TD_params;tdPosL];
    
    tspan = linspace(tf,tf+tSim,N);
    msStanceOptions = odeset('Events',@msStanceEvent,'AbsTol',1e-12,'RelTol',1e-12);
    [T STATE tf statef] = ode45(@leftStanceDynamics,tspan,statef,msStanceOptions);
    
    foot=[ones(size(STATE,1),1)*tdPosL,STATE(:,1:3)];
    en=0;
    
    figure(1)
    subplot(411)
    hold on
    plot(T,STATE(:,3));
    
    subplot(412)
    hold on
    plot(T,STATE(:,2));
%     relPos = STATE(1:3) - tdPosL;
%     plot(T,sqrt(transpose(sum(transpose(relPos.^2))))','r')
    subplot(413)
    hold on
    plot(T,STATE(:,4));
    
    subplot(414)
    hold on
    plot(T,-k*(sum((STATE(:,1:3)-tdPosL).^2,2)-1)*cos(tdParams(1)));
    
    figure(3)
    hold on
    plot3(STATE(:,1),STATE(:,2),STATE(:,3));
    
    T_out = [T_out ; T];
    STATE_out = [STATE_out ; STATE];
    FOOT_out = [ FOOT_out ; foot];
    EN_out = [EN_out ; en];
    
    %% 滑模控制趋近律控制
    stepNumbers = stepNumbers + 1;
    stateDistance=(statef-[tdPosL 0 0 0])'-(A^stepNumbers)*xbest;
    stateDistance;
    pinv(Ce*Ju);
    stateDistance=(A^stepNumbers)*stateDistance;
    Sn=Ce*stateDistance
    ueq=-pinv(Ce*Ju)*Ce*(Jx-eye(6,6))*stateDistance;
    [FDn,s1,s2,s3]=FD(stateDistance);
    
    if 2*Sn*Ce*Ju*FDn*stateDistance+(Ce*Ju*FDn*stateDistance)^2>0
        'break after 3'
    end
    ueq+FDn*stateDistance;
    Un=ubest+ueq+FDn*stateDistance;
    
    Snn=Ce*Jx*stateDistance+Ce*Ju*(ueq+FDn*stateDistance);
    tdParams=Un(1:2);
    k=Un(3);
    
    
    
    end

%% 二分之一周期仿真
    function [T_out STATE_out FOOT_out tf statef EN_out lowestState stepLengthNow stepWidthNow stepTimenow] = halfSimulatePeriod(t0,stateInit)
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
    stepLengthNow = 0.25;
    stepWidthNow = 0.25;
    stepTimenow=0;
    
    if length(tf) == 0
            'break after 1'
            return;
    end
    
    tdPosR = [STATE(end,1)+cos(tdParams(2))*sin(tdParams(1)) STATE(end,2)-sin(tdParams(2))*sin(tdParams(1)) 0];
    stepLengthNow =abs(tdPosR(1)-tdPosL(1));
    stepWidthNow = abs(tdPosR(2)-tdPosL(2));
    
    subplot(312)
    hold on
    plot(tdPosR(1),tdPosR(2),'ro');
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
    stepTimenow = tf;
    
    if length(tf) == 0
            'break after 3'
            return;
    end
    
    T_out = [T_out ; T];
    STATE_out = [STATE_out ; STATE];
    FOOT_out = [ FOOT_out ; foot];
    EN_out = [EN_out ; en];
    end


%% 寻找最优解
    function output = periodResidual(params)
     
     stateInit = zeros(6,1);     
     tdParams = params(1:2);
     heightThreshold = L0*cos(tdParams(1));
     k = params(3);
     stateInit(2) = 0.05;
%      stateInit(3) = sqrt(params(4)^2-(tdPosL(2)-stateInit(2))^2);
     stateInit(3) = params(4);
     stateInit(4) = vx;
     t0=0;
     
     [T_out STATE_out FOOT_out tf statef EN_out lowestState stepLengthNow stepWidthNow stepTimenow] =  halfSimulatePeriod(t0,stateInit);
     output = zeros(8,1);
     output(1) = STATE_out(end,4)-vx;
     output(2) = STATE_out(end,5);
     output(3) = (STATE_out(end,2)+stateInit(2))/2-tdPosL(2)/2-tdPosR(2)/2;
     output(4) = (STATE_out(end,1)+stateInit(1))/2-tdPosL(1)/2-tdPosR(1)/2;
     output(5) = STATE_out(end,3)-stateInit(3);
     output(6) = STATE_out(end,6);
     output(7) = lowestState(2)-tdPosL(2)/2-tdPosR(2)/2;
     output(8) = lowestState(1)-tdPosL(1)/2-tdPosR(1)/2;
%      output(9) = stepLengthNow - stepLength;
%      output(9) = stepWidthNow - stepWidth;
    end

    function [Jx,Ju]=sovleJocobi(X,U,tdPos)
        Xbest = X;
        Ubest = U;
        h = 0.00001;
    %     slipFunction(Xn,Un)
        Jx = [];
        Ju = [];
        for j=1:length(Xbest)
           xPositiveH = Xbest;
           xPositive2H = Xbest;
           xNegativeH = Xbest;
           xNegative2H = Xbest;
           xPositiveH(j) = xPositiveH(j) +h;
           xPositive2H(j) = xPositive2H(j) +2*h;
           xNegativeH(j) = xNegativeH(j) - h;
           xNegative2H(j) = xNegative2H(j) - 2*h;
           s1 =(-slipFunction(xPositive2H ,Ubest,tdPos)+8*slipFunction(xPositiveH ,Ubest,tdPos)-8*slipFunction(xNegativeH,Ubest,tdPos)+slipFunction(xNegative2H ,Ubest,tdPos))/(12*h);
           
           xPositiveH = Xbest;
           xPositive2H = Xbest;
           xNegativeH = Xbest;
           xNegative2H = Xbest;
           xPositiveH(j) = xPositiveH(j) +2*h;
           xPositive2H(j) = xPositive2H(j) +2*2*h;
           xNegativeH(j) = xNegativeH(j) - 2*h;
           xNegative2H(j) = xNegative2H(j) - 2*2*h;
           s2 =(-slipFunction(xPositive2H ,Ubest,tdPos)+8*slipFunction(xPositiveH ,Ubest,tdPos)-8*slipFunction(xNegativeH,Ubest,tdPos)+slipFunction(xNegative2H ,Ubest,tdPos))/(24*h);
           s=(16*s1-s2)/15;
          
           Jx = [Jx  s];
        end
        
        for j=1:length(Ubest)
           UPositiveH = Ubest;
           UPositive2H = Ubest;
           UNegativeH = Ubest;
           UNegative2H = Ubest;
           UPositiveH(j) = UPositiveH(j) +h;
           UPositive2H(j) = UPositive2H(j) +2*h;
           UNegativeH(j) = UNegativeH(j) - h;
           UNegative2H(j) = UNegative2H(j) - 2*h;
           s1=(-slipFunction(Xbest ,UPositive2H,tdPos) + 8*slipFunction(Xbest ,UPositiveH,tdPos) - 8*slipFunction(Xbest,UNegativeH,tdPos) + slipFunction(Xbest,UNegative2H,tdPos))/(12*h);
           
           UPositiveH = Ubest;
           UPositive2H = Ubest;
           UNegativeH = Ubest;
           UNegative2H = Ubest;
           UPositiveH(j) = UPositiveH(j) +2*h;
           UPositive2H(j) = UPositive2H(j) +2*2*h;
           UNegativeH(j) = UNegativeH(j) - 2*h;
           UNegative2H(j) = UNegative2H(j) - 2*2*h;
           s2 =(-slipFunction(Xbest ,UPositive2H,tdPos) + 8*slipFunction(Xbest ,UPositiveH,tdPos) - 8*slipFunction(Xbest,UNegativeH,tdPos) + slipFunction(Xbest,UNegative2H,tdPos))/(24*h);
           
           s=(16*s1-s2)/15;
           
           Ju = [Ju  s];
        end
        % s=(-slipFunction(1+2*h,Un)+8*slipFunction(1+h,Un)-8*slipFunction(1-h,Un)+slipFunction(1-2*h,Un))/12*h
    end  %%在最优解附近求解雅克比矩阵

    function [F,s1,s2,s3]=FD(x)
       s=Ce*x;
       matrixCeB=Ce*Ju;
       totalAbsXk=0;
       absTotalXk=0;
       totalXk=0;
       T0=0;
       
       for i=1:6
           totalAbsXk=abs(x(i))+totalAbsXk;
       end
       
       for i=1:6
           absTotalXk=x(i)+absTotalXk;
       end
       absTotalXk=abs(absTotalXk);
       
       for i=1:6
           totalXk=totalXk+x(i);
       end
       
       for i=1:6
            T0=abs(matrixCeB(1)*x(i))+abs(matrixCeB(2)*x(i))+abs(matrixCeB(3)*x(i))+T0;
       end
       
       f0=2*abs(s)*absTotalXk/(totalAbsXk*T0*2);
       
       if s*totalXk > 0.5*f0*totalAbsXk*T0
            if matrixCeB(1)>0
                f1=-f0;
                s1= -s*matrixCeB(1)*f1*totalXk-matrixCeB(1)*f1*(0.5*f0*totalAbsXk*T0);
            else
                f1=+f0;
                s1= -s*matrixCeB(1)*f1*totalXk-matrixCeB(1)*f1*(0.5*f0*totalAbsXk*T0);
            end
            
            if matrixCeB(2)>0
                f2=-f0;
                s2= -s*matrixCeB(2)*f2*totalXk-matrixCeB(2)*f2*(0.5*f0*totalAbsXk*T0);
            else
                f2=+f0;
                s2= -s*matrixCeB(2)*f2*totalXk-matrixCeB(2)*f2*(0.5*f0*totalAbsXk*T0);
            end
            
            if matrixCeB(3)>0
                f3=-f0;
                s3= -s*matrixCeB(3)*f3*totalXk-matrixCeB(3)*f3*(0.5*f0*totalAbsXk*T0);
            else
                f3=+f0;
                s3= -s*matrixCeB(3)*f3*totalXk-matrixCeB(3)*f3*(0.5*f0*totalAbsXk*T0);
            end
       elseif s*totalXk < -0.5*f0*totalAbsXk*T0
            if matrixCeB(1)>0
                f1=f0;
                s1= -s*matrixCeB(1)*f1*totalXk-matrixCeB(1)*f1*(0.5*f0*totalAbsXk*T0);
            else
                f1=-f0;
                s1= -s*matrixCeB(1)*f1*totalXk-matrixCeB(1)*f1*(0.5*f0*totalAbsXk*T0);
            end
            
            if matrixCeB(2)>0
                f2=f0;
                s2= -s*matrixCeB(2)*f2*totalXk-matrixCeB(2)*f2*(0.5*f0*totalAbsXk*T0);
            else
                f2=-f0;
                s2= -s*matrixCeB(2)*f2*totalXk-matrixCeB(2)*f2*(0.5*f0*totalAbsXk*T0);
            end
            
            if matrixCeB(3)>0
                f3=f0;
                s3= -s*matrixCeB(3)*f3*totalXk-matrixCeB(3)*f3*(0.5*f0*totalAbsXk*T0);
            else
                f3=-f0;
                s3= -s*matrixCeB(3)*f3*totalXk-matrixCeB(3)*f3*(0.5*f0*totalAbsXk*T0);
            end
       else
           f1=0;
           f2=0;
           f3=0;
           s1=0;
           s2=0;
           s3=0;
       end
       
       if s1<0||s2<0||s3<0
            'break after 1'
       end
       F=[f1 f1 f1 f1 f1 f1;f2 f2 f2 f2 f2 f2;f3 f3 f3 f3 f3 f3];
%        y=[F s1 s2 s3];        
    end


    tdParams = [];
    tdPosL = [0 0.15 0];
    tdPosR = [0 0 0];
    params0 = [0.25 0.3489 20000 0.992 1];
    
    options = optimset('display','iter','TolFun',1e-9);
    [params resnorm res] = lsqnonlin(@periodResidual,params0,[0.139,0.139,10000,0.95],[0.523,0.523,100000,0.9995],options);
    params
    close all

    %%调整初始点使得关于质心的y方向轨迹关于x轴对称
    stateInit = zeros(6,1);     
    tdParams = params(1:2);
    heightThreshold = L0*cos(tdParams(1));
    k = params(3);
    stateInit(2) = 0.05;
    stateInit(3) = params(4);
    stateInit(4) = vx;
    t0=0;
   [T_out STATE_out FOOT_out tf statef EN_out lowestState stepLengthNow stepWidthNow stepTimenow]=halfSimulatePeriod(t0,stateInit);
   comWidth=abs(statef(2)-stateInit(2)) ;
   stepTime = stepTimenow;
   
   close all
   %%多周期仿真初始参数
   stateInit = zeros(6,1);     
   tdParams = params(1:2);
   heightThreshold = L0*cos(tdParams(1));
   k = params(3);
   stateInit(2) = comWidth/2;
   stateInit(3) = params(4);
   stateInit(4) = vx;
   t0=0;
   figure(1)
   subplot(411)
   title('z方向位移')
   xlabel('x (m)'); ylabel('z (m)'); 

   subplot(412)
   title('足位置以及y方向位移');
   xlabel('Time (s)');
   ylabel('y(m)');
   
   subplot(413)
   title('x方向速度');
   xlabel('Time (s)');
   ylabel('Vx(m/s)');
   
   %% 多个周期仿真
   tdPosL = [0 0.1+comWidth/2 0];
   tdPosR = [0 0 0];
   stateInit
   ts=[];
   States=[];
   foots=[];
   td_position=[];
   xbest = stateInit - [tdPosL 0 0 0]';
   ubest = [tdParams k]';
   [Jx,Ju] = sovleJocobi(xbest,ubest,tdPosL);
   
   Jx
   Ju
%    stateInit(4)=1
   for i=1:2
        [T_out STATE_out FOOT_out tf statef EN_out TD_params] = simulatePeriod(t0,stateInit);
        t0 = tf;
        stateInit = statef; 
        ts = [ts ; T_out]; States = [States ; STATE_out];
        foots=[foots;FOOT_out];
        td_position=[td_position;TD_params];
   end
   
   
   %% 动画
   tic
   figure(3)
   title("质心轨迹");
   xlabel('X(m)');
   ylabel('Y(m)');
   zlabel('Z(m)');
   
   
   figure(2); clf; hold on;
   title(sprintf("3D-SLIP Walking Animation for Vx = %.2f m/s",vx));
   xlabel('X(m)');
   ylabel('Y(m)');
   zlabel('Z(m)');
   for i=1:1:size(td_position,1)
      figure(2)
      plot3(td_position(i,1),td_position(i,2),td_position(i,3),'bo');
      hold on;
   end
   hOA = plot3([stateInit(1) 0],[stateInit(2) 0.169],[ stateInit(3) 0],'r');
   hOB = plot3([stateInit(1) stateInit(1)],[stateInit(2) stateInit(2)],[stateInit(3) stateInit(3)],'r');
   hA = plot3(stateInit(1),stateInit(2),stateInit(3),'ro');
   hHist = plot3(stateInit(1),stateInit(2),stateInit(3),'r--');
   
   axis([0 max(States(:,1))  min(foots(:,5)) max(foots(:,2)) 0 max(States(:,3))]); 
   set(hA,'markerfacecolor','k');
   
   view([28 10]);
   startTime = cputime;
   for i = 1:4:length(ts)
        set(hOA,'zdata',[States(i,3) foots(i,3)],'ydata',[States(i,2) foots(i,2)],'xdata',[States(i,1) foots(i,1)]);
        set(hOB,'zdata',[States(i,3) foots(i,6)],'ydata',[States(i,2) foots(i,5)],'xdata',[States(i,1) foots(i,4)]);
        set(hA,'zdata',States(i,3),'ydata',States(i,2),'xdata',States(i,1));
        set(hHist,'zdata',States(1:i,3),'ydata',States(1:i,2),'xdata',States(1:i,1));
        %while (cputime-startTime) < ts(i)
            pause(.03);
        %end
   end
   toc
   ts(end)
%%
end