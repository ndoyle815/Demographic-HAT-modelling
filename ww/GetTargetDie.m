%This code returns the input parameter "p_targetdie" (the max prob of killing tsetse flies per
%contact with the VC targets) in the simple vector model based on VCreduction (percentage of reduction in vector populations),
%ReductionMeasured (VCreduction is measured in ReductionMeasured days) and
%TargetDeploy (yearly frequency of deploying targets)
    

function p_targetdie = GetTargetDie(VCreduction, TargetFreq, TrapCycle, Paras) %VCreduction:0-100(%) TargetFreq:0-infty(times) TrapCycle:0-infty(days)
    % fixed parameters, will be taken out if calling from the main function
    ODEoptions = odeset('NonNegative', 1:3, 'RelTol',1e-8,'AbsTol',1e-8);
    global K N_V
    K = 100 * Paras.mu_V / (Paras.xi_V^2 * Paras.p_survive * (Paras.p_survive * Paras.B_V / Paras.mu_V - 1));
    
    p_targetdie = 0;

    if VCreduction > 0 && TargetFreq > 0
        IC = 100 * [Paras.mu_V / (Paras.xi_V * Paras.p_survive) (Paras.mu_V / (Paras.mu_V + Paras.alpha)) (Paras.alpha/(Paras.mu_V + Paras.alpha))];

        % initial p_targetdie value
        if VCreduction/TargetFreq >= 60 && VCreduction/TargetFreq < 95
            p_targetdie = 0.05;
        elseif VCreduction/TargetFreq >= 95
            p_targetdie = 0.2;
        else
            p_targetdie = 0.01;
        end
    
        % find p_targetdie value to achieve VCreduction after TrapCycle (with 0.01% absolute tolerance) in our model with intervention 
        PercentSurvive = 100 - VCreduction;
        Delta = VCreduction;
    
        while abs(Delta) > 0.01
            p_targetdie = p_targetdie + 0.0001 * Delta;
            [t,pop] = ode45(@TsetseDyn, [0 TrapCycle], IC, [], TargetFreq, p_targetdie, Paras);
            Delta = N_V - PercentSurvive;
        end
    end
end


function dPop = TsetseDyn(t, pop, TargetFreq, p_targetdie, Paras)
    global K N_V
    
    if p_targetdie==0
        f_T = 0;
    else
        f_T = p_targetdie * (1 - 1/(1+exp(-25/365*(mod(t,365/TargetFreq)-0.35*365))));
    end

    dPop = zeros(3,1);
    P_V = pop(1);
    S_V = pop(2);
    G_V = pop(3);
    N_V = S_V + G_V;

    % Pupa
    dPop(1) = Paras.B_V * N_V - (Paras.xi_V + P_V/K) * P_V;

    % Teneral (feed twice as fast as other adults)
    dPop(2) = Paras.xi_V * Paras.p_survive * P_V - 1 * Paras.alpha * S_V - Paras.mu_V * S_V;

    % Non-teneral
    dPop(3)= Paras.alpha * (1 - f_T) * S_V - Paras.alpha * f_T * G_V - Paras.mu_V * G_V;
end
