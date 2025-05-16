function dxdt = wpt_dynamics_pid(t, x, Ltx, Lrx, M, Ctx, Crx, Rload, Kp, Ki, Kd, r)
    I_tx = x(1);
    I_rx = x(2);
    E_int = x(3);  
    V_Ctx = x(4);  
    V_Crx = x(5);  

   
    Vrx = Rload * I_rx;
    error = r(t) - Vrx;

    
    dVrx_dt = 0;

   
    Vtx = Kp * error + Ki * E_int + Kd * (-dVrx_dt);

    
    A = [Ltx, M;
         M, Lrx];
    rhs = [Vtx - V_Ctx;
           -V_Crx - Rload * I_rx];

    dI = A \ rhs;
    dI_tx = dI(1);
    dI_rx = dI(2);

    
    dV_Ctx = (1 / Ctx) * I_tx;
    dV_Crx = (1 / Crx) * I_rx;
    dE_int = error;

   
    dxdt = [dI_tx;
            dI_rx;
            dE_int;
            dV_Ctx;
            dV_Crx];
end