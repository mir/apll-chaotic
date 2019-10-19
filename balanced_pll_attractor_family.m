clear all;

tau_p1 = 1.65*10^(-3); % 1.65ms
tau_p2 = tau_p1;

tau_z1 = 165*10^(-6); % 165us
tau_z2 = tau_z1;

sim_t_step=10^-4;
tspan = (0:sim_t_step:3);
options = odeset('RelTol', 1e-14, 'AbsTol', 1e-14);

F_tf = @(s)((1+s*tau_z1).*(1+s*tau_z2)./((1+s*tau_p1).*(1+s*tau_p2)));

% filter equations
[A,B,C,D]=tf2ss([tau_z1*tau_z2 (tau_z1+tau_z2) 1],[tau_p1*tau_p2 (tau_p1+tau_p2) 1]);
duty_1 = 0.4;
duty_2 = 0.4;
pd_char = @(x)balanced_pd_char(x,duty_1,duty_2);
x_0 = [A^-1*B*0.5; -0.5];

G_0_dbs = [57.1 71.1 78.6 79.3 81.0 81.3 83.8 85.7 87.3 88.6 89.8 89.9 90.8 91.7...
    92.5 93.3 94.0 94.6 95.2 95.8 96.3 96.8 97.3 97.7];
G_0s = db2mag(G_0_dbs);
omega_es = [0 700 1400];
for G_0_db=G_0_dbs
    for omega_e = omega_es
        G_0 = db2mag(G_0_db);
        K_vco = G_0;
        
        rh = @(t,x)([A*x(1:end-1)+B*pd_char(x(end));...
            omega_e-K_vco*(C*x(1:end-1)+D*pd_char(x(end)))]);
        % transients
        [T,X]=ode45(rh,[0,1],x_0,options);
        x_1 =X(end,:);

        [T,X]=ode45(rh,tspan,x_1,options);
        X_1 = X(:,1);
        X_2 = X(:,2);
        Theta = X(:,3);
        %% plot the phase portrait
        fig_1 = figure(1);
        plot3(X_1(1:1:end),...
            X_2(1:1:end),...
            Theta(1:1:end));
        hold on;
        plot3(x_0(1),x_0(2),x_0(3),'*');
        hold off;
        xlabel('x_1');
        ylabel('x_2');
        zlabel('\Theta');
        fig_1_name = sprintf('table_pp/G0_%.0f_omega_%.0f_pp',G_0_db*10,omega_e);
        print(fig_1,fig_1_name,'-dpng');

        fig_2 = figure(2);
        pd_output = pd_char(Theta);
        vco_control_voltage = [X_1,X_2]*C'+D*pd_output;
        plot(pd_output,vco_control_voltage,...
            pd_output(1),vco_control_voltage(1),'*');
        xlabel('PD output');
        ylabel('VCO control voltage');
        fig_2_name = sprintf('table_pv/G0_%.0f_omega_%.0f_pv',G_0_db*10,omega_e);
        print(fig_2,fig_2_name,'-dpng');
    end
end
