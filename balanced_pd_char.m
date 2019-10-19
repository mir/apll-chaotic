function v_e = balanced_pd_char(theta_es,duty_1,duty_2)
v_e = zeros(size(theta_es));
for k=1:length(theta_es)
    theta_e = theta_es(k);
    assert(duty_1<=duty_2);
    v_max = 1-2*abs(duty_1-duty_2);
    v_min = 1-2*min([duty_1+duty_2,2-duty_1-duty_2]);
    theta_delta = pi/2*(v_max-v_min);
        if duty_1<=duty_2
        %     rising
            theta_11 = 2*pi*(duty_1-duty_2) - theta_delta;
            theta_12 = 2*pi*(duty_1-duty_2);
        %     falling
            theta_21 = 0;
            theta_22 = theta_delta;
        else
            %     rising
            theta_11 = - theta_delta;
            theta_12 = 0;
        %     falling
            theta_21 = 2*pi*(duty_1-duty_2);
            theta_22 = theta_21+theta_delta;
        end
    theta = mod(theta_e,2*pi);
    if theta<=theta_22
        v_e(k) = v_max-2/pi*theta;
    elseif theta <= theta_22+(2*pi+theta_11-theta_22)
        v_e(k) = v_min;
    elseif theta <= 2*pi+theta_12
        v_e(k) = v_min+2/pi*(theta-(theta_22+(2*pi+theta_11-theta_22)));
    else
        v_e(k) = v_max;
    end
end
end
