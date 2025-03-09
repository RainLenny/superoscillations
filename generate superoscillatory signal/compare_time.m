syms t w

% input_amps = [44.9368628213046	-181.152733973378	305.156183297628 -263.410775218804	114.366673605787	-18.9916426780696]
input_amps = real(transpose(input_amps));
superoscilation = sum(input_amps .* exp(1i * omega_n * t));

sinc_signal = sinc(t*scaling_factor/pi);

t_span = [-130,130]/scaling_factor;
%% Plot
figure;
hold on;

fplot(real(superoscilation),t_span, 'r', 'LineWidth', 4, 'DisplayName', 'superoscilation');
fplot(sinc_signal,t_span, 'b', 'LineWidth', 4, 'DisplayName', 'sinc');

% Customize plot
xlabel('Time (s)');
ylabel('Voltage (V)');
legend('show','fontweight','bold','fontsize',12);
grid on;
% ylim([0 1])

hold off;


