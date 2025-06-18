syms t w

angular_freqs = [0.4, 0.45, 0.5, 0.55, 0.6, 0.65];
input_amps = [17.551969178838473, -49.93615200929133, 48.91892825922387, -7.222327926051708, -17.23049063370233, 8.876842152413422];
superoscilation = sum(input_amps .* exp(1i * angular_freqs * t));

sinc_signal = sinc(t/pi);

scaling_factor = 1;

%% Fourier transform:
F_superoscilation = fourier(superoscilation, t, w);
F_sinc_signal = fourier(sinc_signal, t, w);

%% sample the Fourier transform because delta cannot be plotted
w_axis = (-2:1/20:2);

sampled_F_superoscilation = double(subs(F_superoscilation,w_axis));
infty_idx = abs(sampled_F_superoscilation) == Inf; % find Inf
sampled_F_superoscilation(infty_idx) = 1; % set Inf to finite valus

sampled_F_sinc = double(subs(F_sinc_signal,w_axis));
%% Plot
figure;
hold on;

plot(w_axis*scaling_factor, sampled_F_superoscilation,'-o','color', 'r', 'LineWidth', 2, 'DisplayName', 'superoscilation');
plot(w_axis*scaling_factor,sampled_F_sinc,'-o','color', 'b', 'LineWidth', 2, 'DisplayName', 'sinc');

% Customize plot
xlabel('Frequency (rad/s)');
ylabel('Amplitude (Normalized)');
legend('show','fontweight','bold','fontsize',12);
grid on;
% ylim([0 1])
ax=gca;
ax.FontWeight = 'bold';
ax.FontSize = 13;
hold off;


