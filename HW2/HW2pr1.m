% HW2 Pr1

T1 = [200:100:1500];
T2 = [50:30:300];
FA = [pi, 2/3 * pi, pi/3];
echo_spacing = 5; % ms

for fa_idx=1:length(FA)
    echo = zeros(64, length(T1), length(T2)); % for each T1,T2 pair, there are 64 echoess
    for T1_idx=1:length(T1)
        for T2_idx=1:length(T2)
            m = [0 0 1]'; % equilibrium state magnetization
            m = EPG_RF(m, pi/2, 0); % 90x excitation pulse
        
            % Although not specified in the question, as we need to count
            % for relaxation, so assuming the 90 pulse and alpha pulse 
            % spacing is half of the echo spacing = 2.5 ms
            
            % Assume crusher gradient for 1 unit cycle twist
            m = EPG_relax(m, T1(T1_idx), T2(T2_idx), echo_spacing/2); % relaxation before crusher
            m = EPG_grad(m, 1); % crusher gradient
            m = EPG_RF(m, pi/2 + FA(fa_idx)/2, pi/2); % 90 + alpha/2 pulse y-axis
            m = EPG_grad(m, 1); % crusher gradient
            m = EPG_relax(m, T1(T1_idx), T2(T2_idx), echo_spacing/2); % relaxation after crusher
            
            echo(1, T1_idx, T2_idx) = abs(m(1,1));
            
            for echo_num=2:64
                m = EPG_relax(m, T1(T1_idx), T2(T2_idx), echo_spacing/2); % relaxation before crusher
                m = EPG_grad(m, 1);
                m = EPG_RF(m, FA(fa_idx), pi/2); % alpha pulse y-axis
                m = EPG_grad(m, 1);
                m = EPG_relax(m, T1(T1_idx), T2(T2_idx), echo_spacing/2); % relaxation after crushe

                echo(echo_num, T1_idx, T2_idx) = abs(m(1,1));
            end
        end % T1
    end % T2


    % Pr1(a) amplitude plots for selected T1, T2

    % [T1,T2] = [(200,170) (800,170) (1500,50) (1500,170) (1500,290)]
    plot_idx = [1 5; 7 5; 14 1; 14 5; 14 9]; % indexing the selected (T1,T2)
    
    %  [T1,T2] = [(200,50) (500,110) (800,170) (1100,230) (1500,290)]
    % plot_idx = [1 1; 4 3; 7 5; 10 7; 14 9]; % indexing the selected (T1,T2)
    
    figure
    hold on; 
    legends = cell(1, 5);
    for plot_num=1:5
        T1_value = T1(plot_idx(plot_num,1));
        T2_value = T2(plot_idx(plot_num,2));
        
        plot(1:64, echo(:, plot_idx(plot_num,1), plot_idx(plot_num,2))', 'LineWidth', 1)
        legends{plot_num} = ['FA = ' num2str(FA(fa_idx)*180/pi) ', T1 = ' num2str(T1_value) ', T2 = ' num2str(T2_value)];
    end
    xlabel('echo number')
    ylabel('echo amplitude')
    legend(legends, 'Location', 'best')
    hold off;
    title(['FA = ' num2str(FA(fa_idx)*180/pi)])

    % Save figure
    plot_name = ['FA = ' num2str(FA(fa_idx)*180/pi)];
    folder = "plot/pr1_a/";
    filename = folder + plot_name + ".png";
    exportgraphics(gcf, filename, "ContentType","image")
    

    % Pr1(b)
    echo_indices = [6 16 32 48];
    
    for echo_idx = echo_indices
        figure
        echo_surface = squeeze(echo(echo_idx, :, :));
        contourf(T1, T2, echo_surface', 30);
        colorbar;
        xlabel('T1 values');
        ylabel('T2 values');
        zlabel('Echo amplitude');
        title(['Contour plot for echo' num2str(echo_idx) ' FA =' num2str(FA(fa_idx)*180/pi)]);
        plot_name = ['Contour plot for echo' num2str(echo_idx) ' FA =' num2str(FA(fa_idx)*180/pi)];
        folder = "plot/pr1_b/";
        filename = folder + plot_name + ".png";
        exportgraphics(gcf, filename, "ContentType","image");
    end

end % FA

