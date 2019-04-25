function channel_labels = homologous_chans()

% Match up homologous channels on the left and right

% Load the magnetometer layout
cfg = [];
cfg.layout = 'neuromag306mag.lay';
lay = ft_prepare_layout(cfg);

% Find which channels are on each side
x = lay.pos;
left_inx = find(x(:,1) < -0.01);
right_inx = find(x(:,1) > 0.01);

% % Make sure that the left and right channels align
% x_left = x(left_inx, :);
% x_right = x(right_inx, :);
% plot(x_left(:,1), x_left(:,2), 'ro');
% hold on
% plot(-x_right(:,1), x_right(:,2), 'g+');
% hold off

% Find each pair of homologous channels
left_homolog_inx = nan(size(left_inx));
for i_chan = 1:length(left_inx)
    l_inx = left_inx(i_chan);
    left_chan_pos = x(l_inx,:);
    text(left_chan_pos(1), left_chan_pos(2), num2str(l_inx))
    right_pos_dist = abs(x - left_chan_pos);
    right_x_inx = abs(-x(:,1) - left_chan_pos(1)) < 0.01;
    right_y_inx = abs(x(:,2) - left_chan_pos(2)) < 0.01;
    right_chan_inx = find(right_x_inx & right_y_inx);
    if length(right_chan_inx) == 1 % It found 1 homologous channel
        left_homolog_inx(i_chan) = right_chan_inx;
    elseif length(right_chan_inx) > 1
        error('Found too many homologous channels')
    elseif length(right_chan_inx) < 0
        error('Didn''t find any homologous channels')
    end
end
clear right_x_inx right_y_inx right_chan_inx right_pos_dist

% Make the list of channel names
channel_labels = cell([0 2]);
for i_chan = 1:length(left_inx)
    % Append the magnetometers
    mag_label_left = lay.label{left_inx(i_chan)};
    mag_label_right = lay.label{left_homolog_inx(i_chan)};
    if ~startsWith(mag_label_left, 'MEG')
        continue
    end
    channel_labels{end + 1, 1} = mag_label_left;
    channel_labels{end, 2} = mag_label_right;
    % Append the gradiometers
    for n_grad = 2:3
        mag2grad = @(s) [s(1:(end-1)) num2str(n_grad)];
        channel_labels{end + 1, 1} = mag2grad(mag_label_left);
        channel_labels{end, 2} = mag2grad(mag_label_right);
    end
end
