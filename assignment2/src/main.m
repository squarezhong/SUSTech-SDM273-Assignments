clear
clc
warning('off','all')

traj_id=0;
tdoa_sigma=1e-4; % TDOA noise sigma
init_sigmas= [0.5 1 2]; % Initial value noise sigma, you can increase it, i.e. 0.5, 1, 2, ...
Mic_num=8; % Mic. Number
lim_t=2; % Time upper limit

for init_sigma = init_sigmas
    su_err = [];
    for epoch = 1:500
        ori_g=gt_generation(tdoa_sigma, traj_id, Mic_num);
        total_t = tic;
        while toc(total_t) < lim_t
            g = init_generation2(ori_g, init_sigma);
            sg = sg_generation(g);
            [sg, norm_dk, value_f] = GN_Solver(sg, 1, lim_t);
            if norm_dk < sg.dk_p || value_f < sg.f_p
                su_err = [su_err; [sg.rec, toc(total_t)]];
                break
            end
        end
    end

    % calculate the mean of the error
    % su_err(6) is the error of the microphone location of each epoch
    disp(size(su_err));
    mean_err = mean(su_err(:,6));
    sprintf('Average Mic Loc err of sigma %f: %f m', init_sigma, mean_err)
end




