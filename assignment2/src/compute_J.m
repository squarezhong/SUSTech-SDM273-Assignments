% Compute Jacobian of non-linear least square of TDOA and odometry
function [J,r] = compute_J(g) % Compute Jacobian
    J = [];
    r = [];
    for j=1:g.K
        J_mic = zeros(g.M-1,g.m*g.M-8);
        J_s = zeros(g.M-1,g.n*g.K);
        s_loc = g.x(g.M+j,1:3); % s_j
        for i=2:g.M
            mic_loc = g.x(i,1:3); % x_i
            dx=mic_loc-s_loc; % x_i-s_j
            off = g.x(i,4); % \tau_i
            dri = g.x(i,5); % \delta_i
            % \Delta t_j=sum(g.dt(1:j))
            % g.cc is sound speed
            % t_{i,j}=g.tdoa(i-1,j), left side of i equals to 2,3,...,N
            if i==2
                J_mic(i-1,1:3)=[dx(1)/norm(dx)/g.cc,1,sum(g.dt(1:j))];
            elseif i==3
                J_mic(i-1,4:7)=[dx(1)/norm(dx)/g.cc,dx(2)/norm(dx)/g.cc,1,sum(g.dt(1:j))];
            else
                Eij = [dx / norm(dx) / g.cc, 1, sum(g.dt(1:j))];
                J_mic(i-1,7+g.m*(i-4)+1:7+g.m*(i-3))=Eij;
            end
            Gij = (-dx / norm(dx) - s_loc / norm(s_loc)) / g.cc;
            J_s(i-1,g.n*(j-1)+1:g.n*j)=Gij;
            r_tdoaij = (norm(dx)-norm(s_loc)) / g.cc + off + sum(g.dt(1:j)) * dri - g.tdoa(i-1,j);
            r = [r;r_tdoaij];
        end
        J = [J;[J_mic,J_s]];
        if j<g.K
            J_odo=zeros(3,g.m*g.M+g.n*g.K-8);
            J_odo(:,g.m*g.M-8+g.n*(j-1)+1:g.m*g.M-8+g.n*j)...
                =-eye(3);
            J_odo(:,g.m*g.M-8+g.n*j+1:g.m*g.M-8+g.n*(j+1))...
                =eye(3);
            s_now = g.x(g.M+j,1:3)'; % s_j
            s_next = g.x(g.M+j+1,1:3)'; % s_{j+1}
            % m_j=g.S(:,j+1)-g.S(:,j)
            r_odoj = s_next - s_now - g.S(:,j+1) + g.S(:,j);
            r = [r;r_odoj];
            J = [J;J_odo];
        end
    end
end