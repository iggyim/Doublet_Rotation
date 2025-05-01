% Rotational doublet on circular adhesive micropattern
% Calina Copos
% Last updated: 3/10/2025

warning('off','all')
%close all;
clear;

plot_flag       = 1;
video_flag      = 0;

% Micropattern surface
r_circle    = 0.5;
th          = linspace(0,2*pi,100); 
bdd         = [r_circle*cos(th);r_circle*sin(th)]';

% Position of the doublet cluster 
pos0        = [0.2,0;-0.2,0];
n           = min(size(pos0));

% Parameters
dt          = 0.001; % be careful! stiff system, need 0.001
Tmax        = 10.0;
mu          = -0.1;
epsilon     = 0.01;%0.01;%0.1; 
dl0         = 0.5;
k           = 1.0;
xi          = [1,1];%[2.9,3.0];%[0.9,1.22];%[1.4,1]; 
Torient     = 0.5; % Suspect unphysical results for Torient below 0.1

% Number of trials
n_trials = 1000;

% Relative strength of macroscopic forces
A = 0.05; % strength of attractive component of cell-cell adh
R = 0.2; % strength of repulsive component of cell-cell adh
c_elastic = 5; % minimum 10
c_wall = 1.5;%3;
k_poly = 100;

% Metrics
ccw_count   = 0;
cw_count    = 0;
switch_count= 0;
angvel      = zeros(n_trials,2); % time average angular velocity (ccw, cw)
rcc_count   = zeros(n_trials,2); % time average cell-cell separation (ccw, cw)

% Run multiple model realizations for a fixed parameter set
n = 2;

% Bookkeeping
Nt = Tmax/dt;
X1 = zeros(Nt,n_trials);
Y1 = zeros(Nt,n_trials);
T1 = zeros(Nt,n_trials); 
L1 = zeros(n_trials,1); % label: switch:1 and nonswitch:0
X2 = zeros(Nt,n_trials);
Y2 = zeros(Nt,n_trials);
T2 = zeros(Nt,n_trials); 
L2 = zeros(n_trials,1); % label: switch:1 and nonswitch:0
for j = 1:n_trials

    % Clear variables
    clear t pos phi v prev_pos th_cc sign_th_cc r_cc dot_cc 
    
    % Initialize positions, velocities, polarities of cells in doublet
    pos     = pos0;
    v       = zeros(n,2);
    th_poly = 2*pi*rand(n,1);
    p       = [cos(th_poly),sin(th_poly)];
    
    phi0    = 2*pi*rand(n,1);
    phi     = phi0;
    
    prev_pos = zeros(n,2); % previous position
    v        = p;          % previous velocity (determined by polarity vector only)
    
    % Post-sim metrics
    th_cc  = zeros(Tmax/dt,1);
    sign_th_cc = zeros(Tmax/dt,1);
    r_cc   = zeros(Tmax/dt,1);
    dot_cc   = zeros(Tmax/dt,1);
    
    % Generate random orientation for polarity vector
    r = rand(1);
    th_rho = 2*pi*rand(n,1); %3/17 %-pi+2*pi*rand(n,1);
    
    if(video_flag)
        vidObj = VideoWriter('test_double_ccw.mp4','MPEG-4');
        open(vidObj);
    end

    %bias = rand(n,1);
    c_poly = 1;
    cond_cil = 0;
    % Evolve system in time
    for t = 1:Nt
        th_vel = atan2(v(:,2),v(:,1));

        % Polymerization forces due to spontaneous polarity AND 
        % active contribution, velocity alignment
        % As in Camley et al. 2014
        % periodic extension of -1/Torient*(th_rho-theta_v)
        theta_v     = th_vel;
        th_rho      = th_rho + dt*((1/Torient)*asin(cos(th_rho).*sin(theta_v)-sin(th_rho).*cos(theta_v)) + mu*rand(n,1));
        F_poly      = [cos(th_rho),sin(th_rho)];

        r_cc(t) = norm(pos(2,:)-pos(1,:),2);
        if t>1000 && r_cc(t)<2*r
            cond_cil = 1;
            c_poly = c_poly - dt*k_poly*c_poly; 
            % d(c_poly)/dt = -k_poly*(c_poly)
        end
       
        % Confinement/steric repulsion force to keep cell(s) 
        % on micropatterned surface
        F_wall = compute_wallinteraction(pos,r_circle,epsilon);

        % Cell-cell adhesion force 
        % v1: linear, elastic spring 
        %F_cc = compute_cellcell(pos,k,dl0);
        % v2: Morse potential 
        r = 0.05;
        a = 0.2;
        F_cc = 0.1*compute_cellcellv2(pos,R,A,r,a) + compute_repulsive(pos,0.2,0.1);

        % Report symmetry of cell-cell adhesion forces
        S1 = sum(F_cc);
        if abs(S1)>1e-13 
           sprintf('TROUBLE WITH SUM OF FORCES!!!\f')
        end
   
        % Determine cell-cell plane (visual effect only)
        p       = polyfit(pos(:,1),pos(:,2),1);
        m       = p(1);
        x1      = mean(pos(:,1)); y1=mean(pos(:,2));
        r_cc(t) = norm(pos(2,:)-pos(1,:),2);
        x_cc    = linspace(-0.5,0.5,1000)';
        y_cc    = (-1/m)*x_cc+(y1+x1/m); %%%%%% !!!!!! (half-domain only?)
        id      = find((sqrt(x_cc.^2+y_cc.^2))<=0.5);
        [~,ii]  = max(y_cc(id));
        r_cc(t) = norm(pos(2,:)-pos(1,:),2);
        if (isempty(ii)==0)
            th_cc(t)   = atan2(y_cc(id(ii)),x_cc(id(ii)));
            sign_th_cc(t) = sign(th_cc(t)-th_cc(t-1));
            %w_cc(t) = abs(th_cc(t)-th_cc(t-1))/dt;
        end

        dot_cc(t)   = dot(v(2,:),v(1,:));
        th_cc(t)    = atan2(pos(2,2)-pos(1,2),pos(2,1)-pos(1,1));
    
        % Track directionality
        if (t==1 || (abs(th_cc(t)-th_cc(t-1))<1e-6))
            c = 'ks';
            sign_th_cc(t) = 0;
        elseif sign(th_cc(t)-th_cc(t-1))<0
            c = 'b*';
            sign_th_cc(t) = -1;
        elseif sign(th_cc(t)-th_cc(t-1))>0
            c = 'ro';
            sign_th_cc(t) = 1;
        end

        % Bookkeeping
        X1(t,j) = pos(1,1); X2(t,j) = pos(2,1);
        Y1(t,j) = pos(1,2); Y2(t,j) = pos(2,2);
        T1(t,j) = t*dt; T2(t,j) = t*dt;
    
        % Plot
        if plot_flag && mod(t*dt,0.1)==0
            figure(4);
            subplot(1,2,1);
            ps = polyshape([bdd(:,1);bdd(1,1)],[bdd(:,2);bdd(1,2)]);
            plot(ps,'facecolor',[228 246 248]/256,'facealpha',1); hold on;
            plot([bdd(:,1);bdd(1,1)],[bdd(:,2);bdd(1,2)],'-b','linewidth',2); 
            plot(x_cc(id,:),y_cc(id,:),'-b','linewidth',2);
            scatter(pos(1,1),pos(1,2),200,'o','markerfacecolor',[135 135 135]/256,'markeredgecolor','k'); 
            scatter(pos(2,1),pos(2,2),200,'o','markerfacecolor',[0 0 0]/256,'markeredgecolor','k'); 
            quiver(pos(:,1),pos(:,2),0.1*v(:,1),0.1*v(:,2),'r','linewidth',2);
            quiver(pos(:,1),pos(:,2),F_wall(:,1),F_wall(:,2),'b','linewidth',1,'autoscale','off');
            quiver(pos(:,1),pos(:,2),F_poly(:,1),F_poly(:,2),'k','linewidth',2,'autoscale','off');
           % quiver(pos(:,1),pos(:,2),F_cc(:,1),F_cc(:,2),'g','linewidth',1,'autoscale','off');
            xlim([-1.0 1.0]); ylim([-1.0 1.0]); box on; grid off;
            set(gca,'plotBoxAspectRatio',[1 1 1]);
            set(gca,'FontSize',20,'fontname','dejavu sans'); %set(gca,'Color','k')
            set(gcf,'color','w'); %set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
            x0=200;y0=500;width=1000;height=400;
            set(gcf,'position',[x0,y0,width,height]);
            currFrame = getframe(gcf);
            pause(0.1)
            %keyboard()
            hold off;
    
            subplot(1,2,2);
            scatter(t*dt,th_cc(t),c); hold on;
            set(gca,'FontSize',20);
            set(gcf,'color','w'); box on;
            xlim([0 Tmax]); ylim([-3.5 3.5]); grid on;
            title('Blue star: CW, Red circle: CCW')
            set(gca,'fontname','dejavu sans');
            xlabel('Time'); ylabel('Rotational angle');
        end
     
        % Newton's 2nd law (overdamped approximation)
        F = c_elastic*F_cc + c_wall*F_wall + c_poly*F_poly;

        prev_pos    = pos;
        v           = F./xi;
        pos         = pos + v*dt;
    
        if(video_flag)
            writeVideo(vidObj,currFrame);
        end
       
    end
    
    % Report mean cell-cell separation metric
    %rcc_count(j,1) = mean(r_cc(500:end).*(sign_th_cc(500:end)>0));
    %rcc_count(j,2) = mean(r_cc(500:end).*(sign_th_cc(500:end)<0));

    % Clean up angular velocity and find a `clean' fit
    ttt = th_cc';
    tf  = ischange(ttt,'linear'); % 1: indicates abrupt change in mean
    [~,s1,s2] = ischange(ttt,'linear','threshold',200);
    idsofchange = find(tf==1);
    segline = s1.*(1:length(tf)) + s2;
    %p = polyfit((idsofchange(end):length(tf))*dt,segline(idsofchange(end):length(tf)),1);

    % Keep track of rotational directionality
    cw_perc = sum(sign_th_cc<0)/(Tmax/dt);
    ccw_perc = sum(sign_th_cc>0)/(Tmax/dt);
    stuck_perc = sum(sign_th_cc==0)/(Tmax/dt);
    dot_perc = sum(abs(dot_cc)<0.2)/(Tmax/dt); % amount/percent parrallel
    if stuck_perc>0.1
        nr = 'yes'; s = 'yes';
    %elseif (dot_perc>0.15) || (round(cw_perc,2)>=0.2 && round(ccw_perc,2)>=0.2) % prior to 4/9
    elseif (round(cw_perc,2)>=0.2 && round(ccw_perc,2)>=0.2)
    %elseif (cw_perc>0.2 && ccw_perc>0.2) || abs(p(1))<0.1 || (cw_perc>0.1 && ccw_perc>0.1 && abs(p(1))<1.5) %(mean(abs(diff(th_cc(5000:end))))<1e-4)
        s = 'yes';
    else 
        s = 'no';
        if cw_perc>0.2
            cw_count = cw_count+1;
        else 
            ccw_count = ccw_count+1;
        end
    end

    % Bookkeeping
    L1(j,1) = 1*strcmp(s,'yes') + 0*strcmp(s,'no');
    L2(j,1) = L1(j,1);

    % Report angular velocity
    if strcmp(s,'no')
        % Use the last switch and smooth out the data
        start_int = idsofchange(end);
        end_int = length(tf);
        %keyboard()
        if cw_perc<0.2
            angvel(j,1) = abs(p(1));
            rcc_count(j,1) = mean(r_cc(start_int:end_int).*(sign_th_cc(start_int:end_int)>0));
        else
            angvel(j,2) = abs(p(1));
            rcc_count(j,2) = mean(r_cc(start_int:end_int).*(sign_th_cc(start_int:end_int)<0));
        end

        % Compute angular velocity as the slope over a window size of 500
        % in a region of no abrupt change (computed by tf)
        % and away from the start of the simulation
        % trial = 0;
        % for i=1:(length(tf)-500)
        %     size_domain = 500; % arbitrarily chosen
        %     if nnz(tf(i:i+size_domain)) == 0 
        %         if trial>500 % arbitrarily chosen to avoid initial domain
        %             if cw_perc<0.2
        %                 angvel(j,1) = mean(nonzeros(w_cc(i:i+size_domain).*(sign_th_cc(i:i+size_domain)>0))); %ccw
        %                 rcc_count(j,1) = mean(r_cc(i:i+size_domain).*(sign_th_cc(i:i+size_domain)>0));
        %             else
        %                 angvel(j,2) = mean(nonzeros(w_cc(i:i+size_domain).*(sign_th_cc(i:i+size_domain)<0))); %cw
        %                 rcc_count(j,2) = mean(r_cc(i:i+size_domain).*(sign_th_cc(i:i+size_domain)<0));
        %             end
        %         end
        %         trial = trial + 1;
        %     end
        % end
    end

    if plot_flag == 1
        formatSpec = "CCW percentage: %.3f, CW percentage: %.3f, switch: %s";
        sprintf(formatSpec,ccw_perc,cw_perc,s)
        keyboard()
    end

    if(video_flag)
        close(vidObj);
    end
end

% Truncate nans from angular velocity and mean separation
tmp_angvel = angvel;
truncated_angvel_ccw = tmp_angvel(~isnan(tmp_angvel(:,1)),1); %ccw
truncated_angvel_cw = tmp_angvel(~isnan(tmp_angvel(:,2)),2);  %cw

tmp_rcc = rcc_count;
truncated_rcc_ccw = tmp_rcc(~isnan(tmp_rcc(:,1)),1); %ccw
truncated_rcc_cw = tmp_rcc(~isnan(tmp_rcc(:,2)),2); %cw

% Metrics
formatSpec = "CCW (out of rotating): %.3f, CW (out of rotating): %.3f, R: %.3f";
sprintf(formatSpec,ccw_count/(ccw_count+cw_count),cw_count/(ccw_count+cw_count),(ccw_count+cw_count)/n_trials)

formatSpec = "CCW: avg rotational speed: %.3f, avg cell-cell separation: %.3f";
sprintf(formatSpec,mean(nonzeros(truncated_angvel_ccw)),mean(nonzeros(truncated_rcc_ccw)))

formatSpec = "CW: avg rotational speed: %.3f, avg cell-cell separation: %.3f";
sprintf(formatSpec,mean(nonzeros(truncated_angvel_cw)),mean(nonzeros(truncated_rcc_cw)))

[h1,p1,~,~] = ttest2(nonzeros(truncated_angvel_ccw),nonzeros(truncated_angvel_cw));
[h2,p2,~,~] = ttest2(nonzeros(truncated_rcc_ccw),nonzeros(truncated_rcc_cw));
formatSpec = "CCW/CW p-value for rotational speed: %.6f, cell-cell separation: %.6f";
sprintf(formatSpec,p1,p2)


%% OPTIONAL for post-processing
figure;
for i=1:100
    subplot(10,10,i); hold on;
    plot(X1(:,i),Y1(:,i),'-b');
    plot(X2(:,i),Y2(:,i),'-r');
    %plot(T1(:,i),X1(:,i),'-b');
    title(L1(i));
end
%save('runs_label','T1','X1','Y1','X2','Y2','L1');
%%

%%%%%%%%%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CIL-like repulsion 
function [F_repulsive] = compute_repulsive(pos,k_rep,dl_crit)
    F_repulsive = zeros(length(pos),2);
    n = length(pos);

    i = 1; j = 2;
    %for i=1:n
    %    for j=1:n
            dl = sqrt( (pos(i,1)-pos(j,1))^2 + (pos(i,2)-pos(j,2))^2 );
            if (dl<dl_crit)%&&(j~=i)
                dl = sqrt( (pos(i,1)-pos(j,1))^2 + (pos(i,2)-pos(j,2))^2 );
                F_repulsive(j,:) = -k_rep*(pos(i,:)-pos(j,:))/dl;
                F_repulsive(i,:) = -k_rep*(pos(j,:)-pos(i,:))/dl;
            end
    %    end
    %end
end 

% Cell-cell adhesion forces (Morse potential)
function [F_cc] = compute_cellcellv2(pos,R,A,r,a)
    n = min(size(pos));
    F_cc = zeros(n,2);

    i = 1; j = 2;
    %for i=1:n
        %for j=1:n
            %if j~=i
                dl = sqrt( (pos(i,1)-pos(j,1))^2 + (pos(i,2)-pos(j,2))^2 );
                
                % fixed sign (4/14)
                F_cc(j,:) = F_cc(j,:) + sign(pos(i,:)-pos(j,:))*(A*exp(-dl/a)-R*exp(-dl/r));
                F_cc(i,:) = F_cc(i,:) + sign(pos(j,:)-pos(i,:))*(A*exp(-dl/a)-R*exp(-dl/r));
            %end
        %end
    %end

end

%Cell-cell adhesion forces (linear elastic spring)
function [F_cc] = compute_cellcell(pos,k,dl0)
    n = min(size(pos));
    F_cc = zeros(n,2);

    i = 1; j = 2;
    %for i=1:n
    %    for j=1:n
    %        if j~=i
                dl = sqrt( (pos(i,1)-pos(j,1))^2 + (pos(i,2)-pos(j,2))^2 );
                F_cc(j,:) = F_cc(j,:) + k*(dl/dl0-1.0)*(pos(i,:)-pos(j,:))/dl;
                F_cc(i,:) = F_cc(i,:) + k*(dl/dl0-1.0)*(pos(j,:)-pos(i,:))/dl;
    %        end
    %    end
    %end
end

% Cell-matrix or confinement forces
function [F_wall] = compute_wallinteraction(pos,r_circle, epsilon)
    r_pos = sqrt(pos(:,1).^2+pos(:,2).^2);
    th_pos = atan2(pos(:,2),pos(:,1));
    F_wall = zeros(min(size(pos)),2);

    % Andreas' confinement force
    F_wall(:,1) = exp((r_pos-r_circle)/epsilon).*(-cos(th_pos));
    F_wall(:,2) = exp((r_pos-r_circle)/epsilon).*(-sin(th_pos));

    % Linear confinement force
    %F_wall(:,1) = 0.5*(r_pos).*(-cos(th_pos));
    %F_wall(:,2) = 0.5*(r_pos).*(-sin(th_pos));
end
