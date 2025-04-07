% Rotational doublet on circular adhesive micropattern
% Calina Copos
% Last updated: 3/10/2025

%close all;
clear;

plot_flag       = 0;
video_flag      = 0;

% Micropattern surface
r_circle    = 0.5;
th          = linspace(0,2*pi,100); 
bdd         = [r_circle*cos(th);r_circle*sin(th)]';

% Position of the doublet cluster 
pos0        = [0.1,0;-0.1,0];
n           = min(size(pos0));

% Parameters
dt          = 0.001; % be careful! stiff system, need 0.001
Tmax        = 10.0;
mu          = -0.1;%-0.1;%-0.1;
%c_poly      = 5;
%c_rep       = 5;
%c_core      = 2;%5
%r_bdd       = 0.05;
%A           = 0.4342;
%r_cutoff    = -log(1/10)*A;
%r_core      = 0.25;
epsilon     = 0.01; %0.01;
%k_rep       = 5.0; % CIL
%dl0         = 0.5;
%k           = 1.0;
%dl_crit     = 0.25; % CIL
xi          = [0.92,1.5];%[0.9,1.22];%[0.88,0.92]; 
Torient     = 0.05; 

% Number of trials
n_trials = 1000;

% Relative strength of macroscopic forces
c_elastic = 0;   %20;
c_wall = 1.5;%20      ;%1.5;%2;
%c_cil = 5;

% Metrics
ccw_count   = 0;
cw_count    = 0;
switch_count= 0;
angvel      = zeros(n_trials,2); % time average angular velocity (ccw, cw)
rcc_count   = zeros(n_trials,2); % time average cell-cell separation (ccw, cw)

% Run multiple model realizations for a fixed parameter set
n = 2;
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

    % Evolve system in time
    for t = 1:Tmax/dt
        %n = length(pos);
        %th_pos = atan2(pos(:,2),pos(:,1)); %com position angle
        th_vel = atan2(v(:,2),v(:,1));
    
        % % Polymerization forces due to Rho & velocity alignment 
        % % (cell-cell interaction) 
        % argpos      = (pos-prev_pos)/norm(pos-prev_pos,2);
        % argterm     = [cos(th_rho),sin(th_rho)]*argpos';
        % % CW or CCW controlled by d(th_rho)/dt 
        % %th_rho      = th_rho + dt*((1/lambda)*(phi-th_rho));% + 1*2*pi*rand); % CCW
        % %th_rho      = th_rho + dt*(lambda*(th_rho-phi)); %CW
        % th_rho      = th_rho + dt*(lambda*(phi-th_rho)); %CW (fixed 10/7 to match VE formulation: th_rho tends to phi)
        % rho_F       = [cos(th_rho),sin(th_rho)];
        % phi         = phi + dt*((1/Torient)*asin(cos(phi).*sin(th_pos)-sin(phi).*cos(th_pos)) + mu*2*pi*rand);
        % attract_F   = [cos(phi),sin(phi)];
        % F_poly      = rho_F+attract_F;
        
        % Polymerization forces due to cubic polarization and velocity alignment 
        % As in Vercruysse et al. 2024; D'Orsogna et al. 2006
        %th_rho = th_rho + dt*( th_rho.*(1-abs(th_rho).^2) + 0.1*atan2(F_repulsive(:,2),F_repulsive(:,1)) + 0.15*th_vel + 0.05*rand );
        %F_poly = [cos(th_rho),sin(th_rho)];         
        %p = p + dt*( p.*(1-sqrt(p(:,1).^2+p(:,2).^2).^2) + 0.1*F_repulsive + 0.15*v + 0.05*rand(2,2) );
        %p = p./sqrt(p(:,1).^2+p(:,2).^2);
        %F_poly = p;
        % p gets huge, need a normalization maybe?    
        % F_poly = ( 1 - 1*sqrt(v(:,1).^2+v(:,2).^2).^2 ).*v + 0.05*rand(2,2);    

        % Polymerization forces due to spontaneous polarity AND 
        % active contribution, velocity alignment
        % As in Camley et al. 2014
        % periodic extension of -1/Torient*(th_rho-theta_v)
        theta_v     = th_vel;
        th_rho      = th_rho + dt*((1/Torient)*asin(cos(th_rho).*sin(theta_v)-sin(th_rho).*cos(theta_v)) + mu*rand(n,1));
        F_poly      = [cos(th_rho),sin(th_rho)];

        % Confinement/steric repulsion force to keep cell(s) 
        % on micropatterned surface
        %F_wall = compute_wallinteraction(pos,A,c_rep,c_core,r_circle,r_core,r);
        F_wall = compute_wallinteraction(pos,r_circle,epsilon);

        % Cell-cell adhesion force 
        % v1: linear, elastic spring 
        %F_cc = compute_cellcell(pos,k,dl0);
        % v2: Morse potential
        R = 0.2;%0.1; 
        A = 0.05;
        r = 0.01;%0.05;
        a = 0.2;
        F_cc = compute_cellcellv2(pos,R,A,r,a);

        % Report symmetry of cell-cell adhesion forces
        S1 = sum(F_cc);
        if abs(S1)>1e-13 
           sprintf('TROUBLE WITH SUM OF FORCES!!!\f')
        end
    
        % CIL-like repulsive force 
        %F_repulsive = compute_repulsive(pos,k_rep,dl_crit);
   
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
    
        % Plot
        if plot_flag && mod(t*dt,0.1)==0
            figure(4);
            subplot(1,2,1);
            ps = polyshape([bdd(:,1);bdd(1,1)],[bdd(:,2);bdd(1,2)]);
            plot(ps,'facecolor',[228 246 248]/256,'facealpha',1); hold on;
            plot([bdd(:,1);bdd(1,1)],[bdd(:,2);bdd(1,2)],'-b','linewidth',2); 
            plot(x_cc(id,:),y_cc(id,:),'-b','linewidth',2);
            scatter(pos(:,1),pos(:,2),200,'o','markerfacecolor',[135 135 135]/256,'markeredgecolor','k'); 
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
        %F = c_elastic*F_cc + c_wall*F_wall + F_poly + c_cil*F_repulsive;
        F = c_elastic*F_cc + c_wall*F_wall + F_poly;

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
    elseif (dot_perc>0.15) || (cw_perc>0.2 && ccw_perc>0.2)
    %if (cw_perc>0.2 && ccw_perc>0.2) || abs(p(1))<0.1 || (cw_perc>0.1 && ccw_perc>0.1 && abs(p(1))<1.5) %(mean(abs(diff(th_cc(5000:end))))<1e-4)
        s = 'yes';
    else 
        s = 'no';
        if cw_perc>0.2
            cw_count = cw_count+1;
        else 
            ccw_count = ccw_count+1;
        end
    end

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
                F_cc(j,:) = F_cc(j,:) + sign(pos(i,:)-pos(j,:))*(R*exp(-dl/r)-A*exp(-dl/a));
                F_cc(i,:) = F_cc(i,:) + sign(pos(j,:)-pos(i,:))*(R*exp(-dl/r)-A*exp(-dl/a));
                %dl_horiz_signed = pos(i,1)-pos(j,1);
                %F_cc(j,:) = F_cc(j,:) +  dl_horiz_signed*(R*exp(-dl/r)-A*exp(-dl/a))*(pos(i,:)-pos(j,:))/dl;
                %F_cc(i,:) = F_cc(i,:) -
                %(R*exp(-dl/r)-A*exp(-dl/a))*(pos(j,:)-pos(i,:))/dl;  %
                %this is completely unnnecessary! (2/13)
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
    th_tang = th_pos;

    % Andreas' confinement force
    F_wall(:,1) = exp((r_pos-r_circle)/epsilon).*(-cos(th_pos));
    F_wall(:,2) = exp((r_pos-r_circle)/epsilon).*(-sin(th_pos));

    % Linear confinement force
    %F_wall(:,1) = 0.5*r_pos.*(-cos(th_pos));
    %F_wall(:,2) = 0.5*r_pos.*(-sin(th_pos));

    % repulsive force like in Viscek (vahabli, vicsek nat comm phys 2023)
    % c_rep = 2; c_core = 2;%5
    % A = 0.4342; r_core = 0.25;
    % if any(r_pos-(r_circle-0.15)>0)
    %     F_wall(:,1) = ((r_pos-(r_circle-0.15))>0).*c_rep.*(-cos(th_pos));
    %     F_wall(:,2) = ((r_pos-(r_circle-0.15))>0).*c_rep.*(-sin(th_pos));
    %     F_wall(:,1) = F_wall(:,1) + 2*((r_pos-(r_circle-0.15))>0).*c_rep.*(-cos(th_tang));
    %     F_wall(:,2) = F_wall(:,2) + 2*((r_pos-(r_circle-0.15))>0).*c_rep.*(-sin(th_tang));
    % elseif any((r_pos-r_core)>0)
    %     F_wall(:,1) = c_core*((r_pos-r_core)>0).*(exp(r_pos-r_circle-1)/A).*(-cos(th_pos));
    %     F_wall(:,2) = c_core*((r_pos-r_core)>0).*(exp(r_pos-r_circle-1)/A).*(-sin(th_pos));
    %     F_wall(:,1) = F_wall(:,1) + 2*c_core*((r_pos-r_core)>0).*(exp(r_pos-r_circle-1)/A).*(-cos(th_tang));
    %     F_wall(:,2) = F_wall(:,2) + 2*c_core*((r_pos-r_core)>0).*(exp(r_pos-r_circle-1)/A).*(-sin(th_tang));
    % end

    % diff_pos = sqrt( (pos(:,1)-flipud(pos(:,1))).^2 + (pos(:,2)-flipud(pos(:,2))).^2 );
    % v_rep = exp((r_bdd + r_bdd - diff_pos)/A);
    % v_repuls_core = r_bdd + r_bdd - diff_pos;
    % v_repuls_core = 0.*(diff_pos>(r_bdd+r_bdd)) + v_repuls_core.*(diff_pos<=(r_bdd+r_bdd));
    % v_repuls = c_rep*sum(v_rep + c_core*v_repuls_core).*((flipud(pos)-pos)./diff_pos);

    % r_pos = sqrt(pos(:,1).^2+pos(:,2).^2);
    % th_pos = atan2(pos(:,2),pos(:,1));
    % v_wall_mag = 0*((r_circle-r_pos)>0) + c_core*((r_circle-r_pos)<=0).*(r_circle-r_pos-0.25).^2;
    % v_wall = -v_wall_mag.*[cos(th_pos),sin(th_pos)];

end
