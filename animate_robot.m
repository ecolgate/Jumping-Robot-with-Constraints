%% animate_robot.m
%
% Description:
%   Animates the robot according to a list of configurations over time.
%   
% Inputs:
%   q_list: 5xN list of configurations q
%   F_list: 3xN list of constraint forces F
%   params: a struct with many elements, generated by calling init_params.m
%
% Outputs:
%   none

function animate_robot(q_list,F_list,params,varargin)

% Parse input arguments
% Note: a simple robot animation function doesn't need this, but I want to
% write extensible code, so I'm using "varargin" which requires input
% parsing. See the reference below:
%
% https://people.umass.edu/whopper/posts/better-matlab-functions-with-the-inputparser-class/

% Step 1: instantiate an inputParser:
p = inputParser;

% Step 2: create the parsing schema:
%   2a: required inputs:
addRequired(p,'robot_config', ...
    @(q) isnumeric(q_list) && size(q_list,1)==5);
addRequired(p,'constraint_forces', ...
    @(q) isnumeric(F_list) && size(F_list,1)==4);
addRequired(p,'robot_params', ...
    @(params) ~isempty(params));
%   2b: optional inputs:
%       optional name-value pairs to trace different parts of the robot:
addParameter(p, 'trace_foot_com', false);
addParameter(p, 'trace_body_com', false);
addParameter(p, 'trace_spine_tip', false);
addParameter(p, 'show_constraint_forces', false);
addParameter(p, 'video', false);

% Step 3: parse the inputs:
parse(p, q_list, F_list, params, varargin{:});

fig_handle = figure('Renderer', 'painters', 'Position', [10 10 900 600]);

% if tracing is true for anything, then trace
if (p.Results.trace_foot_com || p.Results.trace_body_com ...
        || p.Results.trace_spine_tip)
    tracing = true;
else
    tracing = false;
end

if tracing
    foot.curr.com.x = [];
    foot.curr.com.z = [];

    body.curr.com.x = [];
    body.curr.com.z = [];

    spine.curr.tip.x = [];
    spine.curr.tip.z = [];
end

% if video is desired...
if p.Results.video
    v = VideoWriter('pushoff.avi');
    open(v);
end

for i = 1:size(q_list,2)
    plot_robot(q_list(:,i),params,'new_fig',false);
    
    % for convenience, define some terms
    q = q_list(:,i);
    x_f = q(1);
    z_f = q(2);
    theta_f = q(3);
    F = F_list(:,i);

    if tracing
        % append (x,z) location of foot CoM:
        [FK_com_s,FK_com_b] = fk_com(q,params);
        foot.curr.com.x = [foot.curr.com.x, x_f];
        foot.curr.com.z = [foot.curr.com.z, z_f];

        % append (x,z) location of body CoM:
        body.curr.com.x = [body.curr.com.x, FK_com_b(1)];
        body.curr.com.z = [body.curr.com.z ,FK_com_b(2)];

        % append (x,z) location of pendulum tip:
        [FK_pivot,FK_tip] = fk_viz(q,params);
        spine.curr.tip.x = [spine.curr.tip.x, FK_tip(1)];
        spine.curr.tip.z = [spine.curr.tip.z, FK_tip(2)];

        if p.Results.trace_foot_com
            hold on;
            plot(foot.curr.com.x,foot.curr.com.z,'o-',...
                'Color',params.viz.colors.tracers.foot_com,...
                'MarkerSize',3,'LineWidth',2,...
                'MarkerFaceColor',params.viz.colors.tracers.foot_com,...
                'MarkerEdgeColor',params.viz.colors.tracers.foot_com);
            hold off;
        end
        if p.Results.trace_body_com
            hold on;
            plot(body.curr.com.x,body.curr.com.z,'o-',...
                'Color',params.viz.colors.tracers.body_com,...
                'MarkerSize',3,'LineWidth',2,...
                'MarkerFaceColor',params.viz.colors.tracers.body_com,...
                'MarkerEdgeColor',params.viz.colors.tracers.body_com);
            hold off;
        end
        if p.Results.trace_spine_tip
            hold on;
            plot(spine.curr.tip.x,spine.curr.tip.z,'o-',...
                'Color',params.viz.colors.tracers.spine_tip,...
                'MarkerSize',3,'LineWidth',2,...
                'MarkerFaceColor',params.viz.colors.tracers.spine_tip,...
                'MarkerEdgeColor',params.viz.colors.tracers.spine_tip);
            hold off;
        end
    end
    
    if p.Results.show_constraint_forces
        % find the locations of the left and right bottom corners of the
        % foot        
        T_foot = [cos(theta_f),         -sin(theta_f),        x_f
                  sin(theta_f),         cos(theta_f),         z_f
                  0,                    0,                    1   ];
        foot.curr.corners = T_foot*params.foot.home.corners;  % the ones we want are in rows 1&2 of columns 2&3
        % compute the reaction forces at the two corners, assuming that the
        % x-direction force is distributed between the two according to
        % their z-direction (normal) forces
        Fscale = 2*params.model.geom.foot.hbot/params.model.dyn.body.m/params.model.dyn.g;  % scale factor for visualization
        foot.curr.force.left.x = F(1)*Fscale;
        foot.curr.force.left.z = F(2)*Fscale;
        foot.curr.force.right.x = F(3)*Fscale;
        foot.curr.force.right.z = F(4)*Fscale;
        % create vectors of x and z ends for both left and right
        xleft = [foot.curr.corners(1,2),foot.curr.corners(1,2)+foot.curr.force.left.x];
        zleft = [foot.curr.corners(2,2),foot.curr.corners(2,2)+foot.curr.force.left.z];
        xright = [foot.curr.corners(1,3),foot.curr.corners(1,3)+foot.curr.force.right.x];
        zright = [foot.curr.corners(2,3),foot.curr.corners(2,3)+foot.curr.force.right.z];
        % draw the vectors
        hold on;
        line(xleft,zleft,'Color',params.viz.colors.vectors,'LineWidth',2)
        line(xright,zright,'Color',params.viz.colors.vectors,'LineWidth',2)
        hold off;
    end

    if p.Results.video
        M(i) = getframe(fig_handle);
        writeVideo(v,M(i));
    end
end

if p.Results.video
    close(v);
%         movie(gcf,M); % comment this out if you don't want to see replay
end
    
end