%% plot_robot.m
%
% Description:
%   Plots the robot in its current configuration.
%   
% Inputs:
%   q: robot configuration, q = [x_f; z_f; theta_f; theta_s; theta_m];
%   params: a struct with many elements, generated by calling init_params.m
%   varargin: optional name-value pair arguments:
%       'new_fig': (default: false), if true, plot appears on new figure
%       'trace_foot_com': (default false), if true, plots a tracer on the
%           foot's center of mass (CoM)
%       'trace_body_com': (default false), if true, plots a tracer on the
%           body's center of mass (CoM)
%       'trace_spine_tip': (default false), if true, plots a tracer on the
%           spine's tip
%
% Outputs:
%   none
%

function plot_robot(q,params,varargin)
%% Parse input arguments
% Note: a simple robot plotting function doesn't need this, but I want to
% write extensible code, so I'm using "varargin" which requires input
% parsing. See the reference below:
%
% https://people.umass.edu/whopper/posts/better-matlab-functions-with-the-inputparser-class/

% Step 1: instantiate an inputParser:
p = inputParser;

% Step 2: create the parsing schema:
%      2a: required inputs:
addRequired(p,'robot_config', ...
    @(q) isnumeric(q) && size(q,1)==5 && size(q,2)==1);
addRequired(p,'robot_params', ...
    @(params) ~isempty(params));
%      2b: optional inputs:
addParameter(p, 'new_fig', false); % if true, plot will be on a new figure

% Step 3: parse the inputs:
parse(p, q, params, varargin{:});

% Verification: display the results of parsing:
% disp(p.Results)

%% for convenience, define each generalized coordinate
x_f = q(1);
z_f = q(2);
theta_f = q(3);
theta_s = q(4);
theta_m = q(5);

%% Translate and rotate the foot into current configuration
% Use a homogeneous transformation to translate the CoM and to rotate about the CoM:

T_foot = [cos(theta_f), -sin(theta_f), x_f
          sin(theta_f),  cos(theta_f), z_f
          0,             0,            1  ];

foot.curr.corners = T_foot*params.foot.home.corners;

%% Translate and rotate the spine into current configuration
% use forward kinematics to get the proper translation
[FK_com_s,FK_com_b] = fk_com(q,params);
% [FK_pivot,FK_tip] = fk_viz(q,params);

% Use a homogeneous transformation to translate the CoM and to rotate about the CoM:
T_spine = [cos(theta_f+theta_s), -sin(theta_f+theta_s), FK_com_s(1)
           sin(theta_f+theta_s),  cos(theta_f+theta_s), FK_com_s(2)
           0,                     0,                    1          ];

spine.curr.corners = T_spine*params.spine.home.corners;

%% Translate and rotate the body into current configuration
T_body = [cos(theta_f+theta_s), -sin(theta_f+theta_s), FK_com_b(1)
          sin(theta_f+theta_s),  cos(theta_f+theta_s), FK_com_b(2)
          0,                     0,                    1          ];
       
body.curr.corners = T_body*params.body.home.corners;

%% Display the foot, spine, and body
if p.Results.new_fig
    figure;
end

fill(foot.curr.corners(1,:),foot.curr.corners(2,:),params.viz.colors.foot);
hold on;
fill(body.curr.corners(1,:),body.curr.corners(2,:),params.viz.colors.body);
fill(spine.curr.corners(1,:),spine.curr.corners(2,:),params.viz.colors.spine);
plot(x_f,z_f,'o','MarkerSize',6,...
    'MarkerFaceColor',params.viz.colors.com,...
    'MarkerEdgeColor',params.viz.colors.com);
plot(FK_com_b(1),FK_com_b(2),'o','MarkerSize',6,...
    'MarkerFaceColor',params.viz.colors.com,...
    'MarkerEdgeColor',params.viz.colors.com);
hold off;

axis(params.viz.axis_lims);
daspect([1 1 1]) % no distortion

xlabel('$x$');
ylabel('$y$');

end