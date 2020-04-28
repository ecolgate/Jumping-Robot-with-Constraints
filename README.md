# Jumping-Robot-with-Constraints
Matlab code for simulating the jumping robot including computation and visualization of ground reaction forces.
Robot Design Studio 2020
Developed by Ed Colgate borrowing heavily from Dan Lynch's cart-pendulum code.

## Overview
Application entry point is `main.m`, which calls (through wrapper functions) dynamics-related functions that are automatically generated by running `derive_equations_JR.mlx`.
That means you have to run `derive_equations_JR.mlx` before you run `main.m` for the first time.

`derive_equations_JR.mlx` uses symbolic computation to generate the state-space dynamics of the jumping robot, which are then exported as MATLAB functions (e.g., `autogen_constraints.m` and `autogen_mass_matrix.m`).

### Autogenerated functions and wrappers
As in the cart-pendulum example, extensive use is made of autogenerated functions and wrapper functions.

### MATLAB structs and parsers
In addition to physical parameters (mass, length, etc) that anchor the dynamics in reality, there are many other parameters relevant to simulation, such as timestep size, appearance of the robot, etc.
I have grouped the parameters into a struct called `params`, which is generated by calling `init_params.m`.
This struct get passed around throughout the simulation, and saves you from having to worry about the order in which you supply input arguments.

Another way to avoid input argument order-dependency (a trademark of brittle code) is to use Name-Value pairs, which you have probably experienced if you've used MATLAB before.
Name-Value pairs rely on an "input parser"; check out `plot_robot.m` and `animate_robot.m` to see how this works.

## Nested Functions
`ode45` calls a function named `robot_dynamics.m` that is placed inside of (i.e., nested within) `main.m` as well as a function named `robot_events.m` for determining when an event, such as a constraint force going to zero, occurs.  This allows the three functions to share variables (including `params` and `F`) without explicit passing.  This is helpful since the constraint forces `F` are computed within `robot_dynamics.m` but cannot be easily passed back through `ode45`.
