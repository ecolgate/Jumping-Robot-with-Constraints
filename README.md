# Jumping-Robot-with-Constraints / Pumping Robot

----------
NEW (6/4/20): the digital controller now uses Colocated Partial Feedback Linearization.  One important thing to note:  to do this correctly, it is necessary to consider constraint forces.  The dynamics equation is:

M(q)q_ddot + H(q,q_dot) + A'F = Q

When deriving the CPFL controller, you need to use H + A'F instead of just H.  

In principle, this approach makes it easier to design the feedback controller, and ensures that the feedback controller works properly in all configuration of the robot.  However, I found it still quite challenging to pick the proper gains:  be prepared for some iteration.  Additionally, I had to speed up the sampling rate to 500 Hz to get everything to be well behaved.

-----------

PREVIOUS: 

Feedback control has been added and a simple pumping strategy has been implemented.  Many other changes for cleanliness, but one IMPORTANT CHANGE: 

The equation in `robot_dymamics.m` for acceleration has been changed to include a feedback term (see line 340 in `main.m`).  The feedback term is:

- A'*((A*A')\A)*q_dot/params.control.dt

This is a projection of the velocity (q_dot) onto the space of constraint violations ... basically, it is an error since there should be no constrain violation.  Since it is a velocity term that is being added to an equation for the accelerations, it is scaled by 1/dt, where dt is the controller update rate.  Some other scaling could be used, but this seems to work pretty well.  Note that I originally left this term out, but included it in the expression for velocities (where it still is).  The effect of that was for the integral of the velocities to be correct, but for the integral of accelerations (i.e., the velocities) to be incorrect.  This created weird behavior, especially with control.

Other changes: 

Feedback control has been added.  This has a simple structure with the function `digital_controller.m` computing the control commands (u) and the function `analog_controller.m` updating the state.

A state machine has been added for pumping.  It is also quite simple:  when the angle and angular velocity of the foot have the same sign, the body is pushed to a high position; when they have different signs, the body is lowered.  This strategy is based on:

Stephen Wirkus, Richard Rand & Andy Ruina (1998) How to Pump a Swing,
The College Mathematics Journal, 29:4, 266-275, DOI: 10.1080/07468342.1998.11973953

https://doi.org/10.1080/07468342.1998.11973953

Unilateral constraints (and changing constraints) have been removed for cleanliness.  They could be added back in if, for instance, we wanted the robot to jump.

I stopped using the symbolic inverse mass matrix and instead did a numerical inversion of the symbolic mass matrix.  It is hard to see an impact on performance, at least on my computer.

Following Zack Woodruff, I made the AVI video real-time by setting the frame rate equal to the time period between samples (params.viz.dt).  The video is pretty snappy.

I added energy tracking.  The "Net Energy" = Kinetic + Potential - Integral(Motor Power) and should be flat.  It is not.  I think this is just integration error (especially in estimating the motor power), but it may be evidence of some problem...

