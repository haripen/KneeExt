# KneeExt
A simple Matlab function for simulating a concentric leg-press knee extension task.

Harald Penasso & Sigrid Thaller*

INPUT
 - Properties of the system, defining anthropometry, force-velocity relation, activation dynamics, environment and initial conditions
 
OUTPUT
 - Internal and external forces and velocities, position of the leg-press sledge, the geometrical ratio and muscle activation as a function of time
 
Details are given in KneeExt.m

Also available at https://www.mathworks.com/matlabcentral/fileexchange/65975-kneeext

This is a simple Matlab function to solve an equation of motion for simulating the concentric knee extension model published by Sust, Schmalz & Linnenbecker (1997a) and additionally includes muscle activation (Sust et al., 1997b).

Subject specific values of anthropometry, muscle activation, and force-velocity relationship as well as external loads can be manipulated to understand the effect of geometrical relations and muscle properties differing between individuals. 

*The original form of the function was provided by Sigrid Thaller and was revised by Harald Penasso, both, Institute of Sport Science, University of Graz, Graz, Austria. 

References:

- Sust, M., Schmalz, T., & Linnenbecker, S. (1997a). Relationship between distribution of muscle fibres and invariables of motion. Human Movement Science, 16(4), 533–546. https://doi.org/10.1016/S0167-9457(96)00063-2

 - Sust, M., Schmalz, T., Beyer, L., Rost, R., Hansen, E., & Weiss, T. (1997b). Assessment of isometric contractions performed with maximal subjective effort: corresponding results for EEG changes and force measurements. The International Journal of Neuroscience, 92(1–2), 103–118. https://doi.org/10.3109/00207459708986394
