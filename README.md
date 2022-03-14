# EIS
This progam runs a physics-based electrochemical impedance model for proton-exchange-membrane fuel cell (PEMFC) in MATLAB.

Three approaches for modeling the electrochemical impedance response are compared using a case study of a PEMFC cathode with Tafel kinetics. 

The first is a transient simulation of a sine wave of small amplitude and numerical integration to obtain impedance. This approach is modeled after the experimental techniques for EIS. The frequency response of the system can be obtained from the time-domain signals using a Fourier transform. This approach is the most computationally intensive but requires no additional development from the transient electrochemical model. 

The second approach to modeling impedance is the transformation of the time-domain model equations into frequency-domain model equations using Laplace transforms and linearized about the steady-state model solution. This approach is quick and accurate, but is impractical for highly coupled, nonlinear systems of equations (such as a full cell PEMFC model, which would require numerical linearization). 

The third approach builds on the second approach by splitting up the frequency domain model equations into real and imaginary components such that the total number of equations is doubled. Two sets of governing equations are written, one for the real components of each variable and one for the imaginary components of each variable. This approach takes advantage of the Cauchy-Riemann equations, which allows the derivatives of complex variables to be written in terms of their real and imaginary components. Modern programming languages include a complex number data type, therefore the second approach is easier to implement and splitting up equations into real and imaginary components is unnecessary. However, this approach has been used in a number of important papers in the literature.
