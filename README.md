<h1 align="center">
Geometric Savitzky-Golay Filter for Noisy Rotation Measurements on SO(3), with Direct Angular Velocity and Acceleration Estimations
</h1>
<div align="center">
<h3>
<a href="https://research.tue.nl/en/persons/maarten-jongeneel">Maarten Jongeneel</a>,
<a href="https://research.tue.nl/en/persons/alessandro-saccon">Alessandro Saccon</a>
<br>
<br>
IEEE/RSJ Intelligent Robots and Systems conference (IROS), 2022
<br>
<br>
<a href="https://hal.archives-ouvertes.fr/hal-03603826v1">[Early Paper on HAL]</a>
</h3>
</div>


Introduction
============

The content of this repository is associated to the paper "Geometric Savitzky-Golay Filter for Noisy Rotation Measurements on SO(3), with Direct Angular Velocity and Acceleration Estimations". The objective for this paper was to create a version of the Savitzky-Golay filter that can directly be applied on a sequence of noisy rotation matrices, for example obtained from motion capture systems, to estimate the angular velocity and acceleration profiles. We showed in an example how this filter can be applied on a sequence of noisy rotation matrices representing a rotating frame in three-dimensional space. The code provided here can be used to reproduce the results of the paper. The code is publicly available.


Table of content
================
- [Overview](#overview)
- [Installation](#installation)
- [Using the scripts](#using-the-scripts)
- [Contact](#contact)

# Overview
The code of this repository can be used to create random rotations in SO(3) and test the effectiveness of the filter. The function ``functions/sgolayfiltSO3.m`` can be used directly to apply the filter on a sequence of noisy rotation matrices and returns the estimated rotation matrices, angular velocities, and angular accelerations. 
# Installation
The code of this repository is all written in MATLAB and can directly be pulled from this repository. 

# Using the scripts
## Main script to reproduce paper results
The main script of this repository, ``SavitzkyGolaySO3.m``, is used to demonstrate the effectiveness of the filter and to create the figures as shown in the paper. In this script, one can change the following parameters:

```matlab
Fc = 1;            %Signal frequency                  [Hz]
a  = 2;            %Signal amplitude                  [deg]
te = 2;            %Signal length                     [s]
Fs = 1000;         %Sampling frequency fine grid      [Hz]
m  = 5;            %Down-sampling rate                [-]
sigma = 0.06;      %Standard deviation of added noise [rad]
n  = 15;           %Window size SG-filter             [-]
p  = 2;            %Savitzky Golay filter order       [-]
```
Where `Fc`, `a`, and `te` define the original signal frequency, amplitude, and length. `Fs`, `m`, and `sigma` determine the (fine-grid) sample frequency, the downsampling of this signal to obtain the noisy measurement, and the standard deviation of the added noise. Finally, the parameters `n` and `p` are the window-size and polynomial order of the Savitzky-Golay filter.

## Filter script to apply on own data
The `Functions` folder contains all the functions used in the main script. The most important function here is `sgolayfiltSO3.m`. This script is the main script of the filter and can be directly applied on a sequence of noisy rotation matrices. This script uses the following parameters:
```matlab
% INPUTS:    R         :Noisy sequence of rotation matrices, specified as
%                       a 3-by-3-by-N matrix containing N rotation matrices
%            p         :Polynomial order, specified as a positive integer,
%                       greater than the window size, n
%            n         :Window size, specified as a positive integer.
%            freq      :Sample frequency, specified as positive integer.
```
It returns the following outputs:
```matlab
% OUTPUTS:   R_est     :Estimated rotation matrices, specified as a
%                       3-by-3-by-(N-(2n+1)) matrix containing the
%                       estimated rotation matrices.
%            omg_est   :Estimated angular velocity, specified as a
%                       3-by-(N-(2n+1)) vector containing the estimated
%                       angular velocities at each time step.
%            domg_est  :Estimated angular acceleration, specified as a
%                       3-by-(N-(2n+1)) vector containing the estimated
%                       angular accelerations at each time step.
%            tf        :Time vector of the filtered signal
```
# Contact
In case you have questions or if you encountered an error, please contact us through the "Issues" functionality on GIT. Alternatively, you can send an email to <a href="mailto:m.j.jongeneel@tue.nl">Maarten Jongeneel</a>.