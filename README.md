<h1 align="center">
Geometric Savitzky-Golay Filter for Noisy Rotation Measurements on SO(3), with Direct Angular Velocity and Acceleration Estimations
</h1>
<div align="center">
<h3>
<a href="https://research.tue.nl/en/persons/maarten-jongeneel">Maarten Jongeneel</a>,
<a href="https://www.tue.nl/en/research/researchers/alessandro-saccon/">Alessandro Saccon</a>
<br>
<br>
IEEE/RSJ Intelligent Robots and Systems conference (IROS), 2022
<br>
<br>
<a href="https://www.maartenjongeneel.nl">[Early Paper on HAL]</a>
</h3>
</div>


Introduction
============

The content of this repository is associated to the paper "Geometric Savitzky-Golay Filter for Noisy Rotation Measurements on SO(3), with Direct Angular Velocity and Acceleration Estimations". The objective for this paper was to create a version of the Savitzky-Golay filter that can directly be applied on a sequence of noisy rotation matrices, for example obtained from motion capture systems, to estimate the angular velocity and acceleration profiles. We showed in an example how this filter can be applied on a sequence of noisy rotation matrices representing a rotating frame in three-dimensional space. The code provided here can be used to reproduce the results of the paper. The code is publicly available.


Table of content
================
- [Overview](#overview)
- [Installation](#installation)
- [Usage of the scripts](#usage-of-the-scripts)
- [Contact](#contact)

# Overview

# Installation
The code of this repository is all written in MATLAB and can directly be pulled from this repository. 

# Usage of the scripts
The main script of this repository, ``SavitzkyGolaySO3.m``, is used to demonstrate the effectiveness of the filter and to create the figures as shown in the paper.

# Contact
In case you have questions or if you encountered an error, please contact us through the "Issues" functionality on GIT. Alternatively, you can send an email to <a href="mailto:m.j.jongeneel@tue.nl">Maarten Jongeneel</a>.