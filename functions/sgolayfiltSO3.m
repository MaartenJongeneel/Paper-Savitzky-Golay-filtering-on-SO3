function [R_est,omg_est,domg_est,tf] = sgolayfiltSO3(R,p,n,freq)
% This function applies a Savitzky-Golay finite impulse response (FIR)
% smoothing filter of polynomial order p and frame length n to the data in
% the sequence of noisy rotation matrices R.
%
% INPUTS:    R         :Noisy sequence of rotation matrices, specified as
%                       a 3-by-3-by-N matrix containing N rotation matrices
%            p         :Polynomial order, specified as a positive integer,
%                       greater than the window size, n
%            n         :Window size, specified as a positive integer.
%            freq      :Sample frequency, specified as positive integer.
%
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
%
% M.J.Jongeneel, A.Saccon. Created 11-03-2022
%% ---------------- Savitzky-Golay Filtering on SO(3) ----------------- %%
%% Check inputed values
if nargin < 4
    error('Not enough input arguments.');
end
if p < 0 || p >= n || p ~= floor(p)
    error(['The polynomial order, p, must be a positive integer greater '...
        'than the window size, n.']);
end
if n < 0 || n ~= floor(n)
    error('The window size, n, should be a positive intiger.');
end
if freq <0 || freq~=floor(freq)
    error('The frequency, freq, should be a positive intiger.');
end
validateattributes(R,{'single','double'},...
    {'nonempty','real','3d','size',[3 3 NaN]})

%% Computed values
N = length(R(1,1,:));     %Number of samples in the sequence    [-]
dt = 1/freq;              %Time step lower sampled              [s]
te = N*dt;                %Total lenght of the sequence         [s]
ts = (0:dt:te);           %Signal time vector                   [s]
w = -n:n;                 %Window for Golay                     [-]
I = eye(3);               %Short hand notation                  [-]
tf = ts((n+1):(N-(n+1))); %Time vector filtered signal          [s]

%% Preallocate memory
R_est = NaN(3,3,N-length(w));
omg_est = NaN(3,N-length(w));
domg_est = NaN(3,N-length(w));

%% Savitzky-Golay
%For each time step (where we can apply the window)
cnt = 1;
for ii = (n+1):(N-(n+1))
    %Build matrix A and vector b based on window size w
    row = 1;
    for jj = 1:length(w)
        %Time difference between 0^th element and w(jj)^th element
        Dt = (ts(ii+w(jj))-ts(ii));
        %Determine row of A matrix
        Ajj = I;
        for kk = 1:p
            Ajj = cat(2,Ajj,(1/kk)*Dt^kk*I); %Concatenation based on order n
        end
        A(row:row+length(I)-1,:) = Ajj;
        b(row:row+length(I)-1,:) = vee(logm(R(:,:,ii+w(jj))/R(:,:,ii)));
        row = row+length(I); %Update to next row
    end
    %Solve the LS problem
    rho = (A'*A)\A'*b;
    
    %Obtain the coefficients of rho
    rho0 = rho(1:3);  rho1 = rho(4:6);  rho2 = rho(7:9);
    
    %Compute analytically the rotation matrices, ang. vel., and ang. acc.
    R_est(:,:,cnt) = expSO3(rho0)*R(:,:,ii);
    omg_est(:,cnt) = dexpSO3(rho0)*rho1;
    domg_est(:,cnt) = DdexpSO3(rho0,rho1)*rho1 +  dexpSO3(rho0)*rho2;
    
    %Update the index counter
    cnt = cnt+1;
end
end
