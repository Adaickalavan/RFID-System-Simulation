function [rpOut] = channel(rp)
%% This function computes a seriesRLC-M-parallelRLC antenna channel model 
% and plots its frequency response.

%% Symbolically derive transfer function of resonator channel
%Define symbolic values for circuit elements
syms s R1 R2 C1 C2 L1 L2 M;

%Compute circuit impedances symbolically in Matlab
Z1 = 1/(1/R2 + s*C2);
Z2 = Z1 + s*(L2-M);
Z3 = 1/(1/Z2 + 1/(s*M));
Z4 = Z3 + 1/(s*C1) + s*(L1-M) + R1;

%Compute the transfer function symbolically using voltage divider rule
Vz3_Vpcd = Z3/Z4;
Vpicc_Vz3 = Z1/Z2;
Vpicc_Vpcd = Vpicc_Vz3*Vz3_Vpcd;
Vpicc_Vpcd = simplifyFraction(Vpicc_Vpcd); %simplify symbolic fraction

%print the symbolic channel transfer function 
%{
fprintf('Symbolic channel transfer function =\n');
Vpicc_Vpcd
pretty(Vpicc_Vpcd) 
%}

%% Enter actual circuit values, and compute resonant frequency, and Q factor
%Enter actual values for circuit elements
Cr = 2.86402e-10; %pre-set capacitance values
Cc = 5.58015e-11; %pre-set capacitance values
k = rp.k; %0<=k<=1; k is usually comprised between 0.03 and 0.3

%Compute circuit parameters of series RLC circuit
%Resonant frequency fc
%2*pi*fres = 1/sqrt(LC)
%1/(2*pi*fres) = sqrt(LC)
%1/(2*pi*fres)^2 = LC
%1/(C*(2*pi*fres)^2) = L
Lr = 1/(Cr*(2*pi*rp.fresr)^2);
%Q factor
%Q = 1/(2*pi*fres*R*C)
%1/Q = 2*pi*fres*R*C
%1/(2*pi*fres*C*Q) = R
Rr = 1/(2*pi*rp.fresr*Cr*rp.Qr);
%transfer function of VL/Vin in series circuit
numHseriesS = [Lr*Cr 0 0];
denHseriesS = [Lr*Cr Rr*Cr 1];
HseriesS = tf(numHseriesS,denHseriesS); %Series RLC transfer function in s domain

%Compute circuit parameters of parallel RLC circuit
%Resonant frequency fc
%2*pi*fres = 1/sqrt(LC)
%1/(2*pi*fres) = sqrt(LC)
%1/(2*pi*fres)^2 = LC
%1/(C*(2*pi*fres)^2) = L
Lc = 1/(Cc*(2*pi*rp.fresc)^2);
%Q factor
%Q = 2*pi*fres*C*R
%Q/(2*pi*fres*C) = R
Rc = rp.Qc/(2*pi*rp.fresc*Cc);
%transfer function of IR/Iin in parallel circuit
numHparallelS = [Lc 0];
denHparallelS = [Rc*Lc*Cc Lc Rc];
HparallelS = tf(numHparallelS,denHparallelS); %Parallel RLC transfer function in s domain

%Print resonant frequency and Q factor of reader and card
fprintf('Reader: fresr=%5.2e, Qr=%3.1f\n', 1/(2*pi*sqrt(Lr*Cr)), 1/(2*pi*rp.fresr*Cr*Rr));
fprintf('Card  : fresc=%5.2e, Qc=%3.1f\n', 1/(2*pi*sqrt(Lc*Cc)), 2*pi*rp.fresc*Cc*Rc);
fprintf('Coupling factor: k=%2.2f\n', k)

%% Compute actual discrete time transfer function
%Substitue symbolic expressions with real values 
Mval = k*sqrt(Lr*Lc);
Vpicc_Vpcd_sub = subs(Vpicc_Vpcd, [R1,R2,C1,C2,L1,L2,M], [Rr,Rc,Cr,Cc,Lr,Lc,Mval]);

%Convert the symbolic expression to Matlab transfer function by obtaining 
%the numerator and denominator coefficients of the symbolic expression
[num, den] = numden(Vpicc_Vpcd_sub);
numChannelS = sym2poly(num);
denChannelS = sym2poly(den);
HchannelS = tf(numChannelS,denChannelS); %Channel transfer function in s domain

%Discretize the numerator and denominator coefficients of channel
HchannelZ = c2d(HchannelS,rp.Ts,'matched');
[numChannelZ, denChannelZ] = tfdata(HchannelZ,'V');

%% Map bandpass channel transfer function into its lowpass equivalent
% Method 1 ----------------------------------------------------------------
%Follows the method described on page 143 of "Simulation of Communication
%Systems" by Jeruchim, 2nd ed, 2000

%Decompose transfer function using partial fraction
[r,p,k] = residue(numChannelS,denChannelS);
if ~isempty(k)
    error('Error in partial fraction expansion of Hchannel');
end

%Discard the poles that are located on the negative-frequency half-plane
r_pos = r(imag(p)>=0);
p_pos = p(imag(p)>=0);

%Shift poles located on the positive-frequency half-plane to the zero axis
%by substituting s -> s+jw
w0 = 2*pi*rp.fc;
p_pos_baseband = p_pos - 1j*w0;

%The lowpass equivalent transfer function
[numChannelBbS,denChannelBbS] = residue(r_pos,p_pos_baseband,k);
HchannelBbS = tf(numChannelBbS,denChannelBbS);

%Discretize the numerator and denominator coefficients of lowpass equivalent channel
HchannelBbZ_Sshift = c2d(HchannelBbS,rp.Ts,'matched');
[numChannelBbZ_Sshift,denChannelBbZ_Sshift] = tfdata(HchannelBbZ_Sshift,'V');
%--------------------------------------------------------------------------

% Method 2 ----------------------------------------------------------------
%Follows the method described on page 50 of "Signal Processing techniques 
%for software radios" by Farhang, 2nd ed, 2010. Example of Matlab
%implementation is given on page 364 of this book by line
%"c=c.*exp(-j*2*pi*[0:length(c)-1]’*Ts*fc);" where c is the channel, Ts is
%the sample time, and fc is the carrier frequency.

%Substitute z^x -> z^x*e^(j*2*pi*fc*Ts*x) in the z-domain transfer function 
%of the passband channel. Here, fc is the carrier frequency to donwconvert 
%and Ts is the sampling time
numChannelBbZ_Zshift = numChannelZ.*exp(1j*2*pi*rp.fc*rp.Ts*(length(numChannelZ)-1:-1:0));
denChannelBbZ_Zshift = denChannelZ.*exp(1j*2*pi*rp.fc*rp.Ts*(length(denChannelZ)-1:-1:0));
HchannelBbZ_Zshift = tf(numChannelBbZ_Zshift,denChannelBbZ_Zshift,rp.Ts);
%--------------------------------------------------------------------------
%% Plot the (i) pole-zero plot of bandpass and lowpass equivalent channels, 
% and (ii) frequency response of bandpass and lowpass equivalent channels

if rp.plotTrue == 1
    % Specify position of figure on screen. rect = [left, bottom, width, height]
    scrsz = get(0,'screensize'); %Get the screensize to specify figure size and location later
    figure('OuterPosition',[1 40 scrsz(3) scrsz(4)-40])          
    %----------------------------------------------------------------------
    handle1 = subaxis(2,2,1, 'Spacing', 0.06, 'Padding', 0.01, 'Margin', 0.03);
    ss = 1; legStr = {};
    hold on
    pzmap(HchannelS,'-b'); legStr{ss} = 'HchannelS'; ss=ss+1; 
    pzmap(HchannelBbS,'-g'); legStr{ss} = 'HchannelBbS'; ss=ss+1; 
    hold off
    legend(legStr,'location','northeast','interpreter','none'); 
    %----------------------------------------------------------------------
    handle2 = subaxis(2,2,2, 'Spacing', 0.06, 'Padding', 0.01, 'Margin', 0.03);
    ss = 1; legStr = {};
    hold on
    pzplot(HchannelZ,'b'); legStr{ss} = 'HchannelZ'; ss=ss+1; 
    pzplot(HchannelBbZ_Sshift,'g'); legStr{ss} = 'HchannelBbZ_Sshift'; ss=ss+1; 
    pzplot(HchannelBbZ_Zshift,'r'); legStr{ss} = 'HchannelBbZ_Zshift'; ss=ss+1; 
    hold off
    legend(legStr,'location','southwest','interpreter','none'); 
    %----------------------------------------------------------------------
    handle3 = subaxis(2,2,3, 'Spacing', 0.06, 'Padding', 0.01, 'Margin', 0.03);
    fs = 1/rp.Ts;
    freqAxis = linspace(-fs,fs,200);
    [freqResHseriesS,~] = freqs(numHseriesS,denHseriesS,2*pi*freqAxis);
    [freqResHparallelS,~] = freqs(numHparallelS,denHparallelS,2*pi*freqAxis);
    [freqResChS,~] = freqs(numChannelS,denChannelS,2*pi*freqAxis);
    [freqResChZ,~] = freqz(numChannelZ,denChannelZ,2*pi*freqAxis/fs);
    [freqResChBbS,~] = freqs(numChannelBbS,denChannelBbS,2*pi*freqAxis);
    [freqResChBbZ_Sshift,~] = freqz(numChannelBbZ_Sshift,denChannelBbZ_Sshift,2*pi*freqAxis/fs);
    [freqResChBbZ_Zshift,~] = freqz(numChannelBbZ_Zshift,denChannelBbZ_Zshift,2*pi*freqAxis/fs);
    ss = 1; legStr = {};
    hold on
    plot(freqAxis,20*log10(abs(freqResHseriesS)),'-c'); legStr{ss} = 'HseriesS'; ss=ss+1; 
    plot(freqAxis,20*log10(abs(freqResHparallelS)),'-m'); legStr{ss} = 'HparallelS'; ss=ss+1;
    plot(freqAxis,20*log10(abs(freqResChS)),'-b'); legStr{ss} = 'ChS'; ss=ss+1;
    plot(freqAxis,20*log10(abs(freqResChZ)),'--b'); legStr{ss} = 'ChZ'; ss=ss+1;
    plot(freqAxis,20*log10(abs(freqResChBbS)),'-g'); legStr{ss} = 'ChBbS'; ss=ss+1;
    plot(freqAxis,20*log10(abs(freqResChBbZ_Sshift)),'--g'); legStr{ss} = 'ChBbZ_Sshift'; ss=ss+1;  
    plot(freqAxis,20*log10(abs(freqResChBbZ_Zshift)),'--r'); legStr{ss} = 'ChBbZ_Zshift'; ss=ss+1;  
    hold off
    legend(legStr,'location','southwest','interpreter','none'); 
    xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)'); grid;
    ylim([-140 20]) %[xmin xmax ymin y max] = Set x and y axis limits
    %----------------------------------------------------------------------
end

%% Function output
rpOut.numChannelZ = numChannelZ;
rpOut.denChannelZ = denChannelZ;
rpOut.numChannelBbZ_Sshift = numChannelBbZ_Sshift;
rpOut.denChannelBbZ_Sshift = denChannelBbZ_Sshift;
rpOut.numChannelBbZ_Zshift = numChannelBbZ_Zshift;
rpOut.denChannelBbZ_Zshift = denChannelBbZ_Zshift;

end
