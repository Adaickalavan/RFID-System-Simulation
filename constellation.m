function [S, cp2S, S2cp, preamble] = constellation(M)

%Store constellations for 16PSK
if M == 16
    phaseRange = 60; %unit:degree. Total phase variation range in PSK constellation
    EPI = phaseRange/(M-1); %elementary phase interval
    
    %Physical constellation points
    S(12+1) = exp(1j*8*EPI*pi/180);
    S(13+1) = exp(1j*7*EPI*pi/180);
    S(15+1) = exp(1j*6*EPI*pi/180);
    S(14+1) = exp(1j*5*EPI*pi/180);
    S(10+1) = exp(1j*4*EPI*pi/180);
    S(11+1) = exp(1j*3*EPI*pi/180);
    S(9+1) = exp(1j*2*EPI*pi/180);
    S(8+1) = exp(1j*1*EPI*pi/180);
    S(0+1) = exp(1j*0*EPI*pi/180);
    S(1+1) = exp(1j*-1*EPI*pi/180);
    S(3+1) = exp(1j*-2*EPI*pi/180);
    S(2+1) = exp(1j*-3*EPI*pi/180);
    S(6+1) = exp(1j*-4*EPI*pi/180);
    S(7+1) = exp(1j*-5*EPI*pi/180);
    S(5+1) = exp(1j*-6*EPI*pi/180);
    S(4+1) = exp(1j*-7*EPI*pi/180);   
    S(M+1) = exp(1j*-pi); %Add -180degree for EOC
    
    %Relate constellation point numbers to physical constellation points 
    cp2S(1) = S(1);
    cp2S(2) = S(9);
    cp2S(3) = S(10);
    cp2S(4) = S(12);
    cp2S(5) = S(11);
    cp2S(6) = S(15);
    cp2S(7) = S(16);
    cp2S(8) = S(14);
    cp2S(9) = S(13);
    cp2S(10) = S(5);
    cp2S(11) = S(6);
    cp2S(12) = S(8);
    cp2S(13) = S(7);
    cp2S(14) = S(3);
    cp2S(15) = S(4);
    cp2S(16) = S(2);   
    cp2S(M+1) = S(M+1); %Mapping for EOC
    
    %Relate physical constellation points to constellation point numbers 
    S2cp(1) = 0;
    S2cp(2) = 15;
    S2cp(3) = 13;
    S2cp(4) = 14;
    S2cp(5) = 9;
    S2cp(6) = 10;
    S2cp(7) = 12;
    S2cp(8) = 11;
    S2cp(9) = 1;
    S2cp(10) = 2;
    S2cp(11) = 4;
    S2cp(12) = 3;
    S2cp(13) = 8;
    S2cp(14) = 7;
    S2cp(15) = 5;
    S2cp(16) = 6;   
    S2cp(M+1) = M; %Mapping for EOC
    
    %Form the preamble in terms of constellation point numbers 
    preamble = [7; 7; 9; 9; 7; 7; 9; 9; 7; 7; 9; 9; 7; 7; 9; 9; 7; 7; 9; 9;
                7; 7; 9; 9; 7; 7; 9; 9; 7; 7; 9; 9; 7; 7; 9; 9; 7; 7; 9; 9;
                7; 7; 9; 9; 7; 9; 7; 9; 8; 8; 9; 2; 13; 8; 0; 4; 14; 7; 8; 9;
                15; 11; 13; 7; 4; 11; 10; 6; 13; 11; 5; 1; 4; 14; 12; 12; 4; 7; 11; 9; 
                8; 2; 9; 12; 3; 12; 7; 4; 2; 11; 8; 13; 1; 1; 15; 13; 4; 1; 9; 2; 
                11; 1; 12; 7; 8; 2; 3; 5; 10; 15; 8; 12; 2; 14; 13; 8; 9; 10; 9; 2; 
                9; 8; 9; 8; 9; 1; 6; 4; 0; 5; 8; 1; 13; 3; 11; 6; 7; 10; 9; 8];
    
else
    error('Constellation is unavailable')
end

%{
%------------------------------------------------------------------------------
%Ensure E[|S|^2] = 1 
fprintf('signal_Var = %g\n\n', norm(S)^2/M);

%Get the screensize to specify figure size and location
scrsz = get(0,'ScreenSize'); 
%Plot the constellation diagram
figure('Outerposition',[2*scrsz(3)/4 scrsz(4)/2 scrsz(3)/4 scrsz(4)/2]);
hold on
plot(S,'r*');
circ = exp(1j*(0:pi/50:2*pi)); %Obtain coordinates of unit-radius circle
plot(imag(circ),real(circ),'-.k'); %Plot unit-radius circle for reference
hold off
grid on
axis([-1.5 1.5 -1.5 1.5]) %[xmin xmax ymin y max] = Set x and y axis limits
title('Signal constellation'),
xlabel('In-Phase'),ylabel('Qudrature')  
%------------------------------------------------------------------------------
%}

%{
Constellation map

16PSK
Degree  Gray Bits    S    cp2S
  32      1100      12      8
  28      1101      13      7 
  24      1111      15      6
  20      1110      14      5
  16      1010      10      4 
  12      1011      11      3
   8      1001       9      2
   4      1000       8      1
   0      0000       0      0
  -4      0001       1     15
  -8      0011       3     14
 -12      0010       2     13
 -16      0110       6     12
 -20      0111       7     11
 -24      0101       5     10
 -28      0100       4      9

%}
end
