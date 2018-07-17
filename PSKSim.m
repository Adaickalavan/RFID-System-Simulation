%Code for simulating the modulation and transmission PSK 
clear all
close all
%Print messages to screen
fprintf('***********************\n');

%--------------------------------------------------------------------------
%% Transmitter
%--------------------------------------------------------------------------
%Parameters for us to set
N = 800; %Number of symbols transmitted
M = 16; %M-ary PSK order
fc = 13.56e6; %Carrier frequency, units: Hz
dataRate = fc/4; %Symbol rate or Baud rate
plotTrue = 1; %enable plot
nSampSim = 32; %Number of samples per symbol for simulating transmitter and channel
nSamp = 2; %Number of samples per symbol for simulating receiver
TsSim = 1/(nSampSim*dataRate); %Sample time for simulating transmitter and channel
Ts = 1/(nSamp*dataRate); %Sample time for simulating receiver
%--------------------------------------------------------------------------
%Ensure random number generator is repeatable
rng(1,'twister'); %initialize the random number generator to make the results in this example repeatable
savedRNG = rng; %Save the generator settings in a structure savedRNG
rng(savedRNG); %return the generator to the original state stored in savedRNG
%--------------------------------------------------------------------------
%Generate transmitted signal
[S, cp2S, S2cp, preamble] = constellation(M); %initialize constellation setting
cp = randi([0,M-1],N,1); %Get N constellation points to send
cpDiffpre = 8; %constellation point of 32degree
cpDiff = mod(cpDiffpre + cumsum(cp), M); %Diferential encoding
cpDiffpre = cpDiff(end); %Reset previous sent constellation point

%Form one frame by adding SOC preamble, EOC, silence
silenceStartLen = 1000;
silenceEndLen = 100;
silenceStart = S(1)*ones(silenceStartLen,1); %unmodulated RF carrier with a NPV of 0°.
silenceEnd = S(1)*ones(silenceEndLen,1); %unmodulated RF carrier with a NPV of 0°.
EOC = S(M+1)*ones(8,1); %sequence of 8 NPVs of -180° and silence.
frame = [silenceStart; cp2S(preamble+1).'; cp2S(cpDiff+1).'; EOC; silenceEnd]; %a sequence of 8 NPVs of -180°.

%Pulse shaping
pT = ones(nSampSim,1); %NRZ transmit filter
txSim = conv(pT, upsample(frame,nSampSim)); %Complex baseband pulse-shaping by NRZ filter 
txSim = txSim(1:end-nSampSim+1); %Truncate the unnecessary trailing zero-padded convolution at the end of tx
tx = txSim(1:nSampSim/nSamp:end); %Form downsampled tx for plotting
tSim = (0:TsSim:TsSim*(length(txSim)-1)).'; %Form simulation time vector
t = tSim(1:nSampSim/nSamp:end); %Form downsampled simulation time vector

%Build character(8symbols) dividers for plotting
charDivSilence = [zeros(1,silenceStartLen*nSamp-1), 2];
charDivSOC(8*nSamp:8*nSamp:length(preamble)*nSamp) = 1;
charDivSOC(length(preamble)*nSamp) = 2;
charDivSig(8*nSamp:8*nSamp:N*nSamp) = 1;
charDivSig(N*nSamp) = 2;
charDiv = [charDivSilence, charDivSOC, charDivSig];

%--------------------------------------------------------------------------
%% Channel
%--------------------------------------------------------------------------
fprintf('-----------------------\n');
fprintf('Channel\n')

k = 0.3; % 0<=k<=1; k is usually comprised between 0.03 and 0.3
Qr = 4;
Qc = 3;
fresr = 13.56e6;
fresc = 14e6;

%Obtain channel transfer function
rpChannel.plotTrue = plotTrue;
rpChannel.Ts = TsSim;
rpChannel.fc = fc;
rpChannel.fresr = fresr;
rpChannel.fresc = fresc;
rpChannel.k = k;
rpChannel.Qr = Qr;
rpChannel.Qc = Qc;
rpOutChannel = channel(rpChannel);

%Filter transmitted signal through passband and lowpass equivalent channel
rxPb = filter(rpOutChannel.numChannelZ,rpOutChannel.denChannelZ,real(txSim.*exp(1j*2*pi*fc*tSim)));
rxBb_Sshift = filter(rpOutChannel.numChannelBbZ_Sshift,rpOutChannel.denChannelBbZ_Sshift,txSim);
rxBb_Zshift = filter(rpOutChannel.numChannelBbZ_Zshift,rpOutChannel.denChannelBbZ_Zshift,txSim);

%--------------------------------------------------------------------------
if plotTrue
    %Specify position of figure on screen. rect = [left, bottom, width, height]
    scrsz = get(0,'screensize'); %Get the screensize to specify figure size and location later
    figure('OuterPosition',[1 40 scrsz(3) scrsz(4)-40])          
    %----------------------------------------------------------------------
    ss = 1; legStr = {};
    handle1 = subaxis(2,2,1, 'Spacing', 0.06, 'Padding', 0.01, 'Margin', 0.03);
    hold on
    plot(real(txSim.*exp(1j*2*pi*fc*tSim)),'-k.'); legStr{ss} = 'real(txCarrier)'; ss=ss+1;
    plot(rxPb,'-b.'); legStr{ss} = 'rxPb'; ss=ss+1;
    plot(real(rxBb_Sshift.*exp(1j*2*pi*fc*tSim)),'-g.'); legStr{ss} = 'real(rxBb_Sshift*exp(jwt))'; ss=ss+1;
    plot(real(rxBb_Zshift.*exp(1j*2*pi*fc*tSim)),'-r.'); legStr{ss} = 'real(rxBb_Zshift*exp(jwt))'; ss=ss+1;
    hold off
    legend(legStr,'location','NorthWest','interpreter','none');
    xlabel('Time'); ylabel('Magnitude'); title('Received passband waveshape');grid;
    %----------------------------------------------------------------------
    ss = 1; legStr = {};
    handle2 = subaxis(2,2,2, 'Spacing', 0.06, 'Padding', 0.01, 'Margin', 0.03);
    hold on
    plot(angle(txSim),'-k.'); legStr{ss} = 'angle(txSim)'; ss=ss+1;
    plot(angle(rxBb_Sshift),'-g.'); legStr{ss} = 'angle(rxBb_Sshift)'; ss=ss+1;
    plot(angle(rxBb_Zshift),'-r.'); legStr{ss} = 'angle(rxBb_Zshift)'; ss=ss+1;
    hold off
    legend(legStr,'location','NorthWest','interpreter','none');
    xlabel('Time'); ylabel('Degree (rad)'); title('Received angle (upsampled)'); grid;
    %----------------------------------------------------------------------
    ss = 1; legStr = {};
    handle3 = subaxis(2,2,3, 'Spacing', 0.06, 'Padding', 0.01, 'Margin', 0.03);
    hold on
    plot(angle(tx),'-k.'); legStr{ss} = 'angle(tx)'; ss=ss+1;
    plot(angle(rxBb_Sshift(1:nSampSim/nSamp:end)),'-g.'); legStr{ss} = 'angle(rxBb_Sshift)'; ss=ss+1;
    plot(angle(rxBb_Zshift(1:nSampSim/nSamp:end)),'-r.'); legStr{ss} = 'angle(rxBb_Zshift)'; ss=ss+1;
    hold off
    legend(legStr,'location','NorthWest','interpreter','none');
    xlabel('Time'); ylabel('Degree (rad)'); title('Received angle (downsampled)');grid;
    %----------------------------------------------------------------------
end

%--------------------------------------------------------------------------
%Assign complex baseband signal output from lowpass equivalent channel to
%received signal for baseband signal processing
rxChOut = rxBb_Zshift;
%Downsample signal to desired samples per symbol
rxChOut = rxChOut(1:nSampSim/nSamp:end);
fprintf('-----------------------\n');

%--------------------------------------------------------------------------
%% Phase locked loop
%--------------------------------------------------------------------------
%PLL - Degree-based PLL
fprintf('-----------------------\n');
fprintf('Phase-Locked Loop\n')

rpPLLdegree.phiTx = angle(tx.')*180/pi;
rpPLLdegree.phiRx = angle(rxChOut.')*180/pi;
rpPLLdegree.Ndiv = 8; %frequency divider
rpPLLdegree.Ktdc = (2^9 - 1)/180; %Map 180 degrees to 511 codes. Ktdc = 2.84
rpPLLdegree.Kdco = rpPLLdegree.Ndiv*3e3; %DCO gain, units: s^-1
rpPLLdegree.a = 1/4; %(1/4);
rpPLLdegree.b = (1/2048)/Ts; %(1/2048)/Ts;
rpPLLdegree.fDCOn0 = 108.0e6; %108.48e6; rpPLLdegree.Ndiv*13.435e6; %initial n-multiple frequency of DCO
rpPLLdegree.phiDCO0 = 14; %Initial phase of DCO in degrees within [-180,180)
rpPLLdegree.plotTrue = plotTrue;

rpPLLdegree.quantizationEn = 1;
rpPLLdegree.saturationEn = 1;
rpPLLdegree.phaseWrapEn = 1;
rpPLLdegree.Ts = Ts;
rpPLLdegree.fc = fc;

[tdcOut, rpOutPLLdegree] = ADPLL_degree(rpPLLdegree);
tdcOutDegree = tdcOut/rpOutPLLdegree.Ktdc;

%Received signal's name changed after being processed by PLL
rxPLLOut = tdcOutDegree;
fprintf('-----------------------\n');

%--------------------------------------------------------------------------
%% Digital baseband receiver   
%--------------------------------------------------------------------------
fprintf('-----------------------\n');
fprintf('Baseband receiver\n')

%For equalizer ------------------------------------------------------------
training = 1:(length(preamble))*nSamp;
update = 1:(length(preamble))*nSamp;
% update = 1:(length(preamble)+length(cpDiff)+length(EOC))*nSamp;

trainSeq = [preamble; cpDiff; S2cp(M+1)*ones(8,1)];
trainSeq = trainSeq(:,ones(1,nSamp)).';
trainSeq = trainSeq(:);

kFwd_pos = 0:-1:-2;
kFwd = length(kFwd_pos);
kBwd = 2;

% tapFI = zeros(kFwd,1);
tapBI = zeros(kBwd,1);
equOut = inf(length(rxPLLOut),1);
dcs = inf(length(rxPLLOut),1);
err = zeros(length(rxPLLOut),1);

mu = 0.25; %convergence factor (step)  (0 < mu < 1); 
L = 2; %data reuse factor (L=0 -> NLMS; L=1 -> BNLMS; ...)
gamma = 1; %small positive constant to avoid singularity
bias = 0;
XapMat = zeros(bias+kFwd+kBwd,L+1);
tapW = zeros(bias+kFwd+kBwd,1);
dapVec = zeros(L+1,1);

delay = 1; %Default value of delay = 1, which has no effect on AP equalizer
XapMatTemp = zeros(bias+kFwd+kBwd,delay*L+1);
dapVecTemp = zeros(delay*L+1,1);

fracEn = 1;
decEn = nSamp-1;
errNum = 0;
dcsErr = inf(length(rxPLLOut),1);
dcsErrLoc = inf(length(rxPLLOut),1);

% Note: Simulation works for k=[0.03,0.3,0.9], Q1=[1,4], Q2=[1,4]
% kFwd=3
% kBwd=2
% mu=0.25
% L=2
% gamma=1
% bias=0
% fracEn=1
% delay=1

%--------------------------------------------------------------------------
%For Symbol synchronizer --------------------------------------------------
tapMatchOut = zeros(4*nSamp,1);
tapMatchOutPlot = zeros(length(rxPLLOut),1);
refInd = 2*nSamp + 1;
startEqu = 0;
trainSeqInd = 2*nSamp - 1;
startMatchInd = silenceStartLen*nSamp - 50;
%--------------------------------------------------------------------------

for ii = startMatchInd:length(rxPLLOut)

    evalInd = ii - refInd + 1;

    tapMatchOutPlot(ii) = sum(rxPLLOut(ii:-1:ii-2*nSamp+1));
    tapMatchOut = [tapMatchOutPlot(ii); tapMatchOut(1:end-1)];
    
    %Search for Sync word
    if  tapMatchOut(refInd) > tapMatchOut(refInd+1) && ...
        tapMatchOut(refInd) > tapMatchOut(refInd-1) && ...
        tapMatchOut(refInd) > 20 && ...
        startEqu == 0
        
        startEqu = 1;
        startEquInd = evalInd;
    end
    
    %Start equalizer after Sync word is found 
    if startEqu == 1
        %Form forward long filter 
        tapFI = rxPLLOut(evalInd+kFwd_pos);
        %Form forward and backward filter input 
        tapFilter = [ones(bias,1); tapFI; tapBI];
        %Update decEn counter
        decEn = decEn + 1;
        
        %Increment training sequence index
        trainSeqInd = trainSeqInd+1;
        
        if (fracEn == 1 && decEn == nSamp) || (fracEn == 0)
            %Form forward and backward filter input matrix
            XapMatTemp = [tapFilter, XapMatTemp(:,1:end-1)];   
            XapMat = XapMatTemp(:,1:delay:end);
            %Compute equalizer output vector
            yapVec = XapMat.'*conj(tapW);    
            equOut(evalInd) = yapVec(1);
            %Make a decision
            if ismember(trainSeqInd, training)
                dcs(evalInd) = cp2S(trainSeq(trainSeqInd) + 1);
            else
                [~, ind] = min(abs(exp(1j*equOut(evalInd)*pi/180) - S));
                dcs(evalInd) = S(ind);
            end
            curSymb = angle(dcs(evalInd))*180/pi;
            %Form past decision vector 
            dapVecTemp = [curSymb; dapVecTemp(1:end-1)];
            dapVec = dapVecTemp(1:delay:end);
            %Form error vector
            eapVec = dapVec - yapVec;
            %Affine projection filter coefficient adaptation
            if ismember(trainSeqInd, update)
                tapW = tapW + mu*XapMat*((XapMat'*XapMat + gamma*eye(L+1))\conj(eapVec));
            end
            %Update feedback filter
            if kBwd ~= 0
                tapBI = [curSymb; tapBI(1:end-1)];
            end
            if trainSeqInd <= nSamp*(length(preamble)+length(cpDiff));
                err(evalInd) = eapVec(1);
            end
        end

        %Check for end of frame
        if S(M+1) == dcs(evalInd) || trainSeqInd >= length(trainSeq)- nSamp*7
            fprintf('Counted symbols: %5.2f\n', trainSeqInd/nSamp)
            fprintf('End of frame detected\n')
            break;
        end
        
        if decEn == nSamp            
            if trainSeqInd > length(trainSeq)
                fprintf('@@@@ Erroneous break @@@@\n')
                break;
            end
            %Error check
            if trainSeqInd > training(end)
                if cp2S(trainSeq(trainSeqInd)+1) ~= dcs(evalInd)    
                    errNum = errNum + 1;
                    dcsErr(evalInd) = angle(dcs(evalInd))*180/pi;
                end
            end
            %Reset decEn counter
            decEn = 0;
        end        
    end
end

% decoded_sig = preamble - [0; preamble(1:end-1)]; %Differentially decode the signal      
% est_sig = mod(decoded_sig, M); %Differential decoding    
fprintf('-----------------------\n');

%--------------------------------------------------------------------------
%% Print results
fprintf('-----------------------\n');
fprintf('Error check\n')
fprintf('Num of symb err: %d\n',errNum);
fprintf('-----------------------\n');

%--------------------------------------------------------------------------
%% Plot graphs 
%--------------------------------------------------------------------------
if plotTrue
    % Specify position of figure on screen. rect = [left, bottom, width, height]
    scrsz = get(0,'screensize'); %Get the screensize to specify figure size and location later
    figure('OuterPosition',[1 40 scrsz(3) scrsz(4)-40])          
    %----------------------------------------------------------------------
    ss = 1; legStr = {};
    handle1 = subaxis(4,2,[1:4], 'Spacing', 0.06, 'Padding', 0.01, 'Margin', 0.03);
    hold on
    plot(angle(tx)*180/pi,'-k.'); legStr{ss} = 'tx'; ss = ss + 1;
    plot(angle(rxChOut)*180/pi,'-c.'); legStr{ss} = 'rxChOut'; ss = ss + 1;
    plot(rxPLLOut,'-b.'); legStr{ss} = 'rxPLLOut'; ss = ss + 1;
    plot(tapMatchOutPlot,'-m.'); legStr{ss} = 'matchOut'; ss = ss + 1;
    plot(equOut,'-go'); legStr{ss} = 'equOut'; ss = ss + 1; %TDC output in degrees
    if errNum ~= 0
        plot(dcsErr,'ro','linewidth',3); legStr{ss} = 'dcsErr'; ss = ss + 1;
    end
    stem(startMatchInd,120,'-k*');
    stem(startEquInd,120,'-k*');
    stem(100*charDiv,'-k.','linewidth',2)
    stem(-100*charDiv,'-k.','linewidth',2)
    legend(legStr)
    xlabel('Time'); ylabel('Code levels'); grid;
    hold off
    %----------------------------------------------------------------------
    ss = 1; legStr = {};
    handle2 = subaxis(4,2,[6,8], 'Spacing', 0.06, 'Padding', 0.01, 'Margin', 0.03);
    plot(err.^2,'-b.'); legStr{ss} = 'err^2'; ss = ss + 1;
    hold on
    stem(100*charDiv,'-k.','linewidth',2)
    stem(-100*charDiv,'-k.','linewidth',2)
    legend(legStr)
    xlabel('Time'); ylabel('Magnitude squared error'); grid;
    hold off
    %----------------------------------------------------------------------    
    ss = 1; legStr = {};
    handle3 = subaxis(4,2,[5,7], 'Spacing', 0.06, 'Padding', 0.01, 'Margin', 0.03);
    hold on
    circ = exp(1j*(0:pi/50:2*pi)); %Obtain coordinates of unit-radius circle
    plot(imag(circ),real(circ),'-.k'); legStr{ss} = 'unit radius circle'; ss = ss + 1; %Plot unit-radius circle for reference
    plot(S,'r*'); legStr{ss} = 'S-constellation points'; ss = ss + 1; %Available constellation points
    plot(tx,'ok'); legStr{ss} = 'tx'; ss = ss + 1; %Transmitted complex signal
    plot(rxChOut,'oc'); legStr{ss} = 'rxChOut'; ss = ss + 1; %Received complex signal after the channel
    plot(exp(1j*rxPLLOut(silenceStartLen*nSamp-10:end)*pi/180),'b.'); legStr{ss} = 'rxPLLOut'; ss = ss + 1; %Transmitted complex signal
    plot(exp(1j*equOut(end-800:end)*pi/180),'og'); legStr{ss} = 'equOut'; ss = ss + 1; %Transmitted complex signal
    
    axis([-2 4 -1.5 2.5]) %[xmin xmax ymin y max] = Set x and y axis limits
    legend(legStr,'Location','NorthWest')
    xlabel('In-Phase'); ylabel('Qudrature-Phase'); grid;
    hold off
    %----------------------------------------------------------------------
    linkaxes([handle1 handle2], 'x')
    %----------------------------------------------------------------------
    
end
%--------------------------------------------------------------------------
rpOut.errNum = errNum;
fprintf('***********************\n');