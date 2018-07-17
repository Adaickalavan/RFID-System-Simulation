function [tdcOut, rpOut] = ADPLL_degree(rp)

%ADPLL signals
ntot = length(rp.phiRx);
fDCOn = zeros(ntot,1);                   %ADPLL frequency
fDCO = zeros(ntot,1);                    %ADPLL frequency
phiDCO = [rp.phiDCO0; zeros(ntot,1)];    %ADPLL output phase (degree)  
tdcOut = zeros(ntot,1);                  %TDC output
LFprop = zeros(ntot,1);                  %Loop filter proportional path value
LFint = zeros(ntot,1);                   %Loop filter integral path value
LFOut = zeros(ntot,1);                   %DCO input

%% Start ADPLL
for ii = 1:ntot
    
    %TDC ------------------------------------------------------------------
    phiErr = rp.phiRx(ii) - phiDCO(ii); %Phase detector
    if rp.phaseWrapEn == 1
        phiErrSawTooth = rem(phiErr+sign(phiErr)*180,360) - sign(phiErr)*180; %sawtooth phase detector with range [-180, +180]
    else
        phiErrSawTooth = phiErr;
    end
    tdcOut(ii) = rp.Ktdc*phiErrSawTooth; %Gain of phase detector
    if rp.quantizationEn == 1
        tdcOut(ii) = round(tdcOut(ii)); %round to nearest integer
    end
    if rp.saturationEn == 1
        %TDC output operates on 10 bits
        tdcOut(ii) = max(min(tdcOut(ii),511),-511); %saturation of codes    
    end
    %Loop Filter ----------------------------------------------------------
    LFprop(ii) = rp.a*tdcOut(ii);
    if (ii == 1)
        LFint(ii) = rp.b*rp.Ts*tdcOut(ii); %initial integral path value        
    else
        LFint(ii) = rp.b*rp.Ts*tdcOut(ii) + LFint(ii-1);
    end
    LFOut(ii) = LFprop(ii)+ LFint(ii);
    if rp.quantizationEn == 1
        LFOut(ii) = round(LFOut(ii)); %round to nearest integer
    end
    if rp.saturationEn == 1
        %DCO input operates on 7 bits
        LFOut(ii) = max(min(LFOut(ii),63),-63); %saturation of codes    
    end
    %DCO-------------------------------------------------------------------
    fDCOn(ii) = rp.fDCOn0 + rp.Kdco*LFOut(ii); %n-multiple frequency of DCO
    fDCO(ii) = fDCOn(ii)/rp.Ndiv; %divided frequency output of DCO   
    phiDCO(ii+1) = 360*fDCO(ii)*rp.Ts + phiDCO(ii); %compute phase of DCO output
    if rp.phaseWrapEn == 1
        phiDCO(ii+1) = rem(phiDCO(ii+1)+sign(phiDCO(ii+1))*180,360) - sign(phiDCO(ii+1))*180; %ensure phase remains within [-180, 180) 
    end
    %----------------------------------------------------------------------

end

%Return the Ktdc value in case LUT was used
rpOut.Ktdc = rp.Ktdc;

%% Plot graphs
if rp.plotTrue
    scrsz = get(0,'screensize'); %Get the screensize to specify figure size and location later
    figTitle = 'Degree-based PLL';          
    figure('OuterPosition',[1 40 scrsz(3) scrsz(4)-40],'Name',figTitle)          

    ss = 1; legStr = {};
    handle1 = subaxis(3,1,1, 'Spacing', 0.06, 'Padding', 0.01, 'Margin', 0.03);
    hold on
    plot([1,ntot],[rp.fc, rp.fc],'--k'); legStr{ss} = 'fc'; ss = ss + 1; %TDC output in degrees
    plot(fDCO,'-b.'); legStr{ss} = 'fDCO'; ss = ss + 1; %TDC output in degrees
    legend(legStr); xlabel('Time'); ylabel('Frequency (Hz)'); title('Frequency of DCO output'); grid;
    hold off

    ss = 1; legStr = {};
    handle2 = subaxis(3,1,2, 'Spacing', 0.06, 'Padding', 0.01, 'Margin', 0.03);
    hold on
    plot(rp.phiTx,'-c.'); legStr{ss} = 'phiTx'; ss = ss + 1; %TDC output in degrees
    plot(rp.phiRx,'-k.'); legStr{ss} = 'phiRx'; ss = ss + 1; %TDC input in degrees
    plot(phiDCO,'-b.'); legStr{ss} = 'phiDCO'; ss = ss + 1; %TDC output in degrees
    plot(tdcOut/rp.Ktdc,'-g.'); legStr{ss} = 'tdcOut-degree'; ss = ss + 1; %TDC output in degrees
    legend(legStr); xlabel('Time'); ylabel('Phase (degree)'); title('Phase of DCO and TDC output'); grid;
    hold off

    ss = 1; legStr = {};
    handle3 = subaxis(3,1,3, 'Spacing', 0.06, 'Padding', 0.01, 'Margin', 0.03);
    hold on
    plot(tdcOut,'-k.'); legStr{ss} = 'tdcOut-codes'; ss = ss + 1; %TDC output in codes
    plot(LFint,'-c.'); legStr{ss} = 'LFint'; ss = ss + 1;
    plot(LFprop,'-b.'); legStr{ss} = 'LFprop'; ss = ss + 1;
    plot(LFOut,'-m.'); legStr{ss} = 'LFOut'; ss = ss + 1;
    legend(legStr); xlabel('Time'); ylabel('Codes'); title('Code levels of internal signal'); grid;
    hold off

    linkaxes([handle1 handle2 handle3],'x')
end