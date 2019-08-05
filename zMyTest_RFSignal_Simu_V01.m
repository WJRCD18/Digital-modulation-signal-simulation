%% 程序说明
% 【1】这个程序生成“有频偏 + 存在邻信道”的数字调制信号：仿真实际的盲信号；
% 【2】在【S02】【信号频谱搬移和混合】代码段，为仿真信号增加：频偏和邻信道信号
% 【3】信号参数的设置参见程序注释！
% 依赖的程序文件
%       1）qam16Modulator.m        
%% 程序初始参数设置
clear all; clc;
showFIG = 1;
sampleType = 'C';
% 【信号参数】===>>> 
        moduBaudRates = 1000;     % 想定的符号率
        moduSymbols = 1024 * 6;   % 生成信号的bit数，不是采样点数；
        moduSNR = 24;                    % 信噪比
        SamplesPerSymbol = 8;       % 每个符号的采样点数
        moduOrder = 8;                % 调制阶数；
        exeModType = 'fsk';          % 调制类型： fsk psk msk，qam
        deltaFreq = 0.1;              % 频率偏移，与SamplesPerSymbol相关 -->> 影响到信号星座图
        adjSigGain = 1;                    % 邻信道信号增益：0~100；
        
        moduSampleRates = moduBaudRates * SamplesPerSymbol/log2(moduOrder);   % 根据符号率，调制阶数，符号数，每个符号采样点数计算出的采样率
        str = sprintf('调制类型 = %d-%s; 信噪比 = %.0f;', moduOrder, exeModType, moduSNR);
        display(str);
        
%% 【S01】 调制信号生成
    for GEN = 1:1
        moduBits = randi([0 1], moduOrder * moduSymbols,1);
        hAWGN_ChanA = comm.AWGNChannel('NoiseMethod', 'Signal to noise ratio (SNR)','SNR', moduSNR);
        
        if strcmp(exeModType, 'fsk') == 1
            hMod_Msk    =  comm.FSKModulator(moduOrder, 2 * 50/moduOrder, 'SamplesPerSymbol', SamplesPerSymbol, 'BitInput', true, 'ContinuousPhase', false);
            ModSignal    = step(hMod_Msk, moduBits);
            NoisySignalA = step(hAWGN_ChanA, ModSignal);   
        end
        if strcmp(exeModType, 'msk') == 1 
            if moduOrder ~= 2
                beep;
                display('ERROR#29');
                return;
            end
            hMod_Msk    =  comm.MSKModulator('BitInput', true, 'SamplesPerSymbol', SamplesPerSymbol);
            ModSignal    = step(hMod_Msk, moduBits);
            NoisySignalA = step(hAWGN_ChanA, ModSignal);    
        end
        
        % 对于Q和8PSK，输出采样点数不同于2FSK和BPSK，因为要调制的二进制数据位数是相同的；
        if strcmp(exeModType, 'psk') == 1
            hMod_Psk    = comm.PSKModulator(moduOrder, 'BitInput', true);
            ModSignal    = step(hMod_Psk, moduBits); 
            if SamplesPerSymbol > 1
                RcosSignal   = rcosflt(ModSignal,1, SamplesPerSymbol);  
            else
                RcosSignal = ModSignal;
            end
            NoisySignalA = step(hAWGN_ChanA, RcosSignal);
            
            % 升余弦滤波器延迟消除
            if SamplesPerSymbol > 1
                NoisySignalA = NoisySignalA(1 * 3 * SamplesPerSymbol + 0: length(NoisySignalA));
            end
            
            %length(y) = (length(x) + 2 * delay)*Fs/Fd
        end

        % mQAM调制
        if strcmp(exeModType, 'qam') == 1
            if moduOrder ~= 4 && moduOrder ~= 8 && moduOrder ~= 16 && moduOrder ~= 32 && moduOrder ~= 64
                beep;
                display('ERROR#66');
                return;
            end

            hMod_Qam    = comm.RectangularQAMModulator(moduOrder,'BitInput',true);
            uFF = log2(moduOrder);
            deltaLEN = mod(length(moduBits),uFF);
            ModSignal    = step(hMod_Qam, moduBits(1:length(moduBits)-deltaLEN)); 
            if SamplesPerSymbol > 1
                RcosSignal   = rcosflt(ModSignal,1, SamplesPerSymbol);  
            else
                RcosSignal = ModSignal;
            end
            NoisySignalA = step(hAWGN_ChanA, RcosSignal);
            
            % 升余弦滤波器延迟消除
            if SamplesPerSymbol > 1
                NoisySignalA = NoisySignalA(1 * 3 * SamplesPerSymbol + 0: length(NoisySignalA));
            end
            
            %length(y) = (length(x) + 2 * delay)*Fs/Fd
        end
        
        
    end      
    
%% 【S02】===信号频谱搬移和混合=== 邻信道信号 + 主信号
        for SHIFT = 1:1
            moduBitsB = randi([0 1], moduOrder * moduSymbols,1);  %邻信道信号bit流
            hAWGN_ChanB = comm.AWGNChannel('NoiseMethod', 'Signal to noise ratio (SNR)','SNR', moduSNR/4); %邻信道白噪声
            % ==>> 邻信道信号调制：调制器与主信道相同，但是bit流不同
            for SIG = 1:1
                if strcmp(exeModType, 'fsk') == 1
                    ModSignalB    = step(hMod_Msk, moduBits);
                    NoisySignalB = step(hAWGN_ChanB, ModSignalB);   
                end
                if strcmp(exeModType, 'msk') == 1 
                    ModSignalB    = step(hMod_Msk, moduBitsB);
                    NoisySignalB = step(hAWGN_ChanB, ModSignalB);    
                end

                if strcmp(exeModType, 'psk') == 1
                    ModSignalB    = step(hMod_Psk, moduBitsB); 
                    if SamplesPerSymbol > 1
                        RcosSignalB   = rcosflt(ModSignalB,1, SamplesPerSymbol);  
                    else
                        RcosSignalB = ModSignalB;
                    end
                    NoisySignalB = step(hAWGN_ChanB, RcosSignalB);

                    % 升余弦滤波器延迟消除
                    if SamplesPerSymbol > 1
                        NoisySignalB = NoisySignalB(1 * 3 * SamplesPerSymbol + 0: length(NoisySignalB));
                    end
                end
                
                 if strcmp(exeModType, 'qam') == 1
                    uFF = log2(moduOrder);
                    deltaLEN = mod(length(moduBits),uFF);
                    ModSignalB    = step(hMod_Qam, moduBitsB(1:length(moduBits)-deltaLEN)); 
                    if SamplesPerSymbol > 1
                        RcosSignalB   = rcosflt(ModSignalB,1, SamplesPerSymbol);  
                    else
                        RcosSignalB = ModSignalB;
                    end
                    NoisySignalB = step(hAWGN_ChanB, RcosSignalB);
                    curLEN = length(NoisySignalB);

                    % 升余弦滤波器延迟消除
                    if 0 && SamplesPerSymbol > 1
                        NoisySignalB = NoisySignalB(1 * 3 * SamplesPerSymbol + 0: length(NoisySignalB));
                    else
                        NoisySignalB = NoisySignalB(1:curLEN - 3 * SamplesPerSymbol + 1);      %不消除滤波器延时
                    end
                end
           end
             
            t = 1:length(NoisySignalA);
            freqShift = exp(1j*2*pi*(deltaFreq * 0.5 / SamplesPerSymbol)*t);
            freqShift = freqShift(:);
            NoisySignalC = NoisySignalA.* freqShift;        %主信号频谱偏移
            
            freqShiftB = exp(1j*2*pi*(-0.5/2)*t);
            freqShiftB= freqShiftB(:);
            NoisySignalD = NoisySignalB.* freqShiftB;        %邻信道信号频谱偏移
            
            NoisySignalA = 1 * NoisySignalC + adjSigGain * NoisySignalD;     %【仿真实际信号】“有频偏 + 存在临信道”的数字信号：参见频谱
                        
            b = fir1(48, 1.5/SamplesPerSymbol);
            NoisySignalF = filter(b,1,NoisySignalA);                           %【滤波】模拟实际处理中的低通滤波器；
            
            for SHOW = 1:1
                figure(showFIG);    clf;    hold on; grid on;
                showFIG = showFIG + 1;
                fftRES = 1e-3 * moduSampleRates/length(NoisySignalA);   % 单位KHz
                logFFT = 20 * log10(fftshift(abs(fft(NoisySignalA)))/length(NoisySignalA));
                for n = 1:4
                    logFFT2 = logFFT;
                    for k= 1 + n : length(NoisySignalA) - 1
                        if logFFT(k) < logFFT(k - 1) && logFFT(k) < logFFT(k + 1)
                            logFFT2(k) = max(logFFT(k - 1), logFFT(k + 1));
                        end
                    end
                    logFFT = logFFT2;
                end
                
                plot(fftRES * ((1:length(NoisySignalF)) - length(NoisySignalA)/2), 20 * log10(fftshift(abs(fft(NoisySignalF)))/length(NoisySignalF)), '--ro','LineWidth',1,...
                        'MarkerEdgeColor','r',...
                        'MarkerFaceColor','r',...
                        'MarkerSize',1)
                plot( fftRES * ((1:length(NoisySignalA)) - length(NoisySignalA)/2), logFFT, '--bo','LineWidth',1,...
                        'MarkerEdgeColor','b',...
                        'MarkerFaceColor','b',...
                        'MarkerSize',1)
                title('主 信号频谱');    xlabel('Sample Points');    ylabel('Ampl(dB)');
                hold off;
                
                figure(showFIG);    clf;    hold on; grid on;
                showFIG = showFIG + 1;
                fftRES = 1e-3 * moduSampleRates/length(NoisySignalF);   % 单位KHz
                plot(fftRES * ((1:length(NoisySignalF)) - length(NoisySignalA)/2), 20 * log10(fftshift(abs(fft(NoisySignalF)))/length(NoisySignalF)), '--ro','LineWidth',1,...
                        'MarkerEdgeColor','r',...
                        'MarkerFaceColor','r',...
                        'MarkerSize',1)
                title('临 信号频谱');    xlabel('Sample Points');    ylabel('Ampl(dB)');
                hold off;
                
                figure(showFIG);    clf;    hold on; grid on;
                showFIG = showFIG + 1;
                fftRES = 1e-3 * moduSampleRates/length(freqShift);   % 单位KHz
                plot(fftRES * ((1:length(freqShift)) - length(freqShift)/2), 20 * log10(fftshift(abs(fft(freqShift)))/length(freqShift)), '--ro','LineWidth',1,...
                        'MarkerEdgeColor','r',...
                        'MarkerFaceColor','r',...
                        'MarkerSize',1)
                title('偏置 载波频谱');    xlabel('频率，KHz');    ylabel('Ampl(dB)');
                hold off;
                
                
            end
            
        end
           
%% 【S03】===解调误码率测试=== 仅仅用于验证信号频谱搬移和混合算法的正确性！ ==>> 后续过程中不再需要使用
    for BER = 1:0
        %Create an error rate calculator
        hError = comm.ErrorRate;
        
        if strcmp(exeModType, 'fsk') == 1
            hDemod = comm.FSKDemodulator(moduOrder, 'BitOutput', true, 'FrequencySeparation', 90, 'SamplesPerSymbol', SamplesPerSymbol);
            RcvSig = NoisySignalA;
            receivedData = step(hDemod, RcvSig);
            curSyms = min(length(receivedData),  length(moduBits));
            errorStats = step(hError, moduBits(1:curSyms), receivedData(1:curSyms));
        elseif strcmp(exeModType, 'psk') == 1
            hDemod = comm.PSKDemodulator(moduOrder, 'BitOutput', true);
            RcvSig = decimate(NoisySignalA, SamplesPerSymbol);
            receivedData = step(hDemod, RcvSig);
            curSyms = min(length(receivedData),  length(moduBits));
            errorStats = step(hError, moduBits(1:curSyms), receivedData(1:curSyms));
        elseif strcmp(exeModType, 'qam') == 1
            hDemod = comm.RectangularQAMDemodulator(moduOrder, 'BitOutput', true);
            RcvSig = decimate(NoisySignalA, SamplesPerSymbol);
            receivedData = step(hDemod, RcvSig);
            curSyms = min(length(receivedData),  length(moduBits));
            errorStats = step(hError, moduBits(1:curSyms), receivedData(1:curSyms));
        end
        fprintf('Error rate = %f\nNumber of errors = %d\n', errorStats(1), errorStats(2))        
      
        for SHOW = 1:1
            figure(showFIG);    clf;    hold on;
            showFIG = showFIG + 1;
            plot(real(RcvSig), '--ro','LineWidth',1,...
                    'MarkerEdgeColor','r',...
                    'MarkerFaceColor','r',...
                    'MarkerSize',2)
            plot(imag(RcvSig), '--yo','LineWidth',1,...
                    'MarkerEdgeColor','y',...
                    'MarkerFaceColor','y',...
                    'MarkerSize',2)
            plot(real(ModSignal), '--go','LineWidth',2,...
                    'MarkerEdgeColor','g',...
                    'MarkerFaceColor','g',...
                    'MarkerSize',2)
            plot(imag(ModSignal), '--bo','LineWidth',2,...
                    'MarkerEdgeColor','b',...
                    'MarkerFaceColor','b',...
                    'MarkerSize',2)
            title('Simulation IQ Signa');    xlabel('Time at Sample Points');    ylabel('Ampl');
            hold off;
        end
        
        
        
        for SHOW = 1:1
            figure(showFIG);    clf;    hold on;
            showFIG = showFIG + 1;
            plot(real(NoisySignalA), '--ro','LineWidth',1,...
                    'MarkerEdgeColor','r',...
                    'MarkerFaceColor','r',...
                    'MarkerSize',2)
            plot(imag(NoisySignalA), '--bo','LineWidth',1,...
                    'MarkerEdgeColor','b',...
                    'MarkerFaceColor','b',...
                    'MarkerSize',2)
            title('Simulation IQ Signa');    xlabel('Time at Sample Points');    ylabel('Ampl');
            hold off;
        end
    end

%return;

%% 【S04】星座图恢复：必须考虑《采样率/符号率》是非整数的情况  ==>> 出现频偏时，星座图是无法恢复的！
    for CST = 1:1
        p = 62;
        q = 65;     % 这两个参数调整“采样率/符号率”倍数：非整数
        NoisySignalH = resample(NoisySignalF, p, q);     % resamples the sequence in vector x at p/q times the original sampling rate, 
        
        maxQQ = max(abs(imag(NoisySignalH)));
        maxII = max(abs(real(NoisySignalH)));
        vFactor = max(maxQQ, maxII);
        
        vImage=zeros(101,101); %initialize
        for k= 1:length(NoisySignalH)
            x = int32(50 * real(NoisySignalH(k))/vFactor) + 50 + 1;
            y = int32(50 * imag(NoisySignalH(k))/vFactor) + 50 + 1;
            vImage(x, y) = vImage(x, y) + 1;
        end
        maxVV = max(vImage);
        vImage = vImage / max(maxVV);
        
        figure(showFIG);    clf;    hold on;
        showFIG = showFIG + 1;
        imshow(vImage)
        hold off;
    end

    
