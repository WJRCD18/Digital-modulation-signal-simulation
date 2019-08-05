%% ����˵��
% ��1������������ɡ���Ƶƫ + �������ŵ��������ֵ����źţ�����ʵ�ʵ�ä�źţ�
% ��2���ڡ�S02�����ź�Ƶ�װ��ƺͻ�ϡ�����Σ�Ϊ�����ź����ӣ�Ƶƫ�����ŵ��ź�
% ��3���źŲ��������òμ�����ע�ͣ�
% �����ĳ����ļ�
%       1��qam16Modulator.m        
%% �����ʼ��������
clear all; clc;
showFIG = 1;
sampleType = 'C';
% ���źŲ�����===>>> 
        moduBaudRates = 1000;     % �붨�ķ�����
        moduSymbols = 1024 * 6;   % �����źŵ�bit�������ǲ���������
        moduSNR = 24;                    % �����
        SamplesPerSymbol = 8;       % ÿ�����ŵĲ�������
        moduOrder = 8;                % ���ƽ�����
        exeModType = 'fsk';          % �������ͣ� fsk psk msk��qam
        deltaFreq = 0.1;              % Ƶ��ƫ�ƣ���SamplesPerSymbol��� -->> Ӱ�쵽�ź�����ͼ
        adjSigGain = 1;                    % ���ŵ��ź����棺0~100��
        
        moduSampleRates = moduBaudRates * SamplesPerSymbol/log2(moduOrder);   % ���ݷ����ʣ����ƽ�������������ÿ�����Ų�������������Ĳ�����
        str = sprintf('�������� = %d-%s; ����� = %.0f;', moduOrder, exeModType, moduSNR);
        display(str);
        
%% ��S01�� �����ź�����
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
        
        % ����Q��8PSK���������������ͬ��2FSK��BPSK����ΪҪ���ƵĶ���������λ������ͬ�ģ�
        if strcmp(exeModType, 'psk') == 1
            hMod_Psk    = comm.PSKModulator(moduOrder, 'BitInput', true);
            ModSignal    = step(hMod_Psk, moduBits); 
            if SamplesPerSymbol > 1
                RcosSignal   = rcosflt(ModSignal,1, SamplesPerSymbol);  
            else
                RcosSignal = ModSignal;
            end
            NoisySignalA = step(hAWGN_ChanA, RcosSignal);
            
            % �������˲����ӳ�����
            if SamplesPerSymbol > 1
                NoisySignalA = NoisySignalA(1 * 3 * SamplesPerSymbol + 0: length(NoisySignalA));
            end
            
            %length(y) = (length(x) + 2 * delay)*Fs/Fd
        end

        % mQAM����
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
            
            % �������˲����ӳ�����
            if SamplesPerSymbol > 1
                NoisySignalA = NoisySignalA(1 * 3 * SamplesPerSymbol + 0: length(NoisySignalA));
            end
            
            %length(y) = (length(x) + 2 * delay)*Fs/Fd
        end
        
        
    end      
    
%% ��S02��===�ź�Ƶ�װ��ƺͻ��=== ���ŵ��ź� + ���ź�
        for SHIFT = 1:1
            moduBitsB = randi([0 1], moduOrder * moduSymbols,1);  %���ŵ��ź�bit��
            hAWGN_ChanB = comm.AWGNChannel('NoiseMethod', 'Signal to noise ratio (SNR)','SNR', moduSNR/4); %���ŵ�������
            % ==>> ���ŵ��źŵ��ƣ������������ŵ���ͬ������bit����ͬ
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

                    % �������˲����ӳ�����
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

                    % �������˲����ӳ�����
                    if 0 && SamplesPerSymbol > 1
                        NoisySignalB = NoisySignalB(1 * 3 * SamplesPerSymbol + 0: length(NoisySignalB));
                    else
                        NoisySignalB = NoisySignalB(1:curLEN - 3 * SamplesPerSymbol + 1);      %�������˲�����ʱ
                    end
                end
           end
             
            t = 1:length(NoisySignalA);
            freqShift = exp(1j*2*pi*(deltaFreq * 0.5 / SamplesPerSymbol)*t);
            freqShift = freqShift(:);
            NoisySignalC = NoisySignalA.* freqShift;        %���ź�Ƶ��ƫ��
            
            freqShiftB = exp(1j*2*pi*(-0.5/2)*t);
            freqShiftB= freqShiftB(:);
            NoisySignalD = NoisySignalB.* freqShiftB;        %���ŵ��ź�Ƶ��ƫ��
            
            NoisySignalA = 1 * NoisySignalC + adjSigGain * NoisySignalD;     %������ʵ���źš�����Ƶƫ + �������ŵ����������źţ��μ�Ƶ��
                        
            b = fir1(48, 1.5/SamplesPerSymbol);
            NoisySignalF = filter(b,1,NoisySignalA);                           %���˲���ģ��ʵ�ʴ����еĵ�ͨ�˲�����
            
            for SHOW = 1:1
                figure(showFIG);    clf;    hold on; grid on;
                showFIG = showFIG + 1;
                fftRES = 1e-3 * moduSampleRates/length(NoisySignalA);   % ��λKHz
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
                title('�� �ź�Ƶ��');    xlabel('Sample Points');    ylabel('Ampl(dB)');
                hold off;
                
                figure(showFIG);    clf;    hold on; grid on;
                showFIG = showFIG + 1;
                fftRES = 1e-3 * moduSampleRates/length(NoisySignalF);   % ��λKHz
                plot(fftRES * ((1:length(NoisySignalF)) - length(NoisySignalA)/2), 20 * log10(fftshift(abs(fft(NoisySignalF)))/length(NoisySignalF)), '--ro','LineWidth',1,...
                        'MarkerEdgeColor','r',...
                        'MarkerFaceColor','r',...
                        'MarkerSize',1)
                title('�� �ź�Ƶ��');    xlabel('Sample Points');    ylabel('Ampl(dB)');
                hold off;
                
                figure(showFIG);    clf;    hold on; grid on;
                showFIG = showFIG + 1;
                fftRES = 1e-3 * moduSampleRates/length(freqShift);   % ��λKHz
                plot(fftRES * ((1:length(freqShift)) - length(freqShift)/2), 20 * log10(fftshift(abs(fft(freqShift)))/length(freqShift)), '--ro','LineWidth',1,...
                        'MarkerEdgeColor','r',...
                        'MarkerFaceColor','r',...
                        'MarkerSize',1)
                title('ƫ�� �ز�Ƶ��');    xlabel('Ƶ�ʣ�KHz');    ylabel('Ampl(dB)');
                hold off;
                
                
            end
            
        end
           
%% ��S03��===��������ʲ���=== ����������֤�ź�Ƶ�װ��ƺͻ���㷨����ȷ�ԣ� ==>> ���������в�����Ҫʹ��
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

%% ��S04������ͼ�ָ������뿼�ǡ�������/�����ʡ��Ƿ����������  ==>> ����Ƶƫʱ������ͼ���޷��ָ��ģ�
    for CST = 1:1
        p = 62;
        q = 65;     % ����������������������/�����ʡ�������������
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

    
