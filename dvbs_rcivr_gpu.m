%DVB-S Receiver script
clear;clc;
close all;
%% System inputs
symbolRate  = 27500e3;       %% enter symbol rate for this variable
sampleRate  = 36e6;         %% enter sample rate for this variable
Bandwidth   = 32e6;          %% Bandwidth
CR          = '5/6';         %% Code Rate (convolutional decoder)
filePath    = 'G:\khajezadeh\dvb\short_36MHz_27500KS.bin';
fileFormat  = 'short';
OutputFile  = 'MPEG_OUTPUT.ts';
rc_beta     = 0.35;
SPS         = 2;             %% system samples per symbol
plotRequest = true;
%% check inputs
fid  = fopen(filePath,'r');
fidO = fopen(OutputFile, 'w');
if (fid == -1)
    error('File Not Found!');
end

%% initialize constants
[cr(1), cr(2)] = rat(str2num(CR));
rsInLen = 1632;
%% Construt Receiver units
fprintf('Constructing System units ... \n');

% digital downconverter
dwnConv = dsp.DigitalDownConverter('SampleRate', 250e6, 'DecimationFactor', 2,...
    'Bandwidth', Bandwidth, 'CenterFrequency', 70e6);

% AGC
agc = comm.AGC('AdaptationStepSize', 0.01, 'DesiredOutputPower', 0.5,...
    'AveragingLength', 100);

% resampler
src = dsp.SampleRateConverter('Bandwidth',Bandwidth,'InputSampleRate',...
    sampleRate,'OutputSampleRate',SPS*symbolRate);

% raised cosine filter
H = comm.RaisedCosineReceiveFilter('RolloffFactor',rc_beta,...
    'FilterSpanInSymbols',10,'InputSamplesPerSymbol',SPS,...
    'Gain',1,'DecimationFactor',1);
% Coarse Frequency offset
CFC = comm.CoarseFrequencyCompensator('Modulation','QPSK','Algorithm','FFT-based',...
    'FrequencyResolution',SPS*symbolRate/100,'SampleRate',SPS*symbolRate);
% Carrier frequency offset
carrSynch = comm.CarrierSynchronizer('Modulation','QPSK','SamplesPerSymbol',1,...
    'DampingFactor',0.707,'NormalizedLoopBandwidth',0.001);
% Symbol Synchronizer
SS = comm.SymbolSynchronizer('Modulation','PAM/PSK/QAM','SamplesPerSymbol',SPS,...
    'DampingFactor', 1, 'NormalizedLoopBandwidth', 0.01, 'DetectorGain', 2.7);
% demodulator
DeMod = comm.PSKDemodulator(4,pi/4,'BitOutput',true,'SymbolMapping','Custom',...
    'CustomSymbolMapping',[0,2,3,1],'OutputDataType','double');
% convolutional decoder
switch CR
    case '1/2'
        PT = [1,1]';
        tryCases = 1;
        trBackDep = 30;
    case '2/3'
        PT = [1,1,0,1]';
        tryCases = 3;
        trBackDep = 42;
    case '3/4'
        PT = [1,1,0,1,1,0]';
        tryCases = 2;
        trBackDep = 60;
    case '5/6'
        PT = [1,1,0,1,1,0,0,1,1,0]';
        tryCases = 3;
        trBackDep = 90;
    case '7/8'
        PT = [1,1,0,1,0,1,0,1,1,0,0,1,1,0]';
        tryCases = 4;
        trBackDep = 112;
end
vitdec = comm.gpu.ViterbiDecoder(poly2trellis(7, [171 133]), ...
    'InputFormat', 'Hard','PuncturePatternSource','Property',...
    'PuncturePattern',PT,'TracebackDepth', trBackDep, 'OutputDataType','Full precision');
% interleaver
deinterleaver = comm.ConvolutionalDeinterleaver('NumRegisters',12,'RegisterLengthStep',17);
% RS Decoder
RSdec = comm.RSDecoder ('CodewordLength',255,'MessageLength',239,...
    'ShortMessageLength',188, 'GeneratorPolynomialSource','Property',...
    'GeneratorPolynomial',rsgenpoly(255,239,285,0));
% plot
M = 4;
x = 0:(M-1);
refSig = qammod(x,M, 'UnitAveragePower', 1);
constdiag = comm.ConstellationDiagram('ReferenceConstellation', refSig);
%% processing
fprintf('Processing Start. \n');
SYNCH             = false;
processData       = 0.01; % seconds
vitBuff           = zeros(10 * processData * sampleRate, 1, 'logical');
vitBufLength      = 0;
numMPEGforSYNCH   = 20;
PRBSsynch         = 0;
PRBSoffset        = 0;
numMPEGFrames     = 0;
totalProcesed     = 0;

% fseek(fid,30e6,'bof');
rawData           = fread(fid,1e3, fileFormat);
%env              = dwnConv(rawData);release(dwnConv);
env               = complex(rawData(1:2:end), rawData(2:2:end));
env_normalized    = agc(env);
tic

while(true)
     rawData = fread(fid, 2 * sampleRate * processData, fileFormat);  % each time, process 10ms of Data
     if (feof(fid))
         break;
     end
     env            = complex(rawData(1:2:end), rawData(2:2:end));
     %env           = dwnConv(rawData);
     env_normalized = agc(env);
     env_resampled  = src(env_normalized);
     env_rcfiltered = H(env_resampled);
     env_Cfreq      = CFC(env_rcfiltered);
     env_sym        = SS(env_Cfreq);
     IQ             = carrSynch(env_sym);
    
    if (plotRequest)
        constdiag(IQ(1:1e4));      
    end
    if (~SYNCH)
        numMPEGforSYNCH = 20;
        [FLAG_LOCKED, FLAG_INV, INIT_SYNCH, PUNCH_ALIGN, FLAG_90] = synchFinder(IQ, CR, numMPEGforSYNCH);
        if(FLAG_LOCKED)
            fprintf('System Locked to MPEG SYNCH. \n');
            SYNCH                   = true; release(vitdec); release(DeMod);release(deinterleaver);release(RSdec)
            DeMod.PhaseOffset       = pi/4 + pi/2 * FLAG_90 + pi * FLAG_INV;
            vitdec.TracebackDepth   = trBackDep - mod(INIT_SYNCH, cr(1)) + 1;
            vitBufLength            = vitBufLenCalculator(cr(1), cr(2), rsInLen, length(IQ), numMPEGforSYNCH);
            rawBits                 = DeMod(IQ);
            [iniBitsToViterbi, idx] = vitInitBitsCalculator(rawBits, cr(1), cr(2), INIT_SYNCH, PUNCH_ALIGN, vitBufLength);
            vitBuff                 = rawBits(idx:end);
            numMPEGFrames           = vitBufLength * cr(1) / cr(2) / rsInLen;
            RBinDec                 = randomizerBitsCalc(numMPEGFrames);
            iniBitsToViterbi        = gpuArray(iniBitsToViterbi);
            vitdec(iniBitsToViterbi);
        else
            continue
        end
    elseif(SYNCH)
        rawBits = DeMod(IQ);
        vitBuff = [vitBuff; rawBits];
        while(length(vitBuff) >= vitBufLength)
            gpuVitIn    = gpuArray(vitBuff(1:vitBufLength));
            bitsVitDec  = vitdec(gpuVitIn);
            bitsVitDec  = gather(bitsVitDec);
            vitBuff     = vitBuff((vitBufLength + 1):end);
            bitsVitDec  = reshape(bitsVitDec, 8, []);
            DataDecimal = cast(bi2de(bitsVitDec', 'left-msb'), 'uint8');
            Data_Deint  = deinterleaver(DataDecimal);
            [rData, Er] = RSdec(Data_Deint);
            if(sum(Er) < -(length(Er) / 2))
                SYNCH = false;
                PRBSsynch = false;
                fprintf('SYNCH Lost! \n');
                continue
            end
            if (~PRBSsynch)
                PRBSsynch  = true;
                PRBSoffset = find(rData(1:188:end) == 184, 1, 'first');
                PRBSoffset = 8 - mod(PRBSoffset - 1, 8) + 1;
            else
                PRBSoffset = mod(PRBSoffset - 1 + numMPEGFrames, 8) + 1;
            end
            MPEG_FRAMES = bitxor(rData, RBinDec((1:(numMPEGFrames * 188)) + (PRBSoffset - 1) * 188), 'uint8');
            fwrite(fidO, MPEG_FRAMES, 'uint8');
            elapsedTime = toc;
            totalProcesed = totalProcesed + numMPEGFrames;
            clc;
            fprintf('Total Processed Frames = %d \n', totalProcesed);
            fprintf('Rate (MPEG Frames per second) = %0.2f \n', totalProcesed / elapsedTime);
        end
    end
end
fclose(fid);
fclose(fidO);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function L = vitBufLenCalculator(crNum, crDenum,rsInLength, TRlen, numMPEGforSYNCH)
availableBits = 2 * TRlen - (numMPEGforSYNCH + 1) * rsInLength;
x = crNum / gcd(crNum, rsInLength);
Linit = rsInLength * x * crDenum / crNum;
if (availableBits < Linit)
    error("length of available data is not Enough!")
end

mul = 1;
while(Linit * mul < availableBits)
    L = mul * Linit;
    mul = mul + 1;
end
end
function [initBits, firstIdxToVitBuf] = vitInitBitsCalculator(data_Demod,...
    crNum, crDenum, INIT_SYNCH, PUNCH_ALIGN, vitBufLength)
numBits = (INIT_SYNCH - mod(INIT_SYNCH, crNum)) * crDenum / crNum - 2 * PUNCH_ALIGN;
initBits = [zeros(2 * PUNCH_ALIGN, 1);data_Demod(1:numBits)];
initBits = [zeros(vitBufLength - length(initBits), 1);initBits];
initBits = cast(initBits, 'double');
firstIdxToVitBuf = numBits + 1;

end
function despersalBits = randomizerBitsCalc(numMPEGinFrame)
ini     = [1 0 0 1 0 1 0 1 0 0 0 0 0 0 0];
len     = numMPEGinFrame + 8;
rndPRBS = zeros(1,1503*8);
for i = 1 : 1503*8
    rndPRBS(i) = xor(ini(15),ini(14));
    ini    = circshift(ini,1);
    ini(1) = rndPRBS(i);
end
rndPRBS    = [zeros(1,8),rndPRBS];
En         = [zeros(1,8),ones(1,187*8)];
Enable     = repmat(En,1,8);
PRBS1504   = and(Enable,rndPRBS);
decPRBS    = reshape(PRBS1504, 8, []);
decPRBS    = bi2de(decPRBS', 'left-msb');
decPRBS(1) = 255;
despersalBits = cast(decPRBS(mod(0:(len * 188 - 1), numel(decPRBS)) + 1), 'uint8');
end





