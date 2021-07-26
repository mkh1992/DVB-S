%synchronizer
%% inputs
function [FLAG_LOCKED, FLAG_INV, INIT_SYNCH, PUNCH_ALIGN, FLAG_90] = synchFinder(inData, CR, numMPEGFrames)
FLAG_LOCKED  = 0;
FLAG_INV     = 0;
INIT_SYNCH   = 0;
PUNCH_ALIGN  = 0;
FLAG_90      = 0;
%% initialize parameters
[cr(1), cr(2)] = rat(str2num(CR));
reqLength = floor(204 * 8 * numMPEGFrames * cr(2) / cr(1) / 2);
reqLength = cr(2) - rem(reqLength, cr(2)) + reqLength;
if((length(inData)) < (reqLength + 200))
    fprintf("Not Enough Data For synchFinder \n");
    return
end
%% convolutional decoder
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
    otherwise
        return
end
vitdec = comm.ViterbiDecoder(poly2trellis(7, [171 133]), ...
    'InputFormat', 'Unquantized','PuncturePatternSource','Property',...
    'PuncturePattern',PT,'TracebackDepth', trBackDep,'OutputDataType','logical');

% DeMod1 = comm.PSKDemodulator(4,pi/4,'BitOutput',true,'SymbolMapping','Custom',...
%     'CustomSymbolMapping',[0,2,3,1],'OutputDataType','logical');
% DeMod2 = comm.PSKDemodulator(4,pi/4 + pi/2,'BitOutput',true,'SymbolMapping','Custom',...
%     'CustomSymbolMapping',[0,2,3,1],'OutputDataType','logical');

DeMod1 = comm.PSKDemodulator(4, 'PhaseOffset', pi/4, 'SymbolMapping', 'Custom', ...
    'CustomSymbolMapping', [0 2 3 1], 'BitOutput', 1 ,'DecisionMethod', ...
    'Approximate log-likelihood ratio');
DeMod2 = comm.PSKDemodulator(4, 'PhaseOffset', pi/4 + pi/2, 'SymbolMapping', 'Custom', ...
    'CustomSymbolMapping', [0 2 3 1], 'BitOutput', 1 ,'DecisionMethod', ...
    'Approximate log-likelihood ratio');

Hard1 = DeMod1(inData(1:reqLength));
Hard2 = DeMod2(inData(1:reqLength));
for PA = 0:(tryCases - 1)
    reset(vitdec);
    candidate1 = vitdec(circshift(Hard1, 2 * PA));
    [FLAG_LOCKED, FLAG_INV, INIT_SYNCH] = sf(candidate1);
    if (FLAG_LOCKED)
        PUNCH_ALIGN = PA;
        FLAG_90 = 0;
        break
    end
    reset(vitdec);
    candidate2 = vitdec(circshift(Hard2, 2 * PA));
    [FLAG_LOCKED, FLAG_INV, INIT_SYNCH] = sf(candidate2);
    if (FLAG_LOCKED)
        PUNCH_ALIGN = PA;
        FLAG_90 = 1;
        break
    end
    
end

    function [FLAG_LOCKED, FLAG_INV, INIT_SYNCH]= sf(input)
        FLAG_LOCKED = 0;
        FLAG_INV    = 0;
        INIT_SYNCH  = 0;
        syncBits = cast([0;1;0;0;0;1;1;1], 'logical');
        rsInLen = 1632;
        idx = zeros(500, 2);
        jp   = 0;
        jn   = 0;
        for i = 1:(length(input) - 8)
            if isequal(syncBits, input((0:7) + i))
                jp = jp + 1;
                idx(jp, 1) = i;
            elseif(isequal(not(syncBits), input((0:7) + i)))
                jn = jn + 1;
                idx(jn, 2) = i;
            end
            
        end
        points = zeros(length(idx),2);
        for i = 1:length(idx)
            s = idx(i, 1);
            while(s + rsInLen  <= idx(jp, 1))
                points(i, 1) = points(i, 1) + any(idx(:, 1) == (s + rsInLen));
                s = s + rsInLen;
            end
            s = idx(i, 2);
            while(s + rsInLen  <= idx(jn, 2))
                points(i, 2) = points(i, 2) + any(idx(:, 2) == (s + rsInLen));
                s = s + rsInLen;
            end
        end
        [mx, mxIdx] = max(points);
        if mx(1) > (length(input) / rsInLen) * 2/3
            FLAG_LOCKED = 1;
            FLAG_INV = 0;
            INIT_SYNCH = idx(mxIdx(1), 1);
        elseif mx(2) > (length(input) / rsInLen) * 2/3
            FLAG_LOCKED = 1;
            FLAG_INV = 1;
            INIT_SYNCH = idx(mxIdx(2), 2);
        end
    end
end

