%function code=PPM_Arith_enco(seq,p,symbol)
clear;
    %sym=symbol;
    sym=[1 2 3 4];
    %prob=p;
    prob=[0.1 0.4 0.2 0.3];
    seq=[3 1 4 1 3 4 2];

    lows = 0; highs = prob(1);
    for i = 2:length(sym)
        lows(i) = lows(i-1)+ prob(i-1);
        highs(i) = highs(i-1)+ prob(i);
    end
    InitialTable = [prob' lows' highs'];



    LowEnc = 0; HighEnc = 1; Progress = [LowEnc HighEnc]; prob_message = 1;
    for i = 1:length(seq)
        current_symbol = seq(i);
        j = find( sym == current_symbol );
        los = lows(j); his = highs(j);
        current_interval = HighEnc-LowEnc;
        HighEnc = LowEnc + current_interval*his;
        LowEnc = LowEnc + current_interval*los;
        Progress = [Progress; LowEnc HighEnc];
        prob_message = prob_message * prob(j);
    end


    intLowEnc = floor(LowEnc*2^40);
    binLowEnc = bitget( intLowEnc, 40:-1:1);
    intHighEnc = floor(HighEnc*2^40);
    binHighEnc = bitget( intHighEnc, 40:-1:1);

    ind = find(binHighEnc~= binLowEnc,1);
    ind1 = find( binLowEnc(ind:end) == 0,1);
    ind2 = ind+ind1-1;
    binLowEnc(ind2) = 1;

    code=binLowEnc(1:ind2);
    message = binLowEnc(1:ind2)*(2.^(-(1:ind2)'));
    ideal_codelength = - log2(prob_message);

    %MessageDecoded=PPM_Arith_deco(code,length(seq),sym,prob);

    %average=length(code)/length(seq);
%end
