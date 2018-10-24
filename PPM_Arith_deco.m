function MessageDecoded=PPM_Arith_deco(code,len,symbol,p)

    MessageDecoded = [];
    code_str=['0.' num2str(code)];
    code_str=code_str(~isspace(code_str));
    message = f_b2d(code_str);
    crt_message = message;
    
    lows = 0; highs = p(1);
    for i = 2:length(symbol)
        lows(i) = lows(i-1)+ p(i-1);
        highs(i) = highs(i-1)+ p(i);
    end
    
    
    for i = 1:len
        for j = 1:length(symbol)
            los = lows(j); his = highs(j);
            if( (los <= crt_message) && (his > crt_message) )
                break
            end
        end
        current_interval = his-los;
        crt_message = (crt_message - los)/current_interval;
        MessageDecoded = [MessageDecoded symbol(j)];
    end
end