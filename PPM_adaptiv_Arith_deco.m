clear;
load('code.mat');
load('table.mat');
order=2;


MessageDecoded = [];
code_str=['0.' num2str(code)];
code_str=code_str(~isspace(code_str));
message = f_b2d(code_str);
crt_message = message;


n=order;
for o=1:n+1
    k(o)=1;
end
%Arith_enco set InitialTable
for i=1:71
    init_sym(i,1)=i+79;
    init_sym(i,2)=table{1,2}(i,1)/sum(table{1,2});
    init_sym(i,3)=0;
end



len_message=2;
%len=length(block);

for i = 1:len_message

    lows = 0; highs = init_sym(1,2);
    for h = 2:length(init_sym)
        lows(h) = lows(h-1)+ init_sym(h-1,2);
        highs(h) = highs(h-1)+ init_sym(h,2);
    end   


    for j = 1:length(init_sym)
        los = lows(j); his = highs(j);
        if( (los <= crt_message) && (his > crt_message) )
            break
        end
    end

    MessageDecoded = [MessageDecoded init_sym(j,1)]; 
    current_interval = his-los;
    crt_message = (crt_message - los)/current_interval;     
    
    %calculate new probabilities
    %search block(i) in table, if found, encode
    %search the highest order first, if not found, order-1
    table_enco=table;
    P_higher=1;
    countnum=1;
    counted=[];    
    found=0;
    
    if i<n+1        %not enough data to complete order n
        j=1;
    else            %enough data to complete order n
        j=i-n;
    end
    
    
    while j<=i
        seq=MessageDecoded(j:i); 
        seqstr=num2str(seq);
        contextorder=i-j+1;     %actually should -1 
        len=size(table_enco{contextorder});
        
        if found==0
            for L=1:len(1,1)
                if strcmp(seqstr,table_enco{contextorder,1}(L,1))
                    found=1;
                    contextorder_found=contextorder;
                    break
                else
                    found=0;
                end
            end
        end
        
        if found==1
            %P=P*(1/(sum(table_enco{contextorder,2}(:,1))+1));
            table{contextorder_found,2}(L,1)=table{contextorder_found,2}(L,1)+1;        
        else
            table{contextorder,1}(k(contextorder),1)={seqstr};
            table{contextorder,2}(k(contextorder),1)=1;
            k(contextorder)=k(contextorder)+1;              
        end  
        j=j+1;
    end
 
    table_enco=table;
    P_higher=1;
    countnum=1;
    counted=[];    
    found=0;
    
    if i<n+1        %not enough data to complete order n
        j=1;
    else            %enough data to complete order n
        j=i-n;
    end    
    
    while j<=i
        seq=MessageDecoded(j:i); 
        seqstr=num2str(seq);
        contextorder=i-j+1;     %actually should -1 
        len=size(table_enco{contextorder});
        
        if found==0
            for L=1:len(1,1)
                if strcmp(seqstr,table_enco{contextorder,1}(L,1))
                    found=1;
                    contextorder_found=contextorder;
                    break
                else
                    found=0;
                end
            end
        end
        
        %delet all irrelated context and upgrade table
      
        if length(seq)>1
            comp_block=num2str(seq(1:length(seq)-1));
        else
            comp_block=num2str(seq);
        end
        lastele=num2str(MessageDecoded(i));
        %comp_block=comp_block(1:(length(comp_block)-length(lastele)-2));
        for m=1:len(1,1)
            comp=table_enco{contextorder,1}{m,1};
            if contextorder>1
                comp=comp(1:(length(comp)-length(lastele)-2));
            end
            tablelast_all=str2num(table_enco{contextorder,1}{m,1});
            tablelast=tablelast_all(end);     
            
            if contextorder~=1 
                if ~strcmp(comp_block,comp) || ~isempty(find(counted==tablelast, 1))
                    table_enco{contextorder,1}{m,1}=[];
                    table_enco{contextorder,2}(m,1)=0;
                end
            else
                if ~isempty(find(counted==tablelast, 1))
                    table_enco{contextorder,1}{m,1}=[];
                    table_enco{contextorder,2}(m,1)=0;
                end
            end        
        end 
        
        %calculate all probabilities of the current context
        for m=1:len(1,1)
            if ~isempty(table_enco{contextorder,1}{m,1})
                comp=table_enco{contextorder,1}{m,1};
                comp=comp(1:(length(comp)-length(lastele)-2));

                tablelast_all=str2num(table_enco{contextorder,1}{m,1});
                tablelast=tablelast_all(end);              
                
                if contextorder==1
                    P=P_higher*(table_enco{contextorder,2}(m,1)/sum(table_enco{contextorder,2}(:,1)));
                else
                    P=P_higher*(table_enco{contextorder,2}(m,1)/(sum(table_enco{contextorder,2}(:,1))+1));
                end
                init_pos=find(init_sym(:,1)==tablelast);
                if init_sym(init_pos,3)==0
                    init_sym(init_pos,2)=P;
                    init_sym(init_pos,3)=1;
                end
            end
        end
        P_higher=P_higher*(1/(sum(table_enco{contextorder,2}(:,1))+1));
        
        if found==1
            %P=P*(1/(sum(table_enco{contextorder,2}(:,1))+1));
            %table{contextorder_found,2}(L,1)=table{contextorder_found,2}(L,1)+1;
            
            pos=find(table_enco{contextorder_found,2}(:,1)~=0);
            for t=1:numel(pos)
                count=str2num(table_enco{contextorder_found,1}{pos(t),1});
                counted(countnum)=count(end);
                countnum=countnum+1;
            end 
            
        else
%             table{contextorder,1}(k(contextorder),1)={seqstr};
%             table{contextorder,2}(k(contextorder),1)=1;
%             k(contextorder)=k(contextorder)+1;  
            
            pos=find(table_enco{contextorder,2}(:,1)~=0);
            for t=1:numel(pos)
                count=str2num(table_enco{contextorder,1}{pos(t),1});
                counted(countnum)=count(end);
                countnum=countnum+1;
            end                 
        end
        j=j+1;
    end
    
    for j=1:71
        init_sym(j,3)=0;
    end 
end

save('table_deco','table');

