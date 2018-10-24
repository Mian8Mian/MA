%Coding

clear;

%set file path
quantityoffiles=1;
filelist=dir('D:\Masterarbeit\Quelle\ACE_AUTO');
dest='D:\Masterarbeit\matlab\debug';  

%define block size, set block number
rownr=[1 1];        %[startrownr endrownr]
columnsize=10;    %comment if data is the whole matrix
blocknumber=451;      %set 1 if data is the whole matrix
order=2;
format long;

%select random files
qty=quantityoffiles;
randf=randperm(numel(filelist),qty);


% for ii=1:qty
%     while (randf(ii)==1)||(randf(ii)==2)
%         randf(ii)=randperm(numel(filelist),1);
%     end
% end
% 
% %copy the selected files to destination folder
% for kk=1:qty
%     source=fullfile(filelist(randf(kk)).folder,filelist(randf(kk)).name);
%     copyfile(source,dest);
% end

%get filelist of the current destination folder
cfiles=dir(dest);

%initialization for clone struct
clone=struct();

%initialization for blockp struct
blockp=struct();

%cfiles(1) and cfiles(2) are '.' and '..', file name is stored from cfiles(3)
jj=3;
while jj<=qty+2

    %load data file
    load(fullfile(cfiles(jj).folder,cfiles(jj).name));

    %initialization, convert NaN to 0
    clone(jj-2).file=cfiles(jj).name;
    clone(jj-2).mat=mat;

    %next data file
    jj=jj+1;
end

%cfiles(1) and cfiles(2) are '.' and '..', file name is stored from cfiles(3)
jj=3;
while jj<=qty+2
    
    %initialization, convert NaN to 0
    clone(jj-2).mat(isnan(clone(jj-2).mat))=0;

    %define block size and split the matrix {clone} into smaller blocks 
    %nob=number of blocks
    row=rownr;
     
    column=columnsize;  %use this if data is only a part of the matrix
    %buffer=size(mat);   %use this if data is the whole matrix
    %column=buffer(2);   %use this if data is the whole matrix
    nob=fix(length(clone(jj-2).mat)/column);

    if nob~=length(clone(jj-2).mat)/column
        nob=nob+1;
    end
    
    %preallocate array(part) size
    part=zeros(nob,2);    
    
    %determine if there is more than one block
    if nob<=1
        part(1,1)=1;
        part(1,2)=length(clone(jj-2).mat);
    end

    if nob>1
        for i=1:nob
            if i~=nob 
                part(i,1)=(i-1)*column+1;
                part(i,2)=column*i;
            else
                part(i,1)=(i-1)*column+1;
                part(i,2)=length(clone(jj-2).mat);
            end
        end
    end

    %set block number
    b=blocknumber;
    
    %determine if the block number too large or negative
    if (b>nob)
        disp("Maximum Block Number"),disp(nob); 
    return
    end
    if (b<0)
        disp("Block Number Can Not Be Negative");
    return
    end

    %get the required block
    blockog=clone(jj-2).mat(row(1,1):row(1,2),part(b,1):part(b,2));
    rowsize=row(1,2)-row(1,1)+1;    

    %delete unneeded zeros
    blockp(jj-2).file=clone(jj-2).file;
    blockp(jj-2).mat=blockog;
    remainzero=zeros(1,11);

    c=1;
    del=0;
    while c<=numel(blockp(jj-2).mat)
        if ((c+numel(remainzero)-1)<numel(blockp(jj-2).mat))
            if (blockp(jj-2).mat(c:c+numel(remainzero)-1)==remainzero)
                c=c+numel(remainzero);
                del=1;
            end
        end
    while c<=numel(blockp(jj-2).mat)
        if (del==1)&&(blockp(jj-2).mat(c)==0) 
            blockp(jj-2).mat(c)=[];    
        else
            break
        end
    end
    del=0;
    c=c+1;
    end
    
    %next data file
    jj=jj+1;
end

block=[];
%combine the different block
for i=1:qty
    block=horzcat(block,blockp(i).mat);
end


%readjust column size to the new combined block 
column=numel(block);

%find unique symbols of the matrix
% symbol=unique(block);



%%%%%%%%%%%%%%%%%%%%            PPM               %%%%%%%%%%%%%%%%%%%%%%%%%
block=[123,123,104,123,123,104,114,125,130,132,114,125,130,104,125];
%initialising table
n=order;
for o=1:n+1
    k(o)=1;
end

for i=2:n+1
    table{i,1}{1,1}='0';
    table{i,2}(1,1)=0;
end

for i=1:71
    table{1,1}{i,1}=num2str(i+79);
    table{1,2}(i,1)=1;
end

%Arith_enco set InitialTable
for i=1:71
    init_sym(i,1)=i+79;
    init_sym(i,2)=table{1,2}(i,1)/sum(table{1,2});
end
lows = 0; highs = init_sym(1,2);
for h = 2:length(init_sym)
    lows(h) = lows(h-1)+ init_sym(h-1,2);
    highs(h) = highs(h-1)+ init_sym(h,2);
end
InitialTable = [init_sym(:,1) lows' highs'];
LowEnc = 0; HighEnc = 1; Progress = [LowEnc HighEnc]; prob_message = 1;

%encode
for i=1:5
    %search block(i) in table, if found, encode
    %search the highest order first, if not found, order-1
    table_enco=table;
    P=1;
    countnum=1;
    counted=[];    
    
    if i<n+1        %not enough data to complete order n
        j=1;
    else            %enough data to complete order n
        j=i-n;
    end
    
    while j<=i
        seq=block(j:i); 
        seqstr=num2str(seq);
        contextorder=i-j+1;     %actually should -1 
        len=size(table_enco{contextorder});
        found=0;
        for L=1:len(1,1)
            if strcmp(seqstr,table_enco{contextorder,1}(L,1))
                found=1;
                break
            else
                found=0;
            end
        end
        
        %delet all irrelated context and upgrade table
        comp_block=num2str(seq(1:length(seq-1)));
        lastele=num2str(block(i));
        comp_block=comp_block(1:(length(comp_block)-length(lastele)-2));
        for m=1:len(1,1)
            comp=table_enco{contextorder,1}{m,1};
            comp=comp(1:(length(comp)-length(lastele)-2));
            
            tablelast_all=str2num(table_enco{contextorder,1}{m,1});
            tablelast=tablelast_all(end);     
            
            if ~strcmp(comp_block,comp) || ~isempty(find(counted==tablelast, 1))
                table_enco{contextorder,1}{m,1}=[];
                table_enco{contextorder,2}(m,1)=0;
            end
        end        

        if found==1
            P=P*(table_enco{contextorder,2}(L,1)/(sum(table_enco{contextorder,2}(:,1))+1));
            table{contextorder,2}(L,1)=table{contextorder,2}(L,1)+1;
            break
        else
            P=P*(1/(sum(table_enco{contextorder,2}(:,1))+1));
            
            table{contextorder,1}(k(contextorder),1)={seqstr};
            table{contextorder,2}(k(contextorder),1)=1;
            k(contextorder)=k(contextorder)+1;  
            
            pos=find(table_enco{contextorder,2}(:,1)~=0);
            for t=1:numel(pos)
                count=str2num(table_enco{contextorder,1}{pos(t),1});
                counted(countnum)=count(end);
                countnum=countnum+1;
            end     
        end
        j=j+1;
    end
    
    if i==1
        current_symbol = block(i);
        j = find( init_sym(:,1) == current_symbol );
        los = lows(j); his = highs(j);
        current_interval = HighEnc-LowEnc;
        HighEnc = LowEnc + current_interval*his;
        LowEnc = LowEnc + current_interval*los;
        Progress = [Progress; LowEnc HighEnc];
        prob_message = prob_message * init_sym(j,2);  
    else
        current_symbol = block(i);
        
    end
    
    
end





% %table update
% for i=1:length(block)
%     %search block(i) in table, if found, encode
%     %search the highest order first, if not found, order-1
%     
%     if i<n+1        %not enough data to complete order n
%         j=1;
%     else            %enough data to complete order n
%         j=i-n;
%     end
%     
%     while j<=i
%         seq=block(j:i);
%         seqstr=num2str(seq);
%         contextorder=i-j+1;     %actually should -1 
%         len=size(table{contextorder});
%         found=0;
%         for L=1:len(1,1)
%             if strcmp(seqstr,table{contextorder,1}(L,1))
%                 found=1;
%                 break
%             else
%                 found=0;
%             end
%         end
% 
%         %if not found
%         if found==0
%             %update table
%             table{contextorder,1}(k(contextorder),1)={seqstr};
%             table{contextorder,2}(k(contextorder),1)=1;
%             k(contextorder)=k(contextorder)+1;                    
%             j=j+1;
%         %if found
%         else
%             %update table
%             table{contextorder,2}(L,1)=table{contextorder,2}(L,1)+1;
%             j=j+1;
%         end
%     end
% end
% 
% 
% %encode
% 



% 
% intLowEnc = floor(LowEnc*2^40);
% binLowEnc = bitget( intLowEnc, 40:-1:1);
% intHighEnc = floor(HighEnc*2^40);
% binHighEnc = bitget( intHighEnc, 40:-1:1);
% 
% ind = find(binHighEnc~= binLowEnc,1);
% ind1 = find( binLowEnc(ind:end) == 0,1);
% ind2 = ind+ind1-1;
% binLowEnc(ind2) = 1;
% 
% code=binLowEnc(1:ind2);
% message = binLowEnc(1:ind2)*(2.^(-(1:ind2)'));
% ideal_codelength = - log2(prob_message);
% 
% 





