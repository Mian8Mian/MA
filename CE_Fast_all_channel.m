
%Coding

clear;

%set file path
qty=100;      %quantity of files
filelist=dir('/home/muyan/Documents/masterarbeit/ACE_AUTO');
dest='/home/muyan/Documents/CE_Fast/All channels/100f6o';   

%define row number and order
rownr=[1 6];        %[startrownr endrownr]
order=6;
rowsize=rownr(1,2)-rownr(1,1)+1;

%select random files
randf=randperm(numel(filelist),qty);

for i=1:qty
    while (randf(i)==1)||(randf(i)==2)
        randf(i)=randperm(numel(filelist),1);
    end
end

%copy the selected files to destination folder
for i=1:qty
    source=fullfile(filelist(randf(i)).folder,filelist(randf(i)).name);
    copyfile(source,dest);
end

%get filelist of the current destination folder
cfiles=dir(dest);

%initialization for clone struct
clone=struct;

%initialization for blockog and blockp struct
blockog=struct;
blockp=struct;

%cfiles(1) and cfiles(2) are '.' and '..', file name is stored from cfiles(3)
j=3;
while j<=qty+2

    %load data file
    load(fullfile(cfiles(j).folder,cfiles(j).name));
    
    %initialization, convert NaN to 0
    clone(j-2).file=cfiles(j).name;
    clone(j-2).mat=mat;

    %next data file
    j=j+1;
end

%cfiles(1) and cfiles(2) are '.' and '..', file name is stored from cfiles(3)
j=3;
while j<=qty+2
    
    %initialization, convert NaN to 0
    clone(j-2).mat(isnan(clone(j-2).mat))=0;
     
    %get the required block
    blockog(j-2).file=clone(j-2).file;
    blockog(j-2).mat=clone(j-2).mat(rownr(1,1):rownr(1,2),:);   
    blockp(j-2).file=clone(j-2).file;
    blockp(j-2).mat=clone(j-2).mat(rownr(1,1):rownr(1,2),:);    
    %next data file
    j=j+1;
end

%delete unneeded zeros
j=3;
while j<=qty+2
    remainzero=zeros(rowsize,11); 
    c=1;
    del=0;

    while c<=length(blockp(j-2).mat)       
        if ((c+length(remainzero)-1)<length(blockp(j-2).mat))       
            bufferzero=[];
            sizereze=size(remainzero);
            if (blockp(j-2).mat(:,c:c+sizereze(1,2)-1)==remainzero) 
                c=c+length(remainzero);         
                del=1;        
            end     
        end 
        while c<=length(blockp(j-2).mat)    
            if blockp(j-2).mat(:,c)==zeros(rowsize,1)
                if (del==1)
                    blockp(j-2).mat(:,c)=[];   
                else
                    break
                end
            else
                break
            end
        end
        del=0;
        c=c+1;  
    end  

    %next data file
    j=j+1;
end

block=[];
%combine the different block
for i=1:qty
    block=horzcat(block,blockp(i).mat);
end

%readjust column size to the new combined block 
column=length(block);

H_ALL=zeros(order+1,rowsize);
for b=1:rowsize
    
    %find unique symbols
    symbol=unique(block(b,:));
    
    %calculate conditional entropy
    for n=1:order
        
        %find every symbol's position
        position=zeros(1,numel(block(b,:)));
        for i=1:numel(block(b,:))
            position(i)=find(symbol==block(b,i));
        end   
        
        %build symbol combination of rownr(1,2)-rownr(1,1)+1n.order
        combination_n=struct;

        k=1;
        while k<=numel(block(b,:))-n
            combination_n(k).combination=block(b,k:k+n);
            combination_n(k).position=0;
            for i=1:n
                combination_n(k).position=combination_n(k).position+((position(k+i-1)-1)*numel(symbol)^(n-i+1));
            end
            combination_n(k).position=combination_n(k).position+position(k+n);
            k=k+1;
        end
        
        buffermat_n=[combination_n.position];
        [symbol_n,~,pos_n]=unique(buffermat_n);
        occurrences_n=[symbol_n' accumarray(pos_n,1)];
        occurrences_n_spec=occurrences_n;
        
        %build symbol combination of (n-1).order
        combination_nm=struct;

        k=1;
        while k<=numel(block(b,:))-n+1
            combination_nm(k).combination=block(b,k:k+n-1);
            combination_nm(k).position=0;
            for i=1:n-1
                combination_nm(k).position=combination_nm(k).position+((position(k+i-1)-1)*numel(symbol)^(n-i));
            end
            combination_nm(k).position=combination_nm(k).position+position(k+n-1);
            k=k+1;
        end

        buffermat_nm=[combination_nm.position];
        [symbol_nm,~,pos_nm]=unique(buffermat_nm);
        occurrences_nm=[symbol_nm' accumarray(pos_nm,1)];        

        %calculate Hy for (n-1).order
        occnmsize=size(occurrences_nm);
        Hy=zeros(occnmsize(1,1),1);
        if mod(occurrences_n(1,1),numel(symbol))~=0
            k=idivide(uint64(occurrences_n(1,1)),uint64(numel(symbol)))+1;
        else
            k=idivide(uint64(occurrences_n(1,1)),uint64(numel(symbol)));
        end
        j=1;
        %l=1;
        %cp=struct;
        while ~isempty(occurrences_n)
            occnr=find(occurrences_n(:,1)<=numel(symbol)*k,1,'last');
            if ~isempty(occnr)     
                occsum=sum(occurrences_n(1:occnr,2));
                p_n=zeros(occnr,1);
                for i=1:occnr
                    p_n(i)=occurrences_n(i,2)/occsum;
                    %cp(l).position=occurrences_n(i,1);
                    %cp(l).probability=p_n(i);
                    %l=l+1;
                end
                Hy(j)=-dot(p_n,log2(p_n));
                j=j+1;
                restsize=size(occurrences_n);
                occurrences_n=occurrences_n(occnr+1:restsize(1,1),:);
                if ~isempty(occurrences_n)
                    if mod(occurrences_n(1,1),numel(symbol))~=0
                        k=idivide(uint64(occurrences_n(1,1)),uint64(numel(symbol)))+1;
                    else
                        k=idivide(uint64(occurrences_n(1,1)),uint64(numel(symbol)));
                    end
                end
            end
        end

        %calculate probabilities of (n-1).order
        p_nm=zeros(occnmsize(1,1),2);
        for k=1:occnmsize(1,1)
            p_nm(k,1)=occurrences_nm(k,1);
            p_nm(k,2)=occurrences_nm(k,2)/numel(buffermat_nm);
        end
        
        %final calculation of conditional entropy
        if n==1
            H_ALL(1,b)=-dot(p_nm(:,2),log2(p_nm(:,2)));
        end
        H_ALL(n+1,b)=sum(Hy.*p_nm(:,2));    
    end
    
end

lgd=cell(rowsize,1);
hold on;
for i=1:rowsize
    plot(0:order,H_ALL(:,i),'--*');
    lgd{i}=num2str(rownr(1,1)+i-1);
end

% legend(lgd);
% legend('boxoff');
% xticks(0:1:order);  
% title(legend,'Channel NO.');
% title('Conditional Entropy');
% xlabel('n.Order');
% ylabel('bits/symbol');

save('100f6ods1c1to6');
