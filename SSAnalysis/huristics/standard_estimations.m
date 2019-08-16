clear;
clear;
myDir=uigetdir;
myFiles = dir(fullfile(myDir,'*.mat'));
musicDir='~/Documents/MusicLibrary/WAV1/';
no_plots=5;
Hexponent_RS={};
Hexponent_VT={};
Hexponent_periodogram={};
Hexponent_DFA={};
for k=1:length(myFiles)
%     pinkNoise = dsp.ColoredNoise(1,44.1e3,1);
%     rng default;
%     file_STL{2}=pinkNoise();
%     whiteNoise = dsp.ColoredNoise(0,44.1e3,1);
%     rng default;
%     file_STL{1}=whiteNoise();
%     cn = dsp.ColoredNoise('Color','brown','SamplesPerFrame',10048,...
%     'NumChannels',1);
%     file_STL{3}=cn();
   
    %music_filenames(1).name= 'white noise';
    %music_filenames(2).name= 'Pink noise';
    %music_filenames(3).name= 'Bownian Motion';
    fullFileName{k} = fullfile(myDir,myFiles(k).name);
    musicsubDir=strcat(musicDir ,erase( myFiles(k).name, '.mat'));
    music_filenames= dir(fullfile(musicsubDir, '*.wav'));
    fprintf(1, 'Now reading %s\n', fullFileName{k});
    clear file_STL;
    load (fullFileName{k});
    l=1;
    figure(); 
    pts=10:50:1000;

    for ia = 1:no_plots*4
     ax(ia) = subplot(no_plots,4,ia);
    end
    pos = get(ax, 'position');
    dim = cellfun(@(x) x.*[1 1 1 .2], pos, 'uni',0);
    for i=1: size(file_STL,2)
       clear HRS; clear xs; clear R; clear S;clear Xs; clear Yfits; clear index; clear n; 
       clear varExponent; clear x; clear y; clear yfit; clear alpha_DFA; clear DFA;
       [HRS ,xs, R, S, Xs, Yfits, index,n]=RS(file_STL{i}); %if RS method
       Hexponent_RS{k,i}=HRS; %if RS method
       [varExponent, x, y, yfit]= var_analysis(file_STL{i}); % if variance method
       Hexponent_VT{k,i}= varExponent; %if variance method
       [Hexponent_per, xp,yp, X, Yfit] = periodogram_method((file_STL{i}));%if periodogram
       Hexponent_periodogram{k,i}=Hexponent_per;%if periodogram
       [alpha_DFA, DFA]= DFA_fun( file_STL{i}, pts, 2 );% if DFA
       Hexponent_DFA{k,i}=alpha_DFA(1); %if DFA
       if i <= no_plots
           %%%%% Integrated Variance%%%%%
           plot(ax(l), x, y, 'b.'); hold (ax(l), 'on');
           plot(ax(l), x, yfit, 'r--');
           ylabel(ax(l), music_filenames(i).name(6:min (25, length(music_filenames(i).name))), 'rotation',0, 'VerticalAlignment','middle', 'HorizontalAlignment','right', 'fontWeight','bold');
           annotation('textbox',dim{l},'String',  strcat ("H= " , num2str(round(Hexponent_VT{k,i},2)), " , Alpha= ",num2str(round(2*Hexponent_VT{k,i}-2,2))), 'FontSize',18);           %p= polyfit(log10(x(inds)), log10(y(inds)),1);
           set(ax(l),'FontSize',18);
           %%%% Periodogram %%%%%%
           plot(ax(l+1), xp,yp,'b.'); hold (ax(l+1), 'on'); plot(ax(l+1), X,Yfit,'r-','LineWidth',3);
           annotation('textbox',dim{l+1},'String',  strcat ("H= " , num2str(round(Hexponent_periodogram{k,i},2)), " , Beta= ",num2str(round(2*Hexponent_periodogram{k,i}-1,2))), 'FontSize',18);           %p= polyfit(log10(x(inds)), log10(y(inds)),1);
           set(ax(l+1),'FontSize',18);
           %%%%DFA%%%%%
           scatter(ax(l+2),log10(pts),log10(DFA));
           annotation('textbox',dim{l+2},'String',  strcat ("H= " , num2str(round(Hexponent_DFA{k,i},2)), " , Alpha= ",num2str(round(2*Hexponent_DFA{k,i}-2,2))), 'FontSize',18);           %p= polyfit(log10(x(inds)), log10(y(inds)),1);
           set(ax(l+2),'FontSize',18);
           %%%%RS%%%%
           N=length(file_STL{i});
           bound = ceil(log10(N));
           axis([0 bound 0 0.75*bound]);
           temp = (1:n).*index;
           index = temp(index);
           n2=length(Xs);
           for j = 1:n2
              plot(xs(j),log10(R{index(j)}./S{index(j)}),'b.');
           end
           xs2 = linspace(0,bound,10);
           y1s = 0.5*xs2;
           y2s = xs2;
           u1 = plot(ax(l+3), xs2,y1s,'b--','LineWidth',2); hold(ax(l+3),'on')
           u2= plot( ax(l+3), xs2,y2s,'b-.','LineWidth',2); hold(ax(l+3),'on')
           plot(ax(l+3), Xs,Yfits,'r-','LineWidth',3);
           
           annotation('textbox',dim{l+3},'String',  strcat ("H= " , num2str(round(Hexponent_RS{k,i},2)), " , Alpha= ",num2str(round(2*Hexponent_RS{k,i}-2,2))), 'FontSize',18);           %p= polyfit(log10(x(inds)), log10(y(inds)),1);
           set(ax(l+3),'FontSize',18);
           %%%
           l=l+4;
           if i==1
               h1= get(ax(1),'title') ;
               set(h1, 'String','Aggregate Variance');
               h2= get(ax(2),'title');
               set(h2,'String', 'periodogram');
               h3= get(ax(3),'title') ;
               set(h3,'String', 'DFA');
               h4= get(ax(4),'title') ;
               set(h4,'String', 'RS');
               h5= get(ax(4),'legend') ;
               set(h5, [u1,u2],'slope 1/2','slope 1')
           end
           if  i== no_plots
               x1= get(ax(l-4),'xlabel'); 
               set( x1, 'String', 'log10(k)');
               x2= get(ax(l-3),'xlabel');
               set( x2, 'String', 'log10(S_v)');
               
               x3= get(ax(l-2),'xlabel');
               set( x3, 'String', 'log10(n)');
               
               x4= get(ax(l-1),'xlabel');
               set( x4, 'String', 'log10(m)');
               %ax1.ylabel='log(Var(\bar(X)_k))';
               
           end
           set(gca,'FontSize',18,'fontWeight','bold')
           %set(findall(gcf,'type','text'),'FontSize',18,'fontWeight','bold')
     
       else
           stop=1; 
       end
    end
%     mean_HRS(k)=mean(Hexponent_RS);std_HRS(k)= std(Hexponent_RS);
%     mean_HVT(k)=mean(Hexponent_VT);std_HVT(k)= std(Hexponent_VT);
%     mean_Hper(k)=mean(Hexponent_periodogram);std_Hper(k)= std(Hexponent_periodogram);
%     mean_HDFA(k)=mean(Hexponent_DFA);std_HDFA(k)= std(Hexponent_DFA);

  
end

Hexponent_RS_all=[Hexponent_RS{:}]; 
Hexponent_per_all=[Hexponent_periodogram{:}]; 
Hexponent_var_all=[Hexponent_VT{:}]; 
Hexponent_DFA_all=[Hexponent_DFA{:}];

%%% Identify excerpts
composer_name={};
excerpt_name={};
no_files_per_composer=[];
no_file_sums=[];
p=1;
for r=1:size(fullFileName,2)
    clear music_filenames
    tmp=fullFileName{1,r};
    composer_name{r}=tmp(50: length(tmp)-4);
    fullFileName{r} = fullfile(myDir,myFiles(r).name);
    musicsubDir=strcat(musicDir ,erase( myFiles(r).name, '.mat'));
    music_filenames= dir(fullfile(musicsubDir, '*.wav'));
    no_files_per_composer(r)=size(music_filenames,1);
    if r==1
        no_file_sums(r)= no_files_per_composer(r);
    else
        no_file_sums(r)= no_files_per_composer(r)+ no_file_sums(r-1);
    end
    for q=1:size(music_filenames,1)
        excerpt_name{r,q}= music_filenames(q).name(6:min (14, length(music_filenames(q).name)));
        p=p+1;
    end
end

no_file_sums=[0  no_file_sums];
no_clusters=2;
[idx_DFA, C_DFA]=kmeans(Hexponent_DFA_all',no_clusters);

cluster_file_mat=zeros(no_clusters, size(excerpt_name,1), size(excerpt_name,2));
for c=1: no_clusters
    cluster_indices=find(idx_DFA==c);
    for ci=1:length(cluster_indices)
      nearest_ind=max(find (cluster_indices(ci)- no_file_sums' >0 ));
      file_number=cluster_indices(ci)- no_file_sums( nearest_ind);
      cluster_file_mat(c,nearest_ind,file_number)=1;
    end
end

figure;
colors={'b','r','k','g'};
legends={'c=1','c=2', 'c=3', 'c=4'};
for c=1: no_clusters
    [xvector, yvector]=find(reshape(cluster_file_mat(c,:,:),size(excerpt_name,1), size(excerpt_name,2))==1);
    plot(xvector, yvector, 'o', 'color', colors{c}); hold on;
end    
%%%% By composer %%%%
composer_name={};
figure();
errorbar(mean_HRS, std_HRS, 'r--', 'DisplayName','RS Method'); hold on;
errorbar(mean_HVT, std_HVT, 'b--', 'DisplayName','Aggregated Variance'); hold on;
errorbar(mean_Hper, std_Hper, 'k--', 'DisplayName','Periodogram'); hold on;
errorbar(mean_HDFA, std_HDFA, 'g--', 'DisplayName','DFA'); hold on;
for r=1:size(fullFileName,2)
   tmp=fullFileName{1,r};
   composer_name{r}=tmp(50: length(tmp)-4);
end
set(gca,'xtick',[1:20],'xticklabel',composer_name)
xtickangle(45)
set(gca,'FontSize',20,'fontWeight','bold')
ylabel('Hurst Exponent')
legend

fin=1;