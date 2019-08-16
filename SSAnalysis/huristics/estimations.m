clear;
myDir=uigetdir;
myFiles = dir(fullfile(myDir,'*.mat'));
musicDir='~/Documents/MusicLibrary/WAV1/';
HSTL_all={};
no_plots=8;
for k=2:3 %length(myFiles)
    file_STL={};
    Hexponent_RS=[];
    Hexponent_VT=[];
    Hexponent_Voss=[];
    Hexponent_DFA=[];
    fullFileName = fullfile(myDir,myFiles(k).name);
    musicsubDir=strcat(musicDir ,erase( myFiles(k).name, '.mat'));
    music_filenames= dir(fullfile(musicsubDir, '*.wav'));
    fprintf(1, 'Now reading %s\n', fullFileName);
    load (fullFileName);
    l=1;
    figure();
    pts=10:500:5000;
    for i=1:no_plots % size(file_STL,2) 
       Hexponent_RS=[Hexponent_RS RS(file_STL{i})]; %if RS method
       [varExponent, x, y, inds]= var_analysis(file_STL{i}); % if variance method
       Hexponent_VT=[Hexponent_VT varExponent]; %if variance method
       [ xprim , yprim , ~ , alpha , ~ , ~ ] = Voss(file_STL{i});
       Hexponent_Voss=[Hexponent_Voss (1-alpha)/2]; % if spectral density- Voss
       [F, Pxx]=periodogram_analysis(file_STL{i});
       [alpha_DFA, DFA]= DFA_fun( file_STL{i}, pts, 2 );
       Hexponent_DFA=[Hexponent_DFA alpha_DFA(1)]; %if variance method
       if i <= no_plots
           ax1=subplot(no_plots,4,l);
           plot(ax1, log10(x), log10(y), 'b.'); hold on;
           p= polyfit(log10(x(inds)), log10(y(inds)),1);
           f=polyval(p, log10(x(inds)));
           plot(ax1, log10(x(inds)), f, 'r--')
           ylabel(music_filenames(i).name(6:min (16, length(music_filenames(i).name(i)))));
           ax2=subplot(no_plots,4,l+1);
           plot(ax2,xprim, yprim);
           ax3=subplot(no_plots,4,l+2);
           plot(ax3, log10(F),log10(Pxx),'b.');
           ax4=subplot(no_plots, 4,l+3);
           scatter(ax4,log10(pts),log10(DFA));
           l=l+4;
           if i==1
               h1= get(ax1,'title') ;
               set(h1, 'String','Aggregate Variance vs Sample size');
               h2= get(ax2,'title') ;
               set(h2,'String', 'Smoothed Spectral density');
               h3= get(ax3,'title') ;
               set(h3,'String', 'periodogram');
               h4= get(ax4,'title') ;
               set(h4,'String', 'DFA');
           end
           if  i== no_plots
               x1= get(ax1,'xlabel'); 
               set( x1, 'String', 'log(k)');
               x2= get(ax2,'xlabel');
               set( x2, 'String', 'log(f)');
               x3= get(ax3,'xlabel');
               set( x3, 'String', 'log(S_v)');
               x4= get(ax4,'xlabel');
               set( x4, 'String', 'log(n)');
               %ax1.ylabel='log(Var(\bar(X)_k))';
               
           end
     
       else
           stop=1; 
       end
    end
    mean_HRS(k)=mean(Hexponent_RS);std_HRS(k)= std(Hexponent_RS);
    mean_HVT(k)=mean(Hexponent_VT);std_VT(k)= std(Hexponent_VT);
    mean_HVoss(k)=mean(Hexponent_Voss);std_VT(k)= std(Hexponent_Voss);

  
end