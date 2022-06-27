clear all; clc; close all;
%% Input
prot='CypA';
residue=115; % Please give the residue number as per the PDB file

%% Loading data files of the protein
eval(['load Coupling_310K_',prot,'.mat']);
eval(['load BlockDet',prot,'.dat;']); 
eval(['BlockDet= BlockDet',prot,';']);

nres=BlockDet(end,2); % Total number of blocks
T=310; %Temperature
R_const=0.008314; %Gas constant
%% Opening pdb file
pdbfile=strcat(prot,'.pdb');
bufferfile=strcat('buffer_',prot,'.txt');
pdbr=fopen(pdbfile,'rt'); % Open Raw PDB data
pdbc=fopen(bufferfile,'wt'); % Curated PDB data

line=fgetl(pdbr); %Get line
     
     while feof(pdbr)~=1 %Till end of file
        if(length(line)>4) %If line length is greater than 4 characters
            if(strcmp(line(1:4),'ATOM')==1 && strcmp(line(13),'H')==0 && strcmp(line(13),'Q')==0 && strcmp(line(14),'H')==0 && strcmp(line(18),' ')==0)  %selecting lines of atoms that are not H
                % if(strcmp(line(1:4),'ATOM')==1 && strcmp(line(13),'H')==0 && strcmp(line(13),'Q')==0 && strcmp(line(14),'H')==0 && strcmp(line(18),' ')==0 && strcmp(line(13:14),'CA')==1)
                fprintf(pdbc,'%s\n',line); %copying all heavy atoms in buffer.txt
            end
            if(strcmp(line(1:3),'TER')==1) %if end of model or chain
                break;
            end
        end
        line=fgetl(pdbr); %gets line from pdb
     end
      
    fclose(pdbr);
    fclose(pdbc);

%% Getting C-alpha atom coordinates
f1= fopen(bufferfile,'r'); %open text file containing atom details
l=fgetl(f1); %gets line
ca_xyz=[];
while feof(f1)~=1 %till end of file
    if strcmp(deblank(l(14:16)),'CA')==1 %if the C alpha atom
        ca_xyz=[ca_xyz; str2double(l(32:38)) str2double(l(40:46)) str2double(l(48:54))]; %store coordinates
    end
    l=fgetl(f1); %gets line
end
%% Calculating C_alpha distances between residues
ca_dist=zeros(size(ca_xyz,1));
for i=1:size(ca_xyz,1)
    for j=1:size(ca_xyz,1)
        ca_dist(i,j)= sqrt(((ca_xyz(i,1)-ca_xyz(j,1))^2)+((ca_xyz(i,2)-ca_xyz(j,2))^2)+((ca_xyz(i,3)-ca_xyz(j,3))^2));
    end
end

%% Analysing DCI as a function of distance from mutated site
targetres=find(MutatedRes(:,2)==residue);
efitX=0:0.5:40; %high density data points to generate fit curve
param=[];
chiPwt=chipluswt'; %delG+ of WT
ind1 = find(isinf(chiPwt)==1);
chiPwt(ind1)=NaN; %Replacing infinite values with NaN. (Coupling of residue i with the same residue i is infinite)

x=[]; y=[];   x2=[]; y2=[];
x=ca_dist(MutatedRes(targetres,1),:)'; %Storing inter-residue distance from the mutated site
chiP_mut=chiplus(:,:,targetres+1)'; %storing delG+ of the mutant
chiP_mut(ind1)=NaN; %Replacing infinite values with NaN
DDGp_temp = nanmean(chiP_mut)'-nanmean(chiPwt)'; %Taking an average across rows and calculating DCIs of every block of residue
y=DDGp_temp(BlockDet(:,2)); %Calculating DCIs of every residue
    
ind = find(isnan(y)==1); 
y(ind)=[];
x(ind)=[];

for j=1:50 
   ind2=find((x>=(j-1)/1) & (x<j/1)); %Binning values every 1 angstrom
   if isempty(ind2)==0
   x2=[x2;nanmean(x(ind2))]; %averaging binned distance values 
   y2=[y2;nanmean(abs(y(ind2)))]; %averaging binned absolute values of DCIs
   end
end
      
options=fitoptions('a1.*exp(-(1/a2).*x)+a3','Lower',[0,0 0],'Upper',[30, 30, 30],'Startpoint',[2 6 0.1]);
aaa=fit(x2,abs(y2),'a1.*exp(-(1/a2).*x)+a3',options); %fitting the data to the model to obtain initial guesses
beta0=[aaa.a1,aaa.a2 aaa.a3]; %Initial guesses for model parameters
model=@(b1,x1)((b1(1)*exp(-(1/b1(2))*x1))+b1(3));
[beta,resnorm,R,exitflag,output,lambda,J]=lsqcurvefit(model,beta0,x2,abs(y2),[0,0 0],[30 30 30]); %Least square curve fitting
conf_int=nlparci(beta,R,'jacobian',J,'alpha',.32); %Determining 95% confidence intervals
error=abs(conf_int(:,2)-conf_int(:,1))./2'; %standard deviation
[ypred,delta]=nlpredci(model,efitX,beta,R,'Jacobian',full(J),'alpha',.32); %Generate fit curve for high density data points
[ypred_ori,delta_ori]=nlpredci(model,x2,beta,R,'Jacobian',full(J),'alpha',.32); % Generate fit curve for original data points
SStot = sum((abs(y2)-mean(abs(y2))).^2); % Total Sum-Of-Squares
SSres= sum((abs(y2)-ypred_ori).^2); %Sum of squares of residuals
Rsq = 1-(SSres/SStot); %R-squared
    
param=[param;  beta(1) error(1) beta(2) error(2) beta(3) error(3) Rsq]; %Storing parameters

%% Analysing distribution of free energies
Dfe_all=fe_all(:,targetres+1)-fe_all(:,1);
ind=find(fe_all(:,end)>min_str);
Dfe_nw=Dfe_all(ind,:);
[N1,edges1] = histcounts(Dfe_all,'BinWidth',0.05, 'Normalization','count');
edges1 = edges1(2:end) - (edges1(2)-edges1(1))/2;
edges1=[edges1(1)-0.05 edges1 edges1(end)+0.05];
N1=[0 N1 0];
[N2,edges2] = histcounts(Dfe_nw,'BinWidth',0.05, 'Normalization','count');
edges2 = edges2(2:end) - (edges2(2)-edges2(1))/2;
edges2=[edges2(1)-0.05 edges2 edges2(end)+0.05];
N2=[0 N2 0];

%% Generating images
    
 figure(1)
 fe_wt=fesmat(:,1);
 fe_mut=fesmat(:,targetres+1);
 plot(1:nres,fesmat(:,1), 'k', 'LineWidth', 2); axis([0 inf 0 inf]);
 hold on
 plot(1:nres, fesmat(:,targetres+1),'r','LineWidth',2);
 set(gcf, 'Position', [1 1 1000 500]);
 xlabel('# of Structured Blocks');
 ylabel('FE (kJ/mol)')
 title('1D Free Energy Profile')
 xt = 0:10:nres; xt1 = {}; for i=1:length(xt); xt1{i}=num2str(xt(i)); end;
set(gca, 'XTick', xt, 'XTickLabel', xt1);
yt = 0:20:max(fesmat(:,1)); yt1 = {}; for i=1:length(yt); yt1{i}=num2str(yt(i)); end;
set(gca, 'YTick', yt, 'YTickLabel', yt1);
set(gca,'fontweight','bold','TickDir','out', 'LineWidth', 2, 'FontSize', 14)
a = gca;
set(a,'box','off','color','none')
b = axes('Position',get(a,'Position'),'box','on','xtick',[],'ytick',[]);
set(b,'LineWidth', 2)
axes(a)
linkprop([a b],'position');
xlim([1, nres])
set(gcf,'PaperPositionMode','auto');
print('FE Profile', '-dpng', '-r300')

      
 figure(2)
 cmapmax=max(max(max(chiPwt)),max(max(chiP_mut)));
 cmapmin=min(min(min(chiPwt)),min(min(chiP_mut)));
 set(gcf,'Position',[1,3,1500,750])
 subplot(1,3,1)
 pcolor(chiPwt(BlockDet(:,2),BlockDet(:,2))'); shading interp; colormap jet
 ax=gca;
 ax.LineWidth=1;
 ax.XTick=20:20:BlockDet(end,1);
 ax.YTick=20:20:BlockDet(end,1);
 ax.FontWeight='Bold';
 ax.Position=[0.1 0.35 0.21 0.4];
 ax.Box='on';
 cb=colorbar;
 cb.Position=[0.32 0.3510 0.010 0.3990];
 caxis([cmapmin cmapmax]);
 xlabel('Sequence Index'); ylabel('Sequence Index');
 title('delG_+ WT')
    
 subplot(1,3,2)
 pcolor(chiP_mut(BlockDet(:,2),BlockDet(:,2))'); shading interp; colormap jet
 ax=gca;
 ax.LineWidth=1;
 ax.XTick=20:20:BlockDet(end,1);
 ax.YTick=20:20:BlockDet(end,1);
 ax.FontWeight='Bold';
 ax.Position=[0.4 0.35 0.21 0.4];
 ax.Box='on';
 xlabel('Sequence Index'); ylabel('Sequence Index');
 cb1=colorbar;
 cb1.Position=[0.62 0.3510 0.010 0.3990];
caxis([cmapmin cmapmax]);
 title('delG_+ mutant')
    
subplot(1,3,3)
DCM=chiP_mut(BlockDet(:,2),BlockDet(:,2))'-chiPwt(BlockDet(:,2),BlockDet(:,2))';
pcolor(DCM); shading interp; colormap jet
ax=gca;
ax.LineWidth=1;
ax.XTick=20:20:BlockDet(end,1);
ax.YTick=20:20:BlockDet(end,1);
ax.FontWeight='Bold';
ax.Position=[0.7 0.35 0.21 0.4];
ax.Box='on';
xlabel('Sequence Index'); ylabel('Sequence Index');
cb=colorbar;
cb.Position=[0.92 0.3510 0.010 0.3990];
title('Differential Coupling Matrix')
hold off
print('Coupling Matrices', '-dpng', '-r300')

figure(3)
ax=gca;
plot(x2,abs(y2),'bo','LineWidth',1)
hold on
plot(efitX, ypred,'r-','LineWidth',2)
xlabel('Distance from the Mutated Residue')
ylabel('<|DDG_+|>')
txt= {['d_c = ' num2str(round(param(3)*100)/100) ' ± ' num2str(round(param(4)*100)/100) ' Angstrom']};
text(25,max(abs(y2))-0.5,txt)
title('Exponential Decay of DDG_+')
ax.FontWeight='Bold';
hold off
print('Exponential decay of DDGp', '-dpng', '-r300')

figure(4)
ax=gca;
plot(edges1,N1,'b-','LineWidth',1.5)
xlabel('Difference in FE (kJ/mol)')
ylabel('No of Microstates')
title('Free Energy Redistribution of All Microstates')
ax.FontWeight='Bold';
print('FE Redistribution ', '-dpng', '-r300')

figure(5)
ax=gca;
plot(edges2,N2,'r-','LineWidth',1.5)
xlabel('Difference in FE (kJ/mol)')
ylabel('No of Microstates')
title('Free Energy Redistribution of Folded Well Microstates')
ax.FontWeight='Bold';
print('Folded Well FE Distribution ', '-dpng', '-r300')

%% Generating output pdb file
outfile=strcat(prot,'_',char(mutlist(targetres)),'col.pdb');
DCI=y;
aa=fopen(pdbfile,'r');
bb=fopen(outfile,'w');
l=fgetl(aa); flag = 1; k=1;
while l>0
    if (length(l)>=4 && strcmp(l(1:4),'ATOM'))
        resno=str2double(l(24:26));
        if flag == 1; currres = resno; flag = 2; end
        if (currres-resno~=0); k=k+1; currres=resno; end
            v=num2str(round(DCI(k)*10000)/10000);
        for i=1:10-length(v)
            v = cat(2,v,' ');
        end
        l(61:61+length(v)-1)=v;
    end
    fprintf(bb,'%s\n',l);
    l=fgetl(aa);
end
fclose(aa);
fclose(bb);

 ff = fopen('Fit Parameters.txt','wt');
fprintf(ff,'REMARK %% Fit Parameters of %s mutant\n',char(mutlist(targetres)));
fprintf(ff,'\nAmplitude : %2.2f ± %2.2f\n ',param(1),param(2));
fprintf(ff,'Coupling distance : %2.2f ± %2.2f\n',param(3), param(4));
fprintf(ff,'Shift term : %2.2f ± %2.2f\n',param(5), param(6));
fprintf(ff,'R-squared  : %2.2f\n',param(7));
fclose(ff);
