function layer= Regiondetection(year,month,day,starttime,endtime)
load(['Madrigal','_',num2str(day),'_',num2str(month),'_',num2str(year),'.mat']);
load(['Madrigalac','_',num2str(day),'_',num2str(month),'_',num2str(year),'.mat']);
%load('Madrigal19Mar2015.mat');
timestamp=unique(Madrigal(:,18));
timestampac=unique(Madrigalac(:,18));
Maxlp=[];
Maxac=[];
save_maxNe_lp=[];
save_maxNe_ac=[];
scintrows=find(timestamp>=starttime & timestamp<endtime);
scinttime=timestamp(scintrows);
scintrowsac=find(timestampac>=starttime & timestampac<endtime);
scinttimeac=timestampac(scintrowsac);
E=0;
F=0;
I=0;

% LP

for i=1:length(scinttime)
rows=find(Madrigal(:,1)==-154.3&Madrigal(:,14)==77.5&Madrigal(:,18)==scinttime(i)&Madrigal(:,4)>=200);
%rows=find(Madrigal(:,1)==-34.69&Madrigal(:,14)==66.09&Madrigal(:,18)==scinttime(i));% for 09Feb2014
%rows=find(Madrigal(:,18)==scinttime(i));
Maxlp=max(Madrigal(rows,8));
maxrow=find(Madrigal(rows,8)==Maxlp);
altitude=Madrigal(rows(maxrow),4);
save_maxNe_lp=[save_maxNe_lp;[Madrigal(rows(maxrow),18),Madrigal(rows(maxrow),4),Madrigal(rows(maxrow),8)]];
% data=[Madrigal(rows,18),Madrigal(rows,4),Madrigal(rows,8)];
% MaxNe=[MaxNe;data];
% for j=1:length(altitude)
%     if(altitude(j)<=200)
%         E=E+1;
%     elseif(altitude(j)>200)
%         F=F+1;
%    
%     end
% end
end

%AC

for i=1:length(scinttimeac)
rowsac=find(Madrigalac(:,1)==-154.3&Madrigalac(:,14)==77.5&Madrigalac(:,18)==scinttimeac(i)&Madrigalac(:,4)<200);
Maxac=max(Madrigalac(rowsac,8));
maxrowac=find(Madrigalac(rowsac,8)==Maxac);
altitudeac=Madrigalac(rowsac(maxrowac),4);
save_maxNe_ac=[save_maxNe_ac;[Madrigalac(rowsac(maxrowac),18),Madrigalac(rowsac(maxrowac),4),Madrigalac(rowsac(maxrowac),8)]];
end

%Comparing LP and AC

if (size(save_maxNe_ac,1)<=size(save_maxNe_lp,1))
    for i=1:size(save_maxNe_ac,1) %change scinttime to length of save_maxNe_ac
        for j=1:size(save_maxNe_lp,1)
            diff1=abs(save_maxNe_ac(i,1)-save_maxNe_lp(j,1));
            if (j<size(save_maxNe_lp,1))
            diff2=abs(save_maxNe_ac(i,1)-save_maxNe_lp((j+1),1));
            else
                diff2=diff1;
                end
            if(diff2<diff1)
                continue;
            else
            if (diff1==0 || diff1<=0.0167) %Also check the diff with next timestamp and check which one is the least
                if (save_maxNe_ac(i,3)>save_maxNe_lp(j,3)&&save_maxNe_ac(i,2)<150)
                    E=E+1;
                    break;
                elseif (save_maxNe_ac(i,3)<save_maxNe_lp(j,3))
                    F=F+1;
                    break;
                elseif(save_maxNe_ac(i,3)>save_maxNe_lp(j,3)&&save_maxNe_ac(i,2)>150&&save_maxNe_ac(i,2)<200)
                    I=I+1;
                    break;
                end
            else
                continue;
            end
            end
                    
                    
            
        end
    end
    
else
    for i=1:size(save_maxNe_lp,1)
        for j=1:size(save_maxNe_ac,1)
            diff1=abs(save_maxNe_lp(i,1)-save_maxNe_ac(j,1));
             if (j<size(save_maxNe_ac,1))
            diff2=abs(save_maxNe_lp(i,1)-save_maxNe_ac((j+1),1));
             else
                 diff2=diff1;
             end
            %diff=abs(scinttimeac(i)-scinttime(j));
            if(diff2<diff1)
                continue;
            else
            if (diff1==0 || diff1<=0.0167)
                if (save_maxNe_ac(j,3)>save_maxNe_lp(i,3)&&save_maxNe_ac(j,2)<150)
                    E=E+1;
                    break;
                elseif (save_maxNe_ac(j,3)<save_maxNe_lp(i,3))
                    F=F+1;
                    break;
                elseif(save_maxNe_ac(j,3)>save_maxNe_lp(i,3)&&save_maxNe_ac(j,2)>150&&save_maxNe_ac(j,2)<200)
                    I=I+1;
                    break;
                end
            else
                continue;
            end
            end
                    
                    
            
        end
    end
    
end

%Detecting region for a scintillation event on a day
if(E==0&&F==0&&I==0)
    disp('No PFISR data available'); % N- No PFISR
    layer='N';
    
elseif(E>F&&E>I)
 disp('E region scintillation');  % E- E region
 layer='E';
 
elseif(F>E&&F>I)
    disp('F region scintillation'); % F- F region
    layer='F';
    
elseif(I>E&&I>F)
    disp('Transition region scintillation'); % F- F region
    layer='T';
    
else
    disp('region cannot be determined'); % I- Inconclusive
    layer='I';
end
end