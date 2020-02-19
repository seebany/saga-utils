function [combo1, combo2, combo3, combo1region, combo2region, combo3region]=L1events(year,doy,freq)
load(['L',num2str(freq),'_S4_',num2str(year,'%04i'),'_',num2str(doy,'%03i'),'.mat']);
S4_hr=TSP_hrv0;
load(['L',num2str(freq),'_SP_',num2str(year,'%04i'),'_',num2str(doy,'%03i'),'.mat']);
SP_hr=TSP_hrv0;
S4prn=[];
SPprn=[];
if ~isempty(S4_hr)
S4prn=S4_hr(:,1);
end
if ~isempty(SP_hr)
SPprn=SP_hr(:,1);
end
RINEX.year=year;
RINEX.day=doy;
[gps,abs,utc]=convert_rinex(RINEX,0); %leap =0 as 2014 and 2015 are not leap years
month=utc.mon;
day=utc.day;

%____________________________SP and S4_____________________________________
S4SP=[];
if (~isempty(SP_hr)&&~isempty(S4_hr))
     S4SPprnlist=intersect(S4prn,SPprn);
    
   for i=1:length(S4SPprnlist)
       S4row=find(S4prn==S4SPprnlist(i));
       SProw=find(SPprn==S4SPprnlist(i));
       for j=1:length(S4row)
           
           S4startdatenum=datenum(year,month,day,S4_hr(S4row(j),2),S4_hr(S4row(j),3),0);
           S4enddatenum=datenum(year,month,day,S4_hr(S4row(j),4),S4_hr(S4row(j),5),0);
           
           for k=1:length(SProw)
               SPstartdatenum=datenum(year,month,day,SP_hr(SProw(k),2),SP_hr(SProw(k),3),0);
               SPenddatenum=datenum(year,month,day,SP_hr(SProw(k),4),SP_hr(SProw(k),5),0); 
               %Common time interval between SP and S4 for same prn 
               maximum_startdatenum=max(S4startdatenum,SPstartdatenum);
               minimum_enddatenum=min(S4enddatenum,SPenddatenum);
               if(maximum_startdatenum<minimum_enddatenum)
                   common_startdatenum=maximum_startdatenum;
                   common_enddatenum=minimum_enddatenum;
                   common_startdatevec=datevec(common_startdatenum);
                   common_enddatevec=datevec(common_enddatenum);
                   S4SP=[S4SP;[year doy freq S4SPprnlist(i) common_startdatevec(4) common_startdatevec(5) common_enddatevec(4) common_enddatevec(5)]]; 
               else
                   disp('No time intervals common to both S4 and SP acintillation');
               end
           end      
       end     
   end 
end

%______________________ only S4 ______________________________________________
S4=[];
if isempty(S4_hr)
    disp(['No Amplitude scintillation on the day',num2str(doy,'%03i'),' of the year',num2str(year,'%04i')]);
else

onlyS4prnlist=setdiff(S4prn,SPprn);

for i=1:length(onlyS4prnlist)
    row=find(S4_hr(:,1)==onlyS4prnlist(i));
for j=1:length(row)
S4=[S4;[year doy freq S4_hr(row(j),1) S4_hr(row(j),2) S4_hr(row(j),3) S4_hr(row(j),4) S4_hr(row(j),5)]]
end
end
if (~isempty(SP_hr)&&~isempty(S4_hr))
for k=1:length(S4SPprnlist)
    S4row=find(S4prn==S4SPprnlist(k));
    SProw=find(SPprn==S4SPprnlist(k));
    for l=1:length(S4row)
           n=0
           S4startdatenum=datenum(year,month,day,S4_hr(S4row(l),2),S4_hr(S4row(l),3),0);
           S4enddatenum=datenum(year,month,day,S4_hr(S4row(l),4),S4_hr(S4row(l),5),0);
           
           for m=1:length(SProw)
               SPstartdatenum=datenum(year,month,day,SP_hr(SProw(m),2),SP_hr(SProw(m),3),0);
               SPenddatenum=datenum(year,month,day,SP_hr(SProw(m),4),SP_hr(SProw(m),5),0); 
               %Finding out disjoint time interval between SP and S4 for same prn 
               maximum_startdatenum=max(S4startdatenum,SPstartdatenum);
               minimum_enddatenum=min(S4enddatenum,SPenddatenum);
               if(maximum_startdatenum<minimum_enddatenum)
                   n=n+1;
               end
           end
               if (n==0)
                   S4=[S4;[year doy freq S4SPprnlist(k) S4_hr(S4row(l),2) S4_hr(S4row(l),3) S4_hr(S4row(l),4) S4_hr(S4row(l),5)]];
               end
               
           end
end
end
end

%______________________ only Sp ____________________________________________
SP=[];
if isempty(SP_hr)
    disp(['No Phase scintillation on the day',num2str(doy,'%03i'),' of the year',num2str(year,'%04i')]);
else

onlySPprnlist=setdiff(SPprn,S4prn);

for i=1:length(onlySPprnlist)
    row=find(SP_hr(:,1)==onlySPprnlist(i));
for j=1:length(row)
SP=[SP;[year doy freq SP_hr(row(j),1) SP_hr(row(j),2) SP_hr(row(j),3) SP_hr(row(j),4) SP_hr(row(j),5)]]
end
end
if (~isempty(SP_hr)&&~isempty(S4_hr))
for k=1:length(S4SPprnlist)
    S4row=find(S4prn==S4SPprnlist(k));
    SProw=find(SPprn==S4SPprnlist(k));
    for l=1:length(SProw)
           n=0
           SPstartdatenum=datenum(year,month,day,SP_hr(SProw(l),2),SP_hr(SProw(l),3),0);
           SPenddatenum=datenum(year,month,day,SP_hr(SProw(l),4),SP_hr(SProw(l),5),0);
           for m=1:length(S4row)
               S4startdatenum=datenum(year,month,day,S4_hr(S4row(m),2),S4_hr(S4row(m),3),0);
               S4enddatenum=datenum(year,month,day,S4_hr(S4row(m),4),S4_hr(S4row(m),5),0); 
               %Finding out disjoint time interval between SP and S4 for same prn 
               maximum_startdatenum=max(S4startdatenum,SPstartdatenum);
               minimum_enddatenum=min(S4enddatenum,SPenddatenum);
               if(maximum_startdatenum<minimum_enddatenum)
                   n=n+1;
               end
           end
           if (n==0)
                   SP=[SP;[year doy freq S4SPprnlist(k) SP_hr(SProw(l),2) SP_hr(SProw(l),3) SP_hr(SProw(l),4) SP_hr(SProw(l),5)]];
           end
    end
end
end


end

  combo1=S4SP;
  combo2=S4;
  combo3=SP;

  %--------------------------------------------------------------------------------------------------------------
  %Extracting the PFISR data for the given day
%   if ~exist(['Madrigal','_',num2str(day),'_',num2str(month),'_',num2str(year),'.mat'])
%   M=Madrigal_File_for_a_day(year,month,day);
%   end
%   
%   if ~exist(['Madrigalac','_',num2str(day),'_',num2str(month),'_',num2str(year),'.mat'])
%   Mac=Madrigal_File_for_a_day_ac(year,month,day);
%   end
  %-------------------------------------------------------------------------------------------------------------------
  
 %Identifying Ionospheric Region for only phase events
 SPregion=[];
 for i=1:size(SP)
     starttime=SP(i,5)+(SP(i,6)/60);
     endtime=SP(i,7)+(SP(i,8)/60);
    region = test_regiondetection_ac_lp(year,month,day,starttime,endtime);
    SPregion=[SPregion;region];
     
 end
 S4region=[];
 %Identifying Ionospheric Region for only amplitude events
 for i=1:size(S4)
     starttime=S4(i,5)+(S4(i,6)/60);
     endtime=S4(i,7)+(S4(i,8)/60);
     region = test_regiondetection_ac_lp(year,month,day,starttime,endtime);
     S4region=[S4region;region];
     
 end
 S4SPregion=[];
 %Identifying Ionospheric Region for both-amplitude and phase events
 for i=1:size(S4SP,1)
     starttime=S4SP(i,5)+(S4SP(i,6)/60);
     endtime=S4SP(i,7)+(S4SP(i,8)/60);
     region = test_regiondetection_ac_lp(year,month,day,starttime,endtime);
     S4SPregion=[S4SPregion;region];
     
 end
 combo1region=S4SPregion;
 combo2region=S4region;
 combo3region=SPregion;
end
  