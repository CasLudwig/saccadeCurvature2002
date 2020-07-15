function [metrics]=CurveMetrics(samples)

% FUNCTION WRITTEN BY CJH LUDWIG, JANUARY 2002
% COMMENTS AND QUESTIONS TO C.LUDWIG@BRISTOL.AC.UK
% This function has a matrix as input that contains four columns: blocknumber, trialnumber, x and y starting positions. 
% It goes through this matrix, normalising each movement and fitting a quadratic, and cubic polynomial on the samples.
% Several metrics of curvature are calculated: Initial deviation (Van Gisbergen et al. 1987), average initial deviation (Sheliga
% et al., 1995), maximum raw deviation (Smit et al., 1990; Doyle and Walker, 2001,2002), and an area-based measure (Doyle and Walker,
% pers. communication). In addition to these existing metrics, we calculate two metrics derived from the curve fitting procedure 

metrics=[];
i=1; %matrix row index
NRrows=length(samples);
while (i<=NRrows) %while not end of matrix
   x=[]; %initialize to-be-filled vectors
   y=[];
   xnorm=[];
   ynorm=[];
   n=i; %keeps track of the number of samples on a given trial
   blockNR=samples(i,1);
   trialNR=samples(i,2);
   while (samples(i,2)==trialNR) %fill x and y coordinate vectors
    x=[x; samples(i,3)];
   	y=[y; samples(i,4)];
    i=i+1;
    if(i>NRrows)
      break
    end
   end
   
   if(x(length(x))>x(1))
       direction=1; %rightward saccade
   else
       direction=0; %left saccade
   end;
      
   NRsamples=i-n; %number of samples on this trial

   hordisplacement=x(NRsamples)-x(1);
   vertdisplacement=y(1)-y(NRsamples);
   Hstraight=sqrt((hordisplacement^2)+(vertdisplacement^2));
   SacAngle=atan2(vertdisplacement,hordisplacement)*(180/pi); %calculate the angle of the entire movement
   
   %build up the normalised vectors
   xnorm=[xnorm; 0]; %each movement is normalised so that the starting position coincides with this origin (0,0)
   ynorm=[ynorm; 0];
   xres=[];
   
   for SampleIndex=2:(NRsamples-1) %first and last samples never deviate from the straight trajectory!
    hordisplacement= x(SampleIndex)-x(1);
    vertdisplacement= y(1)-y(SampleIndex);
    Hsample=sqrt((hordisplacement^2)+(vertdisplacement^2)); 
    SamAngle=atan2(vertdisplacement,hordisplacement)*180/pi;
    if(SacAngle>SamAngle)
        devdir=1; %clockwise deviation
        DevAngle=SacAngle-SamAngle;
    elseif(SacAngle<SamAngle)
        devdir=-1; %anti-clockwise deviation
        DevAngle=SamAngle-SacAngle;
    else
        dev=0; %no deviation
    end
    Deviation=sin(DevAngle*(pi/180))*Hsample;
    Deviation=Deviation*devdir;
    xtrue=sqrt((Hsample^2)-(Deviation^2)); %true x-coordinate along the straight line path
    xnorm=[xnorm; xtrue];
    ynorm=[ynorm; Deviation];
   end    
   xnorm=[xnorm; Hstraight];
   ynorm=[ynorm; 0];
   %rescale the x-coordinates so that xstart=-1 and xend=1
   for SampleIndex=1:NRsamples
       res=-1+((xnorm(SampleIndex)/xnorm(NRsamples))*2);
       xres=[xres;res];
   end     
   
   %now calculate the various established curvature metrics
   IniDev=atan2(ynorm(3),xnorm(3))*180/pi; %initial deviation metric of Van Gisbergen et al. (1987)
   IniAD=mean(ynorm(2:3)); %initial average deviation of Sheliga et al. (1995)
   
   [MaxDev,MaxIndex]=max(ynorm); %maximum raw deviation (Smit et al., 1990; Doyle and Walker, 2001,2002)
   [MinDev,MinIndex]=min(ynorm);
   if(abs(MaxDev)>abs(MinDev))
       RawDev=MaxDev;
       DevIndex=MaxIndex;
   elseif(abs(MaxDev)<abs(MinDev))
       RawDev=MinDev;
       DevIndex=MinIndex;
   else
       if(MaxIndex<MinIndex)
           RawDev=MaxDev;
           DevIndex=MaxIndex;
       else 
           RawDev=MinDev;
           DevIndex=MinIndex;
       end
   end
   RawDev=(RawDev/xnorm(NRsamples))*100;
   RawPOC=(xnorm(DevIndex)/xnorm(length(xnorm)))*100; %raw point of curvature
   
   AreaVector=[]; %area based measure (Doyle and Walker, personal communication)
   AreaVector=[AreaVector;0];
   for(AreaIndex=2:length(xnorm))
       area=(xnorm(AreaIndex)-xnorm(AreaIndex-1))*ynorm(AreaIndex);
       AreaVector=[AreaVector;area];
   end
   CurveArea=(sum(AreaVector)/xnorm(NRsamples))*100;
     
   %fit the quadratic function and determine the direction of curvature
   pol2=polyfit(xres, ynorm, 2);
   ypred2=polyval(pol2,xres);
   if(pol2(1)<0) %if quadratic coefficient is negative (upward curve), curvature is clockwise
       pol2(1)=abs(pol2(1));
   else
       pol2(1)=pol2(1)*-1; %if quadratic coefficient is positive (downward curve), curvature is anti-clockwise
   end
   
   pol3=polyfit(xres, ynorm, 3); %derivative of cubic polynomial
   ypred3=polyval(pol3,xres);  
   vertdisplacement=ypred3(1)-ypred3(length(ypred3));
   Hstraight=sqrt((xnorm(length(xnorm))^2)+(vertdisplacement^2));
   SacAngle=atan2(vertdisplacement,xnorm(length(xnorm)))*(180/pi);   
   der3=polyder(pol3); %derivative of cubic polynomial gives a maximum and minimum
   xder3=((-1*der3(2))-sqrt((der3(2)^2)-(4*der3(1)*der3(3))))/(2*der3(1));
   xder3=[xder3; ((-1*der3(2))+sqrt((der3(2)^2)-(4*der3(1)*der3(3))))/(2*der3(1))];
   if ((xder3(1)<xres(1)) | (xder3(1)>xres(length(xres)))) %check whether first maximum/minimum falls within the range of xres
       curve3=0;
       POC3=0;
   else %if yes, then calculate curvature
       ymax3=polyval(pol3,xder3(1));
       POC3=(xder3(1)*std(xnorm))+mean(xnorm);
       POC3=(POC3/xnorm(length(xnorm)))*100;
       hordisplacement=(xder3(1)*std(xnorm))+mean(xnorm);
       vertdisplacement= ypred3(1)-ymax3;
       Hsample=sqrt((hordisplacement^2)+(vertdisplacement^2));
       SamAngle=atan2(vertdisplacement,hordisplacement)*180/pi;
       if(SacAngle>SamAngle)
            devdir=1; %clockwise deviation
            DevAngle=SacAngle-SamAngle;
       elseif(SacAngle<SamAngle)
            devdir=-1; %anti-clockwise deviation
            DevAngle=SamAngle-SacAngle;
       end
       curve3=sin(DevAngle*(pi/180))*Hsample;
       curve3=((curve3*devdir)/xnorm(length(xnorm)))*100;      
   end
   if ((xder3(2)<xres(1)) | (xder3(2)>xres(length(xres)))) %check whether second maximum/minimum falls within the range of xres
       curve3=[curve3;0];
       POC3=[POC3;0];
   else %if yes, then calculate curvature
       ymax3=polyval(pol3,xder3(2));
       POC=(xder3(2)*std(xnorm))+mean(xnorm);
       POC3=[POC3;(POC/xnorm(length(xnorm)))*100];
       hordisplacement=(xder3(2)*std(xnorm))+mean(xnorm);
       vertdisplacement= ypred3(1)-ymax3;
       Hsample=sqrt((hordisplacement^2)+(vertdisplacement^2));
       SamAngle=atan2(vertdisplacement,hordisplacement)*180/pi;
       if(SacAngle>SamAngle)
            devdir=1; %clockwise deviation
            DevAngle=SacAngle-SamAngle;
       elseif(SacAngle<SamAngle)
            devdir=-1; %anti-clockwise deviation
            DevAngle=SamAngle-SacAngle;
       end
       curve=sin(DevAngle*(pi/180))*Hsample;
       curve3=[curve3;((curve*devdir)/xnorm(length(xnorm)))*100];      
   end
   if max(abs(curve3))>0
       [MaxDev,MaxIndex]=max(abs(curve3));
   else MaxIndex=1;
   end
   metrics=[metrics; blockNR trialNR direction IniDev IniAD RawDev RawPOC CurveArea pol2(1) curve3(1) POC3(1) curve3(2) POC3(2) curve3(MaxIndex) POC3(MaxIndex)];   
end
save('curvature metrics', 'metrics', '-ascii');