function [coef]=ComputeCurve(samples)

% FUNCTION WRITTEN BY CJH LUDWIG, JULY 2002 COMMENTS AND QUESTIONS TO
% C.LUDWIG@BRISTOL.AC.UK This function has a matrix as input that contains
% four columns: blocknumber, trialnumber, x and y positions. It
% goes through this matrix, rotates the saccade so that the straight line
% through start and endpoint coincides with the abcissa. The displacement
% along the horizontal is then rescaled so that the x-coordinates of start
% and endpoint are -1 and 1. A quadratic polynomial is fit on the
% samples. The function outputs the polynomial coefficients, and the mean
% of the predicted y-coordinates of the start and endpoint. We recommend
% using the quadratic coefficient as a direct measure of curvature (in
% pixels). This is certainly the case when the means of the predicted start
% and end y-coordinates average out to zero.

coef=[];
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
      
   NRsamples=i-n; %number of samples on this trial

   hordisplacement=x(NRsamples)-x(1);
   vertdisplacement=y(1)-y(NRsamples);
   Hstraight=sqrt((hordisplacement^2)+(vertdisplacement^2));
   SacAngle=atan2(vertdisplacement,hordisplacement)*(180/pi); %calculate the angle of the entire movement
   
   %build up the rotated vectors so that each saccade starts in (0,0)
   xnorm=[xnorm; 0]; 
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
    xtrue=sqrt((Hsample^2)-(Deviation^2));%true x-coordinate along the straight line path
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

   %fit the quadratic function and determine the direction of curvature
   pol2=polyfit(xres, ynorm, 2);
   ypred2=polyval(pol2,xres);
   rsq=corrcoef(ynorm,ypred2);
   rsq=rsq(1,2)^2;
   if(pol2(1)<0) %if quadratic coefficient is negative (upward curve), curvature is clockwise
       pol2(1)=abs(pol2(1));
   else
       pol2(1)=pol2(1)*-1; %if quadratic coefficient is positive (downward curve), curvature is anti-clockwise
   end

   coef=[coef; blockNR trialNR pol2(1) pol2(2) pol2(3) (ypred2(1)+ypred2(length(ypred2)))/2 rsq];   
end