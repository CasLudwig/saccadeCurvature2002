function [Fit]=normcurve(samples)

% FUNCTION WRITTEN BY CJH LUDWIG, DECEMBER 2001
% COMMENTS AND QUESTIONS TO C.LUDWIG@BRISTOL.AC.UK
% This function has a matrix as input that contains four columns: blocknumber, trialnumber, x and y starting positions. 
% It goes through this matrix, normalising each movement and fitting a linear function, a quadratic, a cubic, and a quadrupic 
% polynomial on the samples of each trial. The squared correlation between the normalised and fitted y-coordinates,
% and the mean squareroot error, are calculated to give indications of the goodness-of-fit. The function outputs a matrix with
% 'blocknumber', 'trialnumber', and the relevant data associated with the four fitted functions (coefficients and 
% goodness-of-fit statistic).

Fit=[];
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
    xtrue=sqrt((Hsample^2)-(Deviation^2));
    
    xnorm=[xnorm; xtrue];
    ynorm=[ynorm; Deviation];
   end    
   xnorm=[xnorm; Hstraight]
   ynorm=[ynorm; 0]
   %rescale the x-coordinates so that xstart=-1 and xend=1
   for SampleIndex=1:NRsamples
       res=-1+((xnorm(SampleIndex)/xnorm(NRsamples))*2);
       xres=[xres;res];
   end     
   
   %working with the normalised vectors fit the various functions (as long as number of data points is larger than the 
   %number of parameters estimated by the 4th order polynomial

      if(NRsamples>=4)
       %linear fit
       pol1=polyfit(xres, ynorm, 1);
       ypred1=polyval(pol1,xres); %y-coordinates if x is put into the polynomial
       error=ynorm-ypred1; %error of the fit: difference between observed and fitted y-coordinates
       Rsq1=corrcoef([ynorm ypred1]);
       Rsq1=Rsq1(2,1)^2;
   
       %2nd order polynomial
       pol2=polyfit(xres, ynorm, 2);
       ypred2=polyval(pol2,xres); %y-coordinates if x is put into the polynomial
       error=ynorm-ypred2; %error of the fit: difference between observed and fitted y-coordinates
       Rsq2=corrcoef([ynorm ypred2]);
       Rsq2=Rsq2(2,1)^2;

       %3rd order polynomial
       pol3=polyfit(xres, ynorm, 3);
       ypred3=polyval(pol3,xres); %y-coordinates if x is put into the polynomial
       error=ynorm-ypred3; %error of the fit: difference between observed and fitted y-coordinates
       Rsq3=corrcoef([ynorm ypred3]);
       Rsq3=Rsq3(2,1)^2;

       %4th order polynomial
       pol4=polyfit(xres, ynorm, 4);
       ypred4=polyval(pol4,xres); %y-coordinates if x is put into the polynomial
       Rsq4=corrcoef([ynorm ypred4]);
       Rsq4=Rsq4(2,1)^2;

       %write all the parameter and goodness-of-fit estimates to a matrix
       Fit=[Fit; blockNR trialNR pol1(1) pol1(2) Rsq1 pol2(1) pol2(2) pol2(3) Rsq2 pol3(1) pol3(2) pol3(3) pol3(4) Rsq3 pol4(1) pol4(2) pol4(3) pol4(4) pol4(5) Rsq4];
   end
end
save('data', 'Fit', '-ascii');   
xaxis=0:0.1:1;
Rsq1=Fit(:,5);
Rsq2=Fit(:,9);
Rsq3=Fit(:,14);
Rsq4=Fit(:,20);
subplot(1,4,1);hist(Rsq1,xaxis');
subplot(1,4,2);hist(Rsq2,xaxis');
subplot(1,4,3);hist(Rsq3,xaxis');
subplot(1,4,4);hist(Rsq4,xaxis');
