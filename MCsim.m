function [DistParam, Meds]=MCsim(samples)

% FUNCTION WRITTEN BY CJH LUDWIG, JANUARY 2002
% COMMENTS AND QUESTIONS TO C.LUDWIG@BRISTOL.AC.UK
% This function has a matrix as input that contains four columns: blocknumber, trialnumber, x and y starting positions. 
% It goes through this matrix, normalising each movement and fitting various functions of increasing complexity. To assess
% to what extent the improvement in the goodness-of-fit of these functions is due to simply accounting for more noise as 
% opposed to accounting for a systematic component in the trajectory, the residuals from each fit are randomly mixed up and
% added to the predicted values. The next function is then fit to this "noise" vector. Once each movement has undergone this analysis,
% a median r-squared is calculated for each of the fits. This process is then repeated a number of times, to yield four distributions
% of median r-squared values. The mean and standard deviation of these distributions, serve to evaluate the obtained r-squared values
% for the fits on the observed data.

Meds=[];
DistParam=[];
for MCindex=1:5000 %set the repeats of the entire process
    Distributions=[];
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

   %working with the normalised vectors fit the various functions (as long as number of data points is larger than the 
   %number of parameters estimated by the 4th order polynomial
        if(NRsamples>=4)
        %linear fit
            pol1=polyfit(xres, ynorm, 1);
            ypred1=polyval(pol1,xres); %y-coordinates if x is put into the polynomial
            Rsq1=corrcoef([ynorm ypred1]);
            Rsq1=Rsq1(2,1)^2;
            yrand=ynorm;
            error1=ynorm-ypred1; %error of the fit: difference between observed and fitted y-coordinates
            randvector=randperm(NRsamples);
            for yindex=1:(NRsamples)
                yrand(yindex)=ypred1(yindex)+error1(randvector(yindex));%the residuals are randomly mixed up and then added to the
                                                                          %predicted values to create "noise" vector
            end
   
        %2nd order polynomial
            pol2=polyfit(xres, yrand, 2); %fit the 2nd order function on the "noise" vector
            ypred2=polyval(pol2,xres); 
            Rsq2=corrcoef([yrand ypred2]);
            Rsq2=Rsq2(2,1)^2; %and calculate how well this vector is predicted
            pol2=polyfit(xres, ynorm, 2); %fit the function on the observed data so that a new "noise" vector can be calculated
            ypred2=polyval(pol2,xres); 
            yrand=ynorm;
            error2=ynorm-ypred2; 
            randvector=randperm(NRsamples);
            for yindex=1:(NRsamples)
                yrand(yindex)=ypred2(yindex)+error2(randvector(yindex));
            end       

        %3rd order polynomial
            pol3=polyfit(xres, yrand, 3);
            ypred3=polyval(pol3,xres); 
            Rsq3=corrcoef([yrand ypred3]);
            Rsq3=Rsq3(2,1)^2;
            pol3=polyfit(xres, ynorm, 3)
            ypred3=polyval(pol3,xres)
            yrand=ynorm;
            error3=ynorm-ypred3 
            randvector=randperm(NRsamples);
            for yindex=1:(NRsamples)
                yrand(yindex)=ypred3(yindex)+error3(randvector(yindex));
            end  
            yrand

        %4th order polynomial
            pol4=polyfit(xres, yrand, 4);
            ypred4=polyval(pol4,xres); 
            Rsq4=corrcoef([yrand ypred4]);
            Rsq4=Rsq4(2,1)^2;
        end %if(NRsamples>=4) 
        Distributions=[Distributions; Rsq1 Rsq2 Rsq3 Rsq4]; %Add the four R-squared values to the matrix
    end %while(i<=NRrows)
    med1=median(Distributions(:,1));
    med2=median(Distributions(:,2));
    med3=median(Distributions(:,3));
    med4=median(Distributions(:,4));
    Meds=[Meds; med1 med2 med3 med4];
end %for(MCindex=1:2)
save('medians MC simulation', 'Meds', '-ascii');   
DistParam=[DistParam; mean(Meds(:,1)) mean(Meds(:,2)) mean(Meds(:,3)) mean(Meds(:,4)); std(Meds(:,1)) std(Meds(:,2)) std(Meds(:,3)) std(Meds(:,4))];