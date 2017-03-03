clear
clc
close all
% This Code is only a very basic outline with a few needed functions to
% help guide you in your image analysis.  It will not function in its current
% form.  There are lots of things missing but they should all be achievable
% in the alotted time.  Good Luck - AJ

%% Add path to functions for calculating laminar and turbulent flame speeds.
% NOTE: if you bwfilter.m is not in the same directory that you run this in
% you will need to add a path for it.
Datapath = 'C:\Users\cowand\ME 552 Data';                           %Adds Data path as callable string
bg_path = {'RE_5000_Juan.tiff','RE_10000_Juan.tiff'};                         %Background File Path
Cal_path = 'Spacial_Calibration_Group_Juan.tiff';                       %Calibration File Path

numsamples = 1;                                                                 %Number of sample's
numframes = 360*numsamples;                                                     %Set number of frames from camera settings
home = pwd;                                                                     %Sets variable "home" to current directory Callable

fg_path = {'A2_5k_HTI_01.tiff','A2_5k_HTI_02.tiff','A2_5k_HTI_03.tiff',...
    'A2_5k_LTI_01.tiff','A2_5k_LTI_02.tiff','A2_5k_LTI_03.tiff',...
    'A2_10k_HTI_01.tiff','A2_10k_HTI_02.tiff','A2_10k_HTI_03.tiff',...
    'A2_10k_LTI_01.tiff','A2_10k_LTI_02.tiff','A2_10k_LTI_03.tiff'};          %Forground File Path - This is your image of interest
%fg_path = 'C:\Users\cowand\ME 552 Data\A2_10k_LTI_03.tiff';
%fg_path = 'C:\Users\cowand\ME 552 Data\A2_10k_HTI_02.tiff';

%% Input Parameters
m_dot_total =0.001215;               %Total mass flow (kg/s) of reactants for Re=5,000
%m_dot_total=.00243;                 %Total mass flow (kg/s) of reactants for Re=10,000
rho_u       = 0.7871;                     %Unburned Density (kg/m^3) at 200 C
D           = 12/10;                       %Burner Diameter (cm)

%%  Spacial Calibration

%This will be used for image croping and spatial calibration
info_cal= imfinfo(Cal_path);

im_cal=imread(Cal_path,1,'Info',info_cal,'PixelRegion',{[775,921],[300,900]});
figure(15)
imshow(im_cal,[50,230]); axis image %This will show the image to calibrate to
title('Select locations 70cm and 76cm') %Please select locations 70cm and 76cm
points=ginput(2);
close
point_dist=pdist(points,'euclidean');
cal_val=6/point_dist; %Image calibration  (cm/pixel)

%% Image Cropping

center_l=514;                       %axis of symetry
Lcrop= 440;                       % left edge of the Burner
Rcrop= 1024+Lcrop-2*center_l;                       % right edge of the Burner
Tcrop= 400;                       % location for cropping of the top of the image
Bcrop= 1024-921;                       % location for cropping of the bottom of the image

%% Background images

for f=1:length(bg_path) %finding backgound info for both Re
    info_b=imfinfo(bg_path{f}); %info of image
    for i=1:numframes %calculating the background average
        if i==1
            im_avg_b=double(imread(bg_path{f},1,'Info',info_b,'PixelRegion',{[Tcrop,1024-Bcrop],[Lcrop,1024-Rcrop]}));
        else
            image_i=double(imread(bg_path{f},i,'Info',info_b,'PixelRegion',{[Tcrop,1024-Bcrop],[Lcrop,1024-Rcrop]}));
            im_avg_b=((i-1).*im_avg_b+image_i)./i;
        end
    end
    b_avg(:,:,f)=im_avg_b; %This has both Re values in it.
end

%% Building Meshgrid 

siz=size(im_avg_b); %Used to greate mesh grid
x=0:(siz(2)-1); 
x=cal_val*x; %Calibrated x distance
y=0:(siz(1)-1); 
y=cal_val*y; %Calibrated y distance
y=fliplr(y);

[X,Y]=meshgrid(x,y);
%Y=flipud(Y);

%% Calculating averages and plotting images

count=1;
for f=1:length(fg_path) %Plotting all files
    info_f=imfinfo(fg_path{f}); %info of image


    for i=1:numframes %calculating the image average
        if i==1
            image_0=double(imread(fg_path{f},1,'Info',info_f,'PixelRegion',{[Tcrop,1024-Bcrop],[Lcrop,1024-Rcrop]}));
            im_avg_f=image_0;
        else
            image_i=double(imread(fg_path{f},i,'Info',info_f,'PixelRegion',{[Tcrop,1024-Bcrop],[Lcrop,1024-Rcrop]}));
            im_avg_f=((i-1).*im_avg_f+image_i)./i;
        end
    end
    
    f_avg(:,:,f)=im_avg_f; %Image average per test. This has both Re values in it.
    max_int=90;
    
    if f<6.5 %The following is calculating the flame minus the background
        flame_avg(:,:,f)=f_avg(:,:,f)-b_avg(:,:,1);
    else
        flame_avg(:,:,f)=f_avg(:,:,f)-b_avg(:,:,2);
    end
    
%    max_int=max(max(flame_avg(:,:,f))); %finding maxinum intensity for per test
    flame_filt(:,:,f)=medfilt2(flame_avg(:,:,f)); %filters
    flame(:,:,f)=0.5*(flame_avg(:,1:end,f) + flame_filt(:, end+1-(1:end), f));   % Axis symetric

    
    if rem(f,3)==0 %This is to find the average with all three tests
        avg_flame(:,:,count)=mean(flame(:,:,f-2:f),3);
        if f<3.5 %Plotting average to compare to snap-shots
            figure(1)
            subplot(1,7,1)
            colormap jet
            surf(X,Y,avg_flame(:,:,count),'EdgeColor','None'); axis image
            caxis([0,65])
            view(2)
            xlim([0,max(x)])
            ylim([0,max(y)])
            title('Axis Symetric Corrected Flame (Re=5,000, High Turbulence Intensity)','HorizontalAlignment', 'left')
            xlabel('X-Distance (cm)')
            ylabel('Y-Distance (cm)')
        elseif f<6.5
            figure(2)
            subplot(1,7,1)
            colormap jet
            surf(X,Y,avg_flame(:,:,count),'EdgeColor','None'); axis image
            caxis([0,65])
            view(2)
            grid off
            xlim([0,max(x)])
            ylim([0,max(y)])
            title('Axis Symetric Corrected Flame (Re=5,000, Low Turbulence Intensity)','HorizontalAlignment', 'left')
            xlabel('X-Distance (cm)')
            ylabel('Y-Distance (cm)')
        elseif f<9.5
            figure(3)
            subplot(1,7,1)
            colormap jet
            surf(X,Y,avg_flame(:,:,count),'EdgeColor','None'); axis image
            caxis([0,65])
            view(2)
            grid off
            xlim([0,max(x)])
            ylim([0,max(y)])
            title('Axis Symetric Corrected Flame (Re=10,000, High Turbulence Intensity)','HorizontalAlignment', 'left')
            xlabel('X-Distance (cm)')
            ylabel('Y-Distance (cm)')
        else
            figure(4)
            subplot(1,7,1)
            colormap jet
            surf(X,Y,avg_flame(:,:,count),'EdgeColor','None'); axis image
            caxis([0,65])
            view(2)
            grid off
            xlim([0,max(x)])
            ylim([0,max(y)])
            title('Axis Symetric Corrected Flame (Re=10,000, Low Turbulence Intensity)','HorizontalAlignment', 'left')
            xlabel('X-Distance (cm)')
            ylabel('Y-Distance (cm)')
        end
        count=count+1;
    end
    
    title_loc=4; %Location of center of title over which sub plot
    
    if f<3.5 %The following if statement plots the 6 snap-shots
        figure(1)
        a=f;
        subplot(1,7,a+1)
        colormap jet
        surf(X,Y,image_0-b_avg(:,:,1),'EdgeColor','None'); axis image %one snap shot per test
        caxis([0,65])
        view(2)
        xlim([0,max(x)])
        ylim([0,max(y)])
        set(gca,'ytick',[])
        
        subplot(1,7,a+4)
        colormap jet
        surf(X,Y,image_i-b_avg(:,:,1),'EdgeColor','None'); axis image %The other snap shot per test
        caxis([0,65])
        view(2)
        xlim([0,max(x)])
        ylim([0,max(y)])
        set(gca,'ytick',[])
        
    elseif f<6.5
        figure(2)
        a=f-3;
        subplot(1,7,a+1)
        colormap jet
        surf(X,Y,image_0-b_avg(:,:,1),'EdgeColor','None'); axis image %one snap shot per test
        caxis([0,65])
        view(2)
        xlim([0,max(x)])
        ylim([0,max(y)])
        set(gca,'ytick',[])
        
        subplot(1,7,a+4)
        surf(X,Y,image_i-b_avg(:,:,1),'EdgeColor','None'); axis image %The other snap shot per test
        caxis([0,65])
        view(2)
        xlim([0,max(x)])
        ylim([0,max(y)])
        set(gca,'ytick',[])
        
    elseif f<9.5
        figure(3)
        a=f-6;
        subplot(1,7,a+1)
        colormap jet
        surf(X,Y,image_0-b_avg(:,:,2),'EdgeColor','None'); axis image %one snap shot per test
        caxis([0,65])
        view(2)
        xlim([0,max(x)])
        ylim([0,max(y)])
        set(gca,'ytick',[])
        
        subplot(1,7,a+4)
        colormap jet
        surf(X,Y,image_i-b_avg(:,:,2),'EdgeColor','None'); axis image %The other snap shot per test
        caxis([0,65])
        view(2)
        xlim([0,max(x)])
        ylim([0,max(y)])
        set(gca,'ytick',[])
        
    else
        figure(4)
        a=f-9;
        subplot(1,7,a+1)
        colormap jet
        surf(X,Y,image_0-b_avg(:,:,2),'EdgeColor','None'); axis image %one snap shot per test
        caxis([0,65])
        view(2)
        xlim([0,max(x)])
        ylim([0,max(y)])
        set(gca,'ytick',[])
        
        subplot(1,7,a+4)
        colormap jet
        surf(X,Y,image_i-b_avg(:,:,2),'EdgeColor','None'); axis image %The other snap shot per test
        caxis([0,65])
        view(2)
        xlim([0,max(x)])
        ylim([0,max(y)])
        set(gca,'ytick',[])
        
    end
    
    if f<3.5 %This plots average per test. Set to comment because not needed
%         figure(1)
%         a=f;
%         subplot(1,3,a)
%         surf(X,Y,flame(:,:,f),'EdgeColor','None')
%         view(2)
%         grid off
%         xlim([0,max(x)])
%         ylim([0,max(y)])
%         %imshow(flame(:,:,f),[0,max_int])
%         if a==title_loc 
%             title('Axis Symetric Corrected Flame (Re=5,000, High Turbulence Intensity)')
%         end
%         %legend('L/e=2e3','L/e=5e3','L/e=1e4','L/e=5e4')
%         xlabel('X-Distance (cm)')
%         ylabel('Y-Distance (cm)')
%     elseif f<6.5
%         figure(2)
%         a=f-3;
%         subplot(1,3,a)
%         surf(X,Y,flame(:,:,f),'EdgeColor','None')
%         view(2)
%         grid off
%         xlim([0,max(x)])
%         ylim([0,max(y)])
%         %imshow(flame(:,:,f),[0,max_int])
%         if a==title_loc
%             title('Axis Symetric Corrected Flame (Re=5,000, Low Turbulence Intensity)')
%         end
%         %legend('L/e=2e3','L/e=5e3','L/e=1e4','L/e=5e4')
%         xlabel('X-Distance (cm)')
%         ylabel('Y-Distance (cm)')
%     elseif f<9.5
%         a=f-6;
%         subplot(1,3,a)
%         surf(X,Y,flame(:,:,f),'EdgeColor','None')
%         view(2)
%         grid off
%         xlim([0,max(x)])
%         ylim([0,max(y)])
%         %imshow(flame(:,:,f),[0,max_int])
%         if a==title_loc
%             title('Axis Symetric Corrected Flame (Re=10,000, High Turbulence Intensity)')
%         end
%         %legend('L/e=2e3','L/e=5e3','L/e=1e4','L/e=5e4')
%         xlabel('X-Distance (cm)')
%         ylabel('Y-Distance (cm)')
%     else
%         figure(4)
%         a=f-9;
%         subplot(1,3,a)
%         surf(X,Y,flame(:,:,f),'EdgeColor','None')
%         view(2)
%         grid off
%         xlim([0,max(x)])
%         ylim([0,max(y)])
%         %imshow(flame(:,:,f),[0,max_int])
%         if a==title_loc
%             title('Axis Symetric Corrected Flame (Re=10,000, Low Turbulence Intensity)')
%         end
%         %legend('L/e=2e3','L/e=5e3','L/e=1e4','L/e=5e4')
%         xlabel('X-Distance (cm)')
%         ylabel('Y-Distance (cm)')
    end

end

%% Finding flame area
center=(siz(2)+1)/2; %center index
%figure(11)
%plot(y,movmean(avg_flame(:,center,2),20))
%plot(h,max1h)

count=zeros(1,4);
for f=1:4
    for i=1:siz(2)
        for j=1:siz(1)
            if avg_flame(j,i,f)>24
                count(f)=count(f)+1;
            end
        end
    end
end
appearent_flame_A=cal_val^2*count;


for i=1:4 %This is old incorrect plotting
    [val,indx]=max(movmean(avg_flame(:,center,i),20)); %Locates max temp along center
    height(i)=y(indx); %Finds height of cone
    
    cone_A(i)=(D/2)*pi()*(height(i)^2+(D/2)^2)^(1/2); %cone area

end


%% Additional Plots

plt_titles={'Re=5k HTI','Re=5k LTI','Re=10k HTI','Re=10k LTI'};
%method_plts=zeros(5,siz(1),siz(2),4);
for f=1:4 %This shows all the flames under all four conditions
    figure(5)
    subplot(1,4,f)
    surf(X,Y,avg_flame(:,:,f),'EdgeColor','None'); axis image %one snap shot per test
    colormap jet
    caxis([0,65])
    view(2)
    xlim([0,max(x)])
    ylim([0,max(y)])
    title(plt_titles{f})
    if f==1
        ylabel('Y-Distance (cm)')
    elseif f==2
        xlabel('X-Distance (cm)')
        set(gca,'ytick',[])
    else
        set(gca,'ytick',[])
    end
    
    a=f*3;
    b=ceil(f/2); %The method_plts will show the method
    method_plts(:,:,:,f)=cat(3,b_avg(:,:,b),f_avg(:,:,a),flame_avg(:,:,a),flame_filt(:,:,a),flame(:,:,a));
    
end

plt_titles={'Re=5,000, High Turbulence Intensity',...
    'Re=5,000, Low Turbulence Intensity',...
    'Re=10,000, High Turbulence Intensity',...
    'Re=10,000, Low Turbulence Intensity'};
plt_xlab={'Background','Average Flame','Background Subtracted','Filtered','Axis Symetric'};
for f=1:4
    figure()
    for i=1:5
        subplot(1,5,i)
        colormap jet
        surf(X,Y,method_plts(:,:,i,f),'EdgeColor','None'); axis image %one snap shot per test
        caxis([0,65])
        view(2)
        xlim([0,max(x)])
        ylim([0,max(y)])
        if i==1
            title(plt_titles{f},'HorizontalAlignment', 'left')
            ylabel('Y-Distance (cm)')
            %axes('Position', [0.05 0.05 0.9 0.9], 'Visible', 'off');
            %colorbar
        else
            set(gca,'ytick',[])
        end
        xlabel(plt_xlab{i});
        set(get(gca,'XLabel'),'Rotation',20,'HorizontalAlignment', 'right')
    end
end
        

%% Color Bar Graph
x_cb=0:65;
y_cb=0:5;
[X_cb,Y_cb]=meshgrid(x_cb,y_cb);
figure()
colormap jet
surf(X_cb,Y_cb,X_cb,'EdgeColor','None'); axis image %one snap shot per test
caxis([0,65])
view(2)
set(gca,'ytick',[])
xlabel('Intensity')




%% Was used to find Flame center

for x=1:6 %This was used to find the crop values for centering image
%     max1h=max(flame_avg(:,:,x));
%     max2h=max(flame_avg(:,:,x+6));
% 
%     max1v=max(transpose(flame_avg(:,:,x)));
%     max2v=max(transpose(flame_avg(:,:,x+6)));
%     h=1:length(max1h);
%     v=1:length(max1v);
%     figure(11)
%     subplot(2,1,1)
%     plot(v,max1v)
%     %plot(h,max1h)
%     hold on
%     subplot(2,1,2)
%     plot(v,max2v)
%     %plot(h,max2h)
%     hold on
end


