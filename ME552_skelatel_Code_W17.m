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
D           = 12;                       %Burner Diameter

%% Image Parameters - Selected From Image
%Use imread, imshow, and ginput functions to complete calibration
%You will need to calibrate to the burner width, the outlet plane of the
%burner, the top of the image, and the number of pixels per unit length 
%(mm or in)
%This will be used for image croping and spatial calibration
info_cal= imfinfo(Cal_path);

center_l=514;                       %axis of symetry
Lcrop= 349;                       % left edge of the Burner
Rcrop= 1024+Lcrop-2*center_l;                       % right edge of the Burner
Tcrop= 200;                       % location for cropping of the top of the image
%Tcrop_LRe= 100;                       % location for cropping of the top of the image
%Tcrop_SRe= 500;                       % location for cropping of the top of the image
Bcrop= 1024-900;                       % location for cropping of the bottom of the image



for f=1:length(bg_path) %finding backgound info for both Re
    info_b=imfinfo(bg_path{f}); %info of image
%     if f==1 %This creates the two different crops for both Re values
%         Tcrop=Tcrop_SRe;
%     else
%         Tcrop=Tcrop_LRe;
%     end
    for i=1:1 %calculating the background average
        if i==1
            im_avg_b=double(imread(bg_path{f},1,'Info',info_b,'PixelRegion',{[Tcrop,1024-Bcrop],[Lcrop,1024-Rcrop]}));
        else
            image_i=double(imread(bg_path{f},i,'Info',info_b,'PixelRegion',{[Tcrop,1024-Bcrop],[Lcrop,1024-Rcrop]}));
            im_avg_b=((i-1).*im_avg_b+image_i)./i;
        end
    end
    b_avg(:,:,f)=im_avg_b; %This has both Re values in it.
end

figure()
for f=1:length(fg_path)
    info_f=imfinfo(fg_path{f}); %info of image
%     if f<6.5 %This creates the two different crops for both Re values
%         Tcrop=Tcrop_SRe;
%     else
%         Tcrop=Tcrop_LRe;
%     end

    for i=1:1 %calculating the background average
        if i==1
            im_avg_f=double(imread(fg_path{f},1,'Info',info_f,'PixelRegion',{[Tcrop,1024-Bcrop],[Lcrop,1024-Rcrop]}));
        else
            image_i=double(imread(fg_path{f},i,'Info',info_f,'PixelRegion',{[Tcrop,1024-Bcrop],[Lcrop,1024-Rcrop]}));
            im_avg_f=((i-1).*im_avg_f+image_i)./i;
        end
    end
    f_avg(:,:,f)=im_avg_f; %This has both Re values in it.
    
    if f<6.5 %The following is calculating the flame minus the background
        flame_avg(:,:,f)=f_avg(:,:,f)-b_avg(:,:,1);
    else
        flame_avg(:,:,f)=f_avg(:,:,f)-b_avg(:,:,2);
    end
    max_int=max(max(flame_avg(:,:,f))); %finding maxinum intensity for plot reference
    axis_sym=0.5*(flame_avg(:,1:end,f) + flame_avg(:, end+1-(1:end), f));   % Axis symetric
    flame(:,:,f)=medfilt2(uint8(axis_sym)); %filters and turns back into image
    
    if f<3.5
        figure(1)
        subplot(1,3,f)
        imshow(flame(:,:,f),[0,max_int])
        if f==2
            title('Axis Symetric Corrected Flame (Re=5,000, High Turbulence Intensity)')
        end
        %legend('L/e=2e3','L/e=5e3','L/e=1e4','L/e=5e4')
        xlabel('X-Distance (cm)')
        ylabel('Y-Distance (cm)')
    elseif f<6.5
        figure(2)
        subplot(1,3,f-3)
        imshow(flame(:,:,f),[0,max_int])
        if f==5
            title('Axis Symetric Corrected Flame (Re=5,000, Low Turbulence Intensity)')
        end
        %legend('L/e=2e3','L/e=5e3','L/e=1e4','L/e=5e4')
        xlabel('X-Distance (cm)')
        ylabel('Y-Distance (cm)')
    elseif f<9.5
        figure(3)
        subplot(1,3,f-6)
        imshow(flame(:,:,f),[0,max_int])
        if f==8
            title('Axis Symetric Corrected Flame (Re=10,000, High Turbulence Intensity)')
        end
        %legend('L/e=2e3','L/e=5e3','L/e=1e4','L/e=5e4')
        xlabel('X-Distance (cm)')
        ylabel('Y-Distance (cm)')
    else
        figure(4)
        subplot(1,3,f-9)
        imshow(flame(:,:,f),[0,max_int])
        if f==11
            title('Axis Symetric Corrected Flame (Re=10,000, Low Turbulence Intensity)')
        end
        %legend('L/e=2e3','L/e=5e3','L/e=1e4','L/e=5e4')
        xlabel('X-Distance (cm)')
        ylabel('Y-Distance (cm)')
    end

end

% for x=1:6 %This was used to find the crop values for centering image
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
% end

% flame_avg=im_avg_f-im_avg_b;
% max_f=max(max(flame_avg(1:6)));
% flame=uint8(flame_avg);
% figure, imshow(flame,[0,max_f])

% 
% Left                    = ;                                                     % location of left line on calibration device
% Right                   = ;                                                     % location of right line on calibration device
% 
% %% Image Parameters - User Selected
% %calibration determination
% Cal_Length              = ;                                                     % Calibration Unit Length
% pixels                  = Right - Left;                                         % # of Pixels
% Calibration             = Cal_Length/pixels;                                    % meters/pixel
% 
% %% Image Read In, Average, and Background Subtract
% %read on .tiff, average images in multipage format for backgrounds
% %subtract background form images in foreground

% 
% %% Axisymmetrize the image by averaging both sides
% Axisym = 0.5*(cropped_Flame(:, 1:end) + cropped_Flame(:, end + 1 - (1:end)));   % Use this equation to correct for asymetry
% [f, g] = size(Axisym(:, :));                                                    % get dimensions of the axisymmetrized image
% 
% %% Perform 2D median filter on images
% ImageFilt(:, :) = ;                                                             %
% 
% %% This is were I leave you
% % you will need to develope a method to reduce your data further to
% % calcutate the stated quantities in the lab procedure.
% % One last note,
% % It's dangerous to go alone, take this https://www.mathworks.com/help/index.html


