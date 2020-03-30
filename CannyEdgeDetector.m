function CannyEdgeDetector()
    close all;  % Close figures
    sigma = 1.4; % Gaussian filter sigma
    highThresholdRatio = 0.275; % High threshold ratio
    lowThresholdRatio = 0.25; % Low threshold ratio OF THE high threshold
    
    if(~isdeployed)
      cd(fileparts(which(mfilename)));
    end
    
    im = imread('blocks.tif');  %Buraya resimin adres yolunu yaz
    figure;
    imshow(im);
    title('Original Image');
    %% 1. Step - Smooth with Gaussian 5x5 filter to reduce noise%%
   
    %Ýlk olarak Low pass Filter uyguluyoruz.
    %Resimlerdeki gürütüleri azaltýp cisim kenalarýný daha temiz bir
    %þekilde bulmamýza yarayacak.
    
    imG = im;
    im = double(imgaussfilt(im,sigma));
    imGaus = imgaussfilt(imG,sigma);
    figure; 
    imshow(imGaus);
    title('Gaussian Filter');
    im = (im);
    %% 2. Step - Sobel ile Kenar Belirleme %%
    
 % Gx ve Gy Sobel Matrisleri vardýr. Bu matrisler ile resimi konvülasyona  sokarýz.
 
 %x ve y eksenindeki yüksek frekans deðiþimlerini yani kenar geçiþlerini
 %bulmamýza yarar.
 
 %Matrisler giriþ görüntüsüne ayrý ayrý uygulanabilir. Böylece her bir yön
 %için pikselin deðeri ayrý ayrý ölçülmüþ olur.
    
 function[A] = SobelFilter(A, filterDirection)
    switch filterDirection
        case 'x' 
            Gx = [-1 0 +1; -2 0 +2; -1 0 +1];
            A = imfilter(A, double(Gx), 'conv', 'replicate');
        case 'y'
            Gy = [-1 -2 -1; 0 0 0; +1 +2 +1];
            A = imfilter(A, double(Gy), 'conv', 'replicate');
        otherwise
            error('Bad filter direction - try inputs ''x'' or ''y''');
    end
end
    
 
 
    Gx = SobelFilter(im, 'x');
    Gy = SobelFilter(im, 'y');
    Gx = imgaussfilt(Gx,sigma);
    Gy = imgaussfilt(Gy,sigma);
    figure; imshow(Gx);
    title('Gx Sobel Filter');
    
    figure; imshow(Gy);
    title('Gy Sobel Filter');
   
%Daha sonra bu deðerler,Gx ve Gy, her bir noktada mutlak büyüklüðü ve yönü bulmak üzere birleþtirilir.
%Resimdeki kenar geçiþleri belirgin hale getirilir.    
    Gmag = sqrt(Gx.^2 + Gy.^2);
    angle = atan2(Gy,Gx)*180/pi;
    figure; imshow(Gmag);
    title('Gmag');
   
    
    
    
    %% 3. Step - Non-maximum suppression

    %Oluþturulan görüntü büyük kalýn kenarlara neden olur. Ýdeal olarak, son görüntünün kenarlarý ince olmalýdýr. 
%Bu nedenle, kenarlarý inceltmek için maksimum olmayan baskýlama yapmalýyýz.
%Kenar yönlerinde maksimum deðeri olan pikselleri bulunur.Diðer yan pikseldeðerleri 0'a atýlýr.

   
    [h,w] = size(im);
    X=[-1,0,+1 ;-1,0,+1 ;-1,0,+1];
	Y=[-1,-1,-1 ;0,0,0 ;+1,+1,+1];
    output = zeros(h,w);
    x = [0 1];
    for i=2:h-1 % row
        for j=2:w-1 % col
            if (angle(i,j)>=0 && angle(i,j)<=45) || ...
                    (angle(i,j)<-135 && angle(i,j)>=-180)
                yBot = [Gmag(i,j+1) Gmag(i+1,j+1)];
                yTop = [Gmag(i,j-1) Gmag(i-1,j-1)];
                x_est = abs(Gy(i,j)/Gmag(i,j)); % y
                if (Gmag(i,j) >= ((yBot(2)-yBot(1))*x_est+yBot(1)) && ...
                    Gmag(i,j) >= ((yTop(2)-yTop(1))*x_est+yTop(1))) % interpolation
                    output(i,j)= Gmag(i,j);
                else
                    output(i,j)=0;
                end
            elseif (angle(i,j)>45 && angle(i,j)<=90) || ...
                    (angle(i,j)<-90 && angle(i,j)>=-135)
                yBot = [Gmag(i+1,j) Gmag(i+1,j+1)];
                yTop = [Gmag(i-1,j) Gmag(i-1,j-1)];
                x_est = abs(Gx(i,j)/Gmag(i,j));
                if (Gmag(i,j) >= ((yBot(2)-yBot(1))*x_est+yBot(1)) && ...
                    Gmag(i,j) >= ((yTop(2)-yTop(1))*x_est+yTop(1)))
                    output(i,j)= Gmag(i,j);
                else
                    output(i,j)=0;
                end
            elseif (angle(i,j)>90 && angle(i,j)<=135) || ...
                    (angle(i,j)<-45 && angle(i,j)>=-90)
                yBot = [Gmag(i+1,j) Gmag(i+1,j-1)];
                yTop = [Gmag(i-1,j) Gmag(i-1,j+1)];
                x_est = abs(Gx(i,j)/Gmag(i,j));
                if (Gmag(i,j) >= ((yBot(2)-yBot(1))*x_est+yBot(1)) && ...
                    Gmag(i,j) >= ((yTop(2)-yTop(1))*x_est+yTop(1)))
                    output(i,j)= Gmag(i,j);
                else
                    output(i,j)=0;
                end
            elseif (angle(i,j)>135 && angle(i,j)<=180) || ...
                    (angle(i,j)<0 && angle(i,j)>=-45)
                yBot = [Gmag(i,j-1) Gmag(i+1,j-1)];
                yTop = [Gmag(i,j+1) Gmag(i-1,j+1)];
                x_est = abs(Gx(i,j)/Gmag(i,j));
                if (Gmag(i,j) >= ((yBot(2)-yBot(1))*x_est+yBot(1)) && ...
                    Gmag(i,j) >= ((yTop(2)-yTop(1))*x_est+yTop(1)))
                    output(i,j)= Gmag(i,j);
                else
                    output(i,j)=0;
                end
            end           
        end
    end
    
    Gmag = output;
    figure; imshow(Gmag);
    title('Non Maximum Suppression');
    %% 4. Step - Double Thresholding
 %Çift eþik adýmý, 3 tür pikselin belirlenmesini amaçlamaktadýr: güçlü, zayýf ve konuyla ilgili olmayan.
 %Güçlü pikseller yoðunluðu o kadar yüksek olan piksellerdir; son kenara katkýda bulunduklarýndan eminiz.
 %Zayýf pikseller, yoðunluk deðeri olarak kabul edilmek için yeterli olmayan, 
 %ancak kenar saptama için uygun olmadýðý düþünülebilecek kadar küçük olmayan bir yoðunluk deðerine sahip piksellerdir.
 %Diðer pikseller kenar için uygun deðildir. 
 
 %Bu adýmdaki amaç, kenarlar ile alakasý olmayan pixsel deðerlerini 0'a göndermektir.
 
 
    highThreshold = max(max(Gmag))*highThresholdRatio;
    lowThreshold = highThreshold*lowThresholdRatio;
    strongEdgesRow = zeros(1,h*w); % Keep track of the strong edge row index
    strongEdgesCol = zeros(1,h*w); % Keep track of the strong edge col index
    weakEdgesRow = zeros(1,h*w);  % Keep track of the weak edge row index
    weakEdgesCol = zeros(1,h*w);  % Keep track of the weak edge col index
    strongIndex = 1;
    weakIndex = 1;
    for i=2:h-1 % Satýr
        for j=2:w-1 % Sütun
            if Gmag(i,j) > highThreshold    % Güçlü kenar
                Gmag(i,j) = 1;
                strongEdgesRow(strongIndex) = i;
                strongEdgesCol(strongIndex) = j;
                strongIndex = strongIndex + 1;
            elseif Gmag(i,j) < lowThreshold % Kenarla alakasý olmayan pixel,
                Gmag(i,j) = 0;   % Kenarla alakasý olmayanlar silinir.
            else                            % Zayýf kenar
                weakEdgesRow(weakIndex) = i;
                weakEdgesCol(weakIndex) = j;
                weakIndex = weakIndex + 1;
            end
        end
    end
    figure; imshow(Gmag);
    title('Double Threshold'); 
 
    %% 5. Step - Edge tracking by hysteresis:
%Güçlü kenarlarýn ve zayýf kenarlarýn ne olduðunu belirlediðimize göre, hangi zayýf kenarlarýn 
%gerçek kenarlar olduðunu belirlememiz gerekiyor. 
%Bunu yapmak için bir kenar izleme algoritmasý gerçekleþtiriyoruz. 
%Güçlü kenarlara baðlý zayýf kenarlar gerçek / gerçek kenarlar olacaktýr. 
%Güçlü kenarlara baðlý olmayan zayýf kenarlar kaldýrýlacaktýr. 
%Bu iþlemi hýzlandýrmak için, algoritmam zayýf ve güçlü kenarlarý izler,
%böylece güçlü kenarlarý tekrar tekrar yineleyebilir ve
%görüntüdeki her piksel boyunca yineleme yapmak yerine baðlý zayýf kenarlar olup olmadýðýný görebilirim    
 
    set(0,'RecursionLimit',10000)
    for i=1:strongIndex-1
       
        Gmag = FindConnectedWeakEdges(Gmag, strongEdgesRow(i),...
            strongEdgesCol(i));
    end
    figure; imshow(Gmag);
    title('Edge Tracking Before Clean Up'); 
  
    function[Gmag] = FindConnectedWeakEdges(Gmag, row, col)
    for i = -3:1:3
        for j = -3:1:3
            if (row+i > 0) && (col+j > 0) && (row+i < size(Gmag,1)) && ...
                    (col+j < size(Gmag,2)) % Make sure we are not out of bounds
                if (Gmag(row+i,col+j) > 0) && (Gmag(row+i,col+j) < 1)
                    Gmag(row+i,col+j) = 1;
                    Gmag = FindConnectedWeakEdges(Gmag, row+i, col+j);
                end
            end
        end
    end
    end


    for i=1:weakIndex-1
        if Gmag(weakEdgesRow(i),weakEdgesCol(i)) ~= 1
            Gmag(weakEdgesRow(i),weakEdgesCol(i)) = 0;
        end
    end
    figure; imshow(Gmag);
    title('Edge Tracking After Clean Up'); 
    
    
end




