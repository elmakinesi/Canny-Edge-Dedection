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
   
    %�lk olarak Low pass Filter uyguluyoruz.
    %Resimlerdeki g�r�t�leri azalt�p cisim kenalar�n� daha temiz bir
    %�ekilde bulmam�za yarayacak.
    
    imG = im;
    im = double(imgaussfilt(im,sigma));
    imGaus = imgaussfilt(imG,sigma);
    figure; 
    imshow(imGaus);
    title('Gaussian Filter');
    im = (im);
    %% 2. Step - Sobel ile Kenar Belirleme %%
    
 % Gx ve Gy Sobel Matrisleri vard�r. Bu matrisler ile resimi konv�lasyona  sokar�z.
 
 %x ve y eksenindeki y�ksek frekans de�i�imlerini yani kenar ge�i�lerini
 %bulmam�za yarar.
 
 %Matrisler giri� g�r�nt�s�ne ayr� ayr� uygulanabilir. B�ylece her bir y�n
 %i�in pikselin de�eri ayr� ayr� �l��lm�� olur.
    
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
   
%Daha sonra bu de�erler,Gx ve Gy, her bir noktada mutlak b�y�kl��� ve y�n� bulmak �zere birle�tirilir.
%Resimdeki kenar ge�i�leri belirgin hale getirilir.    
    Gmag = sqrt(Gx.^2 + Gy.^2);
    angle = atan2(Gy,Gx)*180/pi;
    figure; imshow(Gmag);
    title('Gmag');
   
    
    
    
    %% 3. Step - Non-maximum suppression

    %Olu�turulan g�r�nt� b�y�k kal�n kenarlara neden olur. �deal olarak, son g�r�nt�n�n kenarlar� ince olmal�d�r. 
%Bu nedenle, kenarlar� inceltmek i�in maksimum olmayan bask�lama yapmal�y�z.
%Kenar y�nlerinde maksimum de�eri olan pikselleri bulunur.Di�er yan pikselde�erleri 0'a at�l�r.

   
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
 %�ift e�ik ad�m�, 3 t�r pikselin belirlenmesini ama�lamaktad�r: g��l�, zay�f ve konuyla ilgili olmayan.
 %G��l� pikseller yo�unlu�u o kadar y�ksek olan piksellerdir; son kenara katk�da bulunduklar�ndan eminiz.
 %Zay�f pikseller, yo�unluk de�eri olarak kabul edilmek i�in yeterli olmayan, 
 %ancak kenar saptama i�in uygun olmad��� d���n�lebilecek kadar k���k olmayan bir yo�unluk de�erine sahip piksellerdir.
 %Di�er pikseller kenar i�in uygun de�ildir. 
 
 %Bu ad�mdaki ama�, kenarlar ile alakas� olmayan pixsel de�erlerini 0'a g�ndermektir.
 
 
    highThreshold = max(max(Gmag))*highThresholdRatio;
    lowThreshold = highThreshold*lowThresholdRatio;
    strongEdgesRow = zeros(1,h*w); % Keep track of the strong edge row index
    strongEdgesCol = zeros(1,h*w); % Keep track of the strong edge col index
    weakEdgesRow = zeros(1,h*w);  % Keep track of the weak edge row index
    weakEdgesCol = zeros(1,h*w);  % Keep track of the weak edge col index
    strongIndex = 1;
    weakIndex = 1;
    for i=2:h-1 % Sat�r
        for j=2:w-1 % S�tun
            if Gmag(i,j) > highThreshold    % G��l� kenar
                Gmag(i,j) = 1;
                strongEdgesRow(strongIndex) = i;
                strongEdgesCol(strongIndex) = j;
                strongIndex = strongIndex + 1;
            elseif Gmag(i,j) < lowThreshold % Kenarla alakas� olmayan pixel,
                Gmag(i,j) = 0;   % Kenarla alakas� olmayanlar silinir.
            else                            % Zay�f kenar
                weakEdgesRow(weakIndex) = i;
                weakEdgesCol(weakIndex) = j;
                weakIndex = weakIndex + 1;
            end
        end
    end
    figure; imshow(Gmag);
    title('Double Threshold'); 
 
    %% 5. Step - Edge tracking by hysteresis:
%G��l� kenarlar�n ve zay�f kenarlar�n ne oldu�unu belirledi�imize g�re, hangi zay�f kenarlar�n 
%ger�ek kenarlar oldu�unu belirlememiz gerekiyor. 
%Bunu yapmak i�in bir kenar izleme algoritmas� ger�ekle�tiriyoruz. 
%G��l� kenarlara ba�l� zay�f kenarlar ger�ek / ger�ek kenarlar olacakt�r. 
%G��l� kenarlara ba�l� olmayan zay�f kenarlar kald�r�lacakt�r. 
%Bu i�lemi h�zland�rmak i�in, algoritmam zay�f ve g��l� kenarlar� izler,
%b�ylece g��l� kenarlar� tekrar tekrar yineleyebilir ve
%g�r�nt�deki her piksel boyunca yineleme yapmak yerine ba�l� zay�f kenarlar olup olmad���n� g�rebilirim    
 
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




