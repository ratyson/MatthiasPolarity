function main
    doPlot = 1;
    showPrompt = 1;
     
    path = 'C:/Users/rt/Documents/Matthais/DATA/';
    avMicronsActin = 0.15;
    
    
    if(showPrompt)
        path = uigetdir(userpath, 'Folder with quimp data...'); 
        prompt = {'Entre proportion of outline for locating the front (0<x<=1)'};
        def = {'0.15'};
        avMicronsActin = inputdlg(prompt,'Proportion of outline',1, def);
        avMicronsActin = str2num(avMicronsActin{1});  
        if(avMicronsActin <= 0 || avMicronsActin > 1),avMicronsActin = 0.15; end
    end
   
    hFig=figure(1);
    h = waitbar(0,'Reading in data...');
    
    % reads in cells
    cells = readQanalysis(path); 
    
    N_cells = length(cells);
    
    out(N_cells).avMicronsActin = [];
    out(N_cells).cellMask = [];
    out(N_cells).frontPixel = [];
    out(N_cells).centreLine = [];
    out(N_cells).rearPixel = [];  
    out(N_cells).poly = [];
    out(N_cells).frontMask = [];
    out(N_cells).rearMask = [];    
    out(N_cells).centre = [];
    out(N_cells).frontPixels = [];
    out(N_cells).rearPixels = [];
    
    
    
    %cell masks
    for i = 2:N_cells,
       c = cells(i);
        
       waitbar(((i/N_cells)-0.05),h,['proccessing cell ',i,' (', c.name, ')']);
       
       
       o.avMicronsActin = avMicronsActin;
       o.cellMask = m_cellMask(c); %entire cell
       o.frontPixel  = cellFront(c, avMicronsActin);
       o.centreLine = m_centreLine(c, o.frontPixel );
       o.rearPixel = m_cellRear(c, o.frontPixel, o.centreLine);
       o.poly = m_polyMask(c, o.frontPixel, o.rearPixel);
       [o.frontMask,o.rearMask] = m_masks(c, o.cellMask , o.poly );
       o.centre = round([c.stats(1,2), c.stats(1,3)]);
       
       o.frontPixels = zeros( sum(o.frontMask(:))  ,3);
       o.frontPixels(:,1) = m_getPixels( c, 1, o.frontMask);  
       o.frontPixels(:,2) = m_getPixels( c, 2, o.frontMask); 
       o.frontPixels(:,3) = m_getPixels( c, 3, o.frontMask);  
       
       o.rearPixels = zeros( sum(o.rearMask(:))  ,3);
       o.rearPixels(:,1) = m_getPixels( c, 1, o.rearMask);  
       o.rearPixels(:,2) = m_getPixels( c, 2, o.rearMask);  
       o.rearPixels(:,3) = m_getPixels( c, 3, o.rearMask); 
   
       writePixels(c, o);
       saveTiff( o.frontMask.*255, [c.PATH, 'z_' ,c.name, '_frontMask' ]);
       saveTiff( o.rearMask.*255, [c.PATH, 'z_' ,c.name, '_rearMask' ]);

       if(doPlot),
           clf(hFig,'reset');
           plotCell(o);
           th=text(10,15, ['Cell ', num2str(i), ' (', c.name, ')']);
           set(th, 'Color', [1,1,1]);
           saveas(gca, [c.PATH, 'z_' ,c.name, '_plot.png' ], 'png'); 
       end
       out(i) = o;
    
    end
    
    dataPath = [cells(1).PATH, '../' ,'analysisData.mat'];
    waitbar(0.95,h,['Save: ', dataPath ]);
    
    save( dataPath, 'cells','out' ); 
    waitfor(msgbox('FINISHED'));
    
    close(h);
    close all;
    
    fprintf('closing...\n');
    
end

function plotCell(o)
        
        mask = double(o.frontMask);
        mask(o.rearMask(:)~=0) = 2;
        hold off
        imagesc(mask);
        hold on
        plot( o.frontPixel(1), o.frontPixel(2)  , 'wx');
        plot( o.rearPixel(1), o.rearPixel(2)  , 'wo');
        plot( o.centreLine( 1:2), o.centreLine(3:4 ), 'w-'  );
        plot( o.centre(1), o.centre(2)  , 'ro');
        axis off
        axis equal
end

function writePixels(c, o)
    fileName = [c.PATH, 'z_' ,c.name, '_sectors.csv' ];
    
    FID = fopen( fileName, 'wt');
    
    fprintf( FID, '');
    
    fprintf( FID, '%s\n', c.PATH);
    fprintf( FID, '%s\n', c.name);
    fprintf( FID, 'proportion av over:,%.6f\n', o.avMicronsActin );
    
    fprintf( FID, 'front x,fornt y,rear x,rear y,centroid x,centroid x\n');
    fprintf( FID, '%d,%d,%d,%d,%d,%d\n\n', o.frontPixel(1),o.frontPixel(2),...
               o.rearPixel(1),o.rearPixel(2), o.centre(1), o.centre(2) );

     fprintf( FID, 'Front pixels:\n');
     fprintf( FID, 'ch1,ch2,ch3\n');
     for i=1:size(o.frontPixels, 1)
        fprintf( FID, '%d,', o.frontPixels(i, :));
         fprintf( FID, '\n');
     end
     
     fprintf( FID, 'Rear pixels:\n');
     fprintf( FID, 'ch1,ch2,ch3\n');
     for i=1:size(o.rearPixels, 1)
        fprintf( FID, '%d,', o.rearPixels(i, :));
         fprintf( FID, '\n');
     end
     
     fclose(FID);
       
end

function pixels = m_getPixels( c, ch, mask)
    N = sum(mask(:)) ;
    
    if(ch==1),
        if(c.FLUOCH1TIFF == '/'),
            %fprintf('no ch1 data\n');
            pixels = zeros(N,1);
            return;
        end
        im = imread(c.FLUOCH1TIFF);
    elseif(ch==2)
        if(c.FLUOCH2TIFF == '/'),
            %fprintf('no ch2 data\n');
           pixels = zeros(N,1);
            return;
        end
        im = imread(c.FLUOCH2TIFF);
    else
        if(c.FLUOCH3TIFF == '/'),
           % fprintf('no ch3 data\n');
           pixels = zeros(N,1);
            return;
        end
        im = imread(c.FLUOCH3TIFF);
    end
    
    im = double(im(:,:,1)) + 1; % stop zeros being lost
    mask = double(mask);
    
    im = im.* mask; % zero out all but mask
    vals = im~=0;
    pixels = im(vals) -1;

end

function mask = m_cellMask(s)
    %mask from segmentation
    outline = s.outlines{end};
    im = imread(s.FLUOCH1TIFF);
    im = im(:,:,1);
    mask = poly2mask( outline(:,2),outline(:,3), size(im,1), size(im,2)); % inbuilt function does outermask from polygon
end

function [front]=cellFront(s, microns)
    
    profile = s.fluoCh1Map(1,:);
    d = size(profile,2);
    
    %figure(4)
    %hold off
   % plot(profile, 'g');
    
    profile = [profile, profile, profile];
    
    % microns
    %perimeter = s.stats(1,8); %in micron
    %mapPixelsToAv = round((microns*d)/perimeter);
    %if(mapPixelsToAv > d), mapPixelsToAv = d/2; end % bound

    %proportion
    mapPixelsToAv = d*microns;
      
    profile = smooth(profile, mapPixelsToAv);
    profile = profile( (d+1):(d+d+1) );
    [~,i]=max(profile);
    
   % hold on
    %plot(profile, 'b')
    %plot( [i,i] , [0,255], 'r');
    
    front = round([ s.xMap(1,i), s.yMap(1,i)  ]);
end

function line = m_centreLine(s, front)
 %   
    L = sqrt(sum([s.R(2)-s.R(1), s.R(4)-s.R(3)].^2));
    
    line = [front(1), s.stats(1,2)  ,front(2), s.stats(1,3)  ]; % [x,x,y,y]
    
    vec = [line(2)-line(1),line(4)-line(3)] ; % unit vector
    L_vec = sqrt((sum(vec.^2))); % length
    uVec = vec./L_vec;
    
    vec = uVec.*L;

    line = [front(1),front(1)+vec(1)  ,front(2), front(2)+vec(2)   ];
    
end

function rear = m_cellRear(s, front, line)

    outline = s.outlines{1};
    N = size(outline,1);
    
    i_p = [2:N,1]'; % i plus one index, circular
    
    l1 = [ line(1), line(3),line(2), line(4)];

    dist = -1; % square distance to intersect
    
    rear = [0,0];
    
    for i = 1:N,
        
       l2 =  [outline(i,2)', outline(i,3)' , outline(i_p(i),2)', outline(i_p(i),3)'];
        
       [intersect, pI]= fastIntersectM( l1 , l2 );
       
       if(intersect)
          distI =   sum([front(1)-pI(1), front(2)-pI(2)].^2);  
          if(distI > dist),
              dist = distI;
              rear = pI;
          end
       end
        
    end
    
    rear = round(rear);

end

function [doIntersect,  o1 ] = fastIntersectM( s1 , s2 )
    % intersect between two lines
    % x1 y1 x2 y2
    o1 = [0,0];
    
    a1 = s1(4) - s1(2); %  y2 - y1;
    b1 = s1(1) - s1(3); %  x1 - x2;
    c1 = s1(3)*s1(2) - s1(1)*s1(4); % x2 * y1 - x1 * y2;

    r3 = a1 * s2(1) + b1 * s2(2) + c1;
    r4 = a1 * s2(3) + b1 * s2(4) + c1;
    
    if ( r3*r4>=0 ) 
       doIntersect = 0;
       return;
    end

	a2 = s2(4) - s2(2);
    b2 = s2(1) - s2(3);
    c2 = s2(3) * s2(2) - s2(1) * s2(4);

    r1 = a2 * s1(1) + b2 * s1(2) + c2;
    r2 = a2 * s1(3) + b2 * s1(4) + c2;

    if ( r1*r2>=0 )
    	doIntersect = 0;
    	return;
    end
    
    denom = a1 * b2 - a2 * b1;
    if ( denom == 0 ) 
    	doIntersect = 0;
    	return; %co-linear
    end

	o1(1) = (b1 * c2 - b2 * c1) / denom;  
    o1(2) =( a2 * c1 - a1 * c2)/denom;

    doIntersect = 1;
    
    
end

function poly = m_polyMask(s, front, rear)
    
    centre = [s.stats(1,2), s.stats(1,3)];
    vec = rear - front;
    polyDim = sqrt(sum([s.R(2)-s.R(1), s.R(4)-s.R(3)].^2)); % max needed size
    
    L_vec = sqrt((sum(vec.^2))); % length
    uVec = vec./L_vec;

    vec90 = [vec(2), -vec(1)]; % x=y, y = -x cockwise 90 deg
    L_vec90 = sqrt((sum(vec90.^2))); % length
    uVec90 = vec90./L_vec90;

    vec = uVec.*polyDim ;
    vec90 = uVec90.*polyDim ;
    
    
    topL = centre - vec90;
    topR = centre + vec90;
    botL = (centre - vec) - vec90;
    botR = (centre - vec) + vec90;
    
    poly = [topL; topR; botR; botL];

end

function [frontMask, rearMask] = m_masks(s, cellMask , poly )
    im = imread(s.FLUOCH1TIFF);
    im = im(:,:,1);
    
    cellMask = uint8(cellMask);
 
    
    % have to iterate to split the cell evenly down the middle
    N = 20;
    Tol = 0.01; % perc
    shift = [0,0];
    for i = 1:N,
        poly(:,1) = poly(:,1) + shift(1); % x
        poly(:,2) = poly(:,2) + shift(2); % y
        
        polyMask = poly2mask( poly(:,1), poly(:,2), size(im,1), size(im,2));
        polyMask = uint8(polyMask);
        polyMask(polyMask==1) = 2;% polyMask(polyMask) + 10; % 2- background, 3- mask
        tMask = polyMask + cellMask;
        
        frontArea = sum( tMask(:) == 3);
        rearArea = sum( tMask(:) == 1);
        
        diff = abs(frontArea - rearArea)/rearArea;
        if(diff < Tol) break; end
        
        %figure(8);
        %imagesc(tMask);
        break
    end
    
    

    frontMask = tMask==3;
    rearMask = tMask==1;
    
    % figure(9);
    %imagesc(frontMask);
    %  figure(10);
    %imagesc(rearMask);  
        
end

function saveTiff( is, varargin )
%SAVETIFF Summary of this function goes here
%   Detailed explanation goes here

    op = size(varargin,2);

    if(op == 0),
        [filename, pathname] = uiputfile('', 'Save tif images as');
        if isequal(filename,0) || isequal(pathname,0)
            disp('Canceled')
        end
        file = fullfile(pathname, filename);
    elseif(op == 1),
        file = varargin{1,1};
    else,
        error('Too many arguments');
    end
    
    
    for i = 1:size(is,3),     
        path = [file, '_', sprintf('%d',i), '.tif'];
        t = Tiff(path, 'w');
        tagStruct.ImageWidth = size(is, 2);
        tagStruct.ImageLength = size(is, 1);
        tagStruct.Photometric = Tiff.Photometric.MinIsBlack;
        tagStruct.Compression = Tiff.Compression.None; 
        tagStruct.BitsPerSample = 8;
        tagStruct.SampleFormat = Tiff.SampleFormat.UInt;  
        tagStruct.PlanarConfiguration = 1;
        t.setTag(tagStruct); 

        t.write(uint8(is(:, : , i)));
        t.close();
        
       % imwrite( uint8(is(:, : , i)),[file, '_', sprintf('%d',i), '.tif'], 'tiff');
    
    end
end
