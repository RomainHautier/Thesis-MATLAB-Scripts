function [x,y,CIRC] = gammaplotf(speed,conditions,algo,filepath)
    
    %% Defining the azimuthal segmentation of the rotor disk.
    psi = deg2rad(0:10:360);
    np = 36; % Number of aziuthal division
    ns = 32;
    
    %% Defining the filepaths to the standalone rotor cases FWG_b files.
    
    %% Extracting Circulation Data from the FWG_b.DAT files
    for i = 1:length(conditions)
        
        % In case the data to be plotted is that of the standalone rotor
        % cases.
        if conditions{1} == 'MR'
            % Loading the data from the standalone rotor case filepaths.
            FWG_b = load([filepath conditions{i} '\t' num2str(speed) '\FWG_b.dat']);
            
            k=0;
            for p = 1:np
                for s = 1:ns
                    k = k+1;
                    LOCR(s,p) = FWG_b(k,1);
                    CIRC(s,p,i) = FWG_b(k,2);
                end   
            end
    
            r = LOCR(:,1);
            [th,r]=meshgrid(psi,r);
            [x(:,:,i),y(:,:,i)]=pol2cart(th,r);
            CIRC(:,37,i) = CIRC(:,1,i);

        else % If the 2 rotor cases are wanted.
            FWG_b = load([filepath '\d' num2str(speed) '\' conditions{i} '\FWG_b.dat']);
            
            k = 0;
            for r = 1:2
                for p = 1:36
                    for s = 1:32
                       k = k+1;
                       LOCR(s,p) = FWG_b(k,1);
                       CIRC(s,p,r) = FWG_b(k,2);
                    end   
                end
                rad = LOCR(:,1);
                [th,rad]=meshgrid(psi,rad);
                [x(:,:,r),y(:,:,r)]=pol2cart(th,rad);
                CIRC(:,37,r)=CIRC(:,1,r);
            end
        end
        
        % Adding a condition to allow plotting of the gamma plots obtained from
        % HeliUM runs.
        % if algo == 'H' && i == 1
        %     FWG_b = load(['C:\Users\hauti\heliumx2\HeliUM.sample\uh60a_fw_' num2str(speed) 'kt\FWG_B.dat']);
        % end
        
        

    end

    %% Extracting data from the 2 rotor cases.
  
end