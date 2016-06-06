%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [spi1,spi2,spi3,spi4]=firstguess_emi(valori_emi,r,nx,ny)
        
        %compute quadrant emissions
        
        %domain dimensions and other variables
        dimx=nx;
        dimy=ny;
        %prova_emi=valori_emi(:,1);
        emi=reshape(valori_emi,dimy,dimx);
        
        %enlarged emission, with increased zeros all around (used to be able to
        %perform simple products among emi and factor matrix
        emi_enlarged=zeros(size(emi,1)+2*(r-2),size(emi,2)+2*(r-2));
        emi_enlarged(r-2+1:r-2+1+size(emi,1)-1,r-2+1:r-2+1+size(emi,2)-1)=emi;
        [sa sb]=size(emi_enlarged);
        
        %
        %left quadrant
        %factors to perform products of cell emissions
        fattore=tril(ones(r+1,r+1))-diag(repmat(0.5,r+1,1));
        sotto=flipud(fattore);
        fattore=[fattore;sotto(2:end,:)];
        fattore(r+1,r+1)=0.25;
        
        %initialize variable, and set dimensions
        spicchio3=zeros(dimy,dimx);
        
        %down quadrant
        %factors to perform products of cell emissions
        fattoreD=rot90(fattore,3);
        
        %initialize variable, and set dimensions
        spicchio1=zeros(dimy,dimx);
        
        %right quadrant
        %factors to perform products of cell emissions
        fattoreR=rot90(fattore,2);
        
        %initialize variable, and set dimensions
        spicchio4=zeros(dimy,dimx);
        
        %up quadrant
        %factors to perform products of cell emissions
        fattoreU=rot90(fattore,1);
        
        %initialize variable, and set dimensions
        spicchio2=zeros(dimy,dimx);
        
        %loop to aggregate emissions
        for i=r+1:sa-r
            for j=r+1:sb-r
                %left quadrant
                spicchio3(i-r+2,j-r+2)=sum(sum(emi_enlarged(i-r:i+r,j-r:j).*fattore));
                %up quadrant (down as matrix, up geographically speaking)
                spicchio1(i-r+2,j-r+2)=sum(sum(emi_enlarged(i-r:i,j-r:j+r).*fattoreD));
                %right quadrant
                spicchio4(i-r+2,j-r+2)=sum(sum(emi_enlarged(i-r:i+r,j:j+r).*fattoreR));
                %down quadrant (up as matrix, down geographically speaking)
                spicchio2(i-r+2,j-r+2)=sum(sum(emi_enlarged(i:i+r,j-r:j+r).*fattoreU));
            end
        end
        
        spi1=reshape(spicchio1,dimx*dimy,1);
        spi2=reshape(spicchio2,dimx*dimy,1);
        spi3=reshape(spicchio3,dimx*dimy,1);
        spi4=reshape(spicchio4,dimx*dimy,1);
        
    end
